#!/usr/bin/env bash
#SBATCH --job-name=prep_inputs
#SBATCH --output=prep_inputs_%j.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=shared

set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C

# ---------- clean modules & activate micromamba env ----------
module purge || true
if [[ -f "${HOME}/micromamba/etc/profile.d/micromamba.sh" ]]; then
  # shellcheck disable=SC1091
  source "${HOME}/micromamba/etc/profile.d/micromamba.sh"
  micromamba activate rgenesis || true
fi

# prefer micromamba python if available
PY="${HOME}/micromamba/envs/rgenesis/bin/python"
if [[ ! -x "$PY" ]]; then
  PY="$(command -v python || true)"
fi
if [[ -z "${PY:-}" ]]; then
  echo "[err] could not find python" >&2; exit 2
fi

echo "[diag] python: $("$PY" -V 2>&1 | awk '{print $1,$2}')"

# ---------------- PATHS ----------------
ROOT="/gpfs/group/ehlers/mwei"
WORK="${ROOT}/ABCD_2/ABCD/work_provided_assets"
GENO_QC_DIR="${WORK}/imputed_genotype_qc"
BFILE="${GENO_QC_DIR}/imputed_qc_filtered"

COV="${WORK}/final_covariates.tsv"                                 # FID IID age sex site PC1..PC20 motion scanner_*
ICV="${ROOT}/Tables/core/imaging/icv_baseline.tsv"                 # IID  eTIV (z-scored/resid)
PHENO_EXT="${ROOT}/Tables/core/mental-health/ext_longitudinal_phenos.iid_only.tsv"  # IID, EXT_*
PHENO_IMAGING="${ROOT}/Tables/core/imaging/ABCD_Tier1_long.tsv"    # IID + *_baseline/_mean/_slope

OUTDIR="${WORK}/gwas_pipe"
mkdir -p "$OUTDIR"

echo "[diag] ROOT=$ROOT"
echo "[diag] BFILE=$BFILE"
echo "[diag] COV=$COV"
echo "[diag] ICV=$ICV"
echo "[diag] PHENO_EXT=$PHENO_EXT"
echo "[diag] PHENO_IMAGING=$PHENO_IMAGING"
echo "[diag] OUTDIR=$OUTDIR"

# ---------------- sanity checks ----------------
for f in "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam" "$COV" "$ICV" "$PHENO_EXT" "$PHENO_IMAGING"; do
  [[ -s "$f" ]] || { echo "[err] missing or empty file: $f" >&2; exit 2; }
done

# ---------------- get genotype IIDs ----------------
cut -d' ' -f2 "${BFILE}.fam" | sort -u > "${OUTDIR}/ids.geno"

# ---------------- fix covariate header glitches (ONLY the header) ----------------
# Split fused PC tokens like "PC17PC18" -> "PC17<TAB>PC18" in the HEADER ONLY; leave data rows unchanged.
COV_CLEAN="${OUTDIR}/final_covariates.clean.tsv"
"$PY" - <<'PY' "$COV" "$COV_CLEAN"
import sys, re
src, dst = sys.argv[1], sys.argv[2]
with open(src, 'r', encoding='utf-8', errors='ignore') as f:
    lines = f.readlines()
if not lines:
    raise SystemExit("[clean-cov] empty covariate file")

hdr = lines[0].rstrip('\n')

def split_fused_pcs(header: str) -> str:
    if '\t' in header:
        toks = header.split('\t')
        out = []
        for tok in toks:
            cur = tok
            while True:
                new = re.sub(r'(PC\d+)(PC\d+)', r'\1\t\2', cur)
                if new == cur:
                    break
                cur = new
            out.extend(cur.split('\t'))
        return '\t'.join(out)
    else:
        return re.sub(r'(PC\d+)(PC\d+)', r'\1\t\2', header)

hdr_fixed = split_fused_pcs(hdr)

with open(dst, 'w', encoding='utf-8') as w:
    w.write(hdr_fixed + '\n')
    for line in lines[1:]:
        w.write(line)
PY

# ---------------- validate column shapes on cleaned covariates ----------------
"$PY" - <<'PY' "$COV_CLEAN"
import sys, csv
p = sys.argv[1]
with open(p, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    rows = list(r)
hdr = rows[0]
n = len(hdr)
bad = []
for i, row in enumerate(rows[1:201], start=2):  # check first 200 rows
    if len(row) != n:
        bad.append((i, len(row)))
if bad:
    print("[warn] cov header has", n, "columns; mismatched rows (first 200 lines shown):", bad[:10])
else:
    print("[ok] covariate header/data column counts look consistent (first 200 rows)")
PY

# ---------------- intersect IDs across sources ----------------
cut -f2 "$COV_CLEAN" | tail -n +2 | sort -u > "${OUTDIR}/ids.cov"
cut -f1 "$ICV"        | tail -n +2 | sort -u > "${OUTDIR}/ids.icv"
cut -f1 "$PHENO_EXT"  | tail -n +2 | sort -u > "${OUTDIR}/ids.ext"
cut -f1 "$PHENO_IMAGING" | tail -n +2 | sort -u > "${OUTDIR}/ids.img"

# externalizing: genotype ∩ covariates ∩ externalizing
comm -12 "${OUTDIR}/ids.geno" "${OUTDIR}/ids.cov" | comm -12 - "${OUTDIR}/ids.ext" > "${OUTDIR}/keep.ext"
# imaging: genotype ∩ covariates ∩ imaging ∩ ICV
comm -12 "${OUTDIR}/ids.geno" "${OUTDIR}/ids.cov" | comm -12 - "${OUTDIR}/ids.img" | comm -12 - "${OUTDIR}/ids.icv" > "${OUTDIR}/keep.img"

echo "[diag] n_keep.ext = $(wc -l < ${OUTDIR}/keep.ext)"
echo "[diag] n_keep.img = $(wc -l < ${OUTDIR}/keep.img)"

# ---------------- make covariate tables ----------------
# cov_base: for externalizing (FID IID + covs)
COV_BASE="${OUTDIR}/covariates.base.tsv"
{
  head -n1 "$COV_CLEAN"
  awk 'BEGIN{FS=OFS="\t"} FNR==NR{a[$1]; next} ($2 in a){print}' "${OUTDIR}/keep.ext" <(tail -n +2 "$COV_CLEAN")
} > "$COV_BASE"

# cov_base_img: for imaging (FID IID + covs) before adding ICV
COV_BASE_IMG="${OUTDIR}/covariates.base_img.tsv"
{
  head -n1 "$COV_CLEAN"
  awk 'BEGIN{FS=OFS="\t"} FNR==NR{a[$1]; next} ($2 in a){print}' "${OUTDIR}/keep.img" <(tail -n +2 "$COV_CLEAN")
} > "$COV_BASE_IMG"

# merge ICV to cov_base_img -> covariates.withICV.tsv (join by IID; adds eTIV)
COV_WITHICV="${OUTDIR}/covariates.withICV.tsv"
"$PY" - <<'PY' "$COV_BASE_IMG" "$ICV" "$COV_WITHICV"
import sys, csv
cov, icv, outp = sys.argv[1], sys.argv[2], sys.argv[3]

# read ICV: expect header "IID\teTIV"
icv_map = {}
with open(icv, newline='') as f:
    r = csv.DictReader(f, delimiter='\t')
    if 'IID' not in r.fieldnames or 'eTIV' not in r.fieldnames:
        raise SystemExit("[merge-icv] ICV header must contain IID and eTIV")
    for row in r:
        icv_map[row['IID']] = row['eTIV']

with open(cov, newline='') as f:
    rc = csv.reader(f, delimiter='\t')
    rows = list(rc)

hdr = rows[0]
if not (hdr[:2] == ['FID','IID']):
    raise SystemExit("[merge-icv] covariate header must start with FID, IID")

hdr_out = hdr + ['eTIV']
with open(outp, 'w', newline='') as w:
    cw = csv.writer(w, delimiter='\t', lineterminator='\n')
    cw.writerow(hdr_out)
    for row in rows[1:]:
        # IID is column 2
        iid = row[1]
        cw.writerow(row + [icv_map.get(iid, '')])
PY

# ---------------- phenotype tables for REGENIE ----------------
# REGENIE needs FID IID first, then phenos.

# Externalizing phenos
PHENO_EXT_REGENIE="${OUTDIR}/pheno_externalizing.tsv"
"$PY" - <<'PY' "$PHENO_EXT" "${OUTDIR}/keep.ext" "$PHENO_EXT_REGENIE"
import sys, csv
ext, keep_path, outp = sys.argv[1], sys.argv[2], sys.argv[3]
keep = set(x.strip() for x in open(keep_path))
with open(ext, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    hdr = next(r)
    if not hdr or hdr[0] != 'IID':
        raise SystemExit("[pheno-ext] first column must be IID")
    ph_cols = [c for c in hdr[1:] if c.startswith('EXT_')]
    out_hdr = ['FID','IID'] + ph_cols
    rows = []
    for row in r:
        iid = row[0]
        if iid in keep:
            d = dict(zip(hdr, row))
            rows.append([iid, iid] + [d.get(c, '') for c in ph_cols])
with open(outp, 'w', newline='') as w:
    cw = csv.writer(w, delimiter='\t', lineterminator='\n')
    cw.writerow(out_hdr)
    cw.writerows(rows)
PY

# Imaging phenos
PHENO_IMG_REGENIE="${OUTDIR}/pheno_imaging.tsv"
"$PY" - <<'PY' "$PHENO_IMAGING" "${OUTDIR}/keep.img" "$PHENO_IMG_REGENIE"
import sys, csv
img, keep_path, outp = sys.argv[1], sys.argv[2], sys.argv[3]
keep = set(x.strip() for x in open(keep_path))
with open(img, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    hdr = next(r)
    if not hdr or hdr[0] != 'IID':
        raise SystemExit("[pheno-img] first column must be IID")
    ph_cols = [c for c in hdr[1:] if c.endswith(('_baseline','_mean','_slope'))]
    out_hdr = ['FID','IID'] + ph_cols
    rows = []
    for row in r:
        iid = row[0]
        if iid in keep:
            d = dict(zip(hdr, row))
            rows.append([iid, iid] + [d.get(c, '') for c in ph_cols])
with open(outp, 'w', newline='') as w:
    cw = csv.writer(w, delimiter='\t', lineterminator='\n')
    cw.writerow(out_hdr)
    cw.writerows(rows)
PY

# ---------------- phenotype column lists ----------------
# Externalizing:
cut -f3 "$PHENO_EXT_REGENIE" | tr '\t' '\n' > "${OUTDIR}/ext_cols.txt"

# Imaging (lists + SA/TH helpers):
head -n1 "$PHENO_IMG_REGENIE" | cut -f3- | tr '\t' '\n' > "${OUTDIR}/img_all_cols.txt"
grep '_baseline$' "${OUTDIR}/img_all_cols.txt" > "${OUTDIR}/img_cols_baseline.txt" || true
grep '_mean$'     "${OUTDIR}/img_all_cols.txt" > "${OUTDIR}/img_cols_mean.txt"     || true
grep '_slope$'    "${OUTDIR}/img_all_cols.txt" > "${OUTDIR}/img_cols_slope.txt"    || true
grep '_SA_'       "${OUTDIR}/img_all_cols.txt" > "${OUTDIR}/img_cols_SA.txt"       || true
grep '_TH_'       "${OUTDIR}/img_all_cols.txt" > "${OUTDIR}/img_cols_TH.txt"       || true

# ---------------- export handy pointers ----------------
cat > "${OUTDIR}/_paths.env" <<EOF
ROOT=${ROOT}
WORK=${WORK}
BFILE=${BFILE}

# covariates
COV_BASE=${COV_BASE}
COV_BASE_IMG=${COV_BASE_IMG}
COV_WITHICV=${COV_WITHICV}

# phenos for regenie
PHENO_EXT_REGENIE=${PHENO_EXT_REGENIE}
PHENO_IMG_REGENIE=${PHENO_IMG_REGENIE}

# columns
EXT_COLS=${OUTDIR}/ext_cols.txt
IMG_COLS_ALL=${OUTDIR}/img_all_cols.txt
IMG_COLS_BASELINE=${OUTDIR}/img_cols_baseline.txt
IMG_COLS_MEAN=${OUTDIR}/img_cols_mean.txt
IMG_COLS_SLOPE=${OUTDIR}/img_cols_slope.txt
IMG_COLS_SA=${OUTDIR}/img_cols_SA.txt
IMG_COLS_TH=${OUTDIR}/img_cols_TH.txt

# keep lists
KEEP_EXT=${OUTDIR}/keep.ext
KEEP_IMG=${OUTDIR}/keep.img

# out root
OUTDIR=${OUTDIR}
EOF

echo "[done] Inputs prepared in $OUTDIR"

