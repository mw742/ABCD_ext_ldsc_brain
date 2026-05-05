#!/usr/bin/env bash
#SBATCH --job-name=rg_ext_vs_img_pooled
#SBATCH --output=rg_pooled_%A.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=400:00:00
#SBATCH --partition=highmem

set -euo pipefail
IFS=$'\n\t'

echo "[$(date)] Starting rg_ext_vs_img_pooled on host: $(hostname)"

module purge
module load python/2.7.11
module load ldsc/1.0.0

# ---------- Stable Python2/OpenSSL runtime ----------
LIBFIX="/gpfs/home/mwei/libfix"
PY2_BIN="/opt/applications/python/2.7.11/gnu/bin/python2.7"
HASHLIB_SO="/opt/applications/python/2.7.11/gnu/lib/python2.7/lib-dynload/_hashlib.so"

export LD_LIBRARY_PATH="${LIBFIX}:/usr/lib64:/lib64:/opt/applications/python/2.7.11/gnu/lib"
export PYTHONNOUSERSITE=1
unset PYTHONPATH

[[ -x "$PY2_BIN" ]] || { echo "[fatal] python2 binary not found: $PY2_BIN"; exit 1; }
[[ -r "${LIBFIX}/libssl.so.10" ]] || { echo "[fatal] missing ${LIBFIX}/libssl.so.10"; exit 1; }
[[ -r "${LIBFIX}/libcrypto.so.10" ]] || { echo "[fatal] missing ${LIBFIX}/libcrypto.so.10"; exit 1; }

# Wrap every Python2 call with the fixed library path
PY2RUN=(env LD_LIBRARY_PATH="${LIBFIX}:/usr/lib64:/lib64:/opt/applications/python/2.7.11/gnu/lib" "$PY2_BIN" -S)

echo "[diag] hostname=$(hostname)"
echo "[diag] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "[diag] PATH=$PATH"
echo "[diag] which python2=$(command -v python2 || true)"
echo "[diag] PY2_BIN=$PY2_BIN"
echo "[diag] HASHLIB_SO=$HASHLIB_SO"

echo "[diag] checking private OpenSSL compat files"
ls -l "${LIBFIX}/libssl.so.10" \
      "${LIBFIX}/libcrypto.so.10" \
      "${LIBFIX}/libssl.so.1.0.2k" \
      "${LIBFIX}/libcrypto.so.1.0.2k"

echo "[diag] ldd on python2"
ldd "$PY2_BIN" || true

echo "[diag] ldd on _hashlib.so with fixed env"
env LD_LIBRARY_PATH="${LIBFIX}:/usr/lib64:/lib64:/opt/applications/python/2.7.11/gnu/lib" \
    ldd "$HASHLIB_SO" || true

echo "[diag] testing hashlib"
"${PY2RUN[@]}" - <<'PY'
import _hashlib, hashlib
print("OK")
print(hashlib.md5("abc").hexdigest())
PY

# Prefer module-provided binaries
LDSC_BIN="$(command -v ldsc.py || true)"
MUNGE_BIN="$(command -v munge_sumstats.py || true)"
[[ -n "${LDSC_BIN}" ]] || LDSC_BIN="/opt/applications/ldsc/1.0.0/gnu/bin/ldsc.py"
[[ -n "${MUNGE_BIN}" ]] || MUNGE_BIN="/opt/applications/ldsc/1.0.0/gnu/bin/munge_sumstats.py"

# ---------- Clean Python2 env for LDSC ----------
export PYTHONNOUSERSITE=1
unset PYTHONPATH
for p in \
  /opt/applications/ldsc/1.0.0/gnu/lib/python2.7/site-packages \
  /opt/applications/ldsc/1.0.0/gnu/lib64/python2.7/site-packages
do
  [[ -d "$p" ]] && export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$p"
done

echo "[diag] sanity check Python2 deps"
"${PY2RUN[@]}" - <<'PY'
import sys, hashlib
hp = getattr(hashlib, '__file__', '')
if hp and ('site-packages' in hp or 'egg' in hp):
    sys.stderr.write("[fatal] Bad hashlib on sys.path: %s\n" % hp)
    sys.stderr.write("\n".join(sys.path) + "\n")
    raise SystemExit(1)
import numpy, scipy, pandas, ldscore
print("py2 deps OK")
PY

# ---------------- Paths ----------------
ROOT="/gpfs/group/ehlers/mwei/ABCD_2/ABCD"
WORK="${ROOT}/work_provided_assets"

EXT_DIR="${WORK}/gwas_output/externalizing"
IMG_DIR="${WORK}/gwas_output/imaging"

OUTTOP="${WORK}/rg_ldsc_pooled"
MAP_DIR="${OUTTOP}/_maps"
MUNG_EXT="${OUTTOP}/munged_ext"
MUNG_IMG="${OUTTOP}/munged_img"
RG_OUT="${OUTTOP}/rg_results"
mkdir -p "${OUTTOP}" "${MAP_DIR}" "${MUNG_EXT}" "${MUNG_IMG}" "${RG_OUT}"

# RSID map: CHR BP A1 A2 RSID (gz)
RSID_MAP="${MAP_DIR}/chrpos_a1a2_to_rsid.tsv.gz"
[[ -s "${RSID_MAP}" ]] || { echo "[fatal] Missing RSID map: ${RSID_MAP}"; exit 1; }

# ---------------- Use precomputed LD scores ----------------
export REF_LD_CHR="${REF_LD_CHR:-${WORK}/ldsc_ref_pooled_custom/pooled.custom.}"
export W_LD_CHR="${W_LD_CHR:-${WORK}/ldsc_ref_pooled_custom/pooled.custom.}"
for f in "${REF_LD_CHR}1.l2.ldscore.gz" "${W_LD_CHR}1.l2.ldscore.gz"; do
  [[ -s "$f" ]] || { echo "[fatal] Custom LD score missing: $f"; exit 1; }
done

# ---------------- Merge-alleles panel ----------------
PANEL_DIR="${WORK}/ldsc_panel_pooled"
MERGE_ALLELES="${MERGE_ALLELES:-${WORK}/ldsc_ref_pooled_custom/panel.snplist}"
if [[ ! -s "${MERGE_ALLELES}" ]]; then
  echo "[info] Building merge-alleles SNP list from ${PANEL_DIR} -> ${MERGE_ALLELES}"
  mkdir -p "$(dirname "${MERGE_ALLELES}")"
  awk '{print $2}' "${PANEL_DIR}"/pooled_unrel_chr{1..22}.bim 2>/dev/null | sort -u > "${MERGE_ALLELES}" || true
fi
[[ -s "${MERGE_ALLELES}" ]] && echo "[diag] Using --merge-alleles=${MERGE_ALLELES}"

# ---------------- Compression ----------------
if command -v bgzip >/dev/null 2>&1; then
  GZIP_CMD="bgzip -c"
else
  echo "[warn] bgzip not found; using gzip -c"
  GZIP_CMD="gzip -c"
fi

# ---------------- Params ----------------
N_SAMPLES="${N_SAMPLES:-11666}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---------------- Helpers ----------------
prep_for_ldsc() {
  local in="$1" rsmap="$2" out="$3"
  local tmp_dir="${SLURM_TMPDIR:-/tmp}"
  local tmp_norm="${tmp_dir}/norm.$$.$RANDOM.tsv"
  local tmp_map="${tmp_dir}/rsmap.$$.$RANDOM.tsv"

  echo "[prep] $(date) prep_for_ldsc in=${in} out=${out}"

  # rsmap: CHR BP A1 A2 RSID
  zcat "$rsmap" | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > "$tmp_map"

  # Normalize headers and build intermediate with ORIGID, A1, A2, P, BETA, SE, N, CHR, BP, Z
  zcat "$in" \
    | sed '1s/^\xEF\xBB\xBF//' | tr -d '\r' | sed -E 's/[[:space:]]+/\t/g' \
    | awk -vOFS='\t' '
      function toNum(x){return (x==""||x=="NA"||x=="NaN")?"NA":x}
      function pow10(x){return (x=="NA")?"NA":exp(-x*log(10))}
      function zval(b,se){if(b=="NA"||se=="NA"||se==0){return "NA"}; return b/se}
      NR==1{
        for(i=1;i<=NF;i++){
          h=toupper($i)
          gsub(/^CHR$|^CHROMOSOME$|^CHROM$/,"CHROM",h)
          gsub(/^BP$|^POS$|^POSITION$|^GENPOS$/,"GENPOS",h)
          gsub(/^RSID$|^SNPID$|^SNP$/,"ID",h)
          gsub(/^ALLELE1$|^A1_EFF$|^ALT$/,"ALLELE1",h)
          gsub(/^ALLELE0$|^A2_OTHER$|^REF$/,"ALLELE0",h)
          gsub(/^LOG10P$|^LOG10_P$|^LP$/,"LOG10P",h)
          H[i]=h; idx[h]=i
        }
        print "ORIGID","A1","A2","P","BETA","SE","N","CHR","BP","Z"
        next
      }
      {
        id   = idx["ID"]      ? $(idx["ID"])      : "NA"
        a1   = idx["ALLELE1"] ? $(idx["ALLELE1"]) : "NA"
        a2   = idx["ALLELE0"] ? $(idx["ALLELE0"]) : "NA"
        p    = (("P" in idx)  ? $(idx["P"])       : "NA")
        lp   = idx["LOG10P"]  ? $(idx["LOG10P"])  : "NA"
        beta = idx["BETA"]    ? $(idx["BETA"])    : "NA"
        se   = idx["SE"]      ? $(idx["SE"])      : "NA"
        nval = idx["N"]       ? $(idx["N"])       : "NA"
        chr  = idx["CHROM"]   ? $(idx["CHROM"])   : (("CHR" in idx)?$(idx["CHR"]):"NA")
        bp   = idx["GENPOS"]  ? $(idx["GENPOS"])  : (("BP" in idx)?$(idx["BP"]):"NA")
        if(p==""||p=="NA"){ p = pow10(toNum(lp)) }
        if(chr ~ /^chr/){ sub(/^chr/,"",chr) }
        z = zval(toNum(beta), toNum(se))
        print id,a1,a2,p,beta,se,nval,chr,bp,z
      }
    ' > "$tmp_norm"

  # Join rsmap → final LDSC-ready file
  awk -vOFS='\t' '
    BEGIN{
      while((getline < ARGV[1])>0){ key=$1":"$2":"$3":"$4; m[key]=$5 }
      close(ARGV[1]); ARGV[1]=""
      print "SNP","A1","A2","P","BETA","SE","N","CHR","BP","Z"
    }
    NR>1{
      orig=$1; a1=$2; a2=$3; p=$4; b=$5; s=$6; n=$7; c=$8; bp=$9; z=$10
      key=c":"bp":"a1":"a2
      snp = (orig ~ /^rs[0-9]+$/ ? orig : (key in m ? m[key] : key))
      print snp,a1,a2,p,b,s,n,c,bp,z
    }
  ' "$tmp_map" "$tmp_norm" | eval "$GZIP_CMD" > "$out"

  echo -n "[diag] header -> "; zcat "$out" | head -1
  rm -f "$tmp_norm" "$tmp_map"
}

munge_one() {
  local label="$1" in="$2" outdir="$3"
  local prefix="${outdir}/${label}"
  local with="${prefix}.with_rsid.tsv.gz"
  local out="${prefix}.sumstats.gz"

  if [[ -s "$out" ]]; then
    echo "[skip] ${out} already exists; not remunging ${label}"
    return 0
  fi

  echo "[diag] munge ${label} -> ${prefix}.sumstats.gz"
  prep_for_ldsc "$in" "$RSID_MAP" "$with" || { echo "[fatal] prep_for_ldsc failed: $in"; exit 1; }

  if [[ -s "${MERGE_ALLELES}" ]]; then
    "${PY2RUN[@]}" "${MUNGE_BIN}" \
      --signed-sumstats Z,0 \
      --out "${prefix}" \
      --N "${N_SAMPLES}" \
      --a1 A1 --a2 A2 --snp SNP --p P \
      --merge-alleles "${MERGE_ALLELES}" \
      --sumstats "${with}"
  else
    "${PY2RUN[@]}" "${MUNGE_BIN}" \
      --signed-sumstats Z,0 \
      --out "${prefix}" \
      --N "${N_SAMPLES}" \
      --a1 A1 --a2 A2 --snp SNP --p P \
      --sumstats "${with}"
  fi
}

run_rg() {
  local img="$1" ext="$2"
  local out_prefix="${RG_OUT}/${img}__vs__${ext}"
  local out_log="${out_prefix}.log"
  local s1="${MUNG_IMG}/${img}.sumstats.gz"
  local s2="${MUNG_EXT}/${ext}.sumstats.gz"

  if [[ ! -s "$s1" || ! -s "$s2" ]]; then
    echo "[warn] skip rg (missing sumstats): img=${img} ext=${ext}"
    return 0
  fi

  if [[ -s "$out_log" ]]; then
    echo "[skip rg] ${out_log} already exists; img=${img} ext=${ext}"
    return 0
  fi

  echo "[rg] $(date) ${img} vs ${ext}"
  "${PY2RUN[@]}" "${LDSC_BIN}" \
    --rg "${s1},${s2}" \
    --ref-ld-chr "${REF_LD_CHR}" \
    --w-ld-chr "${W_LD_CHR}" \
    --out "${out_prefix}"
}

# ---------------- Discover inputs ----------------
declare -a EXT_BUNCH=()
for b in EXT_baseline EXT_mean EXT_slope; do
  f="${EXT_DIR}/${b}.regenie.gz"
  if [[ -s "$f" ]]; then
    EXT_BUNCH+=("$b")
  else
    echo "[warn] missing EXT file: $f"
  fi
done
((${#EXT_BUNCH[@]})) || { echo "[fatal] No externalizing sumstats found under ${EXT_DIR}"; exit 1; }

# Optionally restrict to a subset of EXT labels via env var EXT_FILTER (comma-separated)
if [[ "${EXT_FILTER:-}" != "" ]]; then
  IFS=',' read -r -a FILTERED <<< "${EXT_FILTER}"
  tmp=()
  for b in "${EXT_BUNCH[@]}"; do
    for f in "${FILTERED[@]}"; do
      if [[ "$b" == "$f" ]]; then tmp+=("$b"); fi
    done
  done
  EXT_BUNCH=("${tmp[@]}")
  echo "[info] Using EXT_FILTER; EXT_BUNCH = ${EXT_BUNCH[*]}"
fi

# Find imaging regenies
mapfile -t IMG_FILES < <(find "${IMG_DIR}" -maxdepth 1 -type f -name "*.regenie.gz" | sort)
((${#IMG_FILES[@]})) || { echo "[fatal] No imaging sumstats found under ${IMG_DIR}"; exit 1; }

# Optional: restrict imaging files via IMG_FILTER (substring match, comma separated)
if [[ "${IMG_FILTER:-}" != "" ]]; then
  IFS=',' read -r -a PATTERNS <<< "${IMG_FILTER}"
  tmp=()
  for path in "${IMG_FILES[@]}"; do
    base="$(basename "$path")"
    for p in "${PATTERNS[@]}"; do
      if [[ "$base" == *"$p"* ]]; then tmp+=("$path"); break; fi
    done
  done
  IMG_FILES=("${tmp[@]}")
  echo "[info] Using IMG_FILTER; number of IMG files = ${#IMG_FILES[@]}"
fi

echo "[info] Externalizing phenotypes: ${EXT_BUNCH[*]}"
echo "[info] Number of imaging files: ${#IMG_FILES[@]}"

# ---------------- MUNGE (EXT) ----------------
echo "[$(date)] Munging externalizing sumstats (resumable)"
for b in "${EXT_BUNCH[@]}"; do
  in="${EXT_DIR}/${b}.regenie.gz"
  out="${MUNG_EXT}/${b}.sumstats.gz"
  if [[ -s "$out" ]]; then
    echo "[skip] EXT ${b}: ${out} already present"
  else
    munge_one "${b}" "${in}" "${MUNG_EXT}"
  fi
done

# ---------------- MUNGE (IMG) ----------------
echo "[$(date)] Munging imaging sumstats (resumable)"
for path in "${IMG_FILES[@]}"; do
  label="$(basename "${path%.regenie.gz}")"
  out="${MUNG_IMG}/${label}.sumstats.gz"
  if [[ -s "$out" ]]; then
    echo "[skip] IMG ${label}: ${out} already present"
  else
    munge_one "${label}" "${path}" "${MUNG_IMG}"
  fi
done

# ---------------- RG: IMG × EXT ----------------
echo "[$(date)] Running LDSC rg for IMG × EXT (resumable)"
for path in "${IMG_FILES[@]}"; do
  img_label="$(basename "${path%.regenie.gz}")"
  for ext_label in "${EXT_BUNCH[@]}"; do
    run_rg "${img_label}" "${ext_label}"
  done
done

echo "[$(date)] [done] All munging + rg outputs -> ${OUTTOP}"

