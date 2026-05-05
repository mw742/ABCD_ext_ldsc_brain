#!/usr/bin/env bash
#SBATCH --job-name=mk_unrel_pooled
#SBATCH --output=mk_unrel_pooled_%j.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=shared

# ----------------------- SAFETY -----------------------
set -euo pipefail
IFS=$'\n\t'
trap 'echo "[fatal] line $LINENO failed" >&2; exit 2' ERR

# ----------------------- USER PATHS -------------------
# Preferred PLINK 2 (has --king-cutoff); fallback PLINK 1.9
PLINK2="/opt/applications/plink/2.00a5.12/bin/plink2"
PLINK1="/opt/applications/plink/1.90b7.2/bin/plink"

ROOT="/gpfs/group/ehlers/mwei/ABCD_2/ABCD"
WORK="${ROOT}/work_provided_assets"
BFILE="${WORK}/imputed_genotype_qc/imputed_qc_filtered"

# One-column lists 
KEEP_ALL_1C="${WORK}/ldsc_panel_pooled/pooled_all_iids.txt"
GENESIS_UNREL_1C="${ROOT}/genomics_sample03/ABCD_GRM/GENESIS/ABCD_202209.updated.nodups.curated_unrelateds.txt"

# Output dir + expected 2-col keeps
OUTDIR="${WORK}/ldsc_panel_pooled"
KEEP_GENESIS_2C="${OUTDIR}/genesis_unrelateds.2col.txt"
KEEP_ALL_2C="${OUTDIR}/pooled_all_iids.2col.txt"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
KING_CUTOFF="0.0884"   # KING kinship ~ 2nd-degree
PIHAT_CUTOFF="0.177"   # ~ 2 * KING kinship as fallback for --rel-cutoff

# ----------------------- HELPERS ----------------------
use_plink2=0
if [[ -x "$PLINK2" ]]; then
  # Try a quick no-op to ensure it runs on this CPU
  if "$PLINK2" --help >/dev/null 2>&1; then
    use_plink2=1
  fi
fi

echo "[diag] Using PLINK2? -> $use_plink2"
if (( use_plink2 )); then
  echo "[diag] PLINK2 at: $PLINK2"
else
  echo "[diag] Falling back to PLINK1 at: $PLINK1"
fi

mkdir -p "${OUTDIR}"

# --- Build/choose KEEP_2COL (FID IID) ---
if [[ -s "$GENESIS_UNREL_1C" && ! -s "$KEEP_GENESIS_2C" ]]; then
  echo "[diag] Converting GENESIS unrelateds (1-col) -> 2-col keep"
  awk 'NR==FNR{iid[$1]=1;next} {if($2 in iid) print $1, $2}' \
    "$GENESIS_UNREL_1C" "${BFILE}.fam" > "$KEEP_GENESIS_2C"
fi
if [[ -s "$KEEP_ALL_1C" && ! -s "$KEEP_ALL_2C" ]]; then
  echo "[diag] Converting pooled_all_iids (1-col) -> 2-col keep"
  awk 'NR==FNR{iid[$1]=1;next} {if($2 in iid) print $1, $2}' \
    "$KEEP_ALL_1C" "${BFILE}.fam" > "$KEEP_ALL_2C"
fi

# Prefer GENESIS unrelateds; else pooled_all_iids; else everyone in fam
if [[ -s "$KEEP_GENESIS_2C" ]]; then
  KEEP_2COL="$KEEP_GENESIS_2C"
  echo "[diag] Using GENESIS curated unrelateds: $KEEP_2COL"
elif [[ -s "$KEEP_ALL_2C" ]]; then
  KEEP_2COL="$KEEP_ALL_2C"
  echo "[diag] Using pooled_all_iids.2col: $KEEP_2COL"
else
  KEEP_2COL="${OUTDIR}/everyone_from_fam.2col.txt"
  echo "[diag] No keep lists provided; using everyone from .fam -> $KEEP_2COL"
  awk '{print $1, $2}' "${BFILE}.fam" > "$KEEP_2COL"
fi

head -3 "$KEEP_2COL" >&2 || true
wc -l  "$KEEP_2COL" >&2  || true

echo "[step] Restricting to keep IIDs"
"$PLINK1" --bfile "$BFILE" \
  --keep "$KEEP_2COL" \
  --threads "$THREADS" \
  --make-bed \
  --out "${OUTDIR}/pooled_keep"

# --- Safety-net unrelated prune ---
# If PLINK2 is available, use KING cutoff; else approximate with PI_HAT cutoff.
if (( use_plink2 )); then
  echo "[step] KING unrelated via PLINK2 (kinship cutoff ${KING_CUTOFF})"
  "$PLINK2" --bfile "${OUTDIR}/pooled_keep" \
    --king-cutoff "${KING_CUTOFF}" \
    --threads "$THREADS" \
    --make-bed \
    --out "${OUTDIR}/pooled_unrel"
else
  echo "[warn] PLINK2 not available on this node; using PLINK1 --rel-cutoff ${PIHAT_CUTOFF} as an approximation."
  "$PLINK1" --bfile "${OUTDIR}/pooled_keep" \
    --rel-cutoff "${PIHAT_CUTOFF}" \
    --threads "$THREADS" \
    --make-bed \
    --out "${OUTDIR}/pooled_unrel"
fi

# --- Split per-chromosome ---
echo "[step] Split per-chromosome (1..22)"
for c in {1..22}; do
  "$PLINK1" --bfile "${OUTDIR}/pooled_unrel" \
    --chr $c \
    --threads "$THREADS" \
    --make-bed \
    --out "${OUTDIR}/pooled_unrel_chr$c"
done

echo "[ok] Pooled unrelated per-chr PLINK ready in ${OUTDIR}"

