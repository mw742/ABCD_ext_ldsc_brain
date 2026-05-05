#!/usr/bin/env bash
#SBATCH --job-name=ldsc_l2_pooled
#SBATCH --output=ldsc_l2_pooled_%A_%a.log
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --partition=shared

# ---------------- SAFETY ----------------
set -euo pipefail
IFS=$'\n\t'

echo "[diag] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# --------------- MODULES ----------------
module purge
module load python/2.7.11
module load ldsc/1.0.0

# Prefer the module-provided wrapper on PATH
LDSC_BIN="$(command -v ldsc.py || true)"
if [[ -z "${LDSC_BIN}" ]]; then
  LDSC_BIN="/opt/applications/ldsc/1.0.0/gnu/bin/ldsc.py"
fi

# --------------- PYTHONPATH FIX ---------
if ! python - <<'PY' >/dev/null 2>&1
import ldscore
print("ok")
PY
then
  for p in \
    /opt/applications/ldsc/1.0.0/gnu/lib/python2.7/site-packages \
    /opt/applications/ldsc/1.0.0/gnu/lib64/python2.7/site-packages
  do
    [[ -d "$p" ]] && export PYTHONPATH="${p}:${PYTHONPATH-}"
  done
fi

# --------------- PATHS ------------------
ROOT="/gpfs/group/ehlers/mwei/ABCD_2/ABCD"
WORK="${ROOT}/work_provided_assets"
PANEL_DIR="${WORK}/ldsc_panel_pooled"
OUT_DIR="${WORK}/ldsc_ref_pooled_custom"
mkdir -p "${OUT_DIR}"

CHR="${SLURM_ARRAY_TASK_ID}"
BFILE_PREFIX="${PANEL_DIR}/pooled_unrel_chr${CHR}"

# --------------- CHECKS -----------------
[[ -s "${BFILE_PREFIX}.bim" && -s "${BFILE_PREFIX}.bed" && -s "${BFILE_PREFIX}.fam" ]] || {
  echo "[err] Missing PLINK files for chr ${CHR}: ${BFILE_PREFIX}.{bed,bim,fam}"
  exit 2
}

echo "[diag] Using Python: $(command -v python)"
python -V
echo "[diag] Using LDSC: ${LDSC_BIN}"
echo "[diag] BIM line count: $(wc -l < "${BFILE_PREFIX}.bim")"

# --------------- WINDOW (KB) ------------
WIND_KB="${WIND_KB:-1000}"   # allow override via env if desired
echo "[diag] Using --ld-wind-kb=${WIND_KB}"

# --------------- RUN --------------------
python "${LDSC_BIN}" \
  --l2 \
  --bfile "${BFILE_PREFIX}" \
  --ld-wind-kb "${WIND_KB}" \
  --out "${OUT_DIR}/pooled.custom.${CHR}"

echo "[ok] ${OUT_DIR}/pooled.custom.${CHR}.l2.ldscore.gz"

