#!/usr/bin/env bash
#SBATCH --job-name=regenie_s2_ext
#SBATCH --output=regenie_s2_ext_%A_%a.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=shared
#SBATCH --chdir=/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe
#SBATCH --array=1-3

set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C
trap 'echo "[fatal] line $LINENO failed" >&2; exit 2' ERR
diag(){ echo "[diag] $*" >&2; }; die(){ echo "[err] $*" >&2; exit 3; }

PHENO_FILE="pheno_externalizing.tsv"
BFILE="/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/imputed_genotype_qc/imputed_qc_filtered"
WORK="/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets"
OUTROOT="${WORK}/gwas_output/externalizing"
S1_OUT="${WORK}/gwas_pipe/regenie_step1_ext"

# prefer fixed v6 if present
COV_BASE_CAND=(
  "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe/covariates.base.fixed.v6.tsv"
  "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe/covariates.base.tsv"
)
COV_BASE=""
for c in "${COV_BASE_CAND[@]}"; do [[ -s "$c" ]] && COV_BASE="$c" && break; done
[[ -n "$COV_BASE" ]] || die "No covariate file found."

module purge   2>/dev/null || true
module load regenie 2>/dev/null || true
module load plink   2>/dev/null || true

[[ -f ${BFILE}.bed && -f ${BFILE}.bim && -f ${BFILE}.fam ]] || die "BFILE trio missing"
[[ -s "$PHENO_FILE"  ]] || die "Missing $PHENO_FILE"
[[ -s "$COV_BASE"    ]] || die "Missing $COV_BASE"

# ---- Step-1 predictions ----
PRED_LIST=""
if   [[ -s "${S1_OUT}/s1_pred.list" ]]; then
  PRED_LIST="${S1_OUT}/s1_pred.list"
elif [[ -s "${S1_OUT}/fit_ext_pred.list" ]]; then
  PRED_LIST="${S1_OUT}/fit_ext_pred.list"
elif [[ -s "${S1_OUT}/s1_pred.list.gz" ]]; then
  PRED_LIST="${S1_OUT}/s1_pred.list.gz"
else
  mapfile -t LOCO < <(ls -1 "${S1_OUT}"/fit_ext_*.loco 2>/dev/null || true)
  (( ${#LOCO[@]} )) || die "Missing Step-1 preds in ${S1_OUT}"
  PRED_LIST="${S1_OUT}/_autogen_pred.list"
  printf "%s\n" "${LOCO[@]}" > "$PRED_LIST"
fi
if [[ "$PRED_LIST" =~ \.gz$ ]]; then
  n_lines=$(zcat "$PRED_LIST" | sed '/^[[:space:]]*$/d' | wc -l)
else
  n_lines=$(sed '/^[[:space:]]*$/d' "$PRED_LIST" | wc -l)
fi
(( n_lines > 0 )) || die "Prediction list $PRED_LIST is empty"

# ---- phenotype list (3 externalizing columns) ----
ALL_COLS="$(head -n1 "$PHENO_FILE" | tr -d '\r' | awk -F'\t' '{for(i=3;i<=NF;i++) print $i}')"
PHENO_LIST="$(printf "%s\n" "$ALL_COLS" | grep -E '^(EXT_baseline|EXT_mean|EXT_slope)$' || true)"
[[ -n "$PHENO_LIST" ]] || die "No matching phenotypes in $PHENO_FILE"
CLEAN_PHENOS="$(mktemp)"; printf "%s\n" "$PHENO_LIST" > "$CLEAN_PHENOS"
N_COLS=$(wc -l < "$CLEAN_PHENOS")

# ---- split covariates: numeric vs categorical ----
declare -A CAT=([site]=1 [scanner_manufacturer]=1 [scanner_device_serial]=1)
NUM_LIST=()
CAT_LIST=()
read -r covhdr < <(head -n1 "$COV_BASE" | tr -d '\r')
IFS=$'\t' read -r -a covarr <<< "$covhdr"
for ((i=2;i<${#covarr[@]};i++)); do
  col="${covarr[$i]}"
  if [[ -n "${CAT[$col]:-}" ]]; then
    CAT_LIST+=("$col")
  else
    NUM_LIST+=("$col")
  fi
done
COVARS="$(IFS=,; echo "${NUM_LIST[*]}")"
CATCOVARS="$(IFS=,; echo "${CAT_LIST[*]}")"
[[ -n "$COVARS" ]] || die "Empty numeric covariate list"
[[ -n "$CATCOVARS" ]] || die "Empty categorical covariate list"

# ---- compute required --maxCatLevels dynamically ----
# allow manual override: export MAX_CAT_LEVELS=50
if [[ -z "${MAX_CAT_LEVELS:-}" ]]; then
  # build index map for categorical columns
  # get 1-based column indices
  get_idx() {
    local name="$1"
    awk -v FS='\t' -v n="$name" 'NR==1{for(i=1;i<=NF;i++) if($i==n) {print i; exit 0}}' "$COV_BASE"
  }
  max_levels=0
  for cname in "${CAT_LIST[@]}"; do
    idx=$(get_idx "$cname" || true)
    [[ -n "$idx" ]] || die "Cannot find column '$cname' in $COV_BASE"
    # count unique non-empty levels
    nlev=$(awk -v FS='\t' -v I="$idx" 'NR>1 && $I!=""{a[$I]=1} END{print length(a)}' "$COV_BASE")
    diag "Levels for ${cname}: ${nlev}"
    (( nlev > max_levels )) && max_levels=$nlev
  done
  # give a little headroom
  MAX_CAT_LEVELS=$(( max_levels + 5 ))
fi
diag "Using --maxCatLevels=${MAX_CAT_LEVELS}"

diag "BFILE=$BFILE"
diag "PHENO_FILE=$PHENO_FILE"
diag "COV_BASE=$COV_BASE"
diag "S1 pred list=$PRED_LIST (n=${n_lines})"
diag "N phenotypes = $N_COLS (expected 3)"
diag "COVARS=$COVARS"
diag "CATCOVARS=$CATCOVARS"

# ---- pick phenotype for this task ----
TASK=${SLURM_ARRAY_TASK_ID:-1}
PHENO_OVERRIDE="${PHENO_OVERRIDE:-}"
if [[ -n "$PHENO_OVERRIDE" ]]; then
  PHENO="$PHENO_OVERRIDE"
else
  (( TASK >= 1 )) || TASK=1
  (( TASK <= N_COLS )) || { echo "[skip] task ${TASK} > N=${N_COLS}"; exit 0; }
  PHENO=$(sed -n "${TASK}p" "$CLEAN_PHENOS")
fi
[[ -n "$PHENO" ]] || die "Empty phenotype"

mkdir -p "$OUTROOT"
OUT_PREFIX="${OUTROOT}/${PHENO}"
if [[ -s "${OUT_PREFIX}.regenie.gz" ]]; then
  diag "Output exists; skipping ${OUT_PREFIX}.regenie.gz"
  echo "[done] ${PHENO} (already present)"
  exit 0
fi

# ---- run REGENIE step 2 ----
diag "PHENO=$PHENO"
diag "Starting REGENIE step 2..."
regenie \
  --step 2 \
  --bed "$BFILE" \
  --phenoFile "$PHENO_FILE" \
  --phenoCol "$PHENO" \
  --covarFile "$COV_BASE" \
  --covarColList "$COVARS" \
  --catCovarList "$CATCOVARS" \
  --maxCatLevels "${MAX_CAT_LEVELS}" \
  --pred "$PRED_LIST" \
  --bsize 400 \
  --minMAC 20 \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --gz \
  --out "$OUT_PREFIX"

diag "Finished REGENIE step 2"
echo "[done] ${PHENO}"

