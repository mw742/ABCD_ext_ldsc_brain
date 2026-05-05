#!/usr/bin/env bash
#SBATCH --job-name=regenie_s2_img
#SBATCH --output=regenie_s2_img_%A_%a.log
#SBATCH --array=1-210
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=shared

# ----------------------- SAFETY -----------------------
set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C
trap 'echo "[fatal] line $LINENO failed" >&2' ERR

diag(){ echo "[diag] $*" >&2; }
die(){  echo "[err]  $*" >&2; exit 2; }

# ----------------------- MODULES -----------------------
module purge || true
module load regenie 2>/dev/null || true
module load plink   2>/dev/null || true

# ----------------------- PATHS.ENV (robust sourcing) -----------------------
SCRIPT_PATH="${BASH_SOURCE[0]-$0}"
SCRIPT_DIR="$(cd -- "$(dirname -- "$SCRIPT_PATH")" >/dev/null 2>&1 && pwd)"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

try_paths=(
  "${SUBMIT_DIR}"
  "${SCRIPT_DIR}"
  "$PWD"
)

PATHS_ENV=""
for d in "${try_paths[@]}"; do
  if [[ -f "${d}/_paths.env" ]]; then
    PATHS_ENV="${d}/_paths.env"
    break
  fi
done
[[ -n "$PATHS_ENV" ]] || die "Cannot find _paths.env in: ${try_paths[*]}"
# shellcheck source=/dev/null
source "$PATHS_ENV"

# ----------------------- REQUIRED VARS + DEFAULTS -----------------------
: "${WORK:=${WORK:-${SUBMIT_DIR}}}"
: "${OUTDIR:=${OUTDIR:-${WORK}/gwas_pipe}}"
: "${BFILE:=${BFILE:-${WORK}/imputed_genotype_qc/imputed_qc_filtered}}"
: "${PHENO_FILE:=${PHENO_FILE:-${WORK}/gwas_pipe/pheno_imaging.tsv}}"
: "${IMG_COLS_ALL:=${IMG_COLS_ALL:-${WORK}/gwas_pipe/img_all_cols.txt}}"
: "${COV_BASE_IMG:=${COV_BASE_IMG:-${WORK}/gwas_pipe/covariates.base.fixed.v6.tsv}}"
: "${COV_WITHICV:=${COV_WITHICV:-${WORK}/gwas_pipe/covariates.base.withICV.fixed.v6.tsv}}"
: "${BSIZE:=${BSIZE:-400}}"
: "${MIN_MAC:=${MIN_MAC:-20}}"
: "${SLURM_CPUS_PER_TASK:=${SLURM_CPUS_PER_TASK:-8}}"
: "${MAX_IMG:=${MAX_IMG:-210}}"

# Fail fast if the essentials are missing on disk
[[ -s "$PHENO_FILE"    ]] || die "PHENO_FILE not found or empty: $PHENO_FILE"
[[ -s "$IMG_COLS_ALL"  ]] || die "IMG_COLS_ALL not found or empty: $IMG_COLS_ALL"
[[ -s "${BFILE}.bed"   ]] || die "BFILE missing: ${BFILE}.bed"
[[ -s "${BFILE}.bim"   ]] || die "BFILE missing: ${BFILE}.bim"
[[ -s "${BFILE}.fam"   ]] || die "BFILE missing: ${BFILE}.fam"

# ----------------------- CHECK STEP1 PRED -----------------------
S1_OUT="${OUTDIR}/regenie_step1_imaging"
if   [[ -s "${S1_OUT}/s1_pred.list.gz" ]]; then PRED="${S1_OUT}/s1_pred.list.gz"
elif [[ -s "${S1_OUT}/s1_pred.list"    ]]; then PRED="${S1_OUT}/s1_pred.list"
else die "Step-1 pred list missing in ${S1_OUT}"
fi
diag "Using PRED=${PRED}"

# ----------------------- ONE-TIME HEADER AUDIT (task 1 only) -----------------------
if [[ "${SLURM_ARRAY_TASK_ID:-1}" -eq 1 ]]; then
  diag "Running header audit (one-time)…"
  HDR="$(mktemp)"; head -n1 "$PHENO_FILE" > "$HDR"
  awk -F'\t' 'NR==1{for(i=3;i<=NF;i++){c[$i]++}} END{for(k in c) if(c[k]>1) printf("[warn] duplicate header: %s (x%d)\n",k,c[k])}' "$PHENO_FILE" || true
  grep -nE 'Posterior_Cingulate_SA_meanPosterior_Cingulate_SA_slope' "$HDR" \
    && echo "[warn] header may be missing a tab between Posterior_Cingulate_SA_mean and Posterior_Cingulate_SA_slope" >&2 || true
  rm -f "$HDR"
fi

# ----------------------- PHENOTYPE SELECTION -----------------------
TMP_COLS="${SLURM_TMPDIR:-/tmp}/img_cols_head_${SLURM_JOB_ID:-$$}.lst"
head -n "${MAX_IMG}" "$IMG_COLS_ALL" > "$TMP_COLS"

N_COLS=$(wc -l < "$TMP_COLS")
[[ "$N_COLS" -gt 0 ]] || die "No phenotype columns found after capping to MAX_IMG=${MAX_IMG}."

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
if [[ "$TASK_ID" -gt "$N_COLS" ]]; then
  echo "[skip] task ${TASK_ID} > N=${N_COLS} (MAX_IMG=${MAX_IMG})"
  exit 0
fi

PHENO="$(sed -n "${TASK_ID}p" "$TMP_COLS")"
[[ -n "$PHENO" ]] || die "Empty phenotype label at task ${TASK_ID}"

# Verify PHENO exists in the pheno header (skip first 2 ID columns)
head -n1 "$PHENO_FILE" | tr -d '\r' | tr '\t' '\n' | awk 'NR>2' | grep -Fxq "$PHENO" \
  || die "PHENO=$PHENO not found in header of $PHENO_FILE"

diag "PHENO=${PHENO} (task ${TASK_ID}/${N_COLS})"

# ----------------------- CHOOSE COVARIATE FILE -----------------------
# SA phenotypes (or GLOBAL_SA_*) get the ICV-adjusted covariates
if [[ "$PHENO" == *_SA_* || "$PHENO" == GLOBAL_SA_* ]]; then
  COVF="${COV_WITHICV}"
else
  COVF="${COV_BASE_IMG}"
fi
[[ -s "$COVF" ]] || die "Covariate file not found: ${COVF}"

# ----------------------- SPLIT NUMERIC vs CATEGORICAL -----------------------
# These are treated as categorical by design
CAT_SET="site,scanner_manufacturer,scanner_device_serial"

hdr=$(head -n1 "$COVF" | tr -d '\r')

COVARS_NUM=$(
  awk -vFS='\t' -v cats="$CAT_SET" '
    BEGIN{ nsplit=split(cats, ctmp, ","); for(i=1;i<=nsplit;i++) cat[ctmp[i]]=1 }
    NR==1{ for(i=3;i<=NF;i++) if(!( $i in cat )) printf("%s%s",$i,(i<NF?",":"")); printf("\n") }
  ' <(printf "%s\n" "$hdr") | sed 's/,,*/,/g; s/,$//'
)
COVARS_CAT=$(
  awk -vFS='\t' -v cats="$CAT_SET" '
    BEGIN{ nsplit=split(cats, ctmp, ","); for(i=1;i<=nsplit;i++) cat[ctmp[i]]=1; first=1 }
    NR==1{ for(i=3;i<=NF;i++) if( $i in cat ){ if(!first) printf(","); printf("%s",$i); first=0 } printf("\n") }
  ' <(printf "%s\n" "$hdr")
)
COVARS_NUM="$(echo -n "$COVARS_NUM" | tr -d '[:space:]')"
COVARS_CAT="$(echo -n "$COVARS_CAT" | tr -d '[:space:]')"
[[ -n "$COVARS_NUM" ]] || die "Empty numeric covariate list derived from ${COVF}"
[[ -n "$COVARS_CAT" ]] && diag "Categorical covars: ${COVARS_CAT}" || diag "No categorical covars found"
diag "Numeric covars: ${COVARS_NUM}"

# ----------------------- CHOOSE maxCatLevels DYNAMICALLY -----------------------
MAX_CAT_LEVELS_OPT=()
if [[ -n "$COVARS_CAT" ]]; then
  # Determine the maximum number of distinct levels among the categorical covariates
  # Assumes first two columns are IDs (FID, IID). We read the actual file (not just header).
  MAX_CAT_LEVELS=$(
    awk -vFS='\t' -v cats="$COVARS_CAT" '
      BEGIN{
        n=split(cats,cc,","); for(i=1;i<=n;i++) want[cc[i]]=0
      }
      NR==1{
        # Map covariate names -> column indices
        for(i=3;i<=NF;i++) if($i in want){ col[$i]=i }
        next
      }
      {
        for (k in col) {
          v=$col[k]
          # treat empty as NA (still a level); change to skip if you prefer
          key=k SUBSEP v
          if(!(key in seen)){ seen[key]=1; cnt[k]++ }
        }
      }
      END{
        mx=0
        for(k in cnt) if(cnt[k]>mx) mx=cnt[k]
        if(mx<1) mx=1
        print mx
      }
    ' "$COVF"
  )
  # Add headroom to avoid edge cases (e.g., 22 -> 32)
  : "${MAX_CAT_HEADROOM:=10}"
  if [[ "$MAX_CAT_LEVELS" =~ ^[0-9]+$ ]]; then
    MAX_CAT_LEVELS=$(( MAX_CAT_LEVELS + MAX_CAT_HEADROOM ))
  else
    MAX_CAT_LEVELS=50
  fi
  # Bound it to something sane to avoid accidental explosions
  if (( MAX_CAT_LEVELS < 25 )); then MAX_CAT_LEVELS=25; fi
  if (( MAX_CAT_LEVELS > 200 )); then MAX_CAT_LEVELS=200; fi
  diag "Setting --maxCatLevels=${MAX_CAT_LEVELS} for categorical covariates: ${COVARS_CAT}"
  MAX_CAT_LEVELS_OPT=( --maxCatLevels "${MAX_CAT_LEVELS}" )
fi

# ----------------------- OUTPUTS -----------------------
OUTROOT="${WORK}/gwas_output/imaging"
mkdir -p "$OUTROOT"

diag "Starting regenie step2…"
diag "BFILE=${BFILE}"
diag "PHENO_FILE=${PHENO_FILE}"
diag "COVF=${COVF}"
diag "OUTROOT=${OUTROOT}"

# ----------------------- RUN REGENIE STEP 2 -----------------------
cmd=( regenie
  --step 2
  --bed "$BFILE"
  --phenoFile "$PHENO_FILE"
  --phenoCol "$PHENO"
  --covarFile "$COVF"
  --covarColList "$COVARS_NUM"
  --pred "$PRED"
  --bsize "${BSIZE}"
  --minMAC "${MIN_MAC}"
  --threads "${SLURM_CPUS_PER_TASK}"
  --gz
  --out "${OUTROOT}/${PHENO}"
)

if [[ -n "$COVARS_CAT" ]]; then
  cmd+=( --catCovarList "$COVARS_CAT" "${MAX_CAT_LEVELS_OPT[@]}" )
fi

printf "[diag] Running: %q " "${cmd[@]}" >&2; echo >&2
"${cmd[@]}"

echo "[done] ${PHENO}"

