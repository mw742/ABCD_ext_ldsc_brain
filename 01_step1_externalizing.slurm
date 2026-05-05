#!/usr/bin/env bash
#SBATCH --job-name=s1_ext
#SBATCH --output=s1_ext_%j.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=shared
#SBATCH --chdir=/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe

# --------------- SAFETY ---------------
set -euo pipefail
IFS=$'\n\t'
export LC_ALL=C LANG=C
trap 'echo "[fatal] line $LINENO failed" >&2; exit 2' ERR
diag(){ echo "[diag] $*" >&2; }
die(){ echo "[err]  $*" >&2; exit 3; }

# --------------- INPUTS ---------------
PHENO_FILE="pheno_externalizing.tsv"   # expects FID IID EXT_baseline EXT_mean EXT_slope
BFILE="/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/imputed_genotype_qc/imputed_qc_filtered"
COV_BASE="/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe/covariates.base.fixed.v6.tsv"

OUTDIR="/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets/gwas_pipe/regenie_step1_ext"
PRUNE_DIR="${OUTDIR}/prune"
LRLD_BED="${PRUNE_DIR}/long_range_ld.bed"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
MEM_MB=$((64*1024))

# --------------- MODULES ---------------
module purge || true
module load plink/1.90b7.2 2>/dev/null || module load plink 2>/dev/null || true
module load regenie 2>/dev/null || true

# --------------- BASIC CHECKS ----------
[[ -f ${BFILE}.bed && -f ${BFILE}.bim && -f ${BFILE}.fam ]] || die "BFILE trio missing at ${BFILE}.{bed,bim,fam}"
[[ -s "$PHENO_FILE" ]] || die "Missing $PHENO_FILE"
[[ -s "$COV_BASE"   ]] || die "Missing $COV_BASE"

mkdir -p "$OUTDIR" "$PRUNE_DIR"

diag "PHENO_FILE -> $(readlink -f "$PHENO_FILE")"
diag "BFILE      -> $(readlink -f "${BFILE}.bed")"
diag "COV_BASE   -> $(readlink -f "$COV_BASE")"
diag "COV_BASE stat -> $(ls -l "$COV_BASE")"
diag "COV_BASE header -> $(head -n1 "$COV_BASE" | sed -E 's/\t/[TAB]/g' | sed $'s/\r/\\r/g')"
diag "COV_BASE header hex -> $(head -n1 "$COV_BASE" | xxd -p | head -n1)"

# --------------- PHENOTYPE COLS --------
PHENOS=$(head -n1 "$PHENO_FILE" | tr -d '\r' | awk -F'\t' '{for(i=3;i<=NF;i++)print $i}' \
  | grep -E '^(EXT_baseline|EXT_mean|EXT_slope)$' | paste -sd, -)
[[ -n "$PHENOS" ]] || die "No matching phenotype cols (EXT_baseline,EXT_mean,EXT_slope) in $PHENO_FILE"
diag "PHENOS=$PHENOS"

# --------------- COVARIATE SPEC --------
COV_CONT="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,motion"
COV_CAT="sex,site,scanner_manufacturer,scanner_device_serial"
diag "COV_CONT=$COV_CONT"
diag "COV_CAT=$COV_CAT"

# Header normalization
cov_hdr_clean=$(head -n1 "$COV_BASE" | sed '1s/^\xEF\xBB\xBF//' | tr -d '\r')
mapfile -t COV_FIELDS < <(awk -vFS='\t' '{for(i=1;i<=NF;i++) print $i}' <<< "$cov_hdr_clean")

has_field() {
  local key="$1"
  for f in "${COV_FIELDS[@]}"; do [[ "$f" == "$key" ]] && return 0; done
  return 1
}

IFS=, read -r -a CONT_ARR <<< "${COV_CONT// /}"
for c in "${CONT_ARR[@]}"; do has_field "$c" || die "Missing continuous covariate '$c' in $COV_BASE"; done
IFS=, read -r -a CAT_ARR <<< "${COV_CAT// /}"
for c in "${CAT_ARR[@]}"; do has_field "$c" || die "Missing categorical covariate '$c' in $COV_BASE"; done

# --------------- NUMERIC AUDIT ---------
NUM_AUDIT="${OUTDIR}/_cov_numeric_audit.tsv"
diag "Auditing numeric covariates -> $NUM_AUDIT (first 50 bad rows if any)"
{
  echo -e "row\tcolumn\tvalue"
  awk -v FS='\t' -v cols="$COV_CONT" '
    BEGIN{ split(cols, C, ","); for(i=1;i<=length(C);i++) want[C[i]]=1 }
    NR==1{ for(i=1;i<=NF;i++) ix[$i]=i; next }
    {
      for (k in want){
        v = $(ix[k])
        if (v == "NA" || v == "") continue
        if (v !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) print NR "\t" k "\t" v
      }
    }
  ' "$COV_BASE"
} | head -n 51 > "$NUM_AUDIT" || true

if (( $(wc -l < "$NUM_AUDIT") > 1 )); then
  diag "Found non-numeric in continuous covariates (showing up to 50 rows):"
  head "$NUM_AUDIT" >&2
  die "Fix numeric issues first."
fi

# --------------- CATEGORICAL LEVELS ----
# Count unique levels per categorical covariate (ignoring empty -> UNK)
LEVEL_REPORT="${OUTDIR}/_cov_cat_levels.tsv"
diag "Counting categorical levels -> $LEVEL_REPORT"
{
  echo -e "covariate\tlevels"
  awk -v FS='\t' -v cats="$COV_CAT" '
    BEGIN{ split(cats, C, ","); for(i=1;i<=length(C);i++) want[C[i]]=1 }
    NR==1{
      for(i=1;i<=NF;i++) ix[$i]=i
      for (k in want) { n[k]=0 }
      next
    }
    {
      for (k in want){
        v=$(ix[k]); if (v=="") v="UNK"
        key=k "\034" v
        if (!(key in seen)) { seen[key]=1; n[k]++ }
      }
    }
    END{
      for (k in want) printf "%s\t%d\n", k, n[k]
    }
  ' "$COV_BASE"
} > "$LEVEL_REPORT"

cat "$LEVEL_REPORT" >&2
OBS_MAX=$(awk 'NR>1{if($2>m)m=$2}END{print (m?m:0)}' "$LEVEL_REPORT")
# set MAXCAT to a safe headroom above observed (obs+5), with a floor 32
MAXCAT=$(( OBS_MAX + 5 ))
(( MAXCAT < 32 )) && MAXCAT=32
diag "Observed max categorical levels = ${OBS_MAX}; using --maxCatLevels ${MAXCAT}"

# Soft warnings for huge categoricals (e.g., device_serial)
awk 'NR>1 && $2>200 {printf("[warn] %s has %d levels; Step-1 may be slower.\n",$1,$2)}' "$LEVEL_REPORT" >&2 || true

# --------------- LONG-RANGE LD BED -----
cat > "$LRLD_BED" <<'BED'
6	25000000	35000000	MHC
8	8000000	12000000	8p23
11	45000000	57000000	11q
17	40000000	65000000	17q
BED

# --------------- PRUNING PARAMS --------
MAF="0.01"; GENO="0.01"; HWE="1e-6"
R2_LIST=("0.10" "0.08" "0.05" "0.02")
WINDOW_KB=1000; STEP=50

# --------------- QC SNPLIST ------------
diag "QC filter -> write snplist"
plink --bfile "$BFILE" --maf $MAF --geno $GENO --hwe $HWE midp \
  --write-snplist --threads $THREADS --memory $MEM_MB \
  --out "${PRUNE_DIR}/qc"

# --------------- EXCLUDE LR-LD ---------
diag "Exclude long-range LD regions"
plink --bfile "$BFILE" \
  --extract "${PRUNE_DIR}/qc.snplist" \
  --exclude range "$LRLD_BED" \
  --write-snplist --threads $THREADS --memory $MEM_MB \
  --out "${PRUNE_DIR}/qc_nolrld"

# --------------- PROGRESSIVE PRUNE -----
PRUNED_IN=""
for R2 in "${R2_LIST[@]}"; do
  diag "LD-prune with r^2=${R2}"
  plink --bfile "$BFILE" \
    --extract "${PRUNE_DIR}/qc_nolrld.snplist" \
    --indep-pairwise ${WINDOW_KB} ${STEP} ${R2} \
    --threads $THREADS --memory $MEM_MB \
    --out "${PRUNE_DIR}/prune_r2_${R2//./}"
  PRUNED_IN="${PRUNE_DIR}/prune_r2_${R2//./}.prune.in"
  N=$(wc -l < "$PRUNED_IN")
  diag "Candidate pruned SNPs at r^2=${R2}: $N"
  if (( N <= 1000000 )); then
    diag "Accepting pruned set with $N SNPs"
    break
  fi
done
[[ -s "$PRUNED_IN" ]] || die "Pruning produced no SNP list"
N_FINAL=$(wc -l < "$PRUNED_IN")
(( N_FINAL <= 1000000 )) || die "Still >1,000,000 SNPs after tightening; consider stronger pruning or lower MAF."

# --------------- MAKE STEP1 BED --------
diag "Creating compact Step-1 bed"
plink --bfile "$BFILE" \
  --extract "$PRUNED_IN" \
  --make-bed --threads $THREADS --memory $MEM_MB \
  --out "${PRUNE_DIR}/step1_pruned"

# --------------- RUN REGENIE S1 --------
diag "Running REGENIE Step-1"
regenie \
  --step 1 \
  --bed "${PRUNE_DIR}/step1_pruned" \
  --phenoFile "$PHENO_FILE" \
  --phenoColList "$PHENOS" \
  --covarFile "$COV_BASE" \
  --covarColList "$COV_CONT" \
  --catCovarList "$COV_CAT" \
  --maxCatLevels "${MAXCAT}" \
  --ref-first \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${OUTDIR}/lm" \
  --threads "$THREADS" \
  --out "${OUTDIR}/fit_ext"

# --------------- PREDICTION LIST -------
find "$OUTDIR" -maxdepth 1 -type f -name "fit_ext_*_pred.list" -print > "${OUTDIR}/s1_pred.list" || true
if [[ ! -s "${OUTDIR}/s1_pred.list" ]]; then
  ls -1 "${OUTDIR}"/fit_ext_*_loco*.gz > "${OUTDIR}/s1_pred.list" || true
fi

if [[ -s "${OUTDIR}/s1_pred.list" ]]; then
  diag "Wrote Step-1 pred list: ${OUTDIR}/s1_pred.list"
else
  diag "Note: verify Step-1 prediction outputs for your REGENIE build."
fi

diag "Done. Step-1 outputs in: $OUTDIR"

