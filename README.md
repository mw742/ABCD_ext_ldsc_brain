# ABCD Externalizing–Imaging GWAS and LDSC Genetic Correlation Pipeline

This repository contains the scripts used to run a GWAS and LDSC workflow for evaluating genetic relationships between longitudinal externalizing behavior phenotypes and brain imaging-derived phenotypes in the ABCD Study.

The workflow has three main goals:

1. Prepare phenotype, covariate, and sample files for REGENIE.
2. Run REGENIE GWAS for externalizing and imaging phenotypes.
3. Build a pooled unrelated LD reference panel and estimate genetic correlations between externalizing and imaging GWAS results using LDSC.

The scripts are written for an HPC environment using SLURM, PLINK, REGENIE, Python, and LDSC.

---

## Repository contents

| Script | Purpose | Main output |
|---|---|---|
| `00_prep_inputs.slurm` | Prepares genotype-matched phenotype and covariate files for REGENIE. It intersects genotype, covariate, externalizing, imaging, and ICV sample IDs; cleans covariate headers; creates phenotype tables; creates phenotype column lists; and writes `_paths.env`. | `gwas_pipe/pheno_externalizing.tsv`, `gwas_pipe/pheno_imaging.tsv`, `gwas_pipe/covariates.*.tsv`, `gwas_pipe/_paths.env`, phenotype column lists |
| `00_prep_img_cols.sh` | Helper script that extracts imaging phenotype columns from the prepared imaging phenotype file. It keeps columns ending in `_baseline`, `_mean`, or `_slope`, grouped by imaging feature base name. | `IMG_COLS_ALL.txt` |
| `01_step1_externalizing.slurm` | Runs REGENIE Step 1 for the three externalizing phenotypes. It performs SNP QC, removes long-range LD regions, LD-prunes variants, creates a compact Step 1 genotype set, and generates LOCO predictions. | `gwas_pipe/regenie_step1_ext/` |
| `01_step1_imaging.slurm` | Runs REGENIE Step 1 for imaging phenotypes. It sources `_paths.env`, optionally limits the number of imaging phenotypes, builds or reuses an LD-pruned training SNP list, and generates LOCO predictions. | `gwas_pipe/regenie_step1_imaging/` |
| `02_step2_externalizing.slurm` | Runs REGENIE Step 2 for externalizing phenotypes as a SLURM array. Each task analyzes one externalizing phenotype. | `gwas_output/externalizing/*.regenie.gz` |
| `02_step2_imaging.slurm` | Runs REGENIE Step 2 for imaging phenotypes as a SLURM array. Each task analyzes one imaging phenotype. Surface-area phenotypes use ICV-adjusted covariates when applicable. | `gwas_output/imaging/*.regenie.gz` |
| `03_mk_unrel_pooled.slurm` | Builds the pooled unrelated genotype panel for LDSC. It uses curated unrelated individuals when available, applies a KING or PI-HAT relatedness filter, and splits the final panel by chromosome. | `ldsc_panel_pooled/pooled_unrel_chr{1..22}.{bed,bim,fam}` |
| `04_ldsc_l2_pooled.slurm` | Computes LD scores for the pooled unrelated reference panel using LDSC. This runs as a chromosome-level SLURM array. | `ldsc_ref_pooled_custom/pooled.custom.{1..22}.l2.ldscore.gz` |
| `05_rg_ext_vs_img_pooled.slurm` | Munges REGENIE summary statistics and runs LDSC genetic correlation between each imaging phenotype and each externalizing phenotype. It is resumable and skips completed outputs. | `rg_ldsc_pooled/munged_ext/`, `rg_ldsc_pooled/munged_img/`, `rg_ldsc_pooled/rg_results/` |


---

## Required software

The scripts assume the following tools are available on the HPC system:

- SLURM
- PLINK 1.9
- PLINK 2.0, preferred for KING relatedness filtering
- REGENIE
- LDSC 1.0.0
- Python 2.7 for LDSC
- Python 3 or a micromamba environment for preprocessing
- `gzip` or `bgzip`

Some scripts contain environment-specific fixes for Python 2 and OpenSSL compatibility. These may need to be removed or modified on a different cluster.

---

## Pipeline overview

```text
Step 0: Prepare REGENIE inputs
        00_prep_inputs.slurm
        00_prep_img_cols.sh

Step 1: Run REGENIE Step 1
        01_step1_externalizing.slurm
        01_step1_imaging.slurm

Step 2: Run REGENIE Step 2
        02_step2_externalizing.slurm
        02_step2_imaging.slurm

Step 3: Build pooled unrelated LDSC panel
        03_mk_unrel_pooled.slurm

Step 4: Compute custom LDSC LD scores
        04_ldsc_l2_pooled.slurm

Step 5: Munge GWAS summary statistics and run LDSC genetic correlation
        05_rg_ext_vs_img_pooled.slurm
```


---

## Citation

If this pipeline is used in a manuscript or repository, describe it as:

Wei, M. and Peng, Q., 2026. An Initial Genetic Correlation Analysis of Externalizing Behavior and Neuroimaging Phenotypes in the ABCD Cohort. medRxiv, pp.2026-07.
