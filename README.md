# Mutation-specific impairment of TET2 and DNMT3A enzymatic activity predicts clonal hematopoiesis disease risk

This repository contains the code and analysis scripts for [Pershad & Zhao et al, medRxiv, 2026](https://www.medrxiv.org/content/10.64898/2026.04.03.26350108v1).

## Table of Contents
1. [Overview](#overview)
2. [CHIP detection and mutation-type classification](#chip-detection-and-mutation-type-classification)
3. [Phenome-wide association studies by mutation type](#phenome-wide-association-studies-by-mutation-type)
4. [EM-seq data processing and beta matrix generation](#em-seq-data-processing-and-beta-matrix-generation)
5. [Epigenome-wide association study (EWAS) for CHIP status](#epigenome-wide-association-study-ewas-for-chip-status)
6. [Sliding window differentially methylated region analysis](#sliding-window-differentially-methylated-region-analysis)
7. [Activity score construction](#activity-score-construction)
8. [Validation in BioVU and CHIVE cohorts](#validation-in-biovu-and-chive-cohorts)
9. [Clinical risk model comparisons](#clinical-risk-model-comparisons)
10. [Data](#data)
11. [Citation](#citation)
12. [License](#license)
13. [Acknowledgements](#acknowledgements)
14. [Contact](#contact)

---

## Overview

We analyzed 1,020,538 individuals across three biobanks (UK Biobank, All of Us, Vanderbilt BioVU) to characterize the mutational landscape, cellular fitness effects, and clinical consequences of TET2 and DNMT3A clonal hematopoiesis of indeterminate potential (CHIP). We found that:

- **TET2 loss-of-function (pLoF)** mutations and **DNMT3A R882** hotspot mutations account for the majority of clinical risk
- These high-risk mutations confer stronger clonal fitness advantage (measured via copy-neutral loss of heterozygosity co-occurrence)
- Peripheral blood DNA methylation provides a direct quantitative readout of enzymatic dysfunction in CHIP patients
- Methylation-based **TET2 and DNMT3A activity scores** outperform and complement existing clinical risk models (CHRS, AHA PREVENT)

This repository is organized by analysis step, from CHIP detection through activity score construction and clinical validation.

---

## CHIP detection and mutation-type classification

Putative somatic SNVs and short indels were called with GATK Mutect2 and filtered using a curated list of 74 canonical CHIP driver genes, as described in [Vlasschaert et al., *Blood*, 2023](https://doi.org/10.1182/blood.2022018825) and applied consistently across BioVU, All of Us, and UK Biobank.

**Variant filters applied:**
- Present in a pre-established list of candidate CHIP variants
- Total sequencing depth ≥ 20
- Alternate allele read depth ≥ 5
- Representation in both sequencing directions (F1R2 ≥ 1 and F2R1 ≥ 1)
- Variant allele fraction (VAF) ≥ 0.02

**TET2 mutation-type classification:**

| Category | Definition |
|---|---|
| pLoF (early) | Nonsense or frameshift variants proximal to amino acid 1129 (preceding catalytic domain) |
| pLoF (late) | Nonsense or frameshift variants within or distal to the catalytic domain |
| Missense | Non-synonymous missense variants meeting CHIP inclusion criteria |

**DNMT3A mutation-type classification:**

| Category | Definition |
|---|---|
| R882 | R882H or R882C hotspot missense substitutions |
| Non-R882 | All other CHIP-qualifying DNMT3A variants (missense, splicing, frameshift, stop-gain) |

Code for mutation-type classification and CN-LOH co-occurrence analysis is available in [notebooks/mutational_profiling_phewas.ipynb](notebooks/mutational_profiling_phewas.ipynb).

---

## Phenome-wide association studies by mutation type

We conducted time-to-event phenome-wide association studies (PheWAS) across 2,180 clinical phenotypes defined by Phecodes derived from ICD-9/ICD-10 codes in electronic health records, performed separately in BioVU, All of Us, and UK Biobank and combined via inverse-variance weighted fixed-effects meta-analysis.

Cox proportional hazards models were adjusted for:
- Age at blood draw (continuous)
- Age² (continuous)
- Genetic sex (categorical)
- Current smoking status (categorical)
- Principal components 1–10 (continuous)

For TET2, models were run separately for pLoF and missense variants. For DNMT3A, models were run separately for R882 and non-R882 variants. Effect size heterogeneity across the phenome was assessed using mixed-effects meta-regression (R package `metafor`).

Code for PheWAS processing and meta-regression is available in [notebooks/mutational_profiling_phewas.ipynb](notebooks/mutational_profiling_phewas.ipynb).

---

## EM-seq data processing and beta matrix generation

Peripheral blood DNA from BioVU and CHIVE cohort participants was profiled using targeted enzymatic methyl-sequencing (EM-seq) via the NEBNext Enzymatic Methyl-seq kit paired with a Twist Biosciences hybrid capture panel targeting ~4 million CpG sites. Libraries were sequenced on an Illumina NovaSeq 6000 (150 bp paired-end) to a mean depth of 30×.

### Step 1: From CX_report files to beta matrix

Raw reads were aligned to GRCh38 and CpG methylation states were extracted into CX_report files using [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). We then constructed per-sample methylation matrices using a parallelized WDL pipeline on Terra.

**Coverage filtering:**
- Minimum coverage per CpG: 10×
- High-confidence threshold: 25×
- Sites retained if ≥ 80% of cohort samples had ≥ 25× coverage

The WDL workflow for matrix generation is available in [wdl/optimized-matrix-generation.wdl](wdl/optimized-matrix-generation.wdl).

The Python/bash notebook for inspecting matrix outputs and merging sample batches is available in [notebooks/BioVU_methylation_processing.ipynb](notebooks/BioVU_methylation_processing.ipynb).

### Step 2: Cell-type deconvolution

Cell-type composition was estimated from methylation data using [MethylCC](https://bioconductor.org/packages/release/bioc/html/methylCC.html) and included as covariates (Gran, CD4T, CD8T, Bcell, Mono, NK) in downstream EWAS models to adjust for differential cell-type composition between CHIP carriers and controls.

---

## Epigenome-wide association study (EWAS) for CHIP status

To identify CpG sites differentially methylated in TET2 and DNMT3A CHIP carriers versus controls, we ran parallel logistic EWAS across all cohort samples. Models were adjusted for age, sex, race/ethnicity, and cell-type proportions.

The WDL workflow for scattered EWAS on Terra is available in [wdl/run-logistic-ewas.wdl](wdl/run-logistic-ewas.wdl). This workflow parallelizes across CpG chunks (default: 50,000 CpGs per chunk) and consolidates output summary statistics.

EWAS results visualization (volcano plots) is available in [notebooks/BioVU_methylation_clean.ipynb](notebooks/BioVU_methylation_clean.ipynb).

---

## Sliding window differentially methylated region analysis

To increase statistical power and biological interpretability, we complemented single-CpG EWAS with a sliding window DMR analysis aggregating methylation over 500 bp non-overlapping windows (as used in the in vitro EWAS) and a sliding window approach (2 kb windows, 1 kb step) for clinical outcome analyses.

The WDL workflow for sliding window DMR EWAS is available in [wdl/run-logistic-dmr-ewas.wdl](wdl/run-logistic-dmr-ewas.wdl).

---

## Activity score construction

To quantify individual-level enzymatic dysfunction from peripheral blood DNA methylation, we constructed TET2 and DNMT3A Activity Scores (AS) using a multi-step regularized regression pipeline.

### In vitro EWAS to identify DMRs

CRISPR-Cas9 was used to introduce TET2 loss-of-function and DNMT3A R882 mutations in primary human CD34+ HSPCs (N = 4–6 donors per group). 5-methylcytosine profiles were generated using the duet evoC assay (Twist Biosciences hybrid capture). An EWAS was performed on 500 bp non-overlapping windows, identifying:
- **TET2**: 27,465 DMRs (25,869 hypermethylated, 1,596 hypomethylated)
- **DNMT3A R882**: 5,288 DMRs (5,151 hypomethylated, 137 hypermethylated)

### Score training in the Framingham Heart Study

The FHS cohort provided peripheral blood Illumina HumanMethylation450 (450K array) data for score training. The workflow for FHS data preparation and cohort metadata is available in [notebooks/FHS_methylation_and_CHIP.ipynb](notebooks/FHS_methylation_and_CHIP.ipynb).

**Activity score construction pipeline** (R code):

```r
# ============================================================
# TET2 and DNMT3A Activity Score Construction
# ============================================================
# Requires:
#   - beta_matrix: CpG x sample matrix of DNA methylation beta values (0–1)
#   - dmr_sites: vector of CpG site IDs within in vitro DMRs
#   - array_cpgs: vector of CpG site IDs present on 450K or EPIC arrays
#   - pheno_df: data frame with columns: sample_id, chip_status (TET2_pLoF/DNMT3A_R882/control),
#               age, sex, and additional covariates
# ============================================================

library(glmnet)
library(Matrix)

# ------------------------------------------------------------------
# Step 1: Restrict to DMR CpGs present on standard methylation arrays
# ------------------------------------------------------------------
# Retain only CpGs within the in vitro DMRs and profiled on 450K/EPIC
# to ensure portability across clinical cohort datasets

activity_score_cpgs <- intersect(dmr_sites, array_cpgs)
beta_dmr <- beta_matrix[activity_score_cpgs, , drop = FALSE]

cat("CpGs in DMRs:", length(dmr_sites), "\n")
cat("CpGs on 450K/EPIC arrays:", length(array_cpgs), "\n")
cat("CpGs retained (intersection):", nrow(beta_dmr), "\n")


# ------------------------------------------------------------------
# Step 2: Correlation-based pruning to reduce multicollinearity
# ------------------------------------------------------------------
# DMRs are ranked by their EWAS association significance (ascending p-value).
# For each pair of CpGs with R^2 > 0.8, the less significant CpG is removed.
# This retains the most statistically robust representative from each
# correlated block.

prune_correlated_cpgs <- function(beta_mat, dmr_pvals, r2_threshold = 0.8) {
  # Order CpGs by significance (most significant first)
  ordered_idx <- order(dmr_pvals)
  beta_ordered <- beta_mat[ordered_idx, , drop = FALSE]
  pvals_ordered <- dmr_pvals[ordered_idx]

  n_cpgs <- nrow(beta_ordered)
  keep <- rep(TRUE, n_cpgs)

  cat("Pruning", n_cpgs, "CpGs at R^2 threshold =", r2_threshold, "...\n")

  for (i in seq_len(n_cpgs - 1)) {
    if (!keep[i]) next
    # Compute correlation of site i with all downstream sites
    for (j in seq(i + 1, n_cpgs)) {
      if (!keep[j]) next
      r2 <- cor(beta_ordered[i, ], beta_ordered[j, ], use = "pairwise.complete.obs")^2
      if (!is.na(r2) && r2 > r2_threshold) {
        keep[j] <- FALSE  # Remove the less significant CpG
      }
    }
  }

  pruned_cpgs <- rownames(beta_ordered)[keep]
  cat("Retained", length(pruned_cpgs), "CpGs after pruning\n")
  return(pruned_cpgs)
}

# Prune correlated CpGs (using EWAS p-values from in vitro study)
pruned_cpgs_tet2  <- prune_correlated_cpgs(beta_dmr_tet2,  dmr_pvals_tet2)
pruned_cpgs_dnmt3a <- prune_correlated_cpgs(beta_dmr_dnmt3a, dmr_pvals_dnmt3a)


# ------------------------------------------------------------------
# Step 3: Elastic Net regularized logistic regression (training)
# ------------------------------------------------------------------
# Models are trained in the Framingham Heart Study (FHS) discovery cohort.
#
# TET2 AS:    TET2 pLoF (N=74) vs. age/sex-matched controls (N=148)
# DNMT3A AS:  DNMT3A R882 (N=20) vs. age/sex-matched controls (N=40)
#
# Elastic Net mixing parameter alpha = 0.5 (equal L1/L2 penalty).
# Lambda is selected by 10-fold cross-validation to minimize deviance.

train_activity_score <- function(beta_mat,
                                  pheno_df,
                                  case_label,
                                  alpha = 0.5,
                                  nfolds = 10,
                                  seed = 42) {
  # Align samples
  shared_samples <- intersect(colnames(beta_mat), pheno_df$sample_id)
  beta_train  <- t(beta_mat[, shared_samples, drop = FALSE])  # samples x CpGs
  pheno_train <- pheno_df[match(shared_samples, pheno_df$sample_id), ]

  y <- as.integer(pheno_train$chip_status == case_label)
  cat("Training", case_label, "Activity Score:\n")
  cat("  Cases:", sum(y == 1), "| Controls:", sum(y == 0), "\n")

  # Scale CpG beta values to zero mean, unit variance
  beta_scaled <- scale(beta_train)
  scale_means <- attr(beta_scaled, "scaled:center")
  scale_sds   <- attr(beta_scaled, "scaled:scale")

  # Fit Elastic Net with cross-validation
  set.seed(seed)
  cv_fit <- cv.glmnet(
    x        = beta_scaled,
    y        = y,
    family   = "binomial",
    alpha    = alpha,
    nfolds   = nfolds,
    type.measure = "deviance",
    standardize = FALSE  # already scaled
  )

  # Extract non-zero coefficients at lambda.min
  coefs <- coef(cv_fit, s = "lambda.min")
  nonzero_coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]
  nonzero_coefs <- nonzero_coefs[-1, , drop = FALSE]  # remove intercept
  cat("  Non-zero CpGs in final model:", nrow(nonzero_coefs), "\n")
  cat("  Optimal lambda:", round(cv_fit$lambda.min, 6), "\n")

  return(list(
    model       = cv_fit,
    coefs       = nonzero_coefs,
    scale_means = scale_means,
    scale_sds   = scale_sds,
    lambda      = cv_fit$lambda.min
  ))
}

# Train TET2 Activity Score
tet2_model <- train_activity_score(
  beta_mat   = beta_dmr_tet2[pruned_cpgs_tet2, ],
  pheno_df   = fhs_pheno_df,
  case_label = "TET2_pLoF"
)
# Final model: 823 CpGs

# Train DNMT3A Activity Score
dnmt3a_model <- train_activity_score(
  beta_mat   = beta_dmr_dnmt3a[pruned_cpgs_dnmt3a, ],
  pheno_df   = fhs_pheno_df,
  case_label = "DNMT3A_R882"
)
# Final model: 542 CpGs


# ------------------------------------------------------------------
# Step 4: Score computation in validation cohorts (BioVU, CHIVE)
# ------------------------------------------------------------------
# Individual activity scores are computed as the dot product of
# standardized beta values and the non-zero model coefficients:
#
#   AS_i = sum_j ( (X_ij - mu_j) / sigma_j ) * w_j
#
# where X_ij is the beta value for sample i at CpG j,
# mu_j and sigma_j are the training-set mean and SD for CpG j,
# and w_j is the elastic net coefficient.
#
# Scores are then standardized to the non-CHIP control distribution
# (mean = 0, SD = 1) for all downstream analyses.

compute_activity_score <- function(beta_mat, model_obj, control_ids) {
  # Subset to CpGs in the model
  model_cpgs <- rownames(model_obj$coefs)
  beta_subset <- beta_mat[model_cpgs, , drop = FALSE]

  # Standardize using training-set parameters
  beta_scaled <- sweep(beta_subset, 1, model_obj$scale_means[model_cpgs], "-")
  beta_scaled <- sweep(beta_scaled, 1, model_obj$scale_sds[model_cpgs],   "/")

  # Dot product with model weights → raw activity score
  weights <- model_obj$coefs[, 1]
  raw_scores <- as.numeric(weights %*% beta_scaled)
  names(raw_scores) <- colnames(beta_mat)

  # Standardize to control population (mean=0, SD=1)
  control_mean <- mean(raw_scores[control_ids], na.rm = TRUE)
  control_sd   <- sd(raw_scores[control_ids],   na.rm = TRUE)
  std_scores   <- (raw_scores - control_mean) / control_sd

  return(std_scores)
}

# Compute TET2 Activity Scores in BioVU
tet2_scores_biovu <- compute_activity_score(
  beta_mat    = biovu_beta_matrix,
  model_obj   = tet2_model,
  control_ids = biovu_control_ids
)

# Compute DNMT3A Activity Scores in BioVU
dnmt3a_scores_biovu <- compute_activity_score(
  beta_mat    = biovu_beta_matrix,
  model_obj   = dnmt3a_model,
  control_ids = biovu_control_ids
)
```

---

## Validation in BioVU and CHIVE cohorts

Activity scores were validated in two independent cohorts with EM-seq data:

- **BioVU**: 356 non-CHIP controls, 207 TET2 CHIP, 371 DNMT3A CHIP
- **CHIVE**: 50 TET2 CHIP, 40 DNMT3A CHIP (independent clinical cohort)

Score distributions were compared across mutation subtypes using pairwise Wilcoxon rank-sum tests with Bonferroni correction.

Code for cohort metadata construction, score application, and in vitro EWAS is available in:
- [notebooks/EM_methylation_data.ipynb](notebooks/EM_methylation_data.ipynb) — In vitro CRISPR model data processing
- [notebooks/BioVU_methylation_processing.ipynb](notebooks/BioVU_methylation_processing.ipynb) — BioVU EM-seq processing and matrix construction
- [notebooks/BioVU_methylation_clean.ipynb](notebooks/BioVU_methylation_clean.ipynb) — Score application, EWAS results, and volcano plots

---

## Clinical risk model comparisons

We evaluated the TET2 and DNMT3A activity scores against:

- **Clonal Hematopoiesis Risk Score (CHRS)** ([Weeks et al., *NEJM Evidence*, 2023](https://evidence.nejm.org/doi/10.1056/EVIDoa2200310)) for prediction of incident persistent cytopenia
- **AHA PREVENT equations** ([Khan et al., *Circulation*, 2024](https://doi.org/10.1161/CIRCULATIONAHA.123.067626)) for prediction of incident atherosclerotic cardiovascular disease (ASCVD)

Model discrimination was quantified using AUROC for 5-year incident risk. Cox proportional hazards models were used for time-to-event analyses adjusted for age, sex, smoking, and the first 10 principal components of ancestry. ASCVD models additionally adjusted for BMI, LDL-C, diabetes, hypertension, and statin use.

CHRS was scored using: mutation count, driver gene, VAF, age, cytopenia, and cell morphology (low-risk: CHRS ≤ 9.5; intermediate: 10–12; high-risk: ≥ 12.5).

PREVENT 5-year risk was calculated using age, sex, systolic BP, BP treatment, smoking, diabetes, total cholesterol, HDL, and eGFR.

---

## Data

This analysis was performed on:
- [UK Biobank DNAnexus Research Analysis Platform](https://ukbiobank.dnanexus.com)
- [All of Us Research Workbench](https://workbench.researchallofus.org/)
- Vanderbilt BioVU Terra environment

Individual-level sequence data and CHIP calls are available to approved researchers through:
- **UK Biobank**: [ukbiobank.ac.uk/register-apply](https://www.ukbiobank.ac.uk/register-apply/)
- **All of Us**: [allofus.nih.gov](https://allofus.nih.gov/)
- **BioVU**: Controlled access via data use agreement with Vanderbilt University Medical Center (contact corresponding author)
- **BioVU methylation data**: dbGaP accession [phs004433.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs004433.v1.p1)
- **Framingham Heart Study methylation data**: dbGaP accession [phs000974.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000974.v1.p1)
- **CHIVE CHIP calls and methylation data**: Available upon request

---

## Citation

If you use this code or find it helpful, please cite:

> Pershad Y\*, Zhao K\*, Van Amburg JC, Corty RW, Parker AC, Silver AJ, Almadani Y, Kishtagari A, Hodges E, Savona MR, Heimlich JB†, Bick AG†. *Mutation-specific impairment of TET2 and DNMT3A enzymatic activity predicts clonal hematopoiesis disease risk.* medRxiv, 2026. \*Equal contribution. †Joint supervision.

---

## License

This project is licensed under the MIT License:

MIT License

Copyright (c) 2026 Yash Pershad

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## Acknowledgements

This work was supported by the following National Institutes of Health grants: UG3 AG097155 (A.G.B., M.R.S.), K08 HL171833 (J.B.H.), R01 AG088657 (A.G.B.), F30 AG099331 (Y.P.), and T32 GM145734 (J.C.V.A.). Additionally, this work was supported by a Burroughs Wellcome Fund Career Award for Medical Scientists (A.G.B.), a Pew-Stewart Scholar for Cancer Research award (A.G.B.), a Hevolution/AFAR New Investigator Award in Aging Biology and Geroscience Research (A.G.B.), and Arthritis National Research Foundation grant 1288083 (R.W.C.). The sequencing of 250,000 WGS individuals from BioVU® was funded by the Alliance for Genomic Discovery, consisting of NashBio, Illumina, and industry partners Amgen, AbbVie, AstraZeneca, Bayer, BMS, GSK, Merck, and Novo Nordisk.

---

## Contact

Yash Pershad, [yash.pershad@vanderbilt.edu](mailto:yash.pershad@vanderbilt.edu)  
Alexander G. Bick, [alexander.bick@vumc.org](mailto:alexander.bick@vumc.org)
