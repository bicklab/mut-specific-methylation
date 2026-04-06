# ============================================================
# Activity Score Construction for TET2 and DNMT3A CHIP
# Pershad & Zhao et al., Nature Medicine, 2025
# ============================================================
#
# This script implements the full pipeline for constructing
# TET2 and DNMT3A methylation-based Activity Scores as described
# in the paper. The pipeline proceeds in four steps:
#
#   Step 1: Restrict to DMR CpGs present on standard methylation arrays
#   Step 2: Correlation-based pruning of redundant CpGs
#   Step 3: Elastic Net logistic regression in the FHS discovery cohort
#   Step 4: Score computation in validation cohorts (BioVU, CHIVE)
#
# Required inputs:
#   - beta_matrix_fhs.rds       : CpG x sample beta matrix (FHS, 450K array)
#   - beta_matrix_biovu.rds     : CpG x sample beta matrix (BioVU, EM-seq)
#   - tet2_dmr_cpgs.txt         : CpG IDs within TET2 LoF in vitro DMRs
#   - dnmt3a_dmr_cpgs.txt       : CpG IDs within DNMT3A R882 in vitro DMRs
#   - illumina_array_cpgs.txt   : CpG IDs on Illumina 450K or EPIC arrays
#   - fhs_pheno.tsv             : FHS sample metadata (sample_id, chip_status, age, sex)
#   - biovu_pheno.tsv           : BioVU sample metadata
#
# Output:
#   - tet2_activity_score_model.rds    : trained elastic net model for TET2 AS
#   - dnmt3a_activity_score_model.rds  : trained elastic net model for DNMT3A AS
#   - biovu_activity_scores.tsv        : per-sample activity scores in BioVU
# ============================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(data.table)
  library(dplyr)
  library(Matrix)
})

set.seed(42)

# ============================================================
# Helper functions
# ============================================================

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}

# ------------------------------------------------------------------
# Step 1: Restrict to DMR CpGs present on standard methylation arrays
# ------------------------------------------------------------------
# We intersect:
#   (a) CpGs within the 500 bp DMRs identified in the in vitro CRISPR model EWAS
#   (b) CpGs profiled on the Illumina HumanMethylation450 or EPIC array
# This ensures the score is portable across clinical datasets that use these arrays.

restrict_to_array_dmr_cpgs <- function(dmr_cpg_file, array_cpg_file, beta_matrix) {
  dmr_cpgs   <- readLines(dmr_cpg_file)
  array_cpgs <- readLines(array_cpg_file)

  # Retain CpGs in both DMRs and array
  candidate_cpgs <- intersect(dmr_cpgs, array_cpgs)
  # Further restrict to CpGs actually present in the beta matrix
  candidate_cpgs <- intersect(candidate_cpgs, rownames(beta_matrix))

  log_step(sprintf(
    "CpGs in DMRs: %d | On array: %d | In beta matrix: %d | Retained: %d",
    length(dmr_cpgs), length(array_cpgs),
    length(intersect(dmr_cpgs, rownames(beta_matrix))),
    length(candidate_cpgs)
  ))

  return(candidate_cpgs)
}


# ------------------------------------------------------------------
# Step 2: Correlation-based pruning of CpG redundancy
# ------------------------------------------------------------------
# DMRs are ranked by EWAS association significance (p-value from in vitro study).
# For any pair of CpGs with Pearson R² > 0.8, the less significant CpG is removed.
# This retains the most statistically robust representative from each correlated block
# and reduces multicollinearity prior to elastic net regression.

prune_correlated_cpgs <- function(beta_mat, cpg_pvals, r2_threshold = 0.8) {
  # Sort by significance (most significant first)
  cpgs_ordered <- names(sort(cpg_pvals))
  cpgs_ordered <- intersect(cpgs_ordered, rownames(beta_mat))
  n_cpgs       <- length(cpgs_ordered)

  log_step(sprintf("Pruning %d CpGs at R² threshold = %.2f...", n_cpgs, r2_threshold))

  keep <- rep(TRUE, n_cpgs)
  names(keep) <- cpgs_ordered

  # Compute correlation matrix in blocks to avoid memory issues
  beta_subset <- t(beta_mat[cpgs_ordered, , drop = FALSE])  # samples x CpGs

  for (i in seq_len(n_cpgs - 1)) {
    if (!keep[cpgs_ordered[i]]) next
    for (j in seq(i + 1, n_cpgs)) {
      if (!keep[cpgs_ordered[j]]) next
      # Pairwise Pearson correlation
      r <- cor(beta_subset[, i], beta_subset[, j], use = "pairwise.complete.obs")
      if (!is.na(r) && r^2 > r2_threshold) {
        keep[cpgs_ordered[j]] <- FALSE
      }
    }
  }

  retained <- names(keep[keep])
  log_step(sprintf("Retained %d CpGs after pruning (removed %d)", length(retained), n_cpgs - length(retained)))
  return(retained)
}


# ------------------------------------------------------------------
# Step 3: Elastic Net regularized logistic regression (training)
# ------------------------------------------------------------------
# Models are trained in the Framingham Heart Study (FHS) discovery cohort.
#
# TET2 Activity Score:
#   Cases: TET2 pLoF (N=74), Controls: age/sex-matched non-CHIP (N=148)
#   Final model: 823 weighted CpGs
#
# DNMT3A Activity Score:
#   Cases: DNMT3A R882 (N=20), Controls: age/sex-matched non-CHIP (N=40)
#   Final model: 542 weighted CpGs
#
# Elastic Net parameters: alpha=0.5 (equal L1/L2 penalty), lambda by 10-fold CV

train_activity_score_model <- function(beta_mat,
                                        pheno_df,
                                        case_label,
                                        sample_id_col = "sample_id",
                                        status_col    = "chip_status",
                                        alpha         = 0.5,
                                        nfolds        = 10) {
  # Align samples between beta matrix and phenotype data
  shared_samples <- intersect(colnames(beta_mat), pheno_df[[sample_id_col]])
  if (length(shared_samples) < 20) {
    stop("Fewer than 20 shared samples between beta matrix and phenotype data.")
  }

  beta_train  <- beta_mat[, shared_samples, drop = FALSE]  # CpGs x samples
  pheno_train <- pheno_df[match(shared_samples, pheno_df[[sample_id_col]]), ]

  # Binary outcome: 1 = case, 0 = control
  y <- as.integer(pheno_train[[status_col]] == case_label)

  log_step(sprintf(
    "Training %s Activity Score | Cases: %d | Controls: %d | CpGs: %d",
    case_label, sum(y == 1), sum(y == 0), nrow(beta_train)
  ))

  # Standardize CpG beta values: (x - mean) / sd
  # Scaling parameters saved for application to validation cohorts
  cpg_means <- rowMeans(beta_train, na.rm = TRUE)
  cpg_sds   <- apply(beta_train, 1, sd, na.rm = TRUE)
  cpg_sds[cpg_sds < 1e-6] <- 1e-6  # avoid division by zero for invariant sites

  beta_scaled <- (beta_train - cpg_means) / cpg_sds
  beta_scaled[is.na(beta_scaled)] <- 0  # impute missing with 0 (= mean after scaling)
  X <- t(beta_scaled)  # samples x CpGs for glmnet

  # 10-fold cross-validated elastic net (alpha=0.5 balances LASSO and Ridge)
  log_step("Running 10-fold cross-validated elastic net...")
  cv_fit <- cv.glmnet(
    x            = X,
    y            = y,
    family       = "binomial",
    alpha        = alpha,
    nfolds       = nfolds,
    type.measure = "deviance",
    standardize  = FALSE  # already standardized
  )

  # Extract non-zero coefficients at the cross-validation-optimal lambda
  coefs_all    <- coef(cv_fit, s = "lambda.min")
  coefs_nonzero <- coefs_all[-1, , drop = FALSE]  # remove intercept
  coefs_nonzero <- coefs_nonzero[coefs_nonzero[, 1] != 0, , drop = FALSE]

  log_step(sprintf(
    "Model trained | Lambda: %.6f | Non-zero CpGs: %d | CV deviance: %.4f",
    cv_fit$lambda.min,
    nrow(coefs_nonzero),
    min(cv_fit$cvm)
  ))

  return(list(
    cv_model    = cv_fit,
    coefficients = coefs_nonzero,
    cpg_means   = cpg_means,
    cpg_sds     = cpg_sds,
    lambda      = cv_fit$lambda.min,
    alpha       = alpha,
    case_label  = case_label,
    n_cases     = sum(y == 1),
    n_controls  = sum(y == 0)
  ))
}


# ------------------------------------------------------------------
# Step 4: Score computation in validation cohorts
# ------------------------------------------------------------------
# Individual activity scores are computed as the dot product of
# standardized CpG beta values and the elastic net model weights:
#
#   AS_i = Σ_j  w_j * (X_ij - μ_j) / σ_j
#
# where:
#   X_ij = beta value for sample i at CpG j
#   μ_j  = training-set mean for CpG j
#   σ_j  = training-set SD for CpG j
#   w_j  = elastic net coefficient
#
# Scores are standardized to the non-CHIP control population (mean=0, SD=1).

compute_activity_scores <- function(beta_mat, model_obj, control_ids) {
  # Subset to CpGs in the model
  model_cpgs   <- rownames(model_obj$coefficients)
  missing_cpgs <- setdiff(model_cpgs, rownames(beta_mat))
  if (length(missing_cpgs) > 0) {
    warning(sprintf("%d model CpGs missing from validation beta matrix; setting to 0.", length(missing_cpgs)))
  }

  # Build CpG x sample matrix for model CpGs (impute missing with 0 = training mean)
  available_cpgs <- intersect(model_cpgs, rownames(beta_mat))
  beta_sub <- matrix(0, nrow = length(model_cpgs), ncol = ncol(beta_mat),
                     dimnames = list(model_cpgs, colnames(beta_mat)))
  beta_sub[available_cpgs, ] <- beta_mat[available_cpgs, , drop = FALSE]

  # Standardize using training-set parameters
  cpg_means <- model_obj$cpg_means[model_cpgs]
  cpg_sds   <- model_obj$cpg_sds[model_cpgs]
  beta_scaled <- (beta_sub - cpg_means) / cpg_sds
  beta_scaled[is.na(beta_scaled)] <- 0

  # Dot product with model weights → raw activity score per sample
  weights     <- model_obj$coefficients[, 1]
  raw_scores  <- as.numeric(weights %*% beta_scaled)
  names(raw_scores) <- colnames(beta_mat)

  # Standardize to control population (mean=0, SD=1)
  control_mean <- mean(raw_scores[control_ids], na.rm = TRUE)
  control_sd   <- sd(raw_scores[control_ids],   na.rm = TRUE)

  if (control_sd < 1e-6) stop("Control SD is near zero; check control_ids.")

  std_scores <- (raw_scores - control_mean) / control_sd

  log_step(sprintf(
    "Scores computed | N=%d samples | Control mean=%.3f SD=%.3f (before standardization)",
    length(std_scores), control_mean, control_sd
  ))

  return(std_scores)
}


# ============================================================
# Main pipeline
# ============================================================

main <- function() {

  # --- Load inputs ---

  log_step("Loading input data...")

  # Beta matrices (CpGs x samples, values in [0,1])
  beta_fhs   <- readRDS("beta_matrix_fhs.rds")     # FHS 450K array
  beta_biovu <- readRDS("beta_matrix_biovu.rds")    # BioVU EM-seq

  # DMR CpG lists (output of in vitro EWAS, restricted to 450K/EPIC probes)
  tet2_dmr_cpgs  <- readLines("tet2_dmr_array_cpgs.txt")
  d3a_dmr_cpgs   <- readLines("dnmt3a_dmr_array_cpgs.txt")
  array_cpgs     <- readLines("illumina_array_cpgs.txt")

  # In vitro EWAS p-values (for pruning; lower = more significant)
  tet2_pvals <- fread("tet2_invitro_ewas_pvals.tsv")   # columns: cpg, pval
  d3a_pvals  <- fread("dnmt3a_invitro_ewas_pvals.tsv")

  tet2_pval_vec <- setNames(tet2_pvals$pval, tet2_pvals$cpg)
  d3a_pval_vec  <- setNames(d3a_pvals$pval,  d3a_pvals$cpg)

  # Phenotype data
  fhs_pheno   <- fread("fhs_tet2_plof_training_metadata.tsv")   # used for TET2 training
  fhs_pheno_d3a <- fread("fhs_d3a_r882_training_metadata.tsv")  # used for DNMT3A training
  biovu_pheno <- fread("biovu_pheno.tsv")

  # BioVU control sample IDs (for standardization)
  biovu_control_ids <- biovu_pheno[chip_status == "Control", sample_id]


  # --- TET2 Activity Score ---

  log_step("=== TET2 Activity Score ===")

  # Step 1: Restrict to array DMR CpGs
  tet2_candidate_cpgs <- restrict_to_array_dmr_cpgs(
    "tet2_dmr_array_cpgs.txt", "illumina_array_cpgs.txt", beta_fhs
  )

  # Step 2: Prune correlated CpGs
  tet2_pruned_cpgs <- prune_correlated_cpgs(
    beta_mat   = beta_fhs[tet2_candidate_cpgs, ],
    cpg_pvals  = tet2_pval_vec[tet2_candidate_cpgs],
    r2_threshold = 0.8
  )

  # Step 3: Train elastic net in FHS
  tet2_model <- train_activity_score_model(
    beta_mat    = beta_fhs[tet2_pruned_cpgs, ],
    pheno_df    = fhs_pheno,
    case_label  = "TET2_pLoF"
  )
  # Expected: 823 CpGs in final model

  saveRDS(tet2_model, "tet2_activity_score_model.rds")
  log_step("TET2 model saved.")


  # --- DNMT3A Activity Score ---

  log_step("=== DNMT3A Activity Score ===")

  # Step 1
  d3a_candidate_cpgs <- restrict_to_array_dmr_cpgs(
    "dnmt3a_dmr_array_cpgs.txt", "illumina_array_cpgs.txt", beta_fhs
  )

  # Step 2
  d3a_pruned_cpgs <- prune_correlated_cpgs(
    beta_mat   = beta_fhs[d3a_candidate_cpgs, ],
    cpg_pvals  = d3a_pval_vec[d3a_candidate_cpgs],
    r2_threshold = 0.8
  )

  # Step 3
  d3a_model <- train_activity_score_model(
    beta_mat    = beta_fhs[d3a_pruned_cpgs, ],
    pheno_df    = fhs_pheno_d3a,
    case_label  = "DNMT3A_R882"
  )
  # Expected: 542 CpGs in final model

  saveRDS(d3a_model, "dnmt3a_activity_score_model.rds")
  log_step("DNMT3A model saved.")


  # --- Apply scores to BioVU validation cohort ---

  log_step("=== Applying scores to BioVU ===")

  tet2_scores  <- compute_activity_scores(beta_biovu, tet2_model,  biovu_control_ids)
  d3a_scores   <- compute_activity_scores(beta_biovu, d3a_model,   biovu_control_ids)

  # Compile output table
  scores_df <- biovu_pheno[, .(sample_id, chip_gene, mutation_class, vaf)] %>%
    as.data.frame() %>%
    mutate(
      tet2_activity_score  = tet2_scores[sample_id],
      dnmt3a_activity_score = d3a_scores[sample_id]
    )

  fwrite(scores_df, "biovu_activity_scores.tsv", sep = "\t")
  log_step(sprintf("Activity scores written for %d BioVU samples.", nrow(scores_df)))


  # --- Summary statistics ---

  log_step("=== Summary ===")
  cat(sprintf("TET2 Activity Score:\n  CpGs: %d\n  Lambda: %.6f\n  N cases: %d\n  N controls: %d\n",
              nrow(tet2_model$coefficients),
              tet2_model$lambda,
              tet2_model$n_cases,
              tet2_model$n_controls))
  cat(sprintf("DNMT3A Activity Score:\n  CpGs: %d\n  Lambda: %.6f\n  N cases: %d\n  N controls: %d\n",
              nrow(d3a_model$coefficients),
              d3a_model$lambda,
              d3a_model$n_cases,
              d3a_model$n_controls))

  log_step("Pipeline complete.")
}

main()
