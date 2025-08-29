# CpG association - overall variability
# Chiara Herzog
# Jul 28 2025

# 0. load libraries -----------
cat("Load Libs ... ")
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EpiDISH))
beta2M <- function(x){
  x[x == 0] <- min(x[x != 0])
  x[x == 1] <- max(x[x != 1])
  log2(x) - log2(1 - x)
}
cat("done.\n")

# 1. Load full methylation dataset ------------------
cat("1.0 | Load data ... ")
load('0-preprocessing/3-output/replication_variability.Rdata')
load("~/Dropbox/data/tirolgesund/pheno/methylation_pheno.Rdata")
load("~/Dropbox/data/tirolgesund/beta_merged.Rdata")
cat("done.\n")

# Match metadata to full beta matrix
cat("1.1 | Match data to pheno... ")
pheno <- pheno[match(colnames(beta_merged), pheno$basename), ]
stopifnot(identical(colnames(beta_merged), pheno$basename))
cat("done.\n")

# Subset by CpGs overlapping with rep set
cat("1.2 | Subset dataset to match replication set... ")
intersect_cpgs <- intersect(rownames(beta_merged), rownames(replicate_sds))
cat(length(intersect_cpgs), 'unique CpGs included')
# intersect_cpgs <- intersect_cpgs[1:3000]  # pilot subset
beta_matrix <- beta_merged[intersect_cpgs, ]
mval_matrix <- beta2M(beta_matrix)
tech_sd <- tech_sd[intersect_cpgs]
cat("done.\n")

# 2. Filter CpGs with zero variance -------------------
cat("2.1 | Compute variability ... \n")

group_variance_mean <- function(mat, group_vector) {
  group_levels <- unique(group_vector)
  group_vars <- matrix(NA_real_, nrow = nrow(mat), ncol = length(group_levels))
  colnames(group_vars) <- group_levels
  
  pb <- txtProgressBar(min = 0, max = length(group_levels), style = 3)
  
  for (i in seq_along(group_levels)) {
    group_name <- group_levels[i]
    idx <- which(group_vector == group_name)
    if (length(idx) > 1) {  # Need >1 sample to compute variance
      group_vars[, i] <- rowSds(mat[, idx, drop = FALSE], na.rm = TRUE)
    } else {
      group_vars[, i] <- NA  # Not enough replicates
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  rowMeans(group_vars, na.rm = TRUE)
}
cat("2.2 | Compute variability by subjectId... \n")
var_by_individual <- group_variance_mean(mval_matrix, pheno$subjectId)
cat("2.2 | Compute variability by tissue \n")
var_by_tissue     <- group_variance_mean(mval_matrix, pheno$sampletype)

overall_sd        <- rowSds(mval_matrix, na.rm = TRUE)
keep_cpgs <- overall_sd > tech_sd | var_by_individual > tech_sd | var_by_tissue > tech_sd

# Track all CpGs for reintegration
all_cpgs <- rownames(beta_matrix)
non_variable_cpgs <- all_cpgs[!keep_cpgs]

cat("done.\n")

cat(sum(keep_cpgs), 'CpGs retained,\n', length(non_variable_cpgs), 'exhibited zero variability, not included in variance partition.')
save(overall_sd, var_by_individual, var_by_tissue, keep_cpgs, non_variable_cpgs,
     file = '1-analyses/variance-analysis-methylation/2-output/variability_overall.Rdata')