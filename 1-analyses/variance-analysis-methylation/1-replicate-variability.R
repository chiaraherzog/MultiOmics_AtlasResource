# Replicate variability
# Chiara Herzog
# Jul 28 2025

# 0. load libraries -----------
cat("Load Libs...")
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
cat("done.\n")

# 1. Load & rep replicate dataset ----------------
cat("Load data ... ")
load("~/Dropbox/data/tirolgesund/beta_merged_replicate.Rdata")
load("~/Dropbox/data/tirolgesund/pheno/methylation_replication_pheno.Rdata")
cat("done.\n")

# Match pheno to beta matrix
cat("Match pheno & data ... ")
pheno_replication <- pheno_replication[match(colnames(beta_merged), pheno_replication$basename),]
stopifnot(identical(pheno_replication$basename, colnames(beta_merged)))

# 2 | Replication in BETA VALUES -------------------
# Identify replicate groups
replicate_groups <- split(colnames(beta_merged), pheno_replication$sampleId)
cat("done.\n")

# Compute CpG-level technical SD from replicate M-values
replicate_sds <- sapply(replicate_groups, function(reps) {
  if (length(reps) > 1) rowSds(as.matrix(beta_merged[, reps]), na.rm = TRUE)
  else rep(NA, nrow(beta_merged))
})

# Compute CpG-level technical SD from replicate M-values
replicate_means <- sapply(replicate_groups, function(reps) {
  if (length(reps) > 1) rowMeans(as.matrix(beta_merged[, reps]), na.rm = TRUE)
  else rep(NA, nrow(beta_merged))
})

df_sd <- data.frame(cpg = rownames(replicate_sds),
                    mean = rowMeans(replicate_means),
                    sd = rowMeans(replicate_sds))

ind_blood <- grepl("BD", colnames(replicate_sds))
ind_buccal <- grepl("BS", colnames(replicate_sds))
ind_cerv <- grepl("CP", colnames(replicate_sds))
stopifnot(identical(colnames(replicate_means), colnames(replicate_sds)))
stopifnot(identical(rownames(replicate_means), rownames(replicate_sds)))

df_sd_blood <- data.frame(cpg = rownames(replicate_sds),
                          mean = rowMeans(replicate_means[,ind_blood]),
                          sd = rowMeans(replicate_sds[,ind_blood]))
df_sd_buccal <- data.frame(cpg = rownames(replicate_sds),
                           mean = rowMeans(replicate_means[,ind_buccal]),
                           sd = rowMeans(replicate_sds[,ind_buccal]))
df_sd_cerv <- data.frame(cpg = rownames(replicate_sds),
                         mean = rowMeans(replicate_means[,ind_cerv]),
                         sd = rowMeans(replicate_sds[,ind_cerv]))

cat("Save outputs (beta) ... ")
save(replicate_means, replicate_sds, df_sd, df_sd_blood, df_sd_buccal, df_sd_cerv, 
     file = '1-analyses/variance-analysis-methylation/1-output/beta-variability.Rdata')
cat("done.\n")
rm(df_sd, replicate_means, replicate_sds, df_sd_blood, df_sd_buccal, df_sd_cerv);gc()


# 3 | Replication in M VALUES -------------------

# Convert replicate beta matrix to M-values (for variance modeling)
cat("Transform beta to M and compute Sds & means ... ")
beta2M <- function(x){
  x[x == 0] <- min(x[x != 0])
  x[x == 1] <- max(x[x != 1])
  log2(x) - log2(1 - x)
}
mval_replicate <- beta2M(beta_merged)

# Identify replicate groups
replicate_groups <- split(colnames(mval_replicate), pheno_replication$sampleId)

# Compute CpG-level technical SD from replicate M-values
replicate_sds <- sapply(replicate_groups, function(reps) {
  if (length(reps) > 1) rowSds(as.matrix(mval_replicate[, reps]), na.rm = TRUE)
  else rep(NA, nrow(mval_replicate))
})

# Compute CpG-level technical SD from replicate M-values
replicate_means <- sapply(replicate_groups, function(reps) {
  if (length(reps) > 1) rowMeans(as.matrix(mval_replicate[, reps]), na.rm = TRUE)
  else rep(NA, nrow(mval_replicate))
})

tech_sd <- rowMeans(replicate_sds, na.rm = TRUE)
cat("done.\n")

cat("Save outputs ... ")
save(tech_sd, replicate_sds, file = '1-analyses/variance-analysis-methylation/1-output/replication_variability.Rdata')
save(replicate_means, file = '1-analyses/variance-analysis-methylation/1-output/replication_mean.Rdata')
cat("done.\n")

# 2. Compute means and SDs per tissue -------------------
cat("By tissue dfs ... ")
ind_blood <- grepl("BD", colnames(replicate_sds))
ind_buccal <- grepl("BS", colnames(replicate_sds))
ind_cerv <- grepl("CP", colnames(replicate_sds))
stopifnot(identical(colnames(replicate_means), colnames(replicate_sds)))
stopifnot(identical(rownames(replicate_means), rownames(replicate_sds)))

df_sd_blood <- data.frame(cpg = rownames(replicate_sds),
                          mean = rowMeans(replicate_means[,ind_blood]),
                          sd = rowMeans(replicate_sds[,ind_blood]))
df_sd_buccal <- data.frame(cpg = rownames(replicate_sds),
                           mean = rowMeans(replicate_means[,ind_buccal]),
                           sd = rowMeans(replicate_sds[,ind_buccal]))
df_sd_cerv <- data.frame(cpg = rownames(replicate_sds),
                         mean = rowMeans(replicate_means[,ind_cerv]),
                         sd = rowMeans(replicate_sds[,ind_cerv]))
cat("done.\n")

cat("Save outputs of by tissue dfs... ")
save(df_sd_blood,
     df_sd_buccal, df_sd_cerv,
     file = '1-analyses/variance-analysis-methylation/1-output/mean_sd_dfs_by_tissue.Rdata')
cat("done.\n")