# Variance partition
# Chiara Herzog
# Jul 28 2025

# 0. load libraries -----------
cat("Load Libs ... ")
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(dplyr))
beta2M <- function(x){
  x[x == 0] <- min(x[x != 0])
  x[x == 1] <- max(x[x != 1])
  log2(x) - log2(1 - x)
}
cat("done.\n")

# 1. Load full methylation dataset ------------------
cat("1.0 | Load data ... ")
load('1-analyses/variance-analysis-methylation/1-output/replication_variability.Rdata') # includes rep SD from 3-replicate-variability.R
load("1-analyses/variance-analysis-methylation/3-output/celltype_association.Rdata") # includes cell type matrix and pheno from 5-celltype-association.R
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
beta_matrix <- beta_merged[intersect_cpgs, ]
mval_matrix <- beta2M(beta_matrix)
tech_sd <- tech_sd[intersect_cpgs]
cat("done.\n")

# 2. Filter CpGs with zero variance -------------------
load('1-analyses/variance-analysis-methylation/2-output/variability_overall.Rdata') # output of 4-variability

# Filter matrices
cat("2 | Filtering matrices ... ")
mval_matrix <- mval_matrix[keep_cpgs, ]
beta_matrix <- beta_matrix[keep_cpgs, ]
tech_sd_keep <- tech_sd[keep_cpgs]
cat("done.\n")

# 3. Variance partitioning ---------------
cat("3.1 | Start variance partition ... ")
pheno$studyarm <- ifelse(pheno$interventionId == 'S', 'S', 'I')

form <- ~ 
  (1 | subjectId) + 
  (1 | sampletype) + 
  (1 | visitId) + 
  (1 | studyarm) + 
  (1 | studyarm:visitId) + 
  (1 | studyarm:sampletype) +
  (1 | sampletype:visitId) + 
  (1 | sampletype:visitId:studyarm)


# Visualise correlation
C <- canCorPairs(form, pheno)
plotCorrMatrix(C)


# Use MulticoreParam 
param <- if (.Platform$OS.type == "unix") {
  MulticoreParam(workers = 8)
} else {
  SnowParam(workers = 8, type = "SOCK")
}

# Chunked variance partitioning
chunk_size <- 20000
chunks <- split(1:nrow(mval_matrix), ceiling(seq_along(1:nrow(mval_matrix)) / chunk_size))

start_time <- Sys.time()
message("Starting variancePartition analysis for ", nrow(mval_matrix), " CpGs in ", length(chunks), " chunks...")

# Clear memory by saving and removing each chunk
pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)
for (i in seq_along(chunks)) {
  idx <- chunks[[i]]
  # message("Running variance partitioning on chunk ", i, "/", length(chunks), ": ", length(idx), " CpGs")
  chunk_mval <- mval_matrix[idx, , drop = FALSE]
  chunk_vp <- fitExtractVarPartModel(chunk_mval, form, pheno, BPPARAM = param)
  saveRDS(chunk_vp, file = paste0("1-analyses/variance-analysis-methylation/4-output/vpExt_results_chunk", i, ".rds"))
  rm(chunk_mval, chunk_vp); gc()
  setTxtProgressBar(pb, i)
}

close(pb)

cat("done.\n")

vp_files <- list.files("1-analyses/variance-analysis-methylation/4-output/", pattern = "vpExt_results_chunk.*\\.rds$", full.names = TRUE)
vp_list <- lapply(vp_files, readRDS)
vp_combined <- do.call(rbind, vp_list)
saveRDS(vp_combined, file = "1-analyses/variance-analysis-methylation/4-output/vpExt_combined_all.rds")