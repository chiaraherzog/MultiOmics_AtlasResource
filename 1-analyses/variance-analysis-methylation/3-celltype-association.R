# CpG association with cell type
# Chiara Herzog
# Jul 28 2025

# 0. load libraries -----------
cat("Load Libs ... ")
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EpiDISH))
suppressPackageStartupMessages(library(pbmcapply))
suppressPackageStartupMessages(library(parallel))
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

# Compute (h)EpiDISH
cat("1.2 | Compute hEpidish ... ")
res2 <- EpiDISH::hepidish(beta_merged, ref1.m = centEpiFibIC.m, ref2.m = cent12CT.m, method = 'RPC',
                          h.CT.idx = 3)
colnames(res2) <- paste0("celltype_", colnames(res2))
pheno <- cbind(pheno, res2)
cat("done.\n")

# Subset by CpGs overlapping with rep set
cat("1.3 | Subset dataset to match replication set... ")
intersect_cpgs <- intersect(rownames(beta_merged), rownames(replicate_sds))
cat(length(intersect_cpgs), 'unique CpGs included')
# intersect_cpgs <- intersect_cpgs[1:3000]  # pilot subset
beta_matrix <- beta_merged[intersect_cpgs, ]
cat("done.\n")

# 5. Cell type classification ---------------------
cat("2.1 | Explore cell type association... ")
n_cores <- parallel::detectCores()/2
n_cpgs <- nrow(beta_matrix)
n <- ncol(beta_matrix)

# Prepare input for each CpG (row index)
all_idxs <- seq_len(n_cpgs)
celltype_matrix <- as.matrix(pheno[, grep("^celltype_", names(pheno))])

# Define function to compute correlation for one CpG
get_best_celltype <- function(i) {
  vals <- beta_matrix[i, ]
  cors <- cor(vals, celltype_matrix, use = "pairwise.complete.obs")
  t_vals <- cors * sqrt((n - 2) / (1 - cors^2))
  p_vals <- 2 * pt(-abs(t_vals), df = n - 2)
  best_idx <- which.max(abs(cors))
  c(cor = cors[best_idx],
    pval = p_vals[best_idx],
    celltype = colnames(celltype_matrix)[best_idx])
}

res_list <- pbmclapply(all_idxs, get_best_celltype, mc.cores = n_cores)

# Convert to data.frame
res_mat <- do.call(rbind, res_list)
celltype_stats <- as.data.frame(res_mat)
celltype_stats$max_cor <- as.numeric(celltype_stats$cor)
celltype_stats$min_p <- as.numeric(celltype_stats$pval)
celltype_stats$best_celltype <- as.character(celltype_stats$celltype)
celltype_stats <- celltype_stats[, c("max_cor", "min_p", "best_celltype")]
rownames(celltype_stats) <- rownames(beta_matrix)

save(celltype_stats, pheno, file = '1-analyses/variance-analysis-methylation/1-output/3-output/celltype_association.Rdata')


# celltype_stats |> 
#   ggplot(aes(x = max_cor,
#              y = -log10(min_p))) +
#   geom_hex()
