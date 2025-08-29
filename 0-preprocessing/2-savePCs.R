# Compute methylation PCs
library(dplyr)
library(MultiAssayExperiment)

# Prep function
computePCs <- function(beta, annofile, n = 50000,
                       threshold = 0.9, 
                       file = ''){
  
  if(file == ''){
    stop("please provide a filename in the file input option")
  }
  
  # 1. Prep annofile
  annofile <- annofile |> dplyr::mutate(primary = paste0(subjectId, visitId))
  
  # 2. Prep beta
  beta <- beta[,annofile$basename]
  if(!identical(colnames(beta), annofile$basename)){
    stop("names not identical - please check")
  }
  
  # 3. Compute variable features
  cat("Getting top", n, "variable CpGs...")
  topVar <- rowSds(beta)
  topVar <- topVar[order(abs(topVar), decreasing = T)]
  topVar <- topVar[1:n]
  cat("done\n")
  
  # 4. Subset beta
  cat("Subsetting beta...")
  beta <- beta[names(topVar),]
  cat("done\n")
  
  # 5. Compute PCs
  cat("Computing PCs...")
  pc <- prcomp(t(beta),
               center = TRUE, scale. = TRUE)
  cat("done\n")
  
  # 6. Get variance explained:
  var_explained <- pc$sdev^2
  prop_var <- var_explained / sum(var_explained)
  cum_prop_var <- cumsum(prop_var) # Cumulative variance
  k <- which(cum_prop_var >= threshold)[1]
  cat("Keeping the first", k, "PCs (", round(cum_prop_var[k] * 100, 1), "% of variance).\n")
  
  # 7. Visualise
  pdf(file = paste0("0-preprocessing/2-output/PCs-", file, ".pdf"),
      width = 4.5, height = 4.5)
  plot(cum_prop_var,
       type = "b", xlab = "Principal Component", ylab = "Cumulative PVE")
  abline(h = threshold, lty = 2) 
  dev.off()
  
  # 8. Format and return results
  cat("Returning results.")
  pcs <- pc$x |> 
    as.data.frame() |> 
    dplyr::select(1:k) |> 
    tibble::rownames_to_column('basename') |> 
    dplyr::left_join(dplyr::select(annofile,
                                   basename, primary),
                     by = 'basename') |> 
    dplyr::relocate("primary", .after = 'basename')
  
  # 9. Get loadings (interpretation)
  pc_weights <- pc$rotation
  
  saveRDS(pcs, file = paste0("0-preprocessing/2-output/PCs-", file, ".Rds"))
  
  saveRDS(pc_weights, file = paste0("0-preprocessing/2-output/PC-loading-", file, ".Rds"))
  
}

# LOAD DATA ----------------
load("~/Dropbox/data/tirolgesund/beta_merged.Rdata")

load("data/data_raw.Rdata")
meta <- as.data.frame(colData(data)) |> 
  dplyr::select(subjectId, visitId,
                basename_blood,
                basename_cervical,
                basename_buccal) |> 
  tidyr::pivot_longer(contains("basename"),
                      names_to = 'sample',
                      values_to = 'basename') |> 
  dplyr::filter(!is.na(basename))



# BASELINE ----------------
base <- meta |> dplyr::filter(visitId == 'M0')
base_blood <- base |> dplyr::filter(grepl("blood", sample))
base_cervical <- base |> dplyr::filter(grepl("cervical", sample))
base_buccal <- base |> dplyr::filter(grepl("buccal", sample))

computePCs(beta = beta_merged,
           annofile = base_blood,
           file = 'base-blood',
           n = 200000,
           threshold = 0.95)

computePCs(beta = beta_merged,
           annofile = base_cervical,
           file = 'base-cervical',
           n = 200000,
           threshold = 0.95)

computePCs(beta = beta_merged,
           annofile = base_buccal,
           file = 'base-buccal',
           n = 200000,
           threshold = 0.95)

# x <- readRDS("0-preprocessing/2-output/PCs-base-blood.Rds")
# x <- readRDS("0-preprocessing/2-output/PCs-base-cervical.Rds")
# x <- readRDS("0-preprocessing/2-output/PCs-base-buccal.Rds")
# plot(x$PC1,
#      x$PC2)

# LONGITUDINAL ----------------
long <- meta |> dplyr::filter(visitId %in% c('M0', "M2", "M4", "M6"))
long_blood <- long |> dplyr::filter(grepl("blood", sample))
long_cervical <- long |> dplyr::filter(grepl("cervical", sample))
long_buccal <- long |> dplyr::filter(grepl("buccal", sample))

computePCs(beta = beta_merged,
           annofile = long_blood,
           file = 'long-blood',
           n = 200000,
           threshold = 0.95)

computePCs(beta = beta_merged,
           annofile = long_cervical,
           file = 'long-cervical',
           n = 200000,
           threshold = 0.95)

computePCs(beta = beta_merged,
           annofile = long_buccal,
           file = 'long-buccal',
           n = 200000,
           threshold = 0.95)