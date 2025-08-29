# CpG Classification by Variance - unifying outputs from multiple scripts
## Bringing together cell type, overall variability, technical noise (rep variability)
# Chiara Herzog
# Jul 28 2025

# 0. load libraries -----------
cat("Load Libs ... ")
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
cat("done.\n")

# 1. Load various data ------------------
cat("1.0 | Load data ... ")
load('1-analyses/variance-analysis-methylation/1-output/replication_variability.Rdata') # includes rep SD from 3-replicate-variability.R
load("1-analyses/variance-analysis-methylation/3-output/celltype_association.Rdata") # includes cell type matrix and pheno from 5-celltype-association.R
load('1-analyses/variance-analysis-methylation/2-output/variability_overall.Rdata') # output of 4-variability
tech_sd_keep <- tech_sd[keep_cpgs]
cat("done.\n")

# 2. Load VP data ---------------
cat("2.1 | Load VP data ... ")
vp_results <- readRDS("1-analyses/variance-analysis-methylation/4-output/vpExt_combined_all.rds")
# Compute total variance fraction
cgs <- rownames(vp_results)

overall_sd <- overall_sd[cgs]
tech_sd <- tech_sd[cgs]
identical(names(overall_sd), names(tech_sd))
tech_var_fraction <- tech_sd^2 / overall_sd

# 3. Classification ---------------------
cat("3.1 | Classify VP data ... ")
vp_results$thresh <- tech_var_fraction[rownames(vp_results)]

classify_vp <- function(vp_results, thresh_col = "thresh",
                                major_cutoff = 0.5) {
  
  # Component names to classify
  comp_cols <- c(
    "sampletype:visitId:studyarm",
    "sampletype:visitId",
    "studyarm:visitId",
    "visitId",
    "studyarm:sampletype",
    "studyarm",
    "subjectId",
    "sampletype"
  )
  
  hierarchy <- c(
    "sampletype:visitId:studyarm" = "tissue × time × intervention (malleable)",
    "sampletype:visitId"          = "tissue × time (malleable)",
    "studyarm:visitId"            = "intervention × time (malleable)",
    "visitId"                     = "time (malleable)",
    "studyarm:sampletype"         = "tissue × intervention (stable)",
    "studyarm"                    = "intervention (stable)",
    "subjectId"                   = "individual (stable)",
    "sampletype"                  = "tissue (stable)"
  )
  
  # Extract the numeric matrix for faster ops
  mat <- as.matrix(vp_results[, comp_cols])
  thresh <- vp_results[[thresh_col]]
  residuals <- vp_results$Residuals
  
  # Compute max component index and value
  max_idx <- max.col(mat, ties.method = "first")
  max_val <- mat[cbind(seq_len(nrow(mat)), max_idx)]
  max_comp <- comp_cols[max_idx]
  
  # Check if any component > threshold per row (vectorized)
  any_above <- rowMaxs(mat) > thresh
  
  # Determine if dominant
  is_dominant <- max_val >= major_cutoff
  
  # Assign categories
  category <- character(nrow(mat))
  category[!any_above] <- "non-variable / residual"
  category[any_above] <- hierarchy[max_comp[any_above]]
  dominant <- character(nrow(mat))
  dominant[any_above & is_dominant] <- 'dominant'
  dominant[any_above & !is_dominant] <- 'non-dominant'
  
  
  # Flag residuals > 0.5
  residual_flag <- residuals > 0.5
  
  # Combine back into a dataframe
  vp_results$category <- category
  vp_results$dominant <- dominant
  vp_results$dominant_component <- max_comp
  vp_results$dominant_fraction <- max_val
  vp_results$residual_flag <- residual_flag
  
  return(vp_results)
}

# Example:
vp_results <- classify_vp(vp_results, thresh_col = "thresh")

# 4. Cell type drivers ---------------------
cat("4.1 | Load cell type data... ")
# celltype states from output of 5-celltype-association

load("1-analyses/variance-analysis-methylation/3-output/celltype_association.Rdata")
celltype_stats <- celltype_stats[match(rownames(vp_results), rownames(celltype_stats)),]
if(!identical(rownames(celltype_stats), rownames(vp_results))){warning("Names not identical, please doublecheck later.")}

vp_results$celltype_corr <- as.numeric(celltype_stats[, "max_cor"])
vp_results$celltype_pval <- as.numeric(celltype_stats[, "min_p"])
vp_results$celltype_type <- celltype_stats[, "best_celltype"]

vp_results$celltype_effect <- case_when(
    vp_results$celltype_corr > 0.3 & vp_results$celltype_pval < 1e-5 ~ "celltype-driven",
    (vp_results$celltype_corr <= 0.3 | vp_results$celltype_pval >= 1e-5) ~ "non-celltype-driven",
  TRUE ~ NA
)

# 5. Dominant intervention -----------------
load("~/Dropbox/data/tirolgesund/beta_merged.Rdata")

# Match metadata to full beta matrix
cat("1.1 | Match data to pheno... ")
pheno <- pheno[match(colnames(beta_merged), pheno$basename), ]
stopifnot(identical(colnames(beta_merged), pheno$basename))

# Subset by CpGs overlapping with rep set
intersect_cpgs <- intersect(rownames(beta_merged), rownames(vp_results))
cat(length(intersect_cpgs), 'unique CpGs included')
beta_matrix <- beta_merged[intersect_cpgs, ]

beta2M <- function(x){
  x[x == 0] <- min(x[x != 0])
  x[x == 1] <- max(x[x != 1])
  log2(x) - log2(1 - x)
}

mval_matrix <- beta2M(beta_matrix)
pheno$studyarm <- ifelse(pheno$interventionId == 'S', 'S', 'I')

compute_dominant_factors <- function(mval_matrix, pheno, vp_results, 
                                     category_col = "category") {
  
  # Output columns
  vp_results$dominant_intervention <- NA_character_
  vp_results$dominant_tissue <- NA_character_
  
  # Factors
  studyarm <- factor(pheno$studyarm)
  visitId  <- factor(pheno$visitId)
  tissue   <- factor(pheno$sampletype)
  
  # Combined groupings
  group_studyarm_visit <- interaction(studyarm, visitId, drop = TRUE)
  group_tissue_visit   <- interaction(tissue, visitId, drop = TRUE)
  
  # Sparse design matrices (Samples × Groups)
  library(Matrix)
  G_iv <- Matrix::sparse.model.matrix(~0 + group_studyarm_visit)  # Intervention × Time
  G_tv <- Matrix::sparse.model.matrix(~0 + group_tissue_visit)    # Tissue × Time
  
  # --- Helper function to compute dominant factor ---
  dominant_factor <- function(mat, G, factor_labels) {
    # mat: CpGs × Samples
    # G: Samples × Groups
    mat <- as.matrix(mat)
    
    # Means per group: (CpGs × Samples) × (Samples × Groups) = CpGs × Groups
    sums <- mat %*% as.matrix(G)
    counts <- (!is.na(mat)) %*% as.matrix(G)
    means <- sums / counts
    
    # Map group columns to top-level factor
    group_names <- colnames(G)
    factor_ids <- sub("\\..*", "", group_names)  # take part before first dot
    unique_factors <- unique(factor_ids)
    
    # Compute max-min delta for each factor
    deltas <- sapply(unique_factors, function(f) {
      cols <- which(factor_ids == f)
      if (length(cols) == 1) {
        # Single timepoint → delta = 0
        rep(0, nrow(means))
      } else {
        apply(means[, cols, drop = FALSE], 1, function(x) {
          if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
        })
      }
    })
    
    if (is.vector(deltas)) deltas <- matrix(deltas, ncol = 1)
    colnames(deltas) <- unique_factors
    
    # Dominant factor = factor with largest delta
    max_idx <- max.col(deltas, ties.method = "first")
    dominant <- unique_factors[max_idx]
    
    # Replace with NA if all deltas are 0 or NA
    dominant[rowSums(is.na(deltas)) == ncol(deltas) | rowSums(deltas, na.rm = TRUE) == 0] <- NA
    
    return(dominant)
  }
  
  # --- Compute dominant intervention for relevant categories ---
  int_cats <- c("intervention × time (malleable)", 
                "tissue × time × intervention (malleable)")
  idx_int <- which(vp_results[[category_col]] %in% int_cats)
  if (length(idx_int) > 0) {
    vp_results$dominant_intervention[idx_int] <- dominant_factor(
      mval_matrix[idx_int, , drop = FALSE], G_iv, levels(studyarm)
    )
  }
  
  # --- Compute dominant tissue for relevant categories ---
  tissue_cats <- c("tissue × time (malleable)", 
                   "tissue × intervention (stable)",
                   "tissue × time × intervention (malleable)")
  idx_tissue <- which(vp_results[[category_col]] %in% tissue_cats)
  if (length(idx_tissue) > 0) {
    vp_results$dominant_tissue[idx_tissue] <- dominant_factor(
      mval_matrix[idx_tissue, , drop = FALSE], G_tv, levels(tissue)
    )
  }
  
  return(vp_results)
}


vp_results <- compute_dominant_factors(mval_matrix, pheno, vp_results)


# --- 5b. Quick screen: which tissue drives arm×time variability? ---

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixStats))

# features to screen: dominated by the triple interaction
triple_lab <- "tissue × time × intervention (malleable)"
idx_triple <- which(vp_results$category == triple_lab & rownames(vp_results) %in% rownames(mval_matrix))
feat_triple <- rownames(vp_results)[idx_triple]

if (length(feat_triple)) {
  
  # helper: row-wise variance of arm×visit group means inside one tissue
  var_by_groups_one_tissue <- function(Y_sub, ph_sub) {
    # Y_sub: features x samples (numeric matrix for one tissue)
    # Build arm×visit group design within this tissue
    grp <- interaction(ph_sub$studyarm, ph_sub$visitId, drop = TRUE)
    X   <- Matrix::sparse.model.matrix(~ 0 + grp)         # S x G
    # sums and counts per group
    sums   <- as.matrix(Y_sub) %*% as.matrix(X)           # F x G
    counts <- (!is.na(Y_sub)) %*% as.matrix(X)            # F x G
    means  <- sums / pmax(counts, 1)                      # avoid divide-by-zero
    # variance across group means per feature
    matrixStats::rowVars(means, na.rm = TRUE)
  }
  
  tissues <- levels(factor(pheno$sampletype))
  out_list <- vector("list", length(tissues))
  
  for (k in seq_along(tissues)) {
    tt  <- tissues[k]
    jdx <- which(pheno$sampletype == tt)
    if (length(jdx) < 2L) {
      out_list[[k]] <- data.frame(feature = feat_triple,
                                  tissue  = tt,
                                  groupmean_var = NA_real_)
      next
    }
    Y_tt <- mval_matrix[feat_triple, jdx, drop = FALSE]
    ph_tt <- droplevels(pheno[jdx, ])
    
    v_tt <- var_by_groups_one_tissue(Y_tt, ph_tt)
    
    out_list[[k]] <- data.frame(
      feature = feat_triple,
      tissue  = tt,
      groupmean_var = as.numeric(v_tt),
      row.names = NULL
    )
  }
  
  screen_tbl <- dplyr::bind_rows(out_list)
  
  # choose the tissue with the largest variance of arm×visit means per feature
  drivers_screen <- screen_tbl |>
    dplyr::group_by(feature) |>
    dplyr::slice_max(groupmean_var, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  # attach to vp_results
  vp_results$driver_tissue_screen <- NA_character_
  vp_results$driver_stat_screen   <- NA_real_
  
  m <- match(rownames(vp_results), drivers_screen$feature)
  hit <- which(!is.na(m))
  vp_results$driver_tissue_screen[hit] <- drivers_screen$tissue[m[hit]]
  vp_results$driver_stat_screen[hit]   <- drivers_screen$groupmean_var[m[hit]]
}


vp_results$driver_stat_screen
table(vp_results$driver_tissue_screen)

x <- vp_results['cg01566127',]
x <- vp_results['cg05000748',]
x <- vp_results['cg05218470',]

vp_results$celltype_allocation <- NA
vp_results <- vp_results |> 
  dplyr::mutate(celltype_allocation = ifelse(category != 'tissue (stable)', NA,
                                             ifelse(celltype_type=='celltype_Epi' & abs(celltype_corr)>0.95,
                                                    "epithelial composition",
                                                    "other")))


# Append non-variable CpGs
load("1-analyses/variance-analysis-methylation/2-output/variability_overall.Rdata")
load('1-analyses/variance-analysis-methylation/1-output/replication_variability.Rdata') # includes rep SD from 3-replicate-variability.R
nonvar <- data.frame(category = rep("non-variable (stable)",length(non_variable_cpgs)),
                     cg = non_variable_cpgs)

vp_results <- vp_results |> 
  tibble::rownames_to_column('cg') |> 
  dplyr::full_join(nonvar) |> 
  dplyr::mutate(category = ifelse(category == 'non-variable / residual', 'residual', category)) |> 
  dplyr::relocate(cg)

overall_sd <- overall_sd[vp_results$cg]
tech_sd <- tech_sd[vp_results$cg]
identical(names(overall_sd), names(tech_sd))
tech_var_fraction <- tech_sd^2 / overall_sd
identical(names(tech_var_fraction), vp_results$cg)

vp_results$thresh <- tech_var_fraction
vp_results$overall_sd <- overall_sd
vp_results$tech_sd <- tech_sd
# cor(vp2$thresh, vp2$thresh2,use = 'complete.obs')

save(vp_results, file = '1-analyses/variance-analysis-methylation/5-output/vp_results.Rdata')
