# Variance partition on data
# Author: Chiara Herzog
# Date: 20 Dec 2024

# Packages
packages <- c('fs', 'variancePartition', 'here', 'dplyr', 'purrr', 'variancePartition',
              'MultiAssayExperiment')

for(p in packages){
  if(!require(p,
              character.only = T)){
    try(install.packages(p),silent = T)
    if(!require(p,
                character.only = T)){
      BiocManager::install(p)
    } 
  }
}

here::i_am('1-analyses/variance-analysis/2-variance-partition.R')


# Load data
load(here('data/data_normalized.Rdata'))

# Keep only M0-M6 data
data <- data[,data$visitId %in% c("M0", "M2", "M4", "M6"),]

# Load variables to use (except methylation)
load(here('data/vars.Rdata')) # Variables were defined previously

vars <- vars |> # Keep ony those assays with >=5 features
  dplyr::group_by(assay) |> 
  dplyr::filter(dplyr::n()>=5) |> 
  dplyr::ungroup()
table(vars$assay)
gc()

# Variance partition by experiment
experiments <- unique(vars$assay)
experiments <- experiments[experiments != 'Body composition']
data$interventionId_bin <- ifelse(data$interventionId=='S', 'S', 'I')

vp <- vector("list", 
             length = length(experiments))

for(e in experiments){
  cat(e, "variance partition starting\n")
  
  if(e != 'Vascular and body sonography'){
    v <- vars[vars$assay==e,]$x
  suppressWarnings(suppressMessages(mat <- data[v,,e][[1]]))
  suppressWarnings(suppressMessages(p <- as.data.frame(colData(data[v,,e]))))
  } else {
    suppressWarnings(suppressMessages(mat <- data[,,e][[1]]))
    suppressWarnings(suppressMessages(p <- as.data.frame(colData(data[,,e]))))
  }
  
  x <- substr(colnames(mat), 1, 6)
  if(identical(rownames(p), x)){colnames(mat) <- rownames(p)} else {stop("check data input")}
  
  ## Formula
  form <- ~ (1 | subjectId) + (1 | interventionId_bin) # binary - smoking or IF
  suppressWarnings(suppressMessages(varPart <- fitExtractVarPartModel(mat, form, p)))
  
  n <- janitor::make_clean_names(e)
  
  vp[[which(experiments==e)]] <- sortCols(varPart)
  names(vp)[which(experiments==e)] <- n
  
  gc()

}

save(vp, file = here('1-analyses/variance-analysis/2-out/vp_w_o_raw_methylation.Rdata'))

# vp_df <- vp |> 
#   purrr::map_dfr(
#     .f = as.data.frame, # Convert each element to a data frame
#     .id = "group"       # Add list names as a column
#   )
# 
# 
# 
# table(vp_df$group)
# 
# filtered <- vp_df |> dplyr::filter(
#   ! group %in% c("saliva_microbiome_as_vs_clr", "saliva_microbiome_families_clr",
#                  "stool_microbiome_as_vs_clr", "stool_microbiome_families_clr")
#   )
# 
# 
# # Barplot
# filtered |> 
#   tibble::rownames_to_column('var') |> 
#   dplyr::arrange(dplyr::desc(subjectId)) |> 
#   tidyr::pivot_longer(subjectId:Residuals,
#                       names_to = 'type',
#                       values_to = 'value') |> 
#   ggplot(aes(x = value,
#              # y = forcats::fct_reorder(var, value),
#              y = var,
#              fill = type)) +
#   geom_col() +
#   theme(axis.text.y = element_blank()) +
#   facet_wrap(~group,scales = 'free')
# 
# 
# 
# # Violin plot
# filtered |>  
#   tidyr::pivot_longer(subjectId:Residuals,
#                       names_to = 'type',
#                       values_to = 'value') |> 
#   ggplot(aes(type,
#              y = value)) +
#   geom_violin() +
#   facet_wrap(~group)
# 
# 
# # ternary plot
# filtered |> 
#   dplyr::filter(group == 'blood_haemogram') |> 
#   ggtern(aes(x=Residuals,
#              y=subjectId,
#              z=interventionId)) +
#   # geom_point() +
#   geom_hex_tern() +
#   labs(title="") +
#   scale_fill_viridis_c() 