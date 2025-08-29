#' plotTopFactorWeights
#' 
#' @name extractTopFactorweights
#' @param mofa mofa object
#' @param factorFilter  which factor to filter - numeric
#' @param top_n how many top features to include
#' @returns a table object
#' 
#' 
extractTopFactorweights <- function(mofa, factorFilter,
                                 top_n = 30){
  
  # Extract weights
  w <- get_weights(mofa,  as.data.frame = T, scale = T) |>
    dplyr::mutate(absval = abs(value)) |> 
    dplyr::filter(factor %in% paste0("Factor", factorFilter))
  
  # Filter & fix labels
  load("data/vars.Rdata")
  
  top_w <- w |> 
    dplyr::arrange(absval) |> 
    dplyr::slice_max(absval, n = top_n) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(feature = gsub("_Composite methylation scores: buccal|_Composite methylation scores: cervical|_Composite methylation scores: blood|_Saliva microbiome: families_clr|_Urine nuclear magnetic resonance: normalized|_Saliva nuclear magnetic resonance: normalized|_Methylation PCs: cervical|_Methylation PCs: blood|_Methylation PCs: buccal|_Stool microbiome: families_clr", "", feature)) |> 
    dplyr::left_join(vars, by = c("feature" = "x",
                                  "view" = "assay")) |> 
    dplyr::mutate(aname2 = case_when(grepl("methylation", view) & grepl("buccal", view) ~ "Buccal methylation scores",
                                     grepl("methylation", view) & grepl("blood", view) ~ "Blood methylation scores",
                                     grepl("methylation", view) & grepl("cervical", view) ~ "Cervical methylation scores",
                                     grepl("Saliva", view, ignore.case = T) & grepl("nuclear", view) ~ "Saliva metabolome",
                                     grepl("Saliva", view, ignore.case = T) & !grepl("metabolome", view) ~ "Saliva microbiome",
                                     grepl("Urine", view, ignore.case = T) & grepl("nuclear", view) ~ "Urine metabolome",
                                     grepl("Stool", view, ignore.case = T) ~ "Faecal microbiome",
                                     TRUE ~ view),
                  label = ifelse(is.na(label), feature, label)) |> 
    dplyr::mutate(feature = label, 
                  view = aname2) |>
    dplyr::rename(absolute_weight = absval,
                  weight = value) |> 
    dplyr::select(factor, view, feature, weight, absolute_weight)
    
  return(top_w)
}
