# Repeated measures correlation across omic features

library(dplyr)
library(here)
library(rmcorr)
library(ggplot2)
library(MultiAssayExperiment)
here::i_am("1-analyses/rmcorr/1-rmcorr.R")

# # Construct data ------------------
# 
# # load in data
# load("data/data_normalized.Rdata")
# 
# # Get variables & experiments
# load("data/vars.Rdata")
# 
# # Format variables
# variables <- vars |>
#   dplyr::filter(!assay %in% c("Cotinine LC/MS Quantitation",
#                               "Immune age: general",
#                               "ImmuneSMK",
#                               "Saliva microbiome: ASVs",
#                               "Saliva microbiome: families",
#                               "Stool microbiome: ASVs",
#                               "Stool microbiome: families")) |>
#   dplyr::filter(x != 'rhr') |>
#   dplyr::filter(!grepl("ASV", assay)) |> 
#   dplyr::group_by(assay) |>
#   dplyr::summarise(features = list(x), .groups = "drop") |>
#   tibble::deframe()
# 
# 
# # Filter data
# df_filtered <- subsetByAssay(data, variables)
# 
# # Get data in long format
# df <- as.data.frame(longForm(df_filtered,
#                                colData = c('subjectId', 'interventionId', 'visitId'))) |>
#   # M0-M6 only
#   dplyr::filter(! visitId %in% c("M12", "M18"))
# 
# # Pivot data to wider, concatenating assay and rowname for unique variables
# dat <- df |>
#   dplyr::mutate(name = paste0(assay, '_', rowname)) |>
#   dplyr::select(-rowname) |>
#   tidyr::pivot_wider(names_from = 'name',
#                      values_from = 'value',
#                      id_cols = c('primary', 'subjectId', 'visitId'))
# 
# dat <- dat |>
#   dplyr::select(where(~!all(is.na(.x))))
# save(dat, file = here("1-analyses/rmcorr/1-out/rmcorr-data.Rdata"))
  
# Run RMCORR ----------------------
load(here("1-analyses/rmcorr/1-out/rmcorr-data.Rdata"))
variables <- colnames(dat)[!colnames(dat) %in% c("visitId", "subjectId", "primary")]

out <- rmcorr::rmcorr_mat(subjectId,
                          variables,
                          dat)

rmcorr <- out$summary |>
  dplyr::mutate(assay1 = stringr::str_split_i(measure1, "_", 1),
                assay2 = stringr::str_split_i(measure2, "_", 1),
                diff_assay = case_when(grepl("cytometry", assay1) & grepl("cytometry", assay2) ~ "no",
                                       assay1 != assay2 ~ 'yes',
                                       TRUE ~ 'no'))

save(rmcorr, file = here("1-analyses/rmcorr/1-out/rmcorr-results.Rdata"))


