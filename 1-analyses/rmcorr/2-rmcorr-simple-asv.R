# Repeated measures correlation across omic features - simple

library(dplyr)
library(here)
library(rmcorr)
library(tidyr)
library(ggplot2)
library(purrr)
library(MultiAssayExperiment)
here::i_am("1-analyses/rmcorr/2-rmcorr-simple-asv.R")

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
#   # dplyr::filter(!grepl("ASV", assay)) |>
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
# save(dat, file = here("1-analyses/rmcorr/1-out/rmcorr-data-asv.Rdata"))

# Run RMCORR ----------------------
load(here("1-analyses/rmcorr/1-out/rmcorr-data-asv.Rdata"))

dat$subjectId <- as.factor(dat$subjectId)
measures1 <- grep("Composite", names(dat), value = TRUE)
measures2 <- setdiff(names(dat), c(measures1, "visitId", "subjectId", "primary"))

pairs <- tidyr::crossing(measure1 = measures1, measure2 = measures2)

run_rmc <- function(m1, m2){
  
  tmp <- rmcorr(participant = subjectId,
                measure1    = m1,
                measure2    = m2, 
                dataset = dat)
  tibble(
    measure1 = m1, measure2 = m2,
    r = tmp$r, p = tmp$p, df = tmp$df,
    ci_lo = tmp$CI[1], ci_hi = tmp$CI[2], ci_lvl = tmp$CI.level
  )
}

# progress = nice console bar/spinner; no pb$tick() needed
res <- purrr::pmap_dfr(
  pairs,
  ~ run_rmc(..1, ..2),
  .progress = "Running repeated-measures correlations"
)


res <- res |>
  dplyr::group_by(measure1) |> 
  dplyr::mutate(padj = p.adjust(p, method = 'fdr')) |> 
  dplyr::ungroup()

save(res, file = here("1-analyses/rmcorr/2-out/rmcorr-simple-asv.Rdata"))

x <- res[res$padj<0.05,]

source("source/rmcorr_plot.R")
