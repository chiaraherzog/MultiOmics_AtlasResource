# Running MEFISTO on multiomic data.
# Author: Chiara Herzog
# Date: Jul 9, 2025

# 0.  Setup -------------------------------------------------

# Packages setup
packages <- c('BiocManager', 'MultiAssayExperiment', 'MOFA2', 'dplyr', 'here')


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


here::i_am('1-analyses/mofa/mefisto.R')
reticulate::use_python("/Users/chiara/opt/anaconda3/bin/python")


# 1. Data preparation -------------------------------------------------


## 1. Loading data
load(here('data/data_normalized_centred.Rdata'))


## 2. Filter relevant variables to be included
load(here('data/vars_if.Rdata'))
vars_if <- vars
load(here('data/vars_smk.Rdata'))
vars_smk <- vars
vars <- vars_if |> dplyr::full_join(vars_smk) |> 
  dplyr::select(assay, x) |> 
  dplyr::filter(!assay %in% c("Cotinine LC/MS Quantitation",
                              "Immune age: general",
                              "ImmuneSMK",
                              "Saliva microbiome: ASVs",
                              "Saliva microbiome: families",
                              "Stool microbiome: ASVs",
                              "Stool microbiome: families",
                              "Saliva microbiome: ASVs_clr",
                              "Stool microbiome: ASVs_clr")) |> 
  dplyr::filter(x != 'rhr') |> 
  dplyr::filter(!x %in% c('ic') & !grepl("epidish_", x)) # keep out = using for later phenotypic associations


# Variables for later correlation:
vars_assoc <- vars |> dplyr::filter(
  assay %in% c("Blood haemogram",
               "Body composition",
               "Skin histology and transepidermal water loss assay",
               "Vascular and body sonography", 
               "Functional sports exam")
)


vars_keep <- vars |> dplyr::filter(
  !assay %in% vars_assoc$assay
)


filter_vars <- vars_keep |> 
  dplyr::group_by(assay) |> 
  dplyr::summarise(features = list(x), .groups = "drop") |>
  tibble::deframe()


## 3. Subset & filter MAE
# 'mae' is the MultiAssayExperiment
#   rows  = features
#   cols  = samples
#   colData(mae) must contain: subjectId, timepoint (visitId), intervention (interventionId)
mae <- subsetByAssay(data, filter_vars)
mae <- mae[ , colData(mae)$visitId %in% c("M0", "M2", "M4", "M6") ]
mae$intervention <- ifelse(mae$interventionId=='S', 'S', 'I') # do not distinguish between I/K




# #### overview of data - centred
# t <- mae[[1]]
# hist(t[1,])
# hist(t[2,])
# hist(t[3,])
# hist(t[10,])
# 
# t <- mae[[4]]
# head(t)
# hist(t[1,])
# hist(t[2,])
# hist(t[3,])
# hist(t[10,])


## 4. Append PCs (bringing them in to right format first - SummarizedExperiment)
scaleDat <- function(mat){
  mat <- as.matrix(mat)
  tmp_scale <- apply(mat, 2, function(i) scale(i,center = T))
  rownames(tmp_scale) <- rownames(mat)
  mat <- tmp_scale
}


PCs_blood <- readRDS(here("0-preprocessing/2-output/PCs-long-blood.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary') |> 
  scaleDat()


PCs_blood <- SummarizedExperiment(
  assays  = list(counts = t(PCs_blood)), # feature x samples
  colData = colData(mae)[rownames(PCs_blood), ] # carry over existing metadata
)


PCs_buccal <- readRDS(here("0-preprocessing/2-output/PCs-long-buccal.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary') |> 
  scaleDat()


PCs_buccal <- SummarizedExperiment(
  assays  = list(counts = t(PCs_buccal)), # feature x samples
  colData = colData(mae)[rownames(PCs_buccal), ] # carry over existing metadata
)


PCs_cervical <- readRDS(here("0-preprocessing/2-output/PCs-long-cervical.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary') |> 
  scaleDat()
PCs_cervical <- SummarizedExperiment(
  assays  = list(counts = t(PCs_cervical)), # feature x samples
  colData = colData(mae)[rownames(PCs_cervical), ] # carry over existing metadata
)


mae <- c(mae, 
         "Methylation PCs: blood" = PCs_blood,
         "Methylation PCs: buccal" = PCs_buccal,
         "Methylation PCs: cervical" = PCs_cervical)


# # 2. MOFA Pipeline (ALL) -------------------------------------------------
# 
# 
# ## 1.  Create the MOFA object
# mofa <- create_mofa_from_MultiAssayExperiment(
#   mae
# )
# 
# 
# ### We want to use actual times, not labels -> extract t
# extract_t <- mae@metadata$`timing of collection` |> 
#   as.data.frame() |> 
#   dplyr::filter(assay == 'Blood haemogram') |> # most complete set
#   dplyr::select(primary, t) |> 
#   dplyr::distinct()
# 
# 
# ### Fix labels
# x <- samples_metadata(mofa) |> 
#   dplyr::left_join(extract_t, by = c("sample" = "primary"))
# 
# 
# ### Append as metadata
# samples_metadata(mofa) <- x
# mofa <- set_covariates(mofa, covariates = c("t")) # day as covariate
# mofa
# get_covariates(mofa, 1)
# rm(x)
# 
# ## Preview
# gg_input <- plot_data_overview(mofa,
#                                show_covariate = TRUE,
#                                show_dimensions = TRUE,
#                                covariate = "t") 
# gg_input
# 
# 
# # 3. Options: we will test a few - as we did for MOFA at baseline :) 
# ###  Modalities -> all continuous, hence gaussian approximation is good for all.
# ### Scaling -> data and views already scaled, so not needed.
# data_opts   = get_default_data_options(mofa)
# 
# mefisto_opts <- get_default_mefisto_options(mofa)
# mefisto_opts$model_groups <- FALSE
# mefisto_opts$warping      <- FALSE   # enable if visits are mis-aligned
# 
# train_opts <- get_default_training_options(mofa)
# train_opts$maxiter <- 1000
# 
# 
# model_opts <- get_default_model_options(mofa)
# model_opts$num_factors # default = 15
# 
# 
# model_opts_n10 <- model_opts
# model_opts_n10$num_factors <- 10
# 
# 
# model_opts_n20 <- model_opts
# model_opts_n20$num_factors <- 20
# # Initial training -> fast convergence.
# 
# 
# model_opts_n25 <- model_opts
# model_opts_n25$num_factors <- 25
# # Initial training -> fast convergence.
# 
# 
# ## 3.  Attach options → prepare_mofa()
# mofa_n15 <- prepare_mofa(mofa,
#                          data_options   = data_opts,
#                          model_options  = model_opts,
#                          training_options = train_opts,
#                          mefisto_options  = mefisto_opts)
# 
# 
# mofa_n20 <- prepare_mofa(mofa,
#                          data_options   = data_opts,
#                          model_options  = model_opts_n20,
#                          training_options = train_opts,
#                          mefisto_options  = mefisto_opts)
# 
# 
# mofa_n10 <- prepare_mofa(mofa,
#                          data_options   = data_opts,
#                          model_options  = model_opts_n10,
#                          training_options = train_opts,
#                          mefisto_options  = mefisto_opts)
# 
# 
# mofa_n25 <- prepare_mofa(mofa,
#                          data_options   = data_opts,
#                          model_options  = model_opts_n25,
#                          training_options = train_opts,
#                          mefisto_options  = mefisto_opts)
# 
# 
# ## 4.  Fit the models → run_mofa()
# fast_n10 <- run_mofa(
#   mofa_n10,
#   outfile = here('1-analyses/mofa/out', "long_fast_n10.hdf5")
# )
# 
# 
# fast_n15 <- run_mofa(
#   mofa_n15,
#   outfile = here('1-analyses/mofa/out', "long_fast_n15.hdf5")
# )
# 
# 
# fast_n20 <- run_mofa(
#   mofa_n20,
#   outfile = here('1-analyses/mofa/out', "long_fast_n20.hdf5")
# )
# 
# 
# fast_n25 <- run_mofa(
#   mofa_n25,
#   outfile = here('1-analyses/mofa/out', "long_fast_n25.hdf5")
# )
# 
# 
# # Diagnostics: let's compare ELBO values.
# model.list <- list(fast_n10, fast_n15, fast_n20, fast_n25)
# elbo <- compare_elbo(model.list)
# elbo$data ## Pretty similar, but model 4 wins (n25)
# 
# 
# plot_variance_explained(fast_n10, x = 'factor',
#                         y = 'view')
# plot_variance_explained(fast_n15, x = 'factor',
#                         y = 'view')
# plot_variance_explained(fast_n20, x = 'factor',
#                         y = 'view')
# plot_variance_explained(fast_n25, x = 'factor',
#                         y = 'view')
# 
# # Let's retrain final model using 25 factors
# train_opts$convergence_mode <- 'medium'
# train_opts$maxiter <- 2000
# train_opts
# 
# mofa_n25_medium <- prepare_mofa(mofa,
#                                 data_options   = data_opts,
#                                 model_options  = model_opts_n25,
#                                 training_options = train_opts,
#                                 mefisto_options  = mefisto_opts)
# 
# 
# mefisto_trained_overall <- run_mofa(
#   mofa_n25_medium,
#   outfile = here('1-analyses/mofa/out', "long_medium.hdf5")
# ) # Final model
# 
# 
# ## 5. Quick checks - Variance explained per view / factor
# plot_variance_explained(mefisto_trained_overall, x = 'factor',
#                         y = 'view')
# plot_factor_cor(mefisto_trained_overall)
# 
# 
# plot_factors_vs_cov(
#   mefisto_trained_overall,
#   covariate = "t"
# )
# 

# 3. Associations - models -------------------------
library(tidyverse)
library(broom)
library(MASS)
outfile <- here('1-analyses/mofa/out', "long_medium.hdf5")
mefisto_trained <- load_model(outfile)

## 1.  Z-score the factors so slopes are comparable
factors <- as.data.frame(get_factors(mefisto_trained, scale = T)[[1]])

## Get features
## 1. Phenotypical variables
feat1 <- colData(mae) |>
  as.data.frame() |> 
  dplyr::mutate(menopause = gsub("[*]", "", mpstatrs),
                studyarm = ifelse(interventionId == 'S', "S", "I"),
                # pseudocode ethanol as numeric
                etohu_curr_pseudonum = as.numeric(as.factor(etohu_curr))) |> 
  dplyr::select(subjectId, studyarm, interventionId, visitId, age_at_consent, time, menopause, bmi_at_consent,
                smoking_ever, smoking_curr, etoh_curr, etohu_curr_pseudonum, 
                diet, preg_ever, ocp_curr,
                compliant, compliance_smkgroup,
                dbmi,
                endomet:risk_autimm,
                comprate_longitudinal,
                dailycig_longitudinal,
                smkstop) |> 
  dplyr::mutate(group_time = paste0(time, studyarm),
                comprate_longitudinal_IFonly = ifelse(interventionId == 'S', NA, comprate_longitudinal)) |> 
  dplyr::mutate(across(everything(), ~ ifelse(. == 'unknown', NA, .))) |>
  dplyr::select(where(~ n_distinct(.,na.rm = T) > 1)) |>
  dplyr::mutate(visitId = factor(visitId, levels = c("M0", "M2", "M4", "M6"), ordered = T),
                interventionId = factor(interventionId, levels= c("I", "K", "S")),
                compliance_smkgroup = factor(compliance_smkgroup, levels = c("non-compliant",
                                                                             "heated tobacco product",
                                                                             "reducer",
                                                                             "complete cessation")),
                studyarm = factor(studyarm),
                group_time = factor(group_time)) |> 
  dplyr::mutate(across(c("menopause", "ocp_curr", "etoh_curr", "endomet", "pcos", "preg_ever", "compliant", "smoking_curr", "smoking_ever", "smkstop"), ~ factor(., levels = c("no", "yes")))) |> 
  dplyr::mutate(across(contains("risk_"), ~ as.numeric(.))) |> 
  dplyr::mutate(diet = factor(diet, levels = c("normal", "pescetarian", "vegetarian", "vegan")))


## 2. Clinical and other features
feat2 <- as.data.frame(wideFormat(data[,data$visitId %in% c('M0', "M2", "M4", "M6"),c(1, 3, 4, 5, 9, 13:15)])) |> 
  tibble::column_to_rownames('primary')
colnames(feat2)

vars <- vars_assoc |> 
  dplyr::add_row(assay = 'Composite methylation scores..buccal', x = 'ic') |> 
  dplyr::add_row(assay = 'Composite methylation scores..cervical', x = 'ic') |> 
  dplyr::add_row(assay = 'Composite methylation scores..blood', x = 'hepidish_Neutro') |> 
  dplyr::mutate(varname = paste0(gsub(" ", ".", assay), "_", x)) |> 
  dplyr::filter(varname %in% colnames(feat2))
feat2 <- feat2[,vars$varname]


## 3. Merge features
feat1 <- feat1[rownames(factors),]
identical(rownames(factors), rownames(feat1))
feat2 <- feat2[rownames(factors),]
identical(rownames(factors), rownames(feat2))
features <- cbind(feat1, feat2)
save(features, file = "1-analyses/mofa/out/features_long.Rdata")


## 2.  All (factor, feature) pairs
pairs_tbl <- tidyr::crossing(
  factor_col  = names(factors),
  feature_col = colnames(features)
) |> 
  dplyr::filter(feature_col != 'subjectId')


### Opt a - spearman/kruskall
kw_eps2 <- function(x, g) {
  ok <- complete.cases(x, g)           # drop NAs the same way the test does
  g  <- as.factor(g[ok])
  x  <- x[ok]
  k  <- nlevels(g)
  n  <- length(x)
  H  <- kruskal.test(x, g)$statistic   # “Kruskal-Wallis chi-squared”
  (H - k + 1) / (n - k)
}


pairs_tbl_spear_krusk <- pairs_tbl |> 
  dplyr::mutate(
    res = purrr::map2(
      factor_col, feature_col,
      ~{
        f_vec <- factors[[.x]]
        x_vec <- features  [[.y]]
        
        if (is.character(x_vec) || is.factor(x_vec)) {
          ## -------- categorical feature  ➜  Kruskal–Wallis
          p  <- kruskal.test(f_vec, as.factor(x_vec))$p.value
          es <- kw_eps2(f_vec, x_vec)              # ε² effect size
          list(p_value = p, effect_size = es, effect_label = "eps2")
        } else {
          ## -------- numeric feature      ➜  Spearman correlation
          ct <- suppressWarnings(cor.test(f_vec, x_vec, method = "spearman"))
          list(p_value     = ct$p.value,
               effect_size = unname(ct$estimate),  # ρ
               effect_label = "rho")
        }
      })
  ) |> 
  tidyr::unnest_wider(res)


pairs_tbl_spear_krusk <- pairs_tbl_spear_krusk |> 
  mutate(q_value = p.adjust(p_value, method = "fdr")) |> 
  arrange(p_value) |> as.data.frame()
# x <- pairs_tbl_spear_krusk[pairs_tbl_spear_krusk$q_value<0.05,]
saveRDS(pairs_tbl_spear_krusk, file = here('1-analyses/mofa/out/pairs_tbl_long_overall_spear_krusk.Rds'))






### Opt b - linear-log models


####Helper that chooses the right model, returns β & p (no correction for subjectId but could include)
assoc_one_nosubject <- function(f_vec, feat_name, features_df) {
  feat_vec <- features_df[[feat_name]]
  subj_vec <- features_df$subjectId
  
  ok <- complete.cases(f_vec, feat_vec, subj_vec)
  df <- data.frame(y = feat_vec[ok],
                   x = f_vec[ok],
                   subjectId = subj_vec[ok])
  
  if (is.numeric(df$y)) {
    fit <- lm(y ~ x, data = df)
    out <- tidy(fit)[2,]
    c(beta = out$estimate, p = out$p.value, model = "linear")
    
  } else if (is.factor(df$y) && !is.ordered(df$y)) {
    fit <- glm(y ~ x, family = binomial, data = df)
    out <- tidy(fit)[2,]
    c(beta = out$estimate, p = out$p.value, model = "logistic")
    
  } else if (is.ordered(df$y)) {
    fit <- polr(y ~ x, data = df, Hess = TRUE)
    co  <- summary(fit)$coefficients[1,]
    z   <- co["t value"]
    p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
    c(beta = unname(co["Value"]), p = p, model = "ordinal")
    
  } else {
    c(beta = NA, p = NA, model = "skipped")
  }
}


pairs_tbl_models <- pairs_tbl |> 
  mutate(res = pmap(list(factor_col, feature_col),
                    ~assoc_one_nosubject(factors[[..1]], ..2, features))) |> 
  unnest_wider(res) |> 
  dplyr::mutate(p = ifelse(is.na(`p.t value`), p, `p.t value`)) |>
  dplyr::mutate(q_value = p.adjust(p, "fdr"))
saveRDS(pairs_tbl_models, file = here('1-analyses/mofa/out/pairs_tbl_long_overall_models.Rds'))


# SAVE output for easy plotting ----------------------
outfile <- here('1-analyses/mofa/out', "long_medium.hdf5")
mofa <- load_model(outfile)
fact <- as.data.frame(get_factors(mofa, scale = T)[[1]]) |> tibble::rownames_to_column('primary')

load(file = "1-analyses/mofa/out/features_long.Rdata") 
features <- features |> tibble::rownames_to_column('primary')

pheno_factors <- features  |> 
  dplyr::left_join(fact, by = "primary")
saveRDS(pheno_factors, file = "1-analyses/mofa/out/pheno_w_mofa_long.Rdata")



# Exploring factors more -------------------
# Exploring F13 more ------------------------

load("data/data_normalized_centred.Rdata")
library(MultiAssayExperiment)
sal <- as.data.frame(longForm(subsetByAssay(data, 'Saliva nuclear magnetic resonance: normalized')))

pheno_w_mofa_F12 <- pheno_w_mofa_long |> 
  dplyr::left_join(sal, by = 'primary')

top_w_f12 <- get_weights(mofa_long, factors = 12, scale = T,
                         as.data.frame = T) |> 
  dplyr::mutate(absval = abs(value)) |>
  dplyr::arrange(desc(absval))

# Any features >0.35 & annotate labels
load("data/vars.Rdata")

enrich_f12 <- top_w_f12 |> dplyr::filter(absval >=0.35) |> 
  dplyr::mutate(feat = gsub("_Saliva nuclear magnetic resonance: normalized", "", feature)) |> 
  dplyr::left_join(vars, by = c("feat" = 'x')) |> 
  dplyr::pull(label) |> unique()

enrich_f20 <- get_weights(mofa_long, factors = 20, scale = T,
                          as.data.frame = T) |> 
  dplyr::mutate(absval = abs(value)) |>
  dplyr::arrange(desc(absval)) |>
  dplyr::filter(absval >=0.35) |> 
  dplyr::mutate(feat = gsub("_Saliva nuclear magnetic resonance: normalized", "", feature)) |> 
  dplyr::left_join(vars, by = c("feat" = 'x')) |> 
  dplyr::pull(label) |> unique()

enrichMetab <- function(metabolites, filename){
  library(MetaboAnalystR)
  mSet<-InitDataObjects("list", "msetora", FALSE, 150)  # Create mSetObj
  mSet<-Setup.MapData(mSet, metabolites); # Set up mSetObj with the list of compounds
  mSet<-CrossReferencing(mSet, "name");
  mSet<-CreateMappingResultTable(mSet)
  mSet<-SetMetabolomeFilter(mSet, F);
  mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
  mSet<-CalculateHyperScore(mSet)
  
  write.table(mSet$analSet$ora.mat, file = paste0('1-analyses/mofa/out/', filename, '.csv'), sep = ",")
  return(mSet)
}

f12 <- enrichMetab(enrich_f12, filename = 'mefisto-f12')
f20 <- enrichMetab(enrich_f20, filename = 'mefisto-f20')
