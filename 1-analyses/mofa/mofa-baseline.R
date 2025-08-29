# Running baseline MOFA on multiomic data.
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


here::i_am('1-analyses/mofa/mofa-baseline.R')


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
mae_base <- mae[ , colData(mae)$visitId == "M0" ]
mae_base$intervention <- ifelse(mae_base$interventionId=='S', 'S', 'I') # do not distinguish between I/K


## 4. Append PCs (bringing them in to right format first - SummarizedExperiment)
scaleDat <- function(mat){
  mat <- as.matrix(mat)
  tmp_scale <- apply(mat, 2, function(i) scale(i,center = T))
  rownames(tmp_scale) <- rownames(mat)
  mat <- tmp_scale
}


PCs_base_blood <- readRDS(here("0-preprocessing/2-output/PCs-base-blood.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary') |> 
  scaleDat()
PCs_base_blood <- SummarizedExperiment(
  assays  = list(counts = t(PCs_base_blood)), # feature x samples
  colData = colData(mae_base)[rownames(PCs_base_blood), ] # carry over existing metadata
)


PCs_base_buccal <- readRDS(here("0-preprocessing/2-output/PCs-base-buccal.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary')|> 
  scaleDat()
PCs_base_buccal <- SummarizedExperiment(
  assays  = list(counts = t(PCs_base_buccal)), # feature x samples
  colData = colData(mae_base)[rownames(PCs_base_buccal), ] # carry over existing metadata
)


PCs_base_cervical <- readRDS(here("0-preprocessing/2-output/PCs-base-cervical.Rds")) |> 
  dplyr::select(-c(basename)) |> 
  tibble::column_to_rownames('primary') |> 
  scaleDat()
PCs_base_cervical <- SummarizedExperiment(
  assays  = list(counts = t(PCs_base_cervical)), # feature x samples
  colData = colData(mae_base)[rownames(PCs_base_cervical), ] # carry over existing metadata
)


mae_base <- c(mae_base, 
              "Methylation PCs: blood" = PCs_base_blood,
              "Methylation PCs: buccal" = PCs_base_buccal,
              "Methylation PCs: cervical" = PCs_base_cervical)


# 3. MOFA Pipeline -------------------------------------------------


## 1.  Create the MOFA object
mofa_base <- create_mofa_from_MultiAssayExperiment(
  mae_base
)


## 1. Set up different options.


### 1.1. Modalities -> all continuous, hence gaussian approximation is good for all.


### 1.2. Scaling -> data and views already scaled, so not needed.
data_opts <- get_default_data_options(mofa_base)
# data_opts$scale_groups
# data_opts$scale_views ## both false already.


### 1.3  Get default option lists
model_opts  <- get_default_model_options(mofa_base)
model_opts$num_factors # default = 15


model_opts_n10 <- model_opts
model_opts_n10$num_factors <- 10


model_opts_n20 <- model_opts
model_opts_n20$num_factors <- 20


train_opts  <- get_default_training_options(mofa_base)
train_opts$maxiter <- 1000


# Initial training -> fast convergence.


## 3.  Attach options → prepare_mofa()
reticulate::use_python("/Users/chiara/opt/anaconda3/bin/python")
mofa_base_n15 <- prepare_mofa(
  object            = mofa_base,
  data_options      = data_opts,
  model_options     = model_opts,
  training_options  = train_opts
)


mofa_base_n10 <- prepare_mofa(
  object            = mofa_base,
  data_options      = data_opts,
  model_options     = model_opts_n10,
  training_options  = train_opts
)


mofa_base_n20 <- prepare_mofa(
  object            = mofa_base,
  data_options      = data_opts,
  model_options     = model_opts_n20,
  training_options  = train_opts
)




## 4.  Fit the models → run_mofa()
fast_n10 <- run_mofa(
  mofa_base_n10,
  outfile = here('1-analyses/mofa/out', "baseline_fast_n10.hdf5")
)


fast_n15 <- run_mofa(
  mofa_base_n15,
  outfile = here('1-analyses/mofa/out', "baseline_fast_n15.hdf5")
)


fast_n20 <- run_mofa(
  mofa_base_n20,
  outfile = here('1-analyses/mofa/out', "baseline_fast_n20.hdf5")
)


# Diagnostics: let's compare ELBO values.
outfile_n10 <- here('1-analyses/mofa/out', "baseline_fast_n10.hdf5")
n10 <- load_model(outfile_n10)

outfile_n20 <- here('1-analyses/mofa/out', "baseline_fast_n20.hdf5")
n20 <- load_model(outfile_n20)

outfile_n15 <- here('1-analyses/mofa/out', "baseline_fast_n15.hdf5")
n15 <- load_model(outfile_n15)

model.list <- list(n10, n15, n20)
elbo <- compare_elbo(model.list)
elbo$data ## Pretty similar, but n10 had the lowest ELBO.

# Let's retrain final model using 10 features.
data_opts <- get_default_data_options(mofa_base)
model_opts <- get_default_model_options(mofa_base)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(mofa_base)
train_opts$convergence_mode <- 'medium'
train_opts$maxiter <- 2000


mofa_baseline <- prepare_mofa(
  object            = mofa_base,
  data_options      = data_opts,
  model_options     = model_opts,
  training_options  = train_opts
)


mofa_baseline_trained <- run_mofa(
  mofa_baseline,
  outfile = here('1-analyses/mofa/out', "baseline_mofa_medium.hdf5")
) # Final model


# # 4. Exploration -------------------------------------------------
# outfile <- here('1-analyses/mofa/out', "baseline_mofa_medium.hdf5")
# mofa_baseline_trained <- load_model(outfile)
# 
# 
# ## Append metadata
# x <- samples_metadata(mofa_baseline_trained)
# y <- colData(mae_base) |> 
#   as.data.frame() |> dplyr::mutate(primary = paste0(subjectId, visitId))
# x <- x |> dplyr::left_join(y,
#                            by = c("sample" = "primary"))
# samples_metadata(mofa_baseline_trained) <- x
# 
# 
# p <- plot_variance_explained(mofa_baseline_trained, plot_total = T)
# p[[1]]
# p[[2]]
# 
# 
# p[[2]]$data
# 
# 
# 
# 
# # Plot factors
# library(ggplot2)
# 
# 
# # AGE
# p <- plot_factor(mofa_baseline_trained, factors = 1:mofa_baseline_trained@dimensions$K,  color_by = "age_at_consent")
# p$layers[[1]]$aes_params$alpha <- 0.7
# p + scale_fill_viridis_c() + theme(legend.position = "none")
# 
# 
# p <- plot_factor(mofa_baseline_trained, factors = 1:mofa_baseline_trained@dimensions$K,  color_by = "subjectId")
# p$layers[[1]]$aes_params$alpha <- 0.7
# p + theme(legend.position = "none")
# 
# 
# p <- plot_factor(mofa_baseline_trained, factors = 1:mofa_baseline_trained@dimensions$K,  color_by = "intervention")
# p$layers[[1]]$aes_params$alpha <- 0.7
# p + theme(legend.position = "none")
# 
# 
# p <- plot_factor(mofa_baseline_trained, factors = 1:mofa_baseline_trained@dimensions$K,  color_by = "bmi_at_consent")
# p$layers[[1]]$aes_params$alpha <- 0.7
# p + scale_fill_viridis_c() + theme(legend.position = "none")
# 
# 
# # Top weights
# weights <- get_weights(mofa_baseline_trained, as.data.frame = T, scale=F, abs=F) 
# 
# 
# w.top50.l <- list()
# for(i in 1:mofa_baseline_trained@dimensions$K){
#   w.top50.l[[i]] <- arrange(weights[weights$factor==paste0("Factor",i),], desc(abs(value)))[1:50,]
#   names(w.top50.l)[i] <- paste0("Factor",i)
# }
# 
# 
# # UMAP
# mofa_baseline_trained <- run_umap(mofa_baseline_trained, n_neighbors = 20)
# plot_dimred(mofa_baseline_trained, method = "UMAP",  color_by = "subjectId", color_name = "subjectId") +
#   scale_fill_viridis_d() +
#   theme(legend.position = "none")
# 
# 
# plot_dimred(mofa_baseline_trained, method = "UMAP",  color_by = "age_at_consent", color_name = "age_at_consent") +
#   scale_fill_viridis_c() +
#   theme(legend.position = "none")
# 
# 
# plot_dimred(mofa_baseline_trained, method = "UMAP",  color_by = "intervention", color_name = "intervention") +
#   scale_fill_viridis_d() +
#   theme(legend.position = "none")
# 
# 
# plot_dimred(mofa_baseline_trained, method = "UMAP",  color_by = "interventionId", color_name = "interventionId") +
#   scale_fill_viridis_d() +
#   theme(legend.position = "none")


# 5. Association analyses with features -------------------------------------------------
outfile <- here('1-analyses/mofa/out', "baseline_mofa_medium.hdf5")
mofa_base_trained <- load_model(outfile)
factors <- as.data.frame(get_factors(mofa_base_trained)[[1]])


## 1. Phenotypical variables
feat1 <- colData(mae_base) |>
  as.data.frame() |> 
  dplyr::mutate(menopause = gsub("[*]", "", mpstatrs),
                studyarm = ifelse(interventionId == 'S', "S", "I"),
                
                # pseudocode ethanol as numeric
                etohu_curr_pseudonum = as.numeric(as.factor(etohu_curr))) |> 
  dplyr::select(studyarm, age_at_consent, menopause, bmi_at_consent,
                smoking_ever, smoking_curr, etoh_curr, etohu_curr_pseudonum, cig_before,
                diet, preg_ever, ocp_curr,
                endomet:risk_autimm) |> 
  dplyr::mutate(cig_before_all = ifelse(is.na(cig_before), 0, cig_before)) |> 
  dplyr::mutate(across(everything(), ~ ifelse(. == 'unknown', NA, .))) |> 
  dplyr::select(where(~ n_distinct(.,na.rm = T) > 1)) |>
  dplyr::mutate(across(c("menopause", "ocp_curr", "etoh_curr", "endomet", "pcos", "preg_ever", "smoking_curr", "smoking_ever"), ~ factor(., levels = c("no", "yes")))) |> 
  dplyr::mutate(across(contains("risk_"), ~ as.numeric(.))) |> 
  dplyr::mutate(diet = factor(diet, levels = c("normal", "pescetarian", "vegetarian", "vegan")),
                studyarm = factor(studyarm))


## 2. Clinical and other features
feat2 <- as.data.frame(wideFormat(data[,data$visitId=='M0',c(1, 3, 4, 5, 9, 13:15)])) |> 
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
identical(rownames(factors), rownames(feat1))
identical(rownames(factors), rownames(feat2))
features <- cbind(feat1, feat2)
saveRDS(features, file = "1-analyses/mofa/out/features_base.Rdata")


## 4.  Prepare a lookup table with every pair of columns to test
pairs_tbl <- tidyr::crossing(
  factor_col  = names(factors),
  feature_col = colnames(features)
)


## 5.  Run the appropriate test for every pair (vapply for vectorisation)
kw_eps2 <- function(x, g) {
  ok <- complete.cases(x, g)           # drop NAs the same way the test does
  g  <- as.factor(g[ok])
  x  <- x[ok]
  k  <- nlevels(g)
  n  <- length(x)
  H  <- kruskal.test(x, g)$statistic   # “Kruskal-Wallis chi-squared”
  (H - k + 1) / (n - k)
}


pairs_tbl <- pairs_tbl |> 
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




## 3.  Multiple-testing correction and nice ordering
pairs_tbl <- pairs_tbl |> 
  mutate(q_value = p.adjust(p_value, method = "fdr")) |> 
  arrange(p_value) |> as.data.frame()
saveRDS(pairs_tbl, file = here('1-analyses/mofa/out/pairs_tbl_baseline.Rds'))




# Alternative analysis: linear/log model -----------------
library(purrr)
library(tidyr)
library(broom)


## 1.  Z-score the factors so slopes are comparable
factors <- as.data.frame(get_factors(mofa_base_trained, scale = T)[[1]])
features <- readRDS("1-analyses/mofa/out/features_base.Rdata")


## 2.  All (factor, feature) pairs
pairs_tbl <- tidyr::crossing(
  factor_col  = names(factors),
  feature_col = colnames(features)
) |> 
  dplyr::filter(feature_col != 'subjectId')


## 3.  Helper that chooses the right model, returns β & p (no correction for subjectId but could include)
assoc <- function(f_vec, feat_name, features_df) {
  feat_vec <- features_df[[feat_name]]
  
  ok <- complete.cases(f_vec, feat_vec)
  df <- data.frame(y = feat_vec[ok],
                   x = f_vec[ok])
  
  if (is.numeric(df$y)) {
    fit <- lm(y ~ x, data = df)
    out <- tidy(fit)[2,]
    c(beta = out$estimate, p = out$p.value, model = "linear")
    
  } else if (!is.numeric(df$y) && !is.ordered(df$y)) {
    fit <- glm(y ~ x, family = binomial, data = df)
    out <- tidy(fit)[2,]
    c(beta = out$estimate, p = out$p.value, model = "logistic")
    
  } else if (!is.numeric(df$y) && is.ordered(df$y)) {
    fit <- polr(y ~ x, data = df, Hess = TRUE)
    co  <- summary(fit)$coefficients[1,]
    z   <- co["t value"]
    p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
    c(beta = unname(co["Value"]), p = p, model = "ordinal")
    
  } else {
    c(beta = NA, p = NA, model = "skipped")
  }
}


## 4.  Vectorised over every pair 
result <- pairs_tbl |> 
  dplyr::mutate(res = pmap(list(factor_col, feature_col),
                           ~assoc(factors[[..1]], ..2, features))) |> 
  tidyr::unnest_wider(res) |> 
  dplyr::mutate(q_value = p.adjust(p, "fdr"))

saveRDS(result, file = "1-analyses/mofa/out/pairs_tbl_regr.Rds")


# Compare linear/log and kruskall/spearman findings
pairs_tbl <- readRDS(here('1-analyses/mofa/out/pairs_tbl_baseline.Rds'))
x <- pairs_tbl[pairs_tbl$q_value < 0.1,]
y <- result[result$q_value < 0.1,]
z <- x |> 
  dplyr::select(feature_col, factor_col, q_value) |> 
  dplyr::inner_join(dplyr::select(y, feature_col, factor_col, q_value),
                    by = c("feature_col", "factor_col")) 
# actually quite consistent results :) 


# SAVE output for easy plotting ----------------------
outfile <- here('1-analyses/mofa/out', "baseline_mofa_medium.hdf5")
mofa_baseline_trained <- load_model(outfile)


fact <- as.data.frame(get_factors(mofa_base_trained, scale = T)[[1]]) |> tibble::rownames_to_column('primary')


f <- readRDS(file = "1-analyses/mofa/out/features_base.Rdata")|> 
  tibble::rownames_to_column("primary")


pheno_factors <- f  |> 
  dplyr::left_join(fact, by = "primary")


saveRDS(pheno_factors, file = "1-analyses/mofa/out/pheno_w_mofa_baseline.Rdata")


# Exploring F9 more ------------------------


load("data/data_normalized_centred.Rdata")
library(MultiAssayExperiment)
sal <- as.data.frame(longForm(data[,,40]))
pheno_w_mofa_F9 <- pheno_w_mofa |> 
  dplyr::left_join(sal, by = 'primary')


top_w_f9 <- get_weights(mofa_base, factors = 9, scale = T,
                        as.data.frame = T) |> 
  dplyr::mutate(absval = abs(value)) |>
  dplyr::arrange(desc(absval))


# Any features >0.3 & annotate labels
load("data/vars.Rdata")


enrich <- top_w_f9 |> dplyr::filter(absval >=0.35) |> 
  dplyr::mutate(feat = gsub("_Saliva nuclear magnetic resonance: normalized", "", feature)) |> 
  dplyr::left_join(vars, by = c("feat" = 'x')) |> 
  dplyr::pull(label) |> unique()




library(MetaboAnalystR)


# Create mSetObj
mSet<-InitDataObjects("list", "msetora", FALSE, 150)
# Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, enrich);
mSet<-CrossReferencing(mSet, "name");


mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
# save in out/download-F9-enrichment


## Pathway Analysis


library(MetaboAnalystR)


# Create mSetObj
mSet<-InitDataObjects("conc", "pathora", FALSE, 150)


# Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, enrich);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
# save in out/download-f9-pathway
