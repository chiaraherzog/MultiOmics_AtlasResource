# Age changes and longitudinal datasame or opposing direction

# libraries ---------------
library(dplyr)
library(tidyr)
library(limma)

# setup 
dt_years <- 0.5 # delta years
eps <- 0.10   # ±10% band for "unchanged"

# functions
cat_from_ratio <- function(r, eps = eps){
  if (is.na(r)) return(NA_character_)
  if (r < 0)             return("age-opposing")  # reversed vs ageing
  if (r <= (1 - eps))    return("attenuated")    # same dir, smaller than expected
  if (r <  (1 + eps))    return("unchanged")     # ~expected
  return("accelerated")                           # same dir, larger than expected
}

# 0: Prep: age-slope & sign per feature ----------

## Load MOWAS data
load("1-analyses/MOWAS/1-out/mowas_nonmethylation_age.Rdata")

## Filter, append sign
total <- total |>
  dplyr::group_by(type) |> 
  dplyr::mutate(padj = p.adjust(p.value, method = 'fdr')) |> 
  dplyr::filter(padj < 0.05) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(sign = sign(estimate))

## Append variable names for filtering & merging of final data
load("data/vars.Rdata")
total_vars <- total |> 
  dplyr::left_join(dplyr::select(vars, assay, x), by = c('var' = 'x'),
                   relationship = 'many-to-many') |> 
  dplyr::filter(!(grepl("buccal", type) & !grepl("buccal", assay)) &
                  !(grepl("blood", type) & !grepl("blood", assay)) & 
                  !(grepl("cervical", type) & !grepl("cervical", assay)) &
                  !(grepl("saliva", type) & !grepl("Saliva", assay)) &
                  !(grepl("urine", type) & !grepl("Urine", assay)) &
                  !(grepl("faecal", type) & !grepl("Stool", assay)) &
                  !(grepl("_stim", type) & grepl("unstimulated", assay))
  ) |> 
  dplyr::rename(rowname = var)

## Extract slope
age_slopes <- total_vars |>
  select(assay, rowname, estimate, sign) |>
  distinct() |>
  rename(beta_age = estimate)

## As list (for subsetByAssay)
varList <- total_vars |> 
  dplyr::select(assay, rowname) |> 
  dplyr::group_by(assay) |> 
  dplyr::summarise(features = list(rowname), .groups = 'drop') |> 
  tibble::deframe()

# 1. Get data ----------
## Load and filter baseline-reference data
### Using normalised data as these were also used for MOWAS
load("data/data_normalized_baseline.Rdata")
library(MultiAssayExperiment)
df_filtered <- subsetByAssay(data, varList)

## Place data into long format
df_delta <- as.data.frame(longForm(df_filtered,
                                   colDataCols = c("age_at_consent", "bmi_at_consent", "comprate", "compliance", "visitId",
                                                   'interventionId',
                                                   'subjectId',
                                                   'compliance_smkgroup'))) |> 
  dplyr::filter(!is.na(value)) |> # must have Δ 
  dplyr::rename(delta = value) |> 
  
  # append slope
  left_join(age_slopes, by = c("assay","rowname")) |>
  filter(!is.na(beta_age), !is.na(sign)) |> 
  
  # Align Δ to ageing direction & expected drift magnitude
  dplyr::mutate(
    delta_aligned = delta * sign,                     # >0 older, <0 younger
    exp_drift     = beta_age * dt_years,              # expected Δ from baseline slope
    exp_mag       = abs(exp_drift)
  )


## 2. Compute changes - IF -------------------------
df_delta_m6_i <- df_delta |> dplyr::filter(interventionId != 'S' & visitId == 'M6') |> 
  
  # binary compliance
  dplyr::mutate(comp_bin = ifelse(comprate > median(comprate), 'high', 'low')) 
# |>
  # dplyr::mutate(comp_bin = ifelse(compliance == 'high', 'high', 'low'))

summary_m6_i <- df_delta_m6_i |> 
  
  dplyr::group_by(assay, rowname) |> 
  dplyr::summarise(
    mean_aligned_hi = mean(delta_aligned[comp_bin=='high'], na.rm = TRUE),
    mean_aligned_lo = mean(delta_aligned[comp_bin!='high'], na.rm = TRUE),
    exp_mag         = dplyr::first(exp_mag),
    beta_age        = dplyr::first(beta_age),
    .groups = "drop"
  ) |> 
  dplyr::mutate(
    ratio_hi = mean_aligned_hi / exp_mag,
    ratio_lo = mean_aligned_lo / exp_mag
  ) |>
  
  # catgorise
  dplyr::rowwise() |> 
  mutate(
    category_hi = cat_from_ratio(ratio_hi, eps),
    category_lo = cat_from_ratio(ratio_lo, eps),
    delta_ratio = ratio_hi - ratio_lo            # <0 => more flattening in high vs low
  ) |> 
  ungroup()

### inference -------------------------
# Build features×subjects matrix of delta_aligned; design uses comp_cont + arm
pheno <- df_delta_m6_i |>  distinct(subjectId, interventionId, comp_bin, comprate)

mat_long <- df_delta_m6_i |> 
  mutate(feat = paste(assay, rowname, sep="|")) |> 
  select(feat, subjectId, delta_aligned)

Y_align <- mat_long |> 
  pivot_wider(names_from = subjectId, values_from = delta_aligned) |> 
  arrange(feat)

feat_ids <- Y_align$feat
Y_align <- as.matrix(Y_align[ , -1, drop = FALSE])
rownames(Y_align) <- feat_ids

# align pheno to columns
pheno <- pheno |> 
  filter(subjectId %in% colnames(Y_align)) |> 
  arrange(factor(subjectId, levels = colnames(Y_align)))
stopifnot(identical(pheno$subjectId, colnames(Y_align)))

design <- model.matrix(~ comprate, data = pheno)
fit <- lmFit(Y_align, design); fit <- eBayes(fit)
tab <- topTable(fit, coef = "comprate", number = Inf, sort.by = "none")
tab$feat <- rownames(tab)
tab_i <- tab |> 
  separate(feat, into = c("assay","rowname"), sep = "\\|", remove = FALSE) |> 
  left_join(summary_m6_i |>  select(assay,rowname,beta_age,ratio_hi,ratio_lo,delta_ratio,category_hi,category_lo),
            by = c("assay","rowname"))

## 3. Compute changes - smk -------------------------
df_delta_m6_s <- df_delta |> dplyr::filter(interventionId == 'S' & visitId == 'M6') |> 
  
  # binary compliance
  dplyr::mutate(comp_bin = ifelse(compliance_smkgroup == 'complete cessation', 'high', 'low'))

summary_m6_s <- df_delta_m6_s |> 
  
  dplyr::group_by(assay, rowname) |> 
  dplyr::summarise(
    mean_aligned_hi = mean(delta_aligned[comp_bin=='high'], na.rm = TRUE),
    mean_aligned_lo = mean(delta_aligned[comp_bin!='high'], na.rm = TRUE),
    exp_mag         = dplyr::first(exp_mag),
    beta_age        = dplyr::first(beta_age),
    .groups = "drop"
  ) |> 
  dplyr::mutate(
    ratio_hi = mean_aligned_hi / exp_mag,
    ratio_lo = mean_aligned_lo / exp_mag
  ) |>
  
  # catgorise
  dplyr::rowwise() |> 
  mutate(
    category_hi = cat_from_ratio(ratio_hi, eps),
    category_lo = cat_from_ratio(ratio_lo, eps),
    delta_ratio = ratio_hi - ratio_lo            # <0 => more flattening in high vs low
  ) |> 
  ungroup()

### inference -------------------------
# Build features×subjects matrix of delta_aligned; design uses comp_cont + arm
pheno <- df_delta_m6_s |>  distinct(subjectId, interventionId, comp_bin, comprate)

mat_long <- df_delta_m6_s |> 
  mutate(feat = paste(assay, rowname, sep="|")) |> 
  select(feat, subjectId, delta_aligned)

Y_align <- mat_long |> 
  pivot_wider(names_from = subjectId, values_from = delta_aligned) |> 
  arrange(feat)

feat_ids <- Y_align$feat
Y_align <- as.matrix(Y_align[ , -1, drop = FALSE])
rownames(Y_align) <- feat_ids

# align pheno to columns
pheno <- pheno |> 
  filter(subjectId %in% colnames(Y_align)) |> 
  arrange(factor(subjectId, levels = colnames(Y_align)))
stopifnot(identical(pheno$subjectId, colnames(Y_align)))

design <- model.matrix(~ comp_bin, data = pheno)
fit <- lmFit(Y_align, design); fit <- eBayes(fit)
tab <- topTable(fit, number = Inf, sort.by = "none")
tab$feat <- rownames(tab)
tab_s <- tab |> 
  separate(feat, into = c("assay","rowname"), sep = "\\|", remove = FALSE) |> 
  left_join(summary_m6_s |>  select(assay,rowname,beta_age,ratio_hi,ratio_lo,delta_ratio,category_hi,category_lo),
            by = c("assay","rowname"))


### Individual change rate --------------------
df_personal <- df_delta |> 
  dplyr::filter(visitId %in% c('M2', 'M4', 'M6')) |> 
  mutate(
    ratio = delta_aligned / exp_mag,
    category = case_when(
      ratio < 0 ~ "opposing",
      ratio <= (1 - eps) ~ "attenuated",
      ratio <  (1 + eps) ~ "unchanged",
      TRUE               ~ "accelerated"
    ),
    val = case_when(
      category == "opposing"    ~  1,
      category == "accelerated" ~ -1,
      TRUE                      ~  0
    )
  ) |> 
  group_by(subjectId, visitId) |> 
  summarise(
    comprate = unique(comprate),
    compliance = unique(compliance),
    compliance_smk = unique(compliance_smkgroup),
    personal_net_index = mean(val, na.rm = TRUE),   # +1 = more opposing, -1 = more accelerated
    .groups = "drop"
  )

df_personal_no_methylation  <- df_delta |> 
  dplyr::filter(visitId %in% c('M2', 'M4', 'M6') & !grepl("methylation", assay)) |> 
  mutate(
    ratio = delta_aligned / exp_mag,
    category = case_when(
      ratio < 0 ~ "opposing",
      ratio <= (1 - eps) ~ "attenuated",
      ratio <  (1 + eps) ~ "unchanged",
      TRUE               ~ "accelerated"
    ),
    val = case_when(
      category == "opposing"    ~  1,
      category == "accelerated" ~ -1,
      TRUE                      ~  0
    )
  ) |> 
  group_by(subjectId, visitId) |> 
  summarise(
    comprate = unique(comprate),
    compliance = unique(compliance),
    compliance_smk = unique(compliance_smkgroup),
    personal_net_index = mean(val, na.rm = TRUE),   # +1 = more opposing, -1 = more accelerated
    .groups = "drop"
  )


# thresholds
p_cut <- 0.05

# 1) pick per-feature effects from the I and S tables
sig_i <- tab_i %>%
  mutate(atten_i = (P.Value < p_cut) & (logFC < 0)) %>%           # <0 = more 'younger' with higher compliance
  select(assay, rowname, logFC_i = logFC, P_i = P.Value, atten_i)

sig_s <- tab_s %>%
  mutate(atten_s = (P.Value < p_cut) & (logFC < 0)) %>%
  select(assay, rowname, logFC_s = logFC, P_s = P.Value, atten_s)

# 2) join and classify overlap
overlap_inf <- sig_i %>%
  full_join(sig_s, by = c("assay","rowname"))%>%
  mutate(altered_by = dplyr::case_when(
   atten_i == T &  atten_s == T ~ "Both",
    atten_i == T & atten_s != T ~ "I only",
   atten_s == T & atten_i != T ~ "S only",
    TRUE                               ~ "Neither"
  ))

# Mark attenuation at high compliance by ratio criterion
rat_i <- summary_m6_i %>% mutate(attn_I = ratio_hi <= (1 - eps))
rat_s <- summary_m6_s %>% mutate(attn_S = ratio_hi <= (1 - eps))

overlap_ratio <- rat_i %>%
  select(assay, rowname, attn_I) %>%
  full_join(rat_s %>% select(assay, rowname, attn_S), by = c("assay","rowname")) %>%
  mutate(group = case_when(
    attn_I & attn_S ~ "Both",
    attn_I & !attn_S ~ "I only",
    !attn_I & attn_S ~ "S only",
    TRUE ~ "Neither"
  ))


save(age_slopes,
     sig_s, sig_i,
     summary_m6_i,
     summary_m6_s,
     tab_s, tab_i,
     overlap_inf,
     overlap_ratio,
     
     df_personal_no_methylation,
     file = '1-analyses/malleability/1-out/ageOut.Rdata')
     
