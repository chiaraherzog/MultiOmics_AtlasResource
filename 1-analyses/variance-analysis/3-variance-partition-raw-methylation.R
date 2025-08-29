# Variance partition of methylation data
# Author: Chiara Herzog
# Date: 20 Dec 2024

# Packages
packages <- c('fs', 'variancePartition', 'here', 'dplyr')

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

here::i_am('1-analyses/variance-analysis/3-variance-partition-raw-methylation.R')

form <- ~ (1 | subjectId) + (1 | interventionId_bin)

## 1 - Blood
cat("Starting blood analysis...\n")
load(here('1-analyses/variance-analysis/1-out/reliable_blood.Rdata'))

mat <- as.matrix(beta_reliable_blood)
pheno_blood$interventionId_bin <- ifelse(pheno_blood$interventionId=='S', 'S', 'I')
p <- pheno_blood
varPart <- fitExtractVarPartModel(mat, form, p)
vp_blood <- sortCols(varPart)
save(vp_blood, file = here('1-analyses/variance-analysis/3-out/vp_blood_methylation.Rdata'))

# vp_df <- as.data.frame(vp_blood) |> 
#   tidyr::pivot_longer(everything(),
#                       names_to = 'type',
#                       values_to = 'value')
# 
# vp_df |> 
#   ggplot(aes(x = type,
#              y = value)) +
#   geom_violin()
# 
# vp_df <- as.data.frame(vp_blood)
# 
# library(ggtern)
# vp_df |> 
#   ggtern(aes(x = subjectId,
#              y = Residuals,
#              z = interventionId)) +
#   geom_hex_tern()

rm(pheno_blood, varPart, vp_blood, beta_reliable_blood, mat, p);gc()

## 2 - Buccal
cat("Starting buccal analysis...\n")
load(here('1-analyses/variance-analysis/1-out/reliable_buccal.Rdata'))

mat <- as.matrix(beta_reliable_buccal)
pheno_buccal$interventionId_bin <- ifelse(pheno_buccal$interventionId=='S', 'S', 'I')
p <- pheno_buccal

varPart <- fitExtractVarPartModel(mat, form, p)
vp_buccal<- sortCols(varPart)
save(vp_buccal, file = here('1-analyses/variance-analysis/3-out/vp_buccal_methylation.Rdata'))

rm(pheno_buccal, varPart, vp_buccal, beta_reliable_buccal, mat, p);gc()

## 2 - cervical
cat("Starting cervical analysis...\n")
load(here('1-analyses/variance-analysis/1-out/reliable_cervical.Rdata'))

mat <- as.matrix(beta_reliable_cervical)
pheno_cervical$interventionId_bin <- ifelse(pheno_cervical$interventionId=='S', 'S', 'I')

p <- pheno_cervical

varPart <- fitExtractVarPartModel(mat, form, p)
vp_cervical<- sortCols(varPart)
save(vp_cervical, file = here('1-analyses/variance-analysis/3-out/vp_cervical_methylation.Rdata'))

rm(pheno_cervical, varPart, vp_cervical, beta_reliable_cervical, mat, p);gc()

# # Some tests 
# library(here)
# 
# 
# load(here('1-analyses/variance-analysis/3-out/vp_buccal_methylation.Rdata'))
# vp_df <- as.data.frame(vp_buccal) |>
#   tidyr::pivot_longer(everything(),
#                       names_to = 'type',
#                       values_to = 'value')
# 
# vp_df |>
#   ggplot(aes(x = type,
#              y = value)) +
#   geom_violin() +
#   geom_boxplot()
# 
# vp_df_buccal <- as.data.frame(vp_buccal)
# 
# library(ggtern)
# buccal <- vp_df_buccal |>
#   ggtern(aes(x = subjectId,
#              y = interventionId,
#              z = Residuals)) +
#   geom_hex_tern(aes(fill = after_stat(log(density))),
#                 bins = 75) +
#   scale_fill_viridis_c(limits = c(-12, -2)) +
#   theme_bw() +
#   labs(x = 'Inter-individual',
#        y = 'Study arm',
#        z = 'Intra-individual')
# 
# vp_df_blood <- as.data.frame(vp_blood)
# 
# library(ggtern)
# blood <- vp_df_blood |>
#   ggtern(aes(x = subjectId,
#              y = interventionId,
#              z = Residuals)) +
#   geom_hex_tern(aes(fill = after_stat(log(density))),
#                 bins = 75) +
#   scale_fill_viridis_c(limits = c(-12, -2)) +
#   theme_bw() +
#   labs(x = 'Inter-individual',
#        y = 'Study arm',
#        z = 'Intra-individual')
# 
# library(patchwork)
# 
# (buccal/blood) & plot_layout(guides = 'collect')

