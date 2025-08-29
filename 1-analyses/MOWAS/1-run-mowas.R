# MOWAS analyses
# Date: 2 Jan 2025
# Author: Chiara Herzog

# Packages
packages <- c('here', 'dplyr', 'janitor', 'broom', 'MultiAssayExperiment')

for(p in packages){
  if(!require(p,
              character.only = T)){
    install.packages(p)
  } 
}

here::i_am('1-analyses/MOWAS/1-run-mowas.R')
source(here('1-analyses/MOWAS/MOWAS.R'))

# # Age phewas
# age_mowas <- MOWAS(phenotype = 'age',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# #
# # # BMI phewas
# bmi <- MOWAS(phenotype = 'bmi',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# #
# # # VO2Max phewas
# vo2max <- MOWAS(phenotype = 'vo2max',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# #
# # # Smoking (current) phewas
# smk <- MOWAS(phenotype = 'smoking_curr',
#                         out.folder = here('1-analyses/MOWAS/1-out/'))
# #
# # # Smoking (ever) phewas
# smk_ever <- MOWAS(phenotype = 'smoking_ever',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# # # Activity
# activity <- MOWAS(phenotype = 'activity',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# # # # Ethanol
# ethanol <- MOWAS(phenotype = 'alcohol',
#                      out.folder = here('1-analyses/MOWAS/1-out/'))
# # Menopause
# meno <- MOWAS(phenotype = 'menopause',
#                    out.folder = here('1-analyses/MOWAS/1-out/'))
# 
ethanol <- MOWAS(phenotype = 'alcohol_units',
                  out.folder = here('1-analyses/MOWAS/1-out/'))