library(dplyr)

# Write supplementary table for variance decomposition
load('1-analyses/variance-analysis-methylation/5-output/vp_results.Rdata')

# Annotate meQTL -----------
load('data/godmc_unique_mqtls.Rdata')
load('data/epic_unique_mqtls.Rdata')

tmp1 <- vp_results |> 
  dplyr::mutate(godmc = cg %in% unique_mqtl,
                epic = cg %in% unique_mqtl_epic,
                meQTL = case_when(category != 'individual (stable)' ~ NA,
                                  godmc == T & epic == T ~ 'both GoDMC & EPIGEN',
                                  godmc == T~ 'GoDMC',
                                  epic == T ~ 'EPIGEN',
                                  TRUE ~ 'new'))
table(tmp1$meQTL, useNA = 'always')


# Format  -------------------

out <- tmp1 |>
  dplyr::rename_with(~paste0("varianceExplained_", .), .cols = c("sampletype":"Residuals")) |> 
  dplyr::rename(noise_threshold = thresh,
                driver_tissue_malleable = driver_tissue_screen) |> 
  dplyr::select(cg, category, dominant, varianceExplained_sampletype:varianceExplained_Residuals, 
                residual_flag,
                noise_threshold,
                meQTL,
                celltype_allocation, celltype_type, celltype_effect, celltype_pval,
                driver_tissue_malleable)

writexl::write_xlsx(out, path = 'out/Suppl_Methylation_Variance.xlsx')  
