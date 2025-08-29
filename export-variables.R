# Export vars table
load("data/vars.Rdata")

vars_out <- vars |> dplyr::filter(!grepl("Immune age|_clr", assay)) |> 
  dplyr::mutate(layer = ifelse(!is.na(assay2), assay2, assay)) |> 
  dplyr::select(layer, label, x)

# For clinical variables, append snomed
load("~/Dropbox/data/tirolgesund/lut.Rdata")

vars_out_lut <- vars_out |> 
  dplyr::left_join(lut, by = c("x" = 'variable')) |> 
  dplyr::mutate(access = ifelse(!is.na(access), access,
                                ifelse(grepl("microbiome|methylation", layer), 'controlled',
                                       'open')),
                intervention = ifelse(is.na(intervention), "both", intervention),
                variable_type = ifelse(is.na(variable_type), "numeric", intervention)) |> 
  dplyr::select(-c(datatype, possible_values)) |> 
  dplyr::rename(feature = label,
                variable_name = x)

write.table(vars_out_lut, file = 'out/Suppl-Variables.csv',row.names = F, sep = ",")
writexl::write_xlsx(vars_out_lut, path = 'out/Suppl-Variables.xlsx')                                
