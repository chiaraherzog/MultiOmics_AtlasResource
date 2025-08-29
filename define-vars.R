# Define vars
library(dplyr)
library(MultiAssayExperiment)
library(here)

load("data/data_raw.Rdata")

# Preparing a list of variables to be included
vars_clin_i <- read.table("data/clinical_variables_i.csv", header = T, sep = ',') |> 
  dplyr::mutate(assay = case_when(grepl("Spiro|test|exercise", assay2) ~ "Functional sports exam",
                                  grepl("vifat|scfat", x) ~ "Vascular and body sonography",
                                  grepl("Vascular", assay2) ~ "Vascular and body sonography",
                                  grepl("composition", assay2) ~ "Body composition",
                                  grepl("Skin", assay2) ~ "Skin histology and transepidermal water loss assay",
                                  TRUE ~ assay2))

vars_clin_s <- read.table("data/clinical_variables_s.csv", header = T, sep = ',') |> 
  dplyr::mutate(assay = case_when(grepl("Spiro|test|exercise", assay2) ~ "Functional sports exam",
                                  grepl("vifat|scfat", x) ~ "Vascular and body sonography",
                                  grepl("Vascular", assay2) ~ "Vascular and body sonography",
                                  grepl("composition", assay2) ~ "Body composition",
                                  grepl("Skin", assay2) ~ "Skin histology and transepidermal water loss assay",
                                  TRUE ~ assay2))
vars_clin <- vars_clin_i |> dplyr::full_join(vars_clin_s);rm(vars_clin_i, vars_clin_s)


vars_cerv_meth_i <- readxl::read_xlsx("data/indices_i.xlsx", sheet = 1) |> 
  dplyr::mutate(assay = 'Composite methylation scores: cervical',
                label = ifelse(x == 'ic', 'immune cell proportion', label))

vars_buccal_meth_i <- readxl::read_xlsx("data/indices_i.xlsx", sheet = 2) |> 
  dplyr::mutate(assay = 'Composite methylation scores: buccal',
                label = ifelse(x == 'ic', 'immune cell proportion', label))

vars_blood_meth_i <- readxl::read_xlsx("data/indices_i.xlsx", sheet = 3)|> 
  dplyr::mutate(assay = 'Composite methylation scores: blood',
                label = ifelse(x == 'ic', 'immune cell proportion', label))


vars_cerv_meth_s <- readxl::read_xlsx("data/indices_s.xlsx", sheet = 1) |> 
  dplyr::mutate(assay = 'Composite methylation scores: cervical',
                label = ifelse(x == 'ic', 'immune cell proportion', label))

vars_buccal_meth_s <- readxl::read_xlsx("data/indices_s.xlsx", sheet = 2) |> 
  dplyr::mutate(assay = 'Composite methylation scores: buccal',
                label = ifelse(x == 'ic', 'immune cell proportion', label))

vars_blood_meth_s <- readxl::read_xlsx("data/indices_s.xlsx", sheet = 3)|> 
  dplyr::mutate(assay = 'Composite methylation scores: blood',
                label = ifelse(x == 'ic', 'immune cell proportion', label))

vars_blood_meth <- vars_blood_meth_i |> dplyr::full_join(vars_blood_meth_s);rm(vars_blood_meth_i, vars_blood_meth_s)
vars_buccal_meth <- vars_buccal_meth_i |> dplyr::full_join(vars_buccal_meth_s);rm(vars_buccal_meth_i, vars_buccal_meth_s)
vars_cerv_meth <- vars_cerv_meth_i |> dplyr::full_join(vars_cerv_meth_s);rm(vars_cerv_meth_i, vars_cerv_meth_s)


load("data/populations_names_annotated.Rdata")

vars_imm <- populations |> 
  dplyr::mutate(x = name,
                label = `population name`,
                assay = ifelse(grepl("wb", staining),
                               "Flow cytometry: white blood cell staining", 
                               ifelse(grepl("stimulation", staining),
                                      "Flow cytometry: T cell cytokines",
                                      "Flow cytometry: T cell staining")))

vars_rest <- as.data.frame(rownames(data)) |> 
  dplyr::filter(grepl("ImmAge_gen", value) | grepl("ASVs|families|normalized", group_name)) |> 
  dplyr::rename(x = value,
                assay = group_name) |> 
  dplyr::mutate(assay2 = case_when(grepl("normalized", assay) & grepl("Urine", assay) ~ "Urine metabolome",
                                   grepl("normalized", assay) & grepl("Saliva", assay) ~ "Saliva metabolome"
  )) |> 
  dplyr::select(-group)

saliva <- readRDS(here("data/RefMet_mapped_saliva.Rds")) |> 
  dplyr::mutate(Input.name = ifelse(grepl("[0-9]", substr(Input.name, 0, 1)), paste0("X", Input.name), Input.name),
                Input.name = ifelse(Input.name == 'Acetate.', "Acetate.mM.", Input.name),
                Input.name = ifelse(Input.name == 'Trimethylamine N-oxide', 'TMA..N.oxide', Input.name))

urine <- readRDS(here("data/RefMet_mapped_urine.Rds")) |> 
  dplyr::mutate(Input.name = ifelse(grepl("[0-9]", substr(Input.name, 0, 1)), paste0("X", Input.name), Input.name),
                Input.name = ifelse(Input.name == 'Acetate.', "Acetate.mM.", Input.name),
                Input.name = ifelse(Input.name == 'Trimethylamine N-oxide', 'TMA..N.oxide', Input.name))
metab <- rbind(saliva, urine) |> dplyr::distinct()


vars_rest <- vars_rest |> 
  dplyr::left_join(dplyr::select(metab, Input.name, Standardized.name), by = c("x" = 'Input.name')) |> 
  dplyr::mutate(label = ifelse(!is.na(Standardized.name), Standardized.name, x))


vars <- plyr::rbind.fill(vars_clin, vars_cerv_meth, vars_blood_meth, vars_buccal_meth, vars_imm, vars_rest)

table(vars$assay)

# Duplicate flow cytometry stainings:
sub <- vars |> dplyr::filter(grepl("cytokine", assay))
sub2 <- sub
sub$assay <- 'Flow cytometry: stimulated T cells'
sub2$assay <- 'Flow cytometry: unstimulated T cells'
vars <- vars|> dplyr::filter(!grepl("cytokine", assay))

vars <- rbind(vars, sub, sub2)
table(vars$assay)
save(vars, file = here("data/vars.Rdata"))
