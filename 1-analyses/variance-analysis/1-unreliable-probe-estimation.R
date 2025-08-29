# Filter methylation data via MI, keeping most reliable probes.
# Author: Chiara Herzog
# Date: 20 Dec 2024

# Packages
packages <- c('fs', 'minfi', 'devtools', 'here')

for(p in packages){
  if(!require(p,
              character.only = T)){
    install.packages(p)
  } 
}

if(!require('epicMI')){
  devtools::install_github("ChVav/epicMI")
}

here::i_am('1-analyses/variance-analysis/1-unreliable-probe-estimation.R')

# Load data - RGset

## Get basenames
db <- fs::path_expand("~/Dropbox/")
load(file.path(db, "data/tirolgesund/pheno/methylation_replication_pheno.Rdata"))


# Compute epicMI on samples
## Separate for each sample type, because they may have different characteristics.
library(epicMI)
sampletype <- c('cervical', 'buccal', 'blood')

for (s in sampletype){
  cat(s, 'sample unreliability assessment starting')
  ids <- pheno_replication[pheno_replication$sampletype==s,]$basename
  
  ## Compile list of files with these basenames
  input <- file.path(db, "data/raw-data/tirolgesund/")
  files <- list.files(input, full.names = T, recursive = T, pattern = ".idat")
  files <- files[grepl(paste0(ids, collapse = "|"), files)]
  
  ## reformat files for metharray read-in
  files <- unique(gsub("_Red.idat|_Grn.idat", "", files))
  
  ## read in this list of files
  assign(paste0("RGset_", s), minfi::read.metharray(basenames = files,
                                 verbose = T,
                                 force = T))
  
  assign(paste0("out_", s),
         eval(parse(text = paste0("unreliability_MI(RGset_", s, ")")))
         )
}

save(out_blood, out_cervical, out_buccal, file = here('1-analyses/variance-analysis/1-out/unreliability_measures.Rdata'))
rm(RGset_blood, RGset_buccal, RGset_cervical);gc()

# Visualise
tmp <- data.frame(blood = out_blood$unreliability,
                  buccal = out_buccal$unreliability,
                  cervical = out_cervical$unreliability)

tmp |> 
  ggplot(aes(x = blood,
             y = buccal)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()

tmp |> 
  ggplot(aes(x = cervical,
             y = buccal)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()

tmp |> 
  ggplot(aes(x = blood,
             y = cervical)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()

tmp <- data.frame(blood = out_blood$MI,
                  buccal = out_buccal$MI,
                  cervical = out_cervical$MI)

tmp |> 
  ggplot(aes(x = blood,
             y = buccal)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()

tmp |> 
  ggplot(aes(x = cervical,
             y = buccal)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()

tmp |> 
  ggplot(aes(x = blood,
             y = cervical)) +
  geom_hex(bins = 500) +
  ggpubr::stat_cor()
# MI correlates better than unreliability across datasets.


# Keeping the top 10% reliable probes that are also in the beta matrix after filtering
load(file.path(db, "data/tirolgesund/beta_merged_replicate.Rdata"))
beta_rep <- beta_merged
load(file.path(db, "data/tirolgesund/beta_merged.Rdata"))

reliable_blood <- out_blood |> 
  dplyr::filter(probe %in% rownames(beta_rep) & probe %in% rownames(beta_merged)) |> 
  dplyr::arrange(dplyr::desc(MI)) |> 
  dplyr::slice(1:(n()*0.1))

reliable_buccal <- out_buccal |> 
  dplyr::filter(probe %in% rownames(beta_rep) & probe %in% rownames(beta_merged)) |> 
  dplyr::arrange(dplyr::desc(MI)) |> 
  dplyr::slice(1:(n()*0.1))

reliable_cervical <- out_cervical |> 
  dplyr::filter(probe %in% rownames(beta_rep) & probe %in% rownames(beta_merged)) |> 
  dplyr::arrange(dplyr::desc(MI)) |> 
  dplyr::slice(1:(n()*0.1))


save(reliable_blood, reliable_buccal, reliable_cervical, file = here('1-analyses/variance-analysis/1-out/reliable_sets.Rdata'))

# Save reliable subsets
load(file.path(db, "data/tirolgesund/pheno/methylation_pheno.Rdata"))
pheno_buccal <- pheno[pheno$sampletype=='buccal',]
beta_reliable_buccal <- beta_merged[reliable_buccal$probe, pheno_buccal$basename]
identical(pheno_buccal$basename, colnames(beta_reliable_buccal))
identical(rownames(beta_reliable_buccal), reliable_buccal$probe)
save(beta_reliable_buccal, pheno_buccal, file = here('1-analyses/variance-analysis/1-out/reliable_buccal.Rdata'))

pheno_blood <- pheno[pheno$sampletype=='blood',]
beta_reliable_blood <- beta_merged[reliable_blood$probe, pheno_blood$basename]
identical(pheno_blood$basename, colnames(beta_reliable_blood))
identical(rownames(beta_reliable_blood), reliable_blood$probe)
save(beta_reliable_blood, pheno_blood, file = here('1-analyses/variance-analysis/1-out/reliable_blood.Rdata'))

pheno_cervical <- pheno[pheno$sampletype=='cervical',]
beta_reliable_cervical <- beta_merged[reliable_cervical$probe, pheno_cervical$basename]
identical(pheno_cervical$basename, colnames(beta_reliable_cervical))
identical(rownames(beta_reliable_cervical), reliable_cervical$probe)
save(beta_reliable_cervical, pheno_cervical, file = here('1-analyses/variance-analysis/1-out/reliable_cervical.Rdata'))

# fin - these sets can be used for computation