# DNA methylation preprocessing pipeline

# Data were preprocessed from IDAT > beta matrix using eutopsQC -------------
if(!require("eutopsQC")){
  devtools::install_github("chiaraherzog/eutopsQC")
}

library(eutopsQC)
preprocessData(input = "~/Dropbox/data/raw-data/tirolgesund/",
               output = "~/Dropbox/data/tirolgesund/",
               report = "~/Dropbox/tg-data-prep/prep-scripts/2-methylation/2-report/",
               array = "EPIC",
               by.dir = T,
               cores = 16,
               path_to_bad_sample_list = "~/Dropbox/tg-data-prep/prep-scripts/2-methylation/1-output/dupes_rm.csv", # identified using previous runs of pipeline on smaller subsets; duplicates with lower signal intensity removed
               pheno = "~/Dropbox/tg-data-prep/prep-scripts/2-methylation/1-output/pheno.Rdata",
               overwrite = T,
               save.rs = T)

# Output from save.rs was used to conduct SNP analysis (not shown)


