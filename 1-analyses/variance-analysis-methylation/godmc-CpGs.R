godmc <- data.table::fread("~/Documents/Work/data.nosync/mQTL_GODMC/assoc_meta_all.csv.gz")
godmc <- godmc[godmc$pval<0.05,]
unique_mqtl <- unique(godmc$cpg) # these are 450k
save(unique_mqtl, file = 'data/godmc_unique_mqtls.Rdata')

# R.utils::gunzip("~/Documents/Work/data.nosync/mQTL_GODMC/meQTL_full.txt.gz")
meqtl_epic <- data.table::fread("~/Documents/Work/data.nosync/mQTL_GODMC/meQTL_full.txt", sep = "\t")
unique_mqtl_epic <- unique(meqtl_epic$CpG)
save(unique_mqtl_epic, file = 'data/epic_unique_mqtls.Rdata')
