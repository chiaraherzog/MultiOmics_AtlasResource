# PID of current job: 3985245
mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
cmpd.vec<-c("Tryptophan","Arginine","Thymine","Hydroxyacetone","Choline","Serine","Alanine","Taurine","Leucine","Methanol","Propionic acid","3-Methylhistidine","Trimethylamine N-oxide","Galactose","Xylose","Creatinine","Threonine","Proline","Fucose","Formic acid","Malonic acid","Desaminotyrosine","Sarcosine","Dimethylamine","3-Methyl-2-oxovaleric acid","Valine","Phenylalanine")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA)
mSet<-SaveTransformedData(mSet)
