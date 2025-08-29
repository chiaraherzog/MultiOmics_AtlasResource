# Extract example CpGs for plotting
load('1-analyses/variance-analysis-methylation/5-output/vp_results.Rdata')
load('1-analyses/variance-analysis-methylation/3-output/celltype_association.Rdata')
load("~/Dropbox/data/tirolgesund/beta_merged.Rdata")

stopifnot(identical(pheno$basename, colnames(beta_merged)))

# Get Top representative cpgs
vp_top_var <- vp_results |> 
  dplyr::filter(category != 'non-variable (stable)') |> 
  dplyr::group_by(category) |> dplyr::arrange(desc(dominant_fraction)) |> dplyr::slice(1:100)

vp_top_nvar <- vp_results |> 
  dplyr::filter(category == 'non-variable (stable)') |> 
  dplyr::arrange(overall_sd) |> dplyr::slice(1:100)

vp_top <- rbind(vp_top_var, vp_top_nvar);rm(vp_top_var, vp_top_nvar)

b <- beta_merged[vp_top$cg,]

pheno_w_cpgs <- cbind(pheno, t(b)) |> 
  tidyr::pivot_longer(any_of(c(vp_top$cg)),
                      names_to = 'cg',
                      values_to = 'beta') |> 
  dplyr::left_join(vp_top, by = 'cg')
pheno_w_cpgs$studyarm <- ifelse(pheno_w_cpgs$interventionId == 'S', 'S', 'I')


save(b, pheno_w_cpgs, vp_top, file = '1-analyses/variance-analysis-methylation/6-output/plotting-outputs.Rdata')

library(ggplot2)
# Individual
indiv <- vp_top[vp_top$category=='individual (stable)',]$cg[1:2]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% indiv) |> 
  ggplot(aes(x = forcats::fct_reorder2(subjectId.x, cg, beta),
             y = beta)) +
  geom_point(aes(colour = sampletype.x)) +
  facet_wrap(~cg)


# sampletype
sample <- vp_top[vp_top$category=='tissue (stable)',]$cg[1:2]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% sample) |> 
  ggplot(aes(x = forcats::fct_reorder2(subjectId.x, cg, beta),
             y = beta)) +
  geom_point(aes(colour = sampletype.x)) +
  facet_wrap(~cg)

# sampletype
intervention <- vp_top[vp_top$category=='intervention (stable)',]$cg[1:2]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% intervention) |> 
  ggplot(aes(x = forcats::fct_reorder2(studyarm, cg, beta),
             y = beta)) +
  geom_point(aes(colour = sampletype.x))



# sampletype
nv <- vp_top[vp_top$category=='non-variable (stable)',]$cg[1:4]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = forcats::fct_reorder2(studyarm, cg, beta),
             y = beta)) +
  geom_point(aes(colour = sampletype.x)) +
  facet_wrap(~cg)


# time
time (malleable)
nv <- vp_top[vp_top$category=='time (malleable)',]$cg[1:4]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = visitId.x,
             y = beta)) +
  geom_point(aes(colour = sampletype.x)) +
  facet_wrap(~cg)

# intervention x time
nv <- vp_top[vp_top$category=='intervention × time (malleable)',]$cg[1:2]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = visitId.x,
             y = beta)) +
  geom_boxplot(aes(fill = studyarm)) +
  facet_wrap(~cg)


# intervention x time
nv <- vp_top[vp_top$category=='tissue × intervention (stable)',]$cg[1:5]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = sampletype.x,
             y = beta)) +
  geom_boxplot(aes(fill = studyarm)) +
  facet_wrap(~cg)

# intervention x time
nv <- vp_top[vp_top$category=='tissue × intervention (stable)',]$cg[1:5]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = sampletype.x,
             y = beta)) +
  geom_boxplot(aes(fill = studyarm)) +
  facet_wrap(~cg)

# intervention x time
nv <- vp_top[vp_top$category=='tissue × time (malleable)',]$cg[1:5]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = sampletype.x,
             y = beta)) +
  geom_boxplot(aes(fill = visitId.x)) +
  facet_wrap(~cg)


# intervention x time
nv <- vp_top[vp_top$category=='tissue × time × intervention (malleable)',]$cg[1:5]

pheno_w_cpgs |> 
  dplyr::filter(cg %in% nv) |> 
  ggplot(aes(x = sampletype.x,
             y = beta)) +
  geom_boxplot(aes(fill = visitId.x)) +
  facet_wrap(studyarm~cg, nrow = 2)




# Plot examples - not really that helpful

sample <- sample(1:nrow(vp_results), size = 50000, replace = F)

vp_sub <- vp_results[sample,]
vp_sub <- na.omit(vp_sub)

m_subsample <- as.matrix(vp_sub[,2:10])
set.seed(1234)
umap <- uwot::umap(m_subsample)
umap

vp_sub_umap <- cbind(vp_sub, umap)

vp_sub_umap$`1`


vp_sub_umap |> 
  ggplot(aes(x = `1`,
             y = `2`)) +
  geom_hex()

vp_sub_umap |> 
  ggplot(aes(x = `1`,
             y = `2`)) +
  geom_point(aes(colour = category),
             size = 0.1, alpha = 0.5)


# PCA -----
sample <- sample(1:nrow(vp_results), size = 30000, replace = F)

vp_sub <- vp_results[sample,]
vp_sub <- na.omit(vp_sub)

m_subsample <- as.matrix(vp_sub[,2:10])
set.seed(1234)
pc <- prcomp(m_subsample, scale. = T, center = T)
vp_sub_pc <- cbind(vp_sub, pc$x[,1:2])

vp_sub_pc |> 
  ggplot(aes(x = PC1,
             y = PC2)) +
  geom_hex()

vp_sub_pc |> 
  ggplot(aes(x = PC1,
             y = PC2)) +
  geom_point(aes(colour = category),
             size = 0.5, alpha = 0.5)


