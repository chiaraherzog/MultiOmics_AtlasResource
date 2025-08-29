# Variance partition on data - expanded
# Author: Chiara Herzog
# Date: 1 Aug 2025

# Packages -----------
cat("Load Libs ... ")
packages <- c('fs', 'variancePartition', 'here', 'dplyr', 'purrr', 'matrixStats', 'MultiAssayExperiment')

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

cat("done.\n")


here::i_am('1-analyses/variance-analysis/2b-variance-partition-extended.R')

# Load data --------
load(here('data/data_normalized.Rdata'))

# Keep only M0-M6 data
data <- data[,data$visitId %in% c("M0", "M2", "M4", "M6"),]

# Load variables to use (except methylation)
load(here('data/vars.Rdata')) # Variables were defined previously for IF analysis.
vars_if <- vars
load(here('data/vars_smk.Rdata')) # Variables were defined previously for SMK analysis.
vars <- rbind(vars_if, vars) |> 
  dplyr::distinct()

vars <- vars |> # Keep ony those assays with >=5 features
  dplyr::group_by(assay) |> 
  dplyr::filter(dplyr::n()>=5) |> 
  dplyr::ungroup()
rm(vars_if);gc()

# Experiment list for variance partition
experiments <- unique(vars$assay)
experiments <- experiments[experiments != 'Body composition']

data$studyarm <- ifelse(data$interventionId=='S', 'S', 'I')

vp <- vector("list", 
             length = length(experiments))

for(e in experiments){
  cat(e, "variance partition starting\n")
  
  suppressWarnings(suppressMessages(mat <- data[,,e][[1]]))
  suppressWarnings(suppressMessages(p <- as.data.frame(colData(data[,,e]))))
  
  x <- substr(colnames(mat), 1, 6)
  if(identical(rownames(p), x)){colnames(mat) <- rownames(p)} else {stop("check data input")}
  
  
  
  ## Formula
  form <- ~ 
    (1 | subjectId) + 
    (1 | visitId) + 
    (1 | studyarm) + 
    (1 | studyarm:visitId)
  
  suppressWarnings(suppressMessages(varPart <- fitExtractVarPartModel(mat, form, p)))
  
  n <- janitor::make_clean_names(e)
  
  vp[[which(experiments==e)]] <- sortCols(varPart)
  names(vp)[which(experiments==e)] <- n
  
  gc()
  
}

save(vp, file = here('1-analyses/variance-analysis/2-out/vp_expanded_w_o_raw_methylation.Rdata'))

vp_df <- vp |> 
  purrr::map_dfr(
    .f = as.data.frame, # Convert each element to a data frame
    .id = "group"       # Add list names as a column
  )

vp_df <- vp_df |> dplyr::filter(
  ! group %in% c("saliva_microbiome_as_vs_clr", "saliva_microbiome_families_clr",
                 "stool_microbiome_as_vs_clr", "stool_microbiome_families_clr"))


# Sort by min -> max intraindividual
order <- vp_df |>
  dplyr::group_by(group) |> 
  dplyr::summarise(x = median(Residuals)) |> 
  dplyr::arrange(dplyr::desc(x))
order <- unique(order$group)

order_vars <- vp_df |> 
  tibble::rownames_to_column('var') |> 
  dplyr::mutate(group = factor(group, levels = order)) |> 
  dplyr::arrange(group, -Residuals)

# Create the group annotations (tile)
load(here('data/vars.Rdata')) # Variables were defined previously for IF analysis.
vars_if <- vars
load(here('data/vars_smk.Rdata')) # Variables were defined previously for SMK analysis.
vars <- rbind(vars_if, vars) |> 
  dplyr::distinct() |> 
  dplyr::mutate(assay3 = janitor::make_clean_names(assay,allow_dupes = T)) |> 
  dplyr::select(assay, assay2, assay3) |> dplyr::distinct()

labels <- vars[match(order, vars$assay3),]
labels <- labels |> dplyr::mutate(label = ifelse(!is.na(assay2), assay2, assay)) |> 
  dplyr::pull(label)

group_annotations <- order_vars |>
  dplyr::group_by(group) |> 
  dplyr::mutate(var_ordered = cur_group_id()) |> 
  dplyr::distinct(var, group, var_ordered) |>
  dplyr::mutate(y = 1.07,  # Fixed y position for the annotation
                height = 0.1) # Height of the tile

# # Set colours
annotation_colours <- colorRampPalette(cols)(length(order))
names(annotation_colours) <- unique(order)

# Create the plot
var_barplot <- order_vars |>
  tidyr::pivot_longer(subjectId:Residuals,
                      names_to = "type",
                      values_to = "value") |>
  dplyr::mutate(group = factor(group, levels = rev(order)),
                var = factor(var, levels = rev(order_vars$var)),
  ) |> 
  ggplot(aes(x = var,
             y = value,
             fill = type)) +
  geom_col() +
  scale_fill_viridis_d(name = 'Variance explained (%)',
                       option = 'G') +
  coord_flip() +
  theme(
    axis.text.y = element_blank()
  ) +
  ggnewscale::new_scale_fill() +  # Reset the fill scale for the group annotation
  coord_cartesian(clip = "off") +
  geom_tile(data = group_annotations,
            aes(x = var,
                y = y,  # Fixed position for the annotation
                fill = group,  # Use group to color tiles
                height = height),
            width = 1, inherit.aes = FALSE) +
  scale_fill_viridis_d(labels = labels)+
  # scale_fill_manual(values = annotation_colours) +
  coord_flip() +
  theme(
    plot.margin = margin(10, 10, 50, 10),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(y = 'Proportion of\nvariance explained')


order_vars <- vp_df |> 
  tibble::rownames_to_column('var') |> 
  dplyr::mutate(group = factor(group, levels = order)) |> 
  dplyr::arrange(group, -Residuals)

# Values for in-text stats:
vals <- order_vars |> 
  dplyr::group_by(group) |> 
  dplyr::reframe(across(c(interventionId_bin, subjectId, Residuals)))

# Checks
table(order_vars$group)
x <- order_vars[order_vars$group=='vascular_and_body_sonography',]
x <- order_vars[order_vars$group=='functional_sports_exam',]

ggplot_labs <- labels
names(ggplot_labs) <- order

# Boxplot
boxplot <- order_vars |>  
  tidyr::pivot_longer(subjectId:Residuals,
                      names_to = 'type',
                      values_to = 'value') |> 
  dplyr::mutate(group = factor(group, levels = order),
                # type = factor(type, levels = c('interventionId_bin', 'Residuals', 'subjectId'))
                ) |> 
  ggplot(aes(x = group,
             y = value)) +
  # geom_violin(aes(fill = type)) +
  geom_boxplot(outlier.shape = NA,
               aes(fill = type)) +
  scale_fill_viridis_d(option = 'G',
                       name = '') +
  theme_bw() +
  labs(x = '') +
  coord_flip() +
  facet_wrap(group ~ 1,
             scales = 'free_y',
             ncol = 1) +
  theme(axis.ticks.y = element_blank(),
        strip.text = element_blank()) +
  labs(y = 'Proportion of variance explained')




# ternary plot
filtered |> 
  dplyr::filter(group == 'blood_haemogram') |> 
  ggtern(aes(x=Residuals,
             y=subjectId,
             z=interventionId)) +
  # geom_point() +
  geom_hex_tern() +
  labs(title="") +
  scale_fill_viridis_c() 