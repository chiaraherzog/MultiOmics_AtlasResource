#' @name plotMOWAS
#' @description
#' Plot function for mowas
#' @param
#' 
#' 

plotMOWAS <- function(mowas = c('age',
                                  'smoking_curr',
                                  'smoking_ever',
                                  'bmi',
                                  'vo2max',
                                  'menopause',
                                  'activity',
                                  'alcohol',
                                  'alcohol_units'),
                       fdr = 0.05,
                       fdr.method = 'bonferroni'){
  
  
  # Set dirs and load files
  load(paste0('1-analyses/MOWAS/1-out/mowas_nonmethylation_', mowas, '.Rdata'))
  load(paste0('1-analyses/MOWAS/1-out/mowas_methylation_blood_', mowas, '.Rdata'))
  load(paste0('1-analyses/MOWAS/1-out/mowas_methylation_buccal_', mowas, '.Rdata'))
  load(paste0('1-analyses/MOWAS/1-out/mowas_methylation_cervical_', mowas, '.Rdata'))
  
  lm_blood$type = 'blood methylation'
  lm_buccal$type <- 'buccal methylation'
  lm_cervical$type <- 'cervical methylation'
  
  df <- rbind(total, lm_blood, lm_buccal, lm_cervical)
  
  filter_terms <- case_when(
    mowas == 'age' ~ 'age',
    mowas == 'vo2max' ~ 'vo2max',
    mowas == 'menopause' ~ c('mpyes:age','mpyes'),
    mowas == 'bmi' ~ 'bmi',
    mowas == 'activity' ~ 'activityyes',
    mowas == 'alcohol' ~ 'alcoholyes',
    mowas == 'smoking_curr' ~ 'smk_curryes,smokingyes',
    mowas == 'smoking_ever' ~ 'smk_everyes,smokingyes',
    mowas == 'alcohol_units' ~ 'alcohol'
  ) |> 
    strsplit(split = ",") |> 
    unlist()
  
  rm(total, lm_blood, lm_buccal, lm_cervical);gc()
  
  # Get labels
  load(here('data/vars.Rdata'))
  
  # colours
  cols <- c("#1b69a1", "#48a0af", "#71b5a9",  "#f39668", "#ec6669", "#bd647d", "#832c9b", "#5f70a8")
  
  # Filter highly correlated variables
  if(mowas == 'vo2max'){
    df <-  df |> dplyr::filter(!var %in% c('relmaxcap', 'relpower60hf', 'relpower80hf', 'vo2max'))
  }
  
  if(mowas == 'bmi'){
    df <-  df |> dplyr::filter(!var %in% c('weight', 'bmi'))
  }
  
  # Group by type and p adjust
  df <- df |> 
    
    # Filter
    dplyr::filter(term %in% filter_terms) |> 
    dplyr::group_by(type) |> 
    dplyr::mutate(padj = p.adjust(p.value, method = fdr.method)) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(sig = ifelse(padj < fdr, 'sig', 'ns')) |> 
    
    # Append labels
    dplyr::left_join(dplyr::select(vars, x, label, assay, assay2),
                     by = c('var' = 'x'),
                     relationship = 'many-to-many') |> 
    dplyr::filter(!((grepl("cervical_comp", type) & !grepl("cervical", assay))) &
                    !((grepl("buccal_comp", type) & !grepl("buccal", assay))) &
                    !((grepl("blood_comp", type) & !grepl("blood", assay))) &
                    !((type == 'faecal_fam' & !grepl("Stool", assay))) &
                    !((type == 'saliva_microb_fam' & !grepl("Saliva", assay))) &
                    !((type == 'flow_stim' & assay != 'Flow cytometry: stimulated T cells')) &
                    !((type == 'flow_unstim' & assay != 'Flow cytometry: unstimulated T cells')) &
                    !(type == 'faecal' & !grepl("clr", assay)) &
                    !(type == 'saliva_microb' & !grepl("clr", assay))  &
                    !(type == 'saliva_metab' & !grepl("Saliva", assay))  &
                    !(type == 'urine_metab' & !grepl("Urine", assay)) 
                  ) |> 
    dplyr::filter(!assay %in% c("Saliva microbiome: families", "Stool microbiome: families"))  |> 
    dplyr::mutate(assay2 = case_when(grepl('methylation', type) ~ stringr::str_to_sentence(type),
                                     type == 'buccal_comp' ~ "Buccal methylation biomarkers",
                                     type == 'blood_comp' ~ "Blood methylation biomarkers",
                                     type == 'cervical_comp' ~ "Cervical methylation biomarkers",
                                     type == 'faecal' ~ "Stool microbiome ASVs",
                                     type == 'faecal_fam' ~ "Stool microbiome families",
                                     type == 'saliva_microb' ~ "Saliva microbiome ASVs",
                                     type == 'saliva_microb_fam' ~ "Saliva microbiome families",
                                     type %in% c('flow_t', 'flow_wb') ~ "Flow cytometry",
                                     type == 'flow_stim' ~ "Flow cytometry: stimulated T cells",
                                     type == 'flow_unstim' ~ "Flow cytometry: unstimulated T cells",
                                     TRUE ~ assay2))
  
  cat(sum(df$sig=='sig'), 'significant features\n')
  
  
  cols_custom <- colorRampPalette(cols)(15)
  cols_custom <- c(cols_custom[1:7],
                   rep(cols_custom[8], 3),
                   cols_custom[9:10],
                   rep(cols_custom[11], 2),
                   rep(cols_custom[12], 2),
                   cols_custom[c(13, 14, 15, 13, 14, 15)])
  names(cols_custom) <- unique(df$assay2)
  
  assayOrder <- c("Blood methylation",
                  "Buccal methylation",
                  "Cervical methylation",
                  "Stool microbiome ASVs",
                  "Saliva microbiome ASVs",
                  "Flow cytometry",
                  "Stool microbiome families",
                  "Saliva microbiome families",
                  "Saliva metabolome",
                  "Urine metabolome",
                  "Blood test",
                  "Cervical methylation biomarkers",
                  "Blood methylation biomarkers",
                  "Buccal methylation biomarkers",
                  "Blood haemogram",
                  "Flow cytometry: unstimulated T cells",
                  "Flow cytometry: stimulated T cells",
                  "Body weight and composition",
                  "Functional exercise capacity",
                  "Skin histology and TEWL",
                  "Spirometry",
                  "Vascular features")
  
  # Plot 1: Pie chart
  df1 <- df |>
    dplyr::filter(sig == 'sig') |> 
    dplyr::mutate(n = n())  |> 
    dplyr::group_by(assay2) |> 
    dplyr::mutate(ntype = n()) |> 
    dplyr::ungroup() |> 
    dplyr::reframe(type = assay2,
                   prop = ntype/n,
                   text = paste0(type, '\nn=', ntype)) |> 
    dplyr::distinct() |> 
    dplyr::arrange(desc(type)) |> 
    dplyr::mutate(text_y = cumsum(prop) - prop/2) |>
    dplyr::ungroup()
  
  pie <- df1 |>
    ggplot(aes(x = '',
               y = prop,
               fill = type)) +
    geom_col(alpha  = 0.8) + 
    geom_text_repel(aes(label = text,
                        y = text_y,
                        colour = type),
                    nudge_x = 0.8,
                    segment.size = 0.2,
                    segment.curvature = -1e-20,
                    size = 3.2) +
    coord_polar(theta = 'y') +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = cols_custom[match(df1$type, names(cols_custom))],
                      aesthetics = c('colour', 'fill')) +
    labs(subtitle = paste0(sum(df$sig=='sig'), ' significantly associated features'))
  
  # Bar
  bar <- df1 |>
    ggplot(aes(x = '',
               y = prop,
               fill = type)) +
    geom_col(alpha  = 0.8) + 
    geom_text_repel(aes(label = text,
                        y = text_y,
                        colour = type),
                    nudge_x = 0.8,
                    segment.size = 0.2,
                    segment.curvature = -1e-20,
                    size = 3.2) +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = cols_custom[match(df1$type, names(cols_custom))],
                      aesthetics = c('colour', 'fill')) +
    labs(subtitle = paste0(sum(df$sig=='sig'), ' significantly associated features'))
  
  
  # plot 2 (proportion of significant features)
  df2 <- df |> 
    dplyr::group_by(assay2) |> 
    dplyr::reframe(nsig = sum(sig=='sig'),
                   nssig = sum(sig != 'sig'),
                   n = n(),
                   prop = nsig/n(),
                   prop_ns = nssig/n()) |> 
    dplyr::ungroup() |> 
    tidyr::pivot_longer(c(prop, prop_ns),
                        names_to = 'prop',
                        values_to = 'value')
  
 prop <- df2 |> 
   dplyr::mutate(assay2 = factor(assay2, levels = rev(assayOrder))) |> 
    ggplot(aes(y = assay2,
               x = value,
               fill = prop,
               text = paste0("signif prop = ", signif(value, 3)))) +
    geom_col(position = 'fill') +
   labs(x = '',
        y = '',
        subtitle = 'proportion of total features',
        title = mowas) +
   theme_bw() +
   theme(axis.ticks.y = element_blank(),
         axis.text.x = element_text(size = 8),
         legend.position = 'bottom',
         panel.grid = element_blank()) +
   scale_x_continuous(breaks = c(0, 0.5, 1)) +
   coord_cartesian(expand = 0) +
   scale_fill_manual(values = c(cols[1], '#ededed'),
                     name = '',
                     labels = c("associated", 'not associated'))
  
 
 
 # Total features
 df_total <- df |>
   dplyr::group_by(assay2) |> 
   dplyr::mutate(ntype = n(),
                 n = sum(sig == 'sig')) |> 
   dplyr::ungroup() |> 
   dplyr::reframe(assay2 = assay2,
                  n = n,
                  n_total = ntype) |> 
   dplyr::distinct()
 
 column <- df_total |> 
   dplyr::mutate(assay2 = factor(assay2, levels = rev(assayOrder))) |> 
   ggplot(aes(y = assay2,
              x = n)) +
   geom_col(aes(fill = assay2)) +
   scale_fill_manual(values = cols_custom[match(df_total$assay2, names(cols_custom))],
                   aesthetics = c('colour', 'fill')) +
   theme_bw() +
   theme(legend.position = 'none',
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         panel.grid = element_blank()) +
   labs(x = 'n', y = '') +
   coord_cartesian(expand = F)
 
  # Manhattan-style plot (excluding methylation)
  df3 <- df |> 
    dplyr::filter(!grepl('methylation', type) & sig == 'sig') |> 
    dplyr::group_by(assay2) |> 
    dplyr::arrange(var) |> 
    dplyr::mutate(id = sequence(n())) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(v = paste0(type, var, id),
                  t = case_when(sig == 'sig' & !grepl('methylation', type) ~ label, 
                                TRUE ~ NA),
                  value = -log10(p.value)*sign(estimate))
  
   manh <- df3 |> 
    dplyr::arrange(v) |> 
    ggplot(aes(x = 1,
               y = value,
               colour = assay2,
               size = -log10(p.value),
               alpha = -log10(p.value),
               text = paste0("assay = ", assay2, "\n",
                              "feature = ", t,
                             "\np = ", p.value))) +
     geom_hline(yintercept = 0,
                colour = 'grey50',
                linetype = 'dashed') +
     geom_point() +
    scale_size(range = c(0.5, 2),
               guide = F) +
    geom_text_repel(aes(label = t),
                    size = 2.8,
                    max.overlaps = 18,
                    segment.curvature = -1e-20,
                    min.segment.length = 0.2) + 
    scale_alpha(limits = c(0.6, 1),
                guide = F) +
     scale_colour_manual(values = cols_custom[match(df3$assay2, names(cols_custom))],
                         name = '',
                         guide = guide_legend(override.aes = list(label = ""))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none') +
     labs(x = '',
          y = '-log10(p) * sign(estimate)') 
  
  
   return(list(pie = pie,
               prop = prop,
               manh = manh,
               column = column,
               bar = bar,
               cols = cols_custom,
               df = df))
  
}
