load('1-analyses/malleability/1-out/ageOut.Rdata')

library(ggplot2)

tab_i |> 
  tidyr::pivot_longer(c(ratio_hi, ratio_lo), names_to = 'ratio', values_to = 'value') |> 
  ggplot(aes(x = value,
             y = forcats::fct_reorder(rowname, value),
             colour = ratio)) +
  geom_line(aes(group = rowname),
            colour = 'grey') +
  geom_vline(xintercept = 1) +
  geom_point() +
  theme_bw() +
  labs(x = 'observed/expected\ndelta (M6)')


tab_s |> 
  tidyr::pivot_longer(c(ratio_hi, ratio_lo), names_to = 'ratio', values_to = 'value') |> 
  ggplot(aes(x = value,
             y = forcats::fct_reorder(rowname, value),
             colour = ratio)) +
  geom_line(aes(group = rowname),
            colour = 'grey') +
  geom_vline(xintercept = 1) +
  geom_point() +
  theme_bw() +
  labs(x = 'observed/expected\ndelta (M6)')


library(ggplot2)
library(dplyr)
p_overlap_bar <- overlap_inf %>%
  count(altered_by) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = altered_by, y = prop, fill = altered_by)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(prop, 1)), vjust = -0.3) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of age-associated features",
       title = "Features attenuated by intervention and/or cessation") +
  theme_bw() +
  theme(legend.position = "none")

p_cut = 0.05
p_effect_map <- overlap_inf %>%
  ggplot(aes(x = P_i, y = P_s)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(alpha = (P_i < p_cut) | (P_s < p_cut),
                 shape = altered_by), size = 1.8) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.3), guide = "none") +
  labs(x = "Compliance effect on aligned Δ (I arm)\n<0 = attenuation",
       y = "Compliance effect on aligned Δ (S arm)\n<0 = attenuation",
       shape = "Altered by",
       title = "Feature-level map of intervention vs cessation effects") +
  theme_bw()


eps = 0.1
# Mark attenuation at high compliance by ratio criterion
rat_i <- summary_m6_i %>% mutate(attn_I = ratio_hi <= (1 - eps))
rat_s <- summary_m6_s %>% mutate(attn_S = ratio_hi <= (1 - eps))

overlap_ratio <- rat_i %>%
  select(assay, rowname, attn_I) %>%
  full_join(rat_s %>% select(assay, rowname, attn_S), by = c("assay","rowname")) %>%
  mutate(group = case_when(
    attn_I & attn_S ~ "Both",
    attn_I & !attn_S ~ "I only",
    !attn_I & attn_S ~ "S only",
    TRUE ~ "Neither"
  ))

ggplot(overlap_ratio %>% count(group),
       aes(x = 1, y = n, fill = group)) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "Features", title = "Attenuated by I vs S (ratio-based)") +
  theme_bw() + theme(legend.position = "none")



ggplot(overlap_ratio %>% count(group), aes(x = group, y = n, fill = group)) +
  geom_col(width = 0.7) +
  labs(x = NULL, y = "Features", title = "Attenuated by I vs S (ratio-based)") +
  theme_bw() + theme(legend.position = "none")

prop_hi <- tab_i |> 
  filter(!is.na(category_hi)) |> 
  count(category_hi, name = "n") |> 
  mutate(prop = n / sum(n)) |> 
  ungroup()



net_index <- tab_i |> 
  filter(!is.na(category_hi)) |> 
  mutate(val = case_when(
    category_hi == "age-opposing" ~  1,
    category_hi == "accelerated"  ~ -1,
    TRUE                          ~  0
  )) |> 
  summarise(net_shift = mean(val), .groups="drop")


net_by_comp <- tab_i |> 
  select(category_hi, category_lo) |> 
  tidyr::pivot_longer(cols = starts_with("category_"),
               names_to = "comp_group", values_to = "category") |> 
  mutate(
    comp_group = ifelse(comp_group == "category_hi", "High", "Low"),
    val = case_when(
      category == "age-opposing" ~  -1,
      category == "accelerated"  ~ 1,
      TRUE                       ~  0
    )
  ) |> 
  group_by(comp_group) |> 
  summarise(net_shift = mean(val, na.rm = TRUE), .groups = "drop")




net_by_comp

plot_props <- prop_hi |> 
  ggplot(aes(x = 1, y = prop, fill = category_hi)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of ageing-associated features",
       fill = "Category",
       title = "High-compliance shift vs expected ageing") +
  theme_bw() +
  theme(legend.position = "right")
plot_props

plot_props + facet_wrap(~ assay, nrow = 1)  # if prop_hi includes an assay column; if not, join it first


plot_net <- net_index |> 
  ggplot(aes(x = 1, y = net_shift)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_segment(aes(xend = 1, y = 0, yend = net_shift), linewidth = 0.6) +
  geom_point(size = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = NULL, y = "Net shift index (opposing − accelerated)",
       title = "Net direction of change among high-compliance features") +
  theme_bw()
plot_net



scatter_df <- summary_m6_i |>
  mutate(
    obs_hi = mean_aligned_hi,
    exp    = exp_mag,  # on aligned magnitude scale
    status = case_when(
      category_hi == "age-opposing" ~ "Opposing (reversal)",
      category_hi == "attenuated"   ~ "Attenuated",
      category_hi == "unchanged"    ~ "Unchanged",
      category_hi == "accelerated"  ~ "Accelerated",
      TRUE                          ~ NA_character_
    )
  )

p_scatter <- ggplot(scatter_df, aes(x = exp, y = obs_hi)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +   # equality line (no attenuation)
  geom_point(alpha = 0.6, size = 1.5, aes(shape = status)) +
  # facet_wrap(~ interventionId) +
  scale_x_continuous(name = "Expected Δ from baseline ageing (|β_age| × Δt)") +
  scale_y_continuous(name = "Observed aligned Δ in high compliance") +
  labs(title = "Expected vs observed change: attenuation below the diagonal") +
  theme_bw() +
  theme(legend.position = "right")
p_scatter



volc <- tab_i |> 
  mutate(
    FDR = p.adjust(P.Value, method = "fdr"),
    sig = FDR < 0.05
  ) |> 
  ggplot(aes(x = logFC, y = -log10(P.Value))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(alpha = sig)) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.3), guide = "none") +
  labs(x = "Effect of compliance on aligned Δ (logFC)\n<0 = more 'younger' with higher compliance",
       y = "-log10(p)",
       title = "Feature-level association with compliance (limma)") +
  theme_bw()
volc


cat_comp <- summary_m6_i |> 
  select(category_hi, category_lo) |> 
  tidyr::pivot_longer(cols = starts_with("category_"),
                      names_to = "comp_group", values_to = "category") %>%
  mutate(comp_group = ifelse(comp_group == "category_hi", "High", "Low"),
         category = factor(category)) %>%
  count(comp_group, category, name = "n") %>%
  group_by(comp_group) %>%
  mutate(prop = n/sum(n)) %>% ungroup()

ggplot(cat_comp, aes(x = comp_group, y = prop, fill = category)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0,0.02))) +
  labs(x = NULL, y = "Proportion of features", fill = "Category",
       title = "High vs Low compliance: category composition") +
  theme_bw()


cat_comp <- summary_m6_s |> 
  select(category_hi, category_lo) |> 
  tidyr::pivot_longer(cols = starts_with("category_"),
                      names_to = "comp_group", values_to = "category") %>%
  mutate(comp_group = ifelse(comp_group == "category_hi", "High", "Low"),
         category = factor(category)) %>%
  count(comp_group, category, name = "n") %>%
  group_by(comp_group) %>%
  mutate(prop = n/sum(n)) %>% ungroup()

ggplot(cat_comp, aes(x = comp_group, y = prop, fill = category)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(c(0,0.02))) +
  labs(x = NULL, y = "Proportion of features", fill = "Category",
       title = "High vs Low compliance: category composition") +
  theme_bw()

