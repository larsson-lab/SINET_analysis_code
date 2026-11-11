source('scripts/analysis/figure2_figure3_read_data.R')

library(patchwork) #v1.3.2
library(dplyr) #v1.2.1
library(ggplot2) #v4.0.3
library(ggpubr) #v1.0.0
library(ggpmisc) #v1.0.0

p1 = ggplot(
  burden_all %>% filter(Clinical_type == 'lymph_met'),
  aes(x = Purity, y = mutation_rate)
) +
  geom_point(
    size = 2, col = "darkblue"
  ) +
  xlim(0,1) +
  ylim(0,60) +
  stat_poly_line(
    method = "lm",
    se = FALSE,
    linewidth = 0.5,
    formula = y ~ x,
    fullrange = TRUE
  ) +
  ggpubr::stat_cor(
    method = "pearson",
    aes(label = paste(
      after_stat(r.label),
      after_stat(p.label),
      sep = "~~~"
    )),
    label.y.npc = 0.95
  ) +
  ggpmisc::stat_poly_eq(
    formula = y ~ x,
    aes(label = after_stat(eq.label)),
    parse = TRUE,
    label.y = "bottom"
  ) +
  ylab("NanoSeq rate (SNVs/year)") +
  xlab("Purity") +
  ggtitle('Lymph node metastases') +
  theme_classic() +
  theme(text = element_text(size = 6)) 

plot_data = burden_all %>% 
  filter(Clinical_type %in% c("primary", "lymph_met", "carcinosis", "liver_met")) %>%
  select(Sample, Clinical_type, mutation_rate, mutation_rate_scaled) %>%
  mutate(
    Clinical_type = factor(
      Clinical_type,
      levels = c("primary", "lymph_met", "carcinosis", "liver_met")
    )
  ) %>%
  pivot_longer(
    cols = c(mutation_rate, mutation_rate_scaled),
    names_to = "Rate",
    values_to = "Value"
  ) 


p2 = ggplot(plot_data, aes(x = Sample, y = Value, fill = Rate)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_grid(
    ~ Clinical_type,
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  scale_fill_manual(
    values = c(
      "mutation_rate" = "darkblue",
      "mutation_rate_scaled" = "lightblue"
    ),
    breaks = c("mutation_rate", "mutation_rate_scaled"),
    labels = c("Uncorrected", "Corrected")
  ) +
  labs(
    x = "",
    y = "NanoSeq rate (SNVs/year)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(angle = 90),
    legend.position = "bottom",
    text = element_text(size = 6)
  )


pdf('figures/supp_figures/Supplementary_Figure_4.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
p1 | p2
dev.off()
