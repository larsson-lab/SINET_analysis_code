# Clean work space and change working directory ####
rm(list = ls())
source('scripts/analysis/figure2_figure3_read_data.R')

# Define palettes and read data ####
palette = c(
  "blood"  = "#d12950",
  "normal" = "#c39c94",
  "liver_met" = "#f5b6d2",
  "primary"  = "#a4d7dd",
  "carcinosis" = "#a985ab",
  "lymph_met" = "#d379b0"
)

# FIGURE 3B (clonal expansion NanoSeq cohort) ####
fractions_scaled = left_join(fractions_scaled, mutation_count %>% select(Sample, recall), by = 'Sample') %>% 
  mutate(fraction_clonal_scaled_correction = fraction_clonal_scaled/recall)

## --- Compare ages ####
age_comparison = inner_join(fractions_scaled, mutation_count, by = c('Sample', 'Patient', 'Tumor', 'Clinical_type', 'Library', 'Age'))
age_comparison = age_comparison %>% select(Sample, fraction_clonal_scaled_correction, Age, Clinical_type, count_corrected, mean, ce_mean)
age_comparison = inner_join(age_comparison, burden_all, by = c('Age', 'Sample', 'Clinical_type'))
age_comparison = age_comparison %>% select(Sample, fraction_clonal_scaled_correction, Age, Clinical_type, count_corrected, mean, ce_mean, mutation_rate_scaled)
age_comparison = age_comparison %>% mutate(label = paste(Sample,round(age_comparison$mutation_rate_scaled), sep = ','))
age_comparison = inner_join(age_comparison, purity %>% select(Purity, Sample), by = c('Sample'))

age_comparison <- age_comparison %>%
  mutate(
    ymin = pmin(count_corrected/mutation_rate_scaled, ce_mean, na.rm = TRUE),
    ymax = pmax(count_corrected/mutation_rate_scaled, ce_mean, na.rm = TRUE)
  )

pdf('figures/supp_figures/age_comparison.pdf', width = 180*0.0393700787, height = 100*0.0393700787)

ggplot(age_comparison, aes(label = Sample)) +
  geom_segment(aes(x = fraction_clonal_scaled_correction*Age, xend = fraction_clonal_scaled_correction*Age, y = ymin, yend = ymax))+
  geom_point(aes(x = fraction_clonal_scaled_correction*Age, y = count_corrected/mutation_rate_scaled), size = 4, shape = 21, colour = "black", stroke = 0.8, fill = 'white') +
  geom_point(aes(x = fraction_clonal_scaled_correction*Age, y = ce_mean, col = Clinical_type), size = 4) +
  theme_bw() +
  scale_colour_manual(values = palette) +
  geom_abline(slope = 1, intercept = 0) +
  theme(legend.position = 'none') +
  xlim(0,130) +
  ylim(0,130) +
  ylab('Estimated age at expansion (rate model)') +
  xlab("Estimated age at expansion (fraction clonal in NanoSeq)") +
  coord_fixed()
dev.off()
