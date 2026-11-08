# Clean work space and change working directory ####
rm(list = ls())

palette = c(
  "blood"  = "#d12950",
  "normal" = "#c39c94",
  "liver_met" = "#f5b6d2",
  "primary"  = "#a4d7dd",
  "carcinosis" = "#a985ab",
  "lymph_met" = "#d379b0",
  "SI-NET" = 'black',
  'Metastasis' = 'grey'
)

palette_subs = c(
  "C>A"  = "#72c7d6",
  "C>G" = "#f07e20",
  "C>T" = "#8fbc8f",
  "T>A"  = "#c8c7c7",
  "T>C" = "#f5adac",
  "T>G" = "#a9a065"
)
source('scripts/analysis/figure2_figure3_read_data.R')

# --- FIGURE 2A (Mutation rate barplot) --- ####
order = mixedsort(burden_all$Sample)
burden_all$Sample = factor(burden_all$Sample, levels = order)

pdf('figures/figure2/figure2_mutationrate.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ggplot(burden_all, aes(x = factor(Clinical_type, levels = c('blood', 'normal', 'primary', 'lymph_met', 'carcinosis', 'liver_met')), y = mutation_rate_scaled, fill = Clinical_type, group = Sample)) +
  geom_bar(stat = 'identity', position = position_dodge2(preserve = 'single', padding = 0), 
           width = 3, 
           aes(group = Sample)) +
  theme_classic() +
  scale_fill_manual(values = palette) +
  geom_errorbar(aes(x = factor(Clinical_type),group = Sample, ymin = mutation_rate_lci_scaled, ymax = mutation_rate_uci_scaled), 
                position = position_dodge2(preserve = 'single'),
                width = 3) +
  geom_text(
    aes(label = Sample),
    y = 6,
    angle = 90,
    hjust = 1,
    position = position_dodge2(preserve = 'single', width = 3, padding = 0)
  ) +
  xlab('') +
  ylab('Mutation rate (subst. per cell per year)') +
  theme(axis.text.x = element_blank())

dev.off()

# --- FIGURE 2A (Regression) --- ####
burden_all = burden_all %>% mutate(group = case_when(Clinical_type %in% c("lymph_met", "carcinosis", 'primary') ~ 'SI-NET',
                                                     Clinical_type == 'blood' ~ 'blood', 
                                                     Clinical_type == 'normal' ~ 'normal')) 

pdf('figures/figure2/figure2_regression.pdf', width = 180*0.0393700787, height = 100*0.0393700787)

burden_mean = burden_all %>% 
  filter(Clinical_type %in% c("lymph_met", "carcinosis", "primary")) %>% 
  group_by(Patient) %>%
  summarise(
    Age = dplyr::first(Age),
    group = dplyr::first(group),
    mean_burden = mean(burden_scaled, na.rm = TRUE),
    SEM = sd(burden_scaled, na.rm = TRUE) / sqrt(sum(!is.na(burden_scaled))),
    .groups = "drop"
  ) 
burden_mean = burden_mean %>% rename(mean_burden = 'burden_scaled')

burden_mean = rbind(burden_mean, burden_all %>% filter(Clinical_type %in% c('blood', 'normal')) %>% select(Patient, Age, group, burden_scaled) %>% mutate(SEM = NA))

ggplot(burden_mean, aes(x = Age, y = burden_scaled, col = group)) +
  geom_point(size = 2) +
  stat_poly_line(method = "lm", se = F, linewidth = 0.5, formula = y ~ x, fullrange = T) +
  xlim(0,85) +
  ylim(0,7000) +
  ggpubr::stat_cor(
    method = "pearson",
    aes(label = paste(after_stat(r.label),
                      after_stat(p.label), sep = "~~~"))
  ) +
  geom_errorbar(
    aes(ymin = burden_scaled - 2*SEM,
        ymax = burden_scaled + 2*SEM),
    width = 0.5
  ) +
  scale_color_manual(values = palette) +
  ylab('Substitutions per cell') +
  xlab('Age (years)') +
  theme_classic() 

dev.off()

# --- FIGURE 2B (Trinucleotide spectrum) --- ####
color = colorRamp2(
  c(0, 0.1),
  c("white", "blue")
)

## --- Plot heatmap --- #
pdf('figures/figure2/figure2_heatmap.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ComplexHeatmap::Heatmap(mat/rowSums(mat), 
                        cluster_columns = F, 
                        cluster_rows = F, 
                        col = color) +
  rowAnnotation(bar = anno_barplot(similarity, bar_width = 0.95), names = profile_all$Clinical_type[1:40])
dev.off()

## --- Plot trinucleotide spectrum clonal vs non-clonal --- #
clonal_nonclonal_profile = clonal_nonclonal_profile %>% mutate(subs = substr(context, start = 3, stop = 5))
order = (clonal_nonclonal_profile %>% arrange(subs, context) %>% select(context))[[1]]
clonal_nonclonal_profile = clonal_nonclonal_profile %>% mutate(context = factor(context, levels = order))

pdf('figures/figure2/figure2_spectrum.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
p1 = ggplot(clonal_nonclonal_profile, aes(x = context, y = clonal, fill = subs)) +
  geom_bar(stat = 'identity') + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  ) +
  ylab('Clonal') +
  xlab('') +
  scale_fill_manual(values = palette_subs) 

p2 = ggplot(clonal_nonclonal_profile, aes(x = context, y = nonclonal, fill = subs)) +
  geom_bar(stat = 'identity') + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  ) +
  ylab('Non-clonal') +
  xlab('') +
  scale_fill_manual(values = palette_subs)
p1 / p2
dev.off()

cosine(clonal_nonclonal_profile$clonal/sum(clonal_nonclonal_profile$clonal), clonal_nonclonal_profile$nonclonal/sum(clonal_nonclonal_profile$nonclonal)) #cosine similarity between profiles
