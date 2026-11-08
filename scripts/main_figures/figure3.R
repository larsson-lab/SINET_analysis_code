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
  "lymph_met" = "#d379b0",
  "SI-NET" = 'black',
  'Metastasis' = 'grey'
)

# FIGURE 3B (clonal expansion NanoSeq cohort) ####
fractions_scaled = left_join(fractions_scaled, mutation_count %>% select(Sample, recall), by = 'Sample') %>% 
  mutate(fraction_clonal_scaled_correction = fraction_clonal_scaled/recall)

fractions_plot = fractions_scaled %>%  mutate(Clinical_type = factor(Clinical_type, levels = c('normal', 'primary', 'lymph_met', 'carcinosis', 'liver_met')))
order = (fractions_plot %>%  arrange(Clinical_type, (fraction_clonal_scaled_correction)*Age) %>% select(Sample))[[1]]
fractions_plot$Sample = factor(fractions_plot$Sample, levels = unique(order))
fractions_plot = fractions_plot %>% mutate(fraction_clonal_scaled_correction = if_else(is.na(fraction_clonal_scaled_correction), true = 0, false = fraction_clonal_scaled_correction))

pdf('figures/figure3/figure3_clonal_fraction.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ggplot(fractions_plot) +
  geom_col(aes(x = as.factor(Sample), y = Age), fill = '#e1a248') +
  geom_col(aes(x = as.factor(Sample), y = (fraction_clonal_scaled_correction*Age), fill = Clinical_type)) +
  scale_fill_manual(values = palette) +
  theme_classic() +
  ylab('Patient age (years)') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

# FIGURE 3C (clonal expansion WGS cohort) ####
order = (mutation_count %>%  
           mutate(Clinical_type = factor (Clinical_type, levels = c('primary', 'lymph_met', 'carcinosis', 'liver_met'))) %>%
           arrange(Clinical_type, count_corrected) %>% 
           select(Sample))[[1]]
mutation_count$Sample = factor(mutation_count$Sample, levels = order)

pdf('figures/figure3/figure3_clonal_expansion.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ggplot(mutation_count %>% filter(Purity > 0.4) %>% filter(!Sample %in% fractions_plot$Sample), aes(x = factor(Sample))) +
  geom_col(
    aes(y = Age, x = factor(Sample)),
    position = position_dodge(preserve = "single", width = 0.8),
    fill = '#e1a248'
  ) +
  geom_col(
    aes(y = ce_mean, x = factor(Sample), fill = Clinical_type),
    position = position_dodge(preserve = "single", width = 0.8)
  ) +
  scale_fill_manual(values = palette) +
  geom_errorbar(aes(x = factor(Sample), ymin = ce_uci_90, ymax = ce_lci_90), position = position_dodge(preserve = 'single', width = 0.8)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  theme(legend.position = 'none') +
  ylab('Patient age (years)') +
  xlab('')
dev.off()

# FIGURE 3D (earliest clonal expansion) ####
mutation_count = mutation_count %>%  mutate(years_before = Age - ce_mean) %>% 
  mutate(years_before_uci90 = Age - ce_uci_90) %>% 
  mutate(years_before_lci90 = Age - ce_lci_90) 

# Change values for NanoSeq samples
for (sample in mutation_count$Sample){
  if (sample %in% fractions_plot$Sample){
    mutation_count[mutation_count$Sample==sample,'years_before'] = mutation_count[mutation_count$Sample==sample,'Age'] - (fractions_plot[fractions_plot$Sample == sample, "fraction_clonal_scaled_correction"] * fractions_plot[fractions_plot$Sample == sample, "Age"])
  }
}

mutation_count = mutation_count %>% mutate(NanoSeq = if_else(Sample %in% fractions_plot$Sample, true = 1, false = 0))

years_before = mutation_count %>%
  filter(Purity > 0.4) %>% 
  group_by(Patient) %>%
  slice_max(years_before, n = 1, with_ties = F) %>% 
  ungroup()
order = (years_before %>%  
           arrange(Patient, Tumor) %>% 
           select(Sample))[[1]]
years_before$Sample = factor(years_before$Sample, levels = order)

pdf('figures/figure3/figure3_earliest_expansion.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ggplot(years_before, aes(x = factor(Sample))) +
  geom_point(
    aes(y = -years_before, x = factor(Sample), col = Clinical_type),
    position = position_dodge(preserve = "single", width = 0.8),
    size = 4
  ) +
  geom_segment(aes(y=0, yend=-Age), linewidth = 0.5, colour = 'grey') +
  scale_colour_manual(values = palette) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  geom_errorbar(data = years_before %>% filter(NanoSeq == 0), aes(x = factor(Sample), ymin = -years_before_lci90, ymax = -years_before_uci90), position = position_dodge(preserve = 'single', width = 0.8)) +
  geom_point(
    aes(y = -Age, x = factor(Sample)),
    position = position_dodge(preserve = "single", width = 0.8),
    size = 2, 
    shape = 22,
    fill = "grey",
    colour = "grey"
  ) +
  ylab('Years before surgery') +
  xlab('') +
  ggtitle('Earliest clonal expansion per patient')
dev.off()

# FIGURE 3E (Phylogenetic trees) ####
# Load data
shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')

## Calculate trunk lengths (based on shared SNVs) and plot phylogenetic trees ####
mutation_len = 0.00095
mutations_yearly = mutation_len*mean 
mutation_len_uci = (mean/(mean+1.64*sd(x)))*mutation_len 
mutation_len_lci = (mean/(mean-1.64*sd(x)))*mutation_len 

### Patient 2, 3, 4, 5, 6, 7, 10, 11, 15, 16 ####
timeline = function(patient, primary, met){
  path = paste('figures/figure3/phylogenies/timeline_', primary, '_', met, '.pdf', sep = "")
  width = 8.27 
  height = 11.69/3 
  # Set margins
  left_margin = 1  
  right_margin = 1  
  correction_primary = mutation_count[mutation_count$Sample == primary, 'correction']
  correction_metastasis = mutation_count[mutation_count$Sample == met, 'correction']
  
  # Calculate plot width and height after considering margins
  plot_width = width - left_margin - right_margin
  plot_height = height 
  pdf(path, width = plot_width, height = plot_height)
  
  if (primary %in% WGD$primary){
    shared_tot = shared_mutations_SNV[[patient]][primary, met]
    shared = ((shared_tot - (WGD[WGD$primary == primary, 'before'])[[1]])/2) + (WGD[WGD$primary == primary, 'before'])[[1]]
    primary_mut = (shared_mutations_SNV[[patient]][primary, primary]-shared_tot)/2
    met_mut = (shared_mutations_SNV[[patient]][met, met]-shared_tot)/2
    
    primary_mut = ((primary_mut+shared)*correction_primary)-shared*correction_primary*correction_metastasis
    met_mut = ((met_mut+shared)*correction_metastasis)-shared*correction_primary*correction_metastasis
    shared = shared*correction_primary*correction_metastasis

  }else{
    shared = shared_mutations_SNV[[patient]][primary, met]
    shared = shared*correction_primary*correction_metastasis
    primary_mut = (shared_mutations_SNV[[patient]][primary, primary]*correction_primary - shared)
    met_mut = (shared_mutations_SNV[[patient]][met, met]*correction_metastasis - shared)

  }
  
  mutation_len = 0.00095 #Mutation size is set to ~0.15 pt (0.002 inches)
  vertical = 1/16
  mutations_yearly = mutation_len*mean
  
  plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), asp = 1, axes = FALSE, xlab = "", ylab = "")
  par(lwd =0.5)
  segments(0, 1, shared*mutation_len, 1)
  segments(shared*mutation_len, 1-vertical, shared*mutation_len, 1+vertical)
  segments(shared*mutation_len, 1-vertical, (shared+met_mut)*mutation_len, 1-vertical)
  segments(shared*mutation_len, 1+vertical, (shared+primary_mut)*mutation_len, 1+vertical)
  #lci
  segments(shared*mutation_len_lci, 1-vertical, shared*mutation_len_lci, 1+vertical, col = 'blue')
  segments((shared+met_mut)*mutation_len_lci, 1-vertical, (shared+met_mut)*mutation_len_lci, 1-2*vertical, col = 'red')
  segments((shared+primary_mut)*mutation_len_lci, 1+vertical, (shared+primary_mut)*mutation_len_lci, 1+2*vertical, col = 'green')
  #uci
  segments(shared*mutation_len_uci, 1-vertical, shared*mutation_len_uci, 1+vertical, col = 'blue')
  segments((shared+met_mut)*mutation_len_uci, 1-vertical, (shared+met_mut)*mutation_len_uci, 1-2*vertical, col = 'red')
  segments((shared+primary_mut)*mutation_len_uci, 1+vertical, (shared+primary_mut)*mutation_len_uci, 1+2*vertical, col = 'green')
  
  text((shared+met_mut)*mutation_len + 0.2, 1-vertical, met, cex = 0.5)
  text((shared+primary_mut)*mutation_len + 0.2, 1+vertical, primary, cex = 0.5)
  segments(metadata_WGS$Age[metadata_WGS$Patient==patient][1]*mutations_yearly, 0.5, metadata_WGS$Age[metadata_WGS$Patient==patient][1]*mutations_yearly, 1.5)
  dev.off()
}

create_timeline = function(patient, primary){
  for (i in 1:ncol(shared_mutations_SNV[[patient]])){
    if (shared_mutations_SNV[[patient]][primary, i] > 10 & colnames(shared_mutations_SNV[[patient]])[i] != primary){
      met = colnames(shared_mutations_SNV[[patient]])[i]
      timeline(patient, primary, met)
    }
  }
}

for (patient in c(2,3,4,5,6,7,10,11,15,16)){
  primaries = (metadata_WGS%>%  filter(Patient == patient & Clinical_type == 'primary') %>% select(Sample))[[1]]
  for (primary in primaries){
    create_timeline(patient, primary) 
  }
}

### Patient 12, 13, 14 ####
create_timeline2 = function(patient){
  names = colnames(shared_mutations_SNV[[patient]])[1:3]
  
  correction_primary = mutation_count[mutation_count$Sample == names[1], 'correction']
  correction_metastasis1 = mutation_count[mutation_count$Sample == names[2], 'correction']
  correction_metastasis2 = mutation_count[mutation_count$Sample == names[3], 'correction']
  
  primaries = (metadata_WGS%>%  filter(Patient == patient & Clinical_type == 'primary') %>% select(Sample))[[1]]
  mrca = (metadata_WGS%>%  filter(Patient == patient & Clinical_type == 'primary') %>% select(Shared))[[1]]
  mrca = mrca*correction_primary*correction_metastasis1*correction_metastasis2
  met_shared = (shared_mutations_SNV[[patient]][2,3]*correction_metastasis1*correction_metastasis2) - mrca
  primary_mut = (shared_mutations_SNV[[patient]][1,1]*correction_primary - mrca)
  met1_mut = (shared_mutations_SNV[[patient]][2,2]*correction_metastasis1 - mrca - met_shared)
  met2_mut = (shared_mutations_SNV[[patient]][3,3]*correction_metastasis2 - mrca - met_shared)

  timeline2(patient, primaries, mrca, met_shared, primary_mut, met1_mut, met2_mut, names)
}

timeline2 = function(patient, primaries, mrca, met_shared, primary_mut, met1_mut, met2_mut, names){
  path = paste('figures/figure3/phylogenies/timeline_', primaries, '.pdf', sep = "")
  width = 8.27  
  height = 11.69/3 
  # Set margins
  left_margin = 1  
  right_margin = 1  
  
  # Calculate plot width and height after considering margins
  plot_width = width - left_margin - right_margin
  plot_height = height 
  pdf(path, width = plot_width, height = plot_height)
  
  mutation_len = 0.00095 #Mutation size is set to ~0.15 pt (0.002 inches)
  vertical = 1/16 # 1/16 inch ~5pt
  mutations_yearly = mutation_len*mean
  
  plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), asp = 1, axes = FALSE, xlab = "", ylab = "")
  par(lwd = 0.5)
  segments(0, 1, mrca*mutation_len, 1) #mrca trunk
  segments(mrca*mutation_len, 1-vertical, mrca*mutation_len, 1+vertical)
  segments(mrca*mutation_len, 1-vertical, (mrca+met_shared)*mutation_len, 1-vertical) #metastases shared
  segments(mrca*mutation_len, 1+vertical, (mrca+primary_mut)*mutation_len, 1+vertical) #primary private
  #lci
  segments(mrca*mutation_len_lci, 1-vertical, mrca*mutation_len_lci, 1+vertical, col = 'blue')
  segments((mrca+met_shared)*mutation_len_lci, 1, (mrca+met_shared)*mutation_len_lci, 1-2*vertical, col = 'red') #metastases shared
  segments((mrca+primary_mut)*mutation_len_lci, 1, (mrca+primary_mut)*mutation_len_lci, 1+2*vertical, col = 'green') #primary private
  #uci
  segments(mrca*mutation_len_uci, 1-vertical, mrca*mutation_len_uci, 1+vertical, col = 'blue')
  segments((mrca+met_shared)*mutation_len_uci, 1, (mrca+met_shared)*mutation_len_uci, 1-2*vertical, col = 'red') #metastases shared
  segments((mrca+primary_mut)*mutation_len_uci, 1, (mrca+primary_mut)*mutation_len_uci, 1+2*vertical, col = 'green') #primary private
  
  text((mrca+primary_mut)*mutation_len + 0.2, 1+vertical, names[1], cex = 0.5) #primary label
  segments((mrca+met_shared)*mutation_len, 1-2*vertical, (mrca+met_shared)*mutation_len, 1)
  
  segments((mrca+met_shared)*mutation_len, 1, (mrca+met_shared+met1_mut)*mutation_len, 1)
  segments((mrca+met_shared)*mutation_len, 1-2*vertical, (mrca+met_shared+met2_mut)*mutation_len, 1-2*vertical)
  
  #lci and uci
  segments((mrca+met_shared+met1_mut)*mutation_len_lci, 1+0.75*vertical, (mrca+met_shared+met1_mut)*mutation_len_lci, 1-1.25*vertical, col = 'pink')
  segments((mrca+met_shared+met2_mut)*mutation_len_uci, 1-1.25*vertical, (mrca+met_shared+met2_mut)*mutation_len_uci, 1-3.25*vertical, col = 'yellow')
  
  segments((mrca+met_shared+met1_mut)*mutation_len_uci, 1+0.75*vertical, (mrca+met_shared+met1_mut)*mutation_len_uci, 1-1.25*vertical, col = 'pink')
  segments((mrca+met_shared+met2_mut)*mutation_len_lci, 1-1.25*vertical, (mrca+met_shared+met2_mut)*mutation_len_lci, 1-3.25*vertical, col = 'yellow')
  
  text((mrca+met_shared+met1_mut)*mutation_len + 0.2, 1, names[2], cex = 0.5) #met1 label
  text((mrca+met_shared+met2_mut)*mutation_len + 0.2, 1-2*vertical, names[3], cex = 0.5) #met2 label
  segments(metadata_WGS$Age[metadata_WGS$Patient==patient][1]*mutations_yearly, 0.5, metadata_WGS$Age[metadata_WGS$Patient==patient][1]*mutations_yearly, 1.5)
  dev.off()
}

for (patient in c(12,13,14)){
  create_timeline2(patient)
}

#Patient 8
primaries = (metadata_WGS%>%  filter(Patient == 8 & Clinical_type == 'primary') %>% select(Sample))[[1]]
correction_8A = mutation_count[mutation_count$Sample == '8A', 'correction']
correction_8C = mutation_count[mutation_count$Sample == '8C', 'correction']
correction_8D = mutation_count[mutation_count$Sample == '8D', 'correction']

mrca = (metadata_WGS%>%  filter(Patient == 8 & Sample == '8A') %>% select(Shared))[[1]]
mrca = mrca*correction_8A*correction_8C*correction_8D
shared_8A_8C= shared_mutations_SNV[[8]][1,3]*correction_8A*correction_8C - mrca
private_8D = (shared_mutations_SNV[[8]][4,4]*correction_8D - mrca)
private_8A = (shared_mutations_SNV[[8]][1,1]*correction_8A - mrca - shared_8A_8C)
private_8C = (shared_mutations_SNV[[8]][3,3]*correction_8C - mrca - shared_8A_8C)

path = 'figures/figure3/phylogenies/timeline_8A.pdf'
width = 8.27  
height = 11.69/3 
# Set margins
left_margin = 1  
right_margin = 1  

# Calculate plot width and height after considering margins
plot_width = width - left_margin - right_margin
plot_height = height 
pdf(path, width = plot_width, height = plot_height)

mutation_len = 0.00095 #Mutation size is set to ~0.15 pt (0.002 inches)
vertical = 1/16 # 1/16 inch 5pt
mutations_yearly = mutation_len*mean

plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), asp = 1, axes = FALSE, xlab = "", ylab = "")
par(lwd = 0.5)
segments(0, 1, mrca*mutation_len, 1) #mrca trunk
segments(mrca*mutation_len, 1-vertical, mrca*mutation_len, 1+vertical)
segments(mrca*mutation_len, 1-vertical, (mrca+private_8D)*mutation_len, 1-vertical) #8D
segments(mrca*mutation_len, 1+vertical, (mrca+shared_8A_8C)*mutation_len, 1+vertical) #shared 8A, 8C
#lci
segments(mrca*mutation_len_lci, 1-vertical, mrca*mutation_len_lci, 1+vertical, col = 'blue')
segments((mrca+private_8D)*mutation_len_lci, 1-1.5*vertical, (mrca+private_8D)*mutation_len_lci, 1-0.5*vertical, col = 'red') #8D
segments((mrca+shared_8A_8C)*mutation_len_lci, 1+1.5*vertical, (mrca+shared_8A_8C)*mutation_len_lci, 1+0.5*vertical, col = 'green') #shared 8A, 8C
#uci
segments(mrca*mutation_len_uci, 1-vertical, mrca*mutation_len_uci, 1+vertical, col = 'blue')
segments((mrca+private_8D)*mutation_len_uci, 1-1.5*vertical, (mrca+private_8D)*mutation_len_uci, 1-0.5*vertical, col = 'red') #8D
segments((mrca+shared_8A_8C)*mutation_len_uci, 1+1.5*vertical, (mrca+shared_8A_8C)*mutation_len_uci, 1+0.5*vertical, col = 'green') #shared 8A, 8C

text((mrca+private_8D)*mutation_len + 0.2, 1-vertical, '8D', cex = 0.5) #8D label
segments((mrca+shared_8A_8C)*mutation_len, 1, (mrca+shared_8A_8C)*mutation_len, 1+2*vertical)

segments((mrca+shared_8A_8C)*mutation_len, 1+2*vertical, (mrca+shared_8A_8C+private_8A)*mutation_len, 1+2*vertical)
segments((mrca+shared_8A_8C)*mutation_len, 1, (mrca+shared_8A_8C+private_8C)*mutation_len, 1)

#lci and uci
segments((mrca+shared_8A_8C+private_8A)*mutation_len_lci, 1+2.5*vertical, (mrca+shared_8A_8C+private_8A)*mutation_len_lci, 1+1.5*vertical, col = 'pink')
segments((mrca+shared_8A_8C+private_8C)*mutation_len_lci, 1+0.5*vertical, (mrca+shared_8A_8C+private_8C)*mutation_len_lci, 1-0.5*vertical, col = 'yellow')

segments((mrca+shared_8A_8C+private_8A)*mutation_len_uci, 1+2.5*vertical, (mrca+shared_8A_8C+private_8A)*mutation_len_uci, 1+1.5*vertical, col = 'pink')
segments((mrca+shared_8A_8C+private_8C)*mutation_len_uci, 1+0.5*vertical, (mrca+shared_8A_8C+private_8C)*mutation_len_uci, 1-0.5*vertical, col = 'yellow')

text((mrca+shared_8A_8C+private_8A)*mutation_len + 0.2, 1+2*vertical, '8A', cex = 0.5) #met1 label
text((mrca+shared_8A_8C+private_8C)*mutation_len + 0.2, 1, '8C', cex = 0.5) #met2 label
segments(metadata_WGS$Age[metadata_WGS$Patient==8][1]*mutations_yearly, 0.5, metadata_WGS$Age[metadata_WGS$Patient==8][1]*mutations_yearly, 1.5)
dev.off()

### Patient 1 ####
correction_1C = mutation_count[mutation_count$Sample == '1C', 'correction']
correction_1G = mutation_count[mutation_count$Sample == '1G', 'correction']
correction_1H = mutation_count[mutation_count$Sample == '1H', 'correction']
correction_1I = mutation_count[mutation_count$Sample == '1I', 'correction']
correction_1J = mutation_count[mutation_count$Sample == '1J', 'correction']
correction_1K = mutation_count[mutation_count$Sample == '1K', 'correction']

shared_1J_1K = inner_join(SNVs %>% filter(sample_id == '3864'), SNVs %>% filter(sample_id == '3865'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1J_1K_1I = inner_join(SNVs %>% filter(sample_id == '3863'), shared_1J_1K, by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1G_1H = inner_join(SNVs %>% filter(sample_id == '3861'), SNVs %>% filter(sample_id == '3862'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1G_1H_1C = inner_join(SNVs %>% filter(sample_id == '3857'), shared_1G_1H, by = c('chr', 'pos', 'ref', 'alt'))

shared_1J_1K_1I_1G_1H_1C = inner_join(shared_1J_1K_1I, shared_1G_1H_1C)

mrca = nrow(shared_1J_1K_1I_1G_1H_1C)*correction_1C*correction_1G*correction_1H*correction_1I*correction_1J*correction_1K
mrca_1I_1J_1K = nrow(shared_1J_1K_1I)*correction_1I*correction_1J*correction_1K - mrca
mrca_1J_1K = nrow(shared_1J_1K)*correction_1J*correction_1K - mrca - mrca_1I_1J_1K

mrca_1C_1G_1H = nrow(shared_1G_1H_1C)*correction_1C*correction_1G*correction_1H - mrca
mrca_1G_1H = nrow(shared_1G_1H)*correction_1G*correction_1H - mrca - mrca_1C_1G_1H

private_1C = (shared_mutations_SNV[[1]]['1C', '1C']*correction_1C - mrca - mrca_1C_1G_1H)
private_1G = (shared_mutations_SNV[[1]]['1G', '1G']*correction_1G - mrca - mrca_1C_1G_1H - mrca_1G_1H)
private_1H = (shared_mutations_SNV[[1]]['1H', '1H']*correction_1H - mrca - mrca_1C_1G_1H - mrca_1G_1H)
private_1I = (shared_mutations_SNV[[1]]['1I', '1I']*correction_1I - mrca - mrca_1I_1J_1K)
private_1J = (shared_mutations_SNV[[1]]['1J', '1J']*correction_1J - mrca - mrca_1I_1J_1K - mrca_1J_1K)
private_1K = (shared_mutations_SNV[[1]]['1K', '1K']*correction_1K - mrca - mrca_1I_1J_1K - mrca_1J_1K)

path = 'figures/figure3/phylogenies/timeline_1C.pdf'
width = 8.27  
height = 11.69/3 
# Set margins
left_margin = 1  
right_margin = 1  

# Calculate plot width and height after considering margins
plot_width = width - left_margin - right_margin
plot_height = height 
pdf(path, width = plot_width, height = plot_height)

mutation_len = 0.00095 #Mutation size is set to ~0.15 pt (0.002 inches)
vertical = 1/16 # 1/16 inch 5pt
mutations_yearly = mutation_len*mean

plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), asp = 1, axes = FALSE, xlab = "", ylab = "")
par(lwd = 0.5)
segments(0, 1, mrca*mutation_len, 1) #mrca trunk

segments(mrca*mutation_len, 1-8*vertical, mrca*mutation_len, 1+8*vertical)
segments(mrca*mutation_len, 1-8*vertical, (mrca+mrca_1I_1J_1K)*mutation_len, 1-8*vertical) #1I,1J,1K
segments(mrca*mutation_len, 1+8*vertical, (mrca+mrca_1C_1G_1H)*mutation_len, 1+8*vertical) #1C,1G,1H
#lci
segments(mrca*mutation_len_lci, 1-0.5*vertical, mrca*mutation_len_lci, 1+0.5*vertical, col = 'blue')
segments((mrca+mrca_1I_1J_1K)*mutation_len_lci, 1-7.5*vertical, (mrca+mrca_1I_1J_1K)*mutation_len_lci, 1-8.5*vertical, col = 'red') #1I,1J,1K
segments((mrca+mrca_1C_1G_1H)*mutation_len_lci, 1+8.5*vertical, (mrca+mrca_1C_1G_1H)*mutation_len_lci, 1+7.5*vertical, col = 'green') #1C,1G,1H
#uci
segments(mrca*mutation_len_uci, 1-0.5*vertical, mrca*mutation_len_uci, 1+0.5*vertical, col = 'blue')
segments((mrca+mrca_1I_1J_1K)*mutation_len_uci, 1-7.5*vertical, (mrca+mrca_1I_1J_1K)*mutation_len_uci, 1-8.5*vertical, col = 'red') #1I,1J,1K
segments((mrca+mrca_1C_1G_1H)*mutation_len_uci, 1+8.5*vertical, (mrca+mrca_1C_1G_1H)*mutation_len_uci, 1+7.5*vertical, col = 'green') #1C,1G,1H

segments((mrca+mrca_1I_1J_1K)*mutation_len, 1-12*vertical, (mrca+mrca_1I_1J_1K)*mutation_len, 1-4*vertical)
segments((mrca+mrca_1C_1G_1H)*mutation_len, 1+4*vertical, (mrca+mrca_1C_1G_1H)*mutation_len, 1+12*vertical)

segments((mrca+mrca_1I_1J_1K)*mutation_len, 1-12*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len, 1-12*vertical)
segments((mrca+mrca_1C_1G_1H)*mutation_len, 1+4*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len, 1+4*vertical)
#lci and uci
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len_lci, 1-11.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len_lci, 1-12.5*vertical, col = 'pink')
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len_lci, 1+4.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len_lci, 1+3.5*vertical, col = 'yellow')
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len_uci, 1-11.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len_uci, 1-12.5*vertical, col = 'pink')
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len_uci, 1+4.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len_uci, 1+3.5*vertical, col = 'yellow')

segments((mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len, 1-14*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len, 1-10*vertical)
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len, 1+2*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len, 1+6*vertical)

segments((mrca+mrca_1I_1J_1K)*mutation_len, 1-4*vertical, (mrca+mrca_1I_1J_1K+private_1I)*mutation_len, 1-4*vertical)
segments((mrca+mrca_1C_1G_1H)*mutation_len, 1+12*vertical, (mrca+mrca_1C_1G_1H+private_1C)*mutation_len, 1+12*vertical)
#lci and uci
segments((mrca+mrca_1I_1J_1K+private_1I)*mutation_len_lci, 1-3.5*vertical, (mrca+mrca_1I_1J_1K+private_1I)*mutation_len_lci, 1-4.5*vertical, col = 'green')
segments((mrca+mrca_1C_1G_1H+private_1C)*mutation_len_lci, 1+12.5*vertical, (mrca+mrca_1C_1G_1H+private_1C)*mutation_len_lci, 1+11.5*vertical, col = 'red')
segments((mrca+mrca_1I_1J_1K+private_1I)*mutation_len_uci, 1-3.5*vertical, (mrca+mrca_1I_1J_1K+private_1I)*mutation_len_uci, 1-4.5*vertical, col = 'green')
segments((mrca+mrca_1C_1G_1H+private_1C)*mutation_len_uci, 1+12.5*vertical, (mrca+mrca_1C_1G_1H+private_1C)*mutation_len_uci, 1+11.5*vertical, col = 'red')

segments((mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len, 1-10*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1J)*mutation_len, 1-10*vertical)
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K)*mutation_len, 1-14*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1K)*mutation_len, 1-14*vertical)
#lci and uci
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1J)*mutation_len_lci, 1-9.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1J)*mutation_len_lci, 1-10.5*vertical, col = 'green')
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1K)*mutation_len_lci, 1-13.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1K)*mutation_len_lci, 1-14.5*vertical, col = 'red')
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1J)*mutation_len_uci, 1-9.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1J)*mutation_len_uci, 1-10.5*vertical, col = 'green')
segments((mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1K)*mutation_len_uci, 1-13.5*vertical, (mrca+mrca_1I_1J_1K+mrca_1J_1K+private_1K)*mutation_len_uci, 1-14.5*vertical, col = 'red')

segments((mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len, 1+6*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1G)*mutation_len, 1+6*vertical)
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H)*mutation_len, 1+2*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1H)*mutation_len, 1+2*vertical)
#lci and uci
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1G)*mutation_len_lci, 1+6.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1G)*mutation_len_lci, 1+5.5*vertical, col = 'green')
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1H)*mutation_len_lci, 1+2.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1H)*mutation_len_lci, 1+1.5*vertical, col = 'red')
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1G)*mutation_len_uci, 1+6.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1G)*mutation_len_uci, 1+5.5*vertical, col = 'green')
segments((mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1H)*mutation_len_uci, 1+2.5*vertical, (mrca+mrca_1C_1G_1H+mrca_1G_1H+private_1H)*mutation_len_uci, 1+1.5*vertical, col = 'red')

segments(metadata_WGS$Age[metadata_WGS$Patient==1][1]*mutations_yearly, 0.5, metadata_WGS$Age[metadata_WGS$Patient==1][1]*mutations_yearly, 1.5)
dev.off()

## Chronological timeline (years) ####
path = 'figures/figure3/phylogenies/chronological_time.pdf'
width = 8.27  # Width of A4 in inches
height = 11.69/3 # Height of A4 in inches
# Set margins
left_margin = 1  # Left margin in inches
right_margin = 1  # Right margin in inches

# Calculate plot width and height after considering margins
plot_width = width - left_margin - right_margin
plot_height = height 
pdf(path, width = plot_width, height = plot_height)

mutation_len = 0.00095 # Mutation size is set to ~0.15 pt (0.002 inches)
vertical = 1/16 # 1/16 inch 5pt
mutations_yearly = mutation_len*mean

plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), asp = 1, axes = FALSE, xlab = "", ylab = "")
par(lwd=0.5)
segments(0, 1, mutations_yearly*100, 1) # Creates 100 years timeline in tree scale
dev.off()

# FIGURE 3F (MRCA) ####
age = metadata_WGS %>% filter(Shared > 0) %>% select(Shared, Age, Sample, Clinical_type, Patient)
age = age %>% mutate(Shared_WGD = Shared) 

age[age$Sample == '3B', 'Shared_WGD'] = (((age[age$Sample == '3B', 'Shared'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]
age[age$Sample == '5E', 'Shared_WGD'] = (((age[age$Sample == '5E', 'Shared'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]
age[age$Sample == '6A', 'Shared_WGD'] = (((age[age$Sample == '6A', 'Shared'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]
age[age$Sample == '6B', 'Shared_WGD'] = (((age[age$Sample == '6B', 'Shared'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]
age[age$Sample == '16A', 'Shared_WGD'] = (((age[age$Sample == '16A', 'Shared'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]

age = age %>% mutate(Shared_corrected = case_when(Sample == '1C' ~ Shared_WGD*correction_1C*correction_1G*correction_1H*correction_1I*correction_1J*correction_1K,
                                                  Sample == '2G' ~ Shared_WGD*mutation_count[mutation_count$Sample == '2G', 'correction']*mutation_count[mutation_count$Sample == '2L', 'correction'],
                                                  Sample == '3B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '3B', 'correction']*mutation_count[mutation_count$Sample == '3D', 'correction'],
                                                  Sample == '4B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '4B', 'correction']*mutation_count[mutation_count$Sample == '4D', 'correction'],
                                                  Sample == '5E' ~ Shared_WGD*mutation_count[mutation_count$Sample == '5E', 'correction']*mutation_count[mutation_count$Sample == '5F', 'correction'],
                                                  Sample == '6A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '6A', 'correction']*mutation_count[mutation_count$Sample == '6C', 'correction'],
                                                  Sample == '6B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '6B', 'correction']*mutation_count[mutation_count$Sample == '6D', 'correction'],
                                                  Sample == '7A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '7A', 'correction']*mutation_count[mutation_count$Sample == '7D', 'correction'],
                                                  Sample == '7B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '7B', 'correction']*mutation_count[mutation_count$Sample == '7C', 'correction'],
                                                  Sample == '8A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '8A', 'correction']*mutation_count[mutation_count$Sample == '8C', 'correction']*mutation_count[mutation_count$Sample == '8D', 'correction'],
                                                  Sample == '10A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '10A', 'correction']*mutation_count[mutation_count$Sample == '10D', 'correction'],
                                                  Sample == '10B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '10B', 'correction']*mutation_count[mutation_count$Sample == '10E', 'correction'],
                                                  Sample == '11B' ~ Shared_WGD*mutation_count[mutation_count$Sample == '11B', 'correction']*mutation_count[mutation_count$Sample == '11D', 'correction'],
                                                  Sample == '12A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '12A', 'correction']*mutation_count[mutation_count$Sample == '12B', 'correction']*mutation_count[mutation_count$Sample == '12C', 'correction'],
                                                  Sample == '13A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '13A', 'correction']*mutation_count[mutation_count$Sample == '13B', 'correction']*mutation_count[mutation_count$Sample == '13C', 'correction'],
                                                  Sample == '14A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '14A', 'correction']*mutation_count[mutation_count$Sample == '14B', 'correction']*mutation_count[mutation_count$Sample == '14C', 'correction'],
                                                  Sample == '15A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '15A', 'correction']*mutation_count[mutation_count$Sample == '15B', 'correction'],
                                                  Sample == '16A' ~ Shared_WGD*mutation_count[mutation_count$Sample == '16A', 'correction']*mutation_count[mutation_count$Sample == '16B', 'correction']
                                              ))
age = age %>% mutate(shared_mean = Shared_corrected/mean) %>% 
  mutate(shared_uci = Shared_corrected/(mean+1.64*sd_mean)) %>% 
  mutate(shared_lci = Shared_corrected/(mean-1.64*sd_mean))

order = (age %>% arrange(shared_mean) %>% select(Sample))[[1]]
age$Sample <- factor(age$Sample, levels = order)
order = (age %>% arrange(Shared_corrected) %>% select(Sample))[[1]]

pdf('figures/figure3/figure3_MRCA.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
ggplot(age %>% filter(!Patient %in% c(3,5,12) ), aes(x = factor(Sample))) + #remove patients where one lesion in MRCA have purity < 0.4
  geom_col(
    aes(y = Age, x = factor(Sample)),
    fill = '#e1a248',
    position = position_dodge(preserve = "single", width = 0.8)
  ) +
  geom_col(
    aes(y = shared_mean, x = factor(Sample)),
    fill = '#9ccb88',
    position = position_dodge(preserve = "single", width = 0.8)
  ) +
  geom_errorbar(aes(x = factor(Sample), ymin = shared_uci, ymax = shared_lci), position = position_dodge(preserve = 'single', width = 0.8)) +
  theme_classic() +
  scale_y_continuous(
    name = "Patient age (years)",
    sec.axis = sec_axis(~ . *mean, name = "SNVs")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  xlab('Primary tumor') +
  ggtitle('Primary/met. MRCA')
dev.off()

# FIGURE 3G (Scatter) ####

# See scripts/supp_figures/Supplementary_Note_1/scatter.R (produces all scatter plots included in supplementary)