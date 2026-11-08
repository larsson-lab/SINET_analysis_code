cosmic_palette = c(
  "other" = "#c4bbc4",
  "SBS1"  = "#9ccb87",
  "SBS3"  = "#ddd9a4",
  "SBS5"  = "#aec7e8",
  "SBS8"  = "#f7c7dd",
  "SBS40" = "#c49cc9",
  "SBS25" = "#e76f51",  
  "SBS89" = "#8e584d",
  "SBS85" = "#f07e23",
  "SBS21" = "#5b9e8a",
  "SBS22" = "#8ecfc9",
  "SBS57" = "#9467bd",  
  "SBS84" = "#88b04b",
  "SBS2"  = "#c6b3e5",
  "SBS13" = "#d65f8d",  
  "SBS39" = "#d9b44a",
  "SBS19" = "#7fa95b",
  "SBS31" = "#4b6584",  
  "SBS92" = "#edc949",
  "SBS12" = "#59a14f",
  "SBS15" = "#b58b54",
  "SBS87" = "#a6761d",  
  "SBS6"  = "#7a7a7a",
  "SBS93" = "#b7aca3",
  "SBS32" = "#17becf",  
  "SBS30" = "#bc5090",  
  "SBS37" = "#76b7b2"
)

metadata_nano = read_delim('input_data/metadata_libraries_SINET.txt', delim = '\t')
source('scripts/analysis/extended_metadata.R')
metadata_WGS = metadata

#---READ SigProfilerAssignment---#
activities_tot = as.matrix(read.table('results/sigprofiler/bootstrap/run1/Assignment_Solution/Activities/Assignment_Solution_Activities.txt', header = TRUE, row.names = 1))

for (i in 1:nrow(activities_tot)){
  activities_tot[i,] = activities_tot[i,]/sum(activities_tot[i,])
}

for (run in 2:500){
  print(run)
  path = paste('results/sigprofiler/bootstrap/run', as.character(run), '/Assignment_Solution/Activities/Assignment_Solution_Activities.txt', sep = '')
  activities = as.matrix(read.table(path, header = TRUE, row.names = 1))
  
  for (i in 1:nrow(activities)){
    activities[i,] = activities[i,]/sum(activities[i,])
  }
  
  activities_tot = activities_tot + activities
}

activities_tot = activities_tot/500

weights_assignment = local({
  weights = as.data.frame(activities_tot) %>%  mutate(Sample = rownames((activities_tot)))
  weights_assignment = pivot_longer(weights, colnames(weights)[-ncol(weights)]) %>% filter(value>0)
  weights_assignment = weights_assignment %>% mutate(name = if_else(value > 0.05, true = name, false = 'other'))
  weights_assignment
}) 

weights_assignment = weights_assignment %>% mutate(Method = if_else(Sample %in% metadata_WGS$Library, true = 'WGS', false = 'NanoSeq'))

weights_assignment = local({
  weights_assignment_WGS = inner_join(weights_assignment %>% rename(Sample = 'Library'), metadata_WGS %>% select(Sample, Library, Patient,SNVs_purity_corrected), by = 'Library')
  weights_assignment_nano = inner_join(weights_assignment %>% rename(Sample = 'sample_name'), metadata_nano %>% select(Sample, sample_name, Patient, Clinical_type), by = 'sample_name')
  weights_assignment_nano = bind_rows(weights_assignment_nano, weights_assignment %>% filter(Sample %in% c('blood', 'normal', 'SI-NET', 'clonal', 'nonclonal')))
  weights_assignment_WGS = weights_assignment_WGS %>% arrange(Patient, SNVs_purity_corrected)
  weights_assignment_nano = weights_assignment_nano %>% 
    mutate(Clinical_type = factor(Clinical_type, levels = c("blood", "normal", "primary","lymph_met", "carcinosis", "liver_met"))) %>% 
    arrange(Clinical_type)
  weights_assignment_WGS = weights_assignment_WGS %>% select(Sample, Patient, name, value, Method, Library)
  weights_assignment_nano = weights_assignment_nano %>% select(Sample, Patient, name, value, Method, sample_name)
  weights_assignment_nano = weights_assignment_nano %>% mutate(Sample = case_when(sample_name == 'blood' ~ 'blood',
                                                                                  sample_name == 'normal' ~ 'normal',
                                                                                  sample_name == 'SI-NET' ~'SI-NET',
                                                                                  sample_name == 'clonal' ~ 'clonal',
                                                                                  sample_name == 'nonclonal' ~ 'nonclonal',
                                                                                  .default = Sample))
  WGS_order = unique(weights_assignment_WGS$Sample)
  nano_order = unique(weights_assignment_nano$Sample)
  weights_assignment_nano = weights_assignment_nano %>% mutate(Sample = factor(Sample, levels = unique(nano_order)))
  weights_assignment_WGS = weights_assignment_WGS %>% mutate(Sample = factor(Sample, levels = unique(WGS_order)))
  colnames(weights_assignment_nano) = c('Sample', 'Patient', 'name', 'value', 'Method', 'Library')
  weights_assignment = rbind(weights_assignment_WGS, weights_assignment_nano)
  weights_assignment = rbind(weights_assignment_WGS, weights_assignment_nano)
  weights_assignment
})

#Sorted by burden SigProfilerAssignment
p1 = ggplot(weights_assignment, aes(x = Sample, y = value, fill = name)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cosmic_palette) +
  labs(x = "", y = "Fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), 
        axis.text.y.left = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 7),
        legend.title = element_blank(), 
        legend.position = 'bottom', 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm')) +
  ggtitle('SigProfilerAssignment') +
  guides(fill = guide_legend(nrow = 2)) +
  facet_grid(~Method, scales = "free_x", space = "free_x")
p1

# Read SigProfilerExtractorData ####
activities_extractor = read_delim('results/sigprofiler/deNovoSig.csv', delim = ',')

activities_extractor = activities_extractor %>% filter(Samples %in% c(metadata_nano$sample_name, metadata_WGS$Library, 'blood', 'normal', 'SI-NET', 'clonal', 'nonclonal'))
activities_extractor = activities_extractor %>% select(Samples, signature, values)
colnames(activities_extractor) = c('Sample', 'name', 'value')

weights_extractor = activities_extractor %>%
  group_by(Sample) %>%
  mutate(value = value / sum(value)) %>%
  ungroup()

weights_extractor = weights_extractor %>% mutate(name = if_else(value > 0.05, true = name, false = 'other'))
weights_extractor = weights_extractor %>% mutate(Method = if_else(Sample %in% metadata_WGS$Library, true = 'WGS', false = 'NanoSeq'))

weights_extractor = local({
  weights_extractor_WGS = inner_join(weights_extractor %>% rename(Sample = 'Library'), metadata_WGS %>% select(Sample, Library, Patient, SNVs_purity_corrected), by = 'Library')
  weights_extractor_nano = inner_join(weights_extractor %>% rename(Sample = 'sample_name'), metadata_nano %>% select(Sample, sample_name, Patient, Clinical_type), by = 'sample_name')
  print(nrow(weights_extractor_WGS))
  print(nrow(weights_extractor_nano))
  weights_extractor_nano = bind_rows(weights_extractor_nano, weights_extractor %>% filter(Sample %in% c('blood', 'normal', 'SI-NET', 'clonal', 'nonclonal')))
  weights_extractor_WGS = weights_extractor_WGS %>% arrange(Patient, SNVs_purity_corrected)
  weights_extractor_nano = weights_extractor_nano %>% 
    mutate(Clinical_type = factor(Clinical_type, levels = c("blood", "normal", "primary","lymph_met", "carcinosis", "liver_met"))) %>% 
    arrange(Clinical_type)
  weights_extractor_WGS = weights_extractor_WGS %>% select(Sample, Patient, name, value, Method, Library)
  weights_extractor_nano = weights_extractor_nano %>% select(Sample, Patient, name, value, Method, sample_name)
  weights_extractor_nano = weights_extractor_nano %>% mutate(Sample = case_when(sample_name == 'blood' ~ 'blood',
                                                                                  sample_name == 'normal' ~ 'normal',
                                                                                  sample_name == 'SI-NET' ~'SI-NET',
                                                                                  sample_name == 'clonal' ~ 'clonal',
                                                                                  sample_name == 'nonclonal' ~ 'nonclonal',
                                                                                  .default = Sample))
  WGS_order = unique(weights_extractor_WGS$Sample)
  nano_order = unique(weights_extractor_nano$Sample)
  weights_extractor_nano = weights_extractor_nano %>% mutate(Sample = factor(Sample, levels = unique(nano_order)))
  weights_extractor_WGS = weights_extractor_WGS %>% mutate(Sample = factor(Sample, levels = unique(WGS_order)))
  colnames(weights_extractor_nano) = c('Sample', 'Patient', 'name', 'value', 'Method', 'Library')
  weights_extractor = rbind(weights_extractor_WGS, weights_extractor_nano)
  weights_extractor = rbind(weights_extractor_WGS, weights_extractor_nano)
  weights_extractor
})

p2 = ggplot(weights_extractor, aes(x =Sample, y = value, fill = name)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cosmic_palette) +
  labs(x = "", y = "Fraction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6), 
        axis.text.y.left = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 7),
        legend.title = element_blank(), 
        legend.position = 'bottom', 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, 'cm')) +
  ggtitle('SigProfilerExtractor') +
  guides(fill = guide_legend(nrow = 2)) +
  facet_grid(~Method, scales = "free_x", space = "free_x")
p2

pdf('figures/supp_figures/signatures.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
p1/p2
dev.off()
