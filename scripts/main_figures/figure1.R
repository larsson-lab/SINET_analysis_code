# Clean work space and change working directory ####
rm(list = ls())

# Read libraries ####
library(dplyr) #v1.2.1
library(readr) #v2.1.5
library(lsa) #v0.73.4
library(gtools) #v3.9.5
library(svglite) #v2.1.3
library(ComplexHeatmap) #v2.18
library(circlize) #v0.4.18
library(deconstructSigs) #v1.9.0
library(BSgenome.Hsapiens.UCSC.hg38) #v1.4.5
library(ggplot2) #v4.0.3
library(patchwork) #v1.3.2
library(tibble) #v3.3.1

# Load needed input data ####
source('scripts/analysis/extended_metadata.R')
SNVs = read.delim('input_data/WGS/SNVs.tsv', sep = '\t')
mutations = read.delim('input_data/WGS/mutations.tsv', sep = '\t')
shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')

# --- FIGURE 1D (Somatic events) --- ####

#---Define functions---#
uniqeness = function(sample_list, patient, count = 1){
  sample = sample_list[[1]]
  names = names(which(shared_mutations_SNV[[patient]][which(colnames(shared_mutations_SNV[[patient]])==sample),]>10))
  
  #Update list 
  sample_list = sample_list[which(!sample_list %in% names)]
  
  #Recursion until empty list
  if(length(sample_list)!=0){
    uniqeness(sample_list, patient, count = count + 1) 
  } else{
    count
  }
}

events = local({
  events = data.frame('All' = numeric(6), 'Unique' = numeric(6))
  rownames(events) = c('CDKN1B', 'chr4', 'chr5', 'chr18', 'chr14', 'chr20')
  
  #---Driver mutations---#
  CDKN1B = mutations %>% filter(SYMBOL == 'CDKN1B')  %>% filter(IMPACT == 'HIGH') 
  colnames(CDKN1B)[1] = 'Library'
  CDKN1B = left_join(CDKN1B, metadata, by = 'Library')
  
  events['CDKN1B', 'All'] = nrow(CDKN1B)
  
  CDKN1B_patients = unique(CDKN1B$Patient)
  unique_events = 0
  for (patient in CDKN1B_patients) {
    samples = (CDKN1B %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['CDKN1B', 'Unique'] = unique_events
  
  #---Chr18---#
  chr18 = metadata %>% filter(chr18_loss %in% c('LOSS', 'CN-LOH'))
  events['chr18', 'All'] = nrow(chr18)
  
  chr18_patients = unique(chr18$Patient)
  unique_events = 0
  for (patient in chr18_patients) {
    samples = (chr18 %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['chr18', 'Unique'] = unique_events + 1 #Patient 8
  
  #---Chr4---#
  chr4 = metadata %>% filter(chr4_gain == 'GAIN')
  events['chr4', 'All'] = nrow(chr4)
  
  chr4_patients = unique(chr4$Patient)
  unique_events = 0
  for (patient in chr4_patients) {
    samples = (chr4 %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['chr4', 'Unique'] = unique_events
  
  #---Chr5---#
  chr5 = metadata %>% filter(chr5_gain == 'GAIN')
  events['chr5', 'All'] = nrow(chr5)
  
  chr5_patients = unique(chr5$Patient)
  unique_events = 0
  for (patient in chr5_patients) {
    samples = (chr5 %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['chr5', 'Unique'] = unique_events
  
  #---Chr14---#
  chr14 = metadata %>% filter(chr14_gain == 'GAIN')
  events['chr14', 'All'] = nrow(chr14)
  
  chr14_patients = unique(chr14$Patient)
  unique_events = 0
  for (patient in chr14_patients) {
    samples = (chr14 %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['chr14', 'Unique'] = unique_events
  
  #---Chr20---#
  chr20 = metadata %>% filter(chr20_gain == 'GAIN') 
  events['chr20', 'All'] = nrow(chr20)
  
  chr20_patients = unique(chr20$Patient)
  unique_events = 0
  for (patient in chr20_patients) {
    samples = (chr20 %>% filter(Patient == patient) %>%  dplyr::select(Sample))[[1]]
    count = uniqeness(samples, patient)
    unique_events = unique_events + count
  }
  
  events['chr20', 'Unique'] = unique_events
  events
})

#---Plot---#
pdf('figures/figure1/figure1_somatic_events.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
par(mfrow = c(1, 1))
barplot(t(events), beside = TRUE, ylim = c(0,50), col = c("#1e78b5", "#aec7e8"))
dev.off()

# --- FIGURE 1E (PCAWG comparison) --- ####
PCAWG = read_csv2('input_data/public_data/PCAWG_supplementary.csv', skip = 2, show_col_types = FALSE)
PCAWG = PCAWG %>% 
  dplyr::select(all.SNVs, age, histology_abbreviation, histology_tier1, histology_tier2, histology_tier3, histology_tier4) %>% 
  mutate(Sample = 'PCAWG')

PCAWG = local({
  tmp = metadata %>% 
    filter(Clinical_type != 'blood') %>% 
    filter(Clinical_type != 'normal') %>% 
    filter(Purity > 0.4) %>% 
    dplyr::select(SNVs_purity_corrected, Age, Sample) %>% 
    mutate(histology_abbreviation = 'SI-NET') %>%  
    mutate(histology_tier1 = 'SI-NET') %>% 
    mutate(histology_tier2 = 'SI-NET') %>% 
    mutate(histology_tier3 = 'SI-NET') %>% 
    mutate(histology_tier4 = 'SI-NET')  
  colnames(tmp)[1:2] = c('all.SNVs', 'age')
  PCAWG = rbind(PCAWG, tmp) %>% 
    mutate(Sample = if_else(Sample == 'PCAWG', true = 'PCAWG', false = 'SI-NET')) %>% 
    mutate(Color = if_else(Sample == 'PCAWG', true = '#C8C7C7', false = '#D72728')) 
  PCAWG
})

group = PCAWG$histology_tier2 # Add group
medians = with(PCAWG, tapply(all.SNVs, group, median)) # Calculate medians for each group
group = factor(group, levels = names(sort(medians))) # Reorder groups based on medians

colors = ifelse(levels(group) == "SI-NET", "#D72728", "#C8C7C7")
pdf('figures/figure1/figure1_PCAWG.pdf', width = 180*0.0393700787, height = 100*0.0393700787)
par(bty='n')
par(lwd = 0.5)
boxplot(all.SNVs~group, 
        data=PCAWG, 
        col = colors, 
        log='y', 
        cex = 0.5,
        cex.axis=5/12, 
        cex.lab = 5/12,
        las = 2,
        ylab = 'SNVs',
        xlab = '',
        lwd = 0.5, 
        ylim = c(0.1,1000000),
        axes = FALSE)
axis(2, lwd = 0.5, lwd.ticks = 0.5, las=2, cex.axis = 0.5)
axis(1, at = 1:length(levels(group)), 
     labels = levels(group), las = 2, cex.axis = 5/12, lwd.ticks = 0.5)

dev.off()

# --- FIGURE 1C and 1F (boxplots) --- ####
metadata = metadata %>%  mutate(Focality = ifelse(Patient %in% c(12,13,14,15,16), 'unifocal', 'multifocal')) %>% relocate(Focality, .before= SNVs)
age_data = metadata %>% dplyr::select(Patient, Age, Focality) %>%  distinct(Patient, .keep_all = TRUE)
SNV_data = metadata %>% filter(!Clinical_type %in% c('blood', 'normal')) %>% filter(Purity > 0.4) %>% dplyr::select(Patient, SNVs_purity_corrected, Focality)

age_data$Focality = factor(
  age_data$Focality,
  levels = c("unifocal", "multifocal")
)
SNV_data$Focality = factor(
  SNV_data$Focality,
  levels = c("unifocal", "multifocal")
)

pdf(
  "figures/figure1/figure1_boxplot.pdf",
  width = 12,
  height = 8
)

p1 = ggplot(age_data, aes(x = Focality, y = Age, fill = Focality)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0,
    size = 2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, 10)
  ) +
  labs(
    x = NULL,
    y = "Age"
  ) +
  scale_fill_manual(
    values = c(
      "unifocal" = '#f29697',
      "multifocal" = '#dbe18d'
    )
  ) +
  theme_classic() +
  theme(legend.position = 'none') 

p2 = ggplot(SNV_data, aes(x = Focality, y = SNVs_purity_corrected, fill = Focality)) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0,
    size = 2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 5000, 1000)
  ) +
  labs(
    x = NULL,
    y = "Mutation burden"
  ) +
  scale_fill_manual(
    values = c(
      "unifocal" = '#f29697',
      "multifocal" = '#dbe18d'
    )
  ) +
  theme_classic() +
  theme(legend.position = 'none') 

p1 / p2

dev.off()

# --- FIGURE 1G (heatmap) --- ####
trinucprof = as.data.frame(t(mut.to.sigs.input(SNVs %>% dplyr::select(sample_id, chr, pos, ref, alt), sample.id = 'sample_id', bsg= BSgenome.Hsapiens.UCSC.hg38)))

## Change from sample_id to sample name ####
for (i in 1:ncol(trinucprof)){
  sample = (metadata %>% filter(Library==colnames(trinucprof)[i]) %>% dplyr::select(Sample))[[1]]
  colnames(trinucprof)[i] = sample
}

## Remove tissue normals and sort by sample name ####
trinucprof = local({
  samples = (metadata %>%  filter(Clinical_type != 'normal') %>%  filter(Clinical_type != 'blood') %>% dplyr::select(Sample))[[1]]
  trinucprof = trinucprof %>% dplyr::select(all_of(samples))
  trinucprof = trinucprof[, mixedsort(colnames(trinucprof))]
})

## Sort by burden ####
trinucprof = local({
  tmp = metadata %>% filter(Clinical_type != 'blood') %>% filter(Clinical_type != 'normal') %>% dplyr::select(Sample, SNVs)
  trinucprof_new = t(trinucprof) %>% as.data.frame() %>%  mutate(Sample = colnames(trinucprof))
  trinucprof_new = trinucprof_new %>% left_join(tmp, by = 'Sample')
  trinucprof_new = trinucprof_new %>%  arrange(SNVs)
  rownames(trinucprof_new) = trinucprof_new$Sample
  trinucprof_new = trinucprof_new[,1:96]
  trinucprof_new = t(trinucprof_new) %>%as.data.frame()
  trinucprof_new[] = lapply(trinucprof_new, as.numeric)
  trinucprof = trinucprof_new
})

## Sort the 96 substitutions and calculate cosine similarity ####
trinucprof = trinucprof %>% mutate(mutation = substr(rownames(trinucprof),2,6)) %>% arrange(mutation)

trinucprof = as.data.frame(t(trinucprof))
trinucprof = trinucprof %>% rownames_to_column()

trinucprof = trinucprof[1:74,] %>%
  column_to_rownames("rowname") %>% 
  as.matrix() 

trinucprof = apply(trinucprof, 2, as.integer)
trinucprof = rbind(trinucprof, sum = colSums(trinucprof))

similarity = cosine(t(trinucprof/rowSums(trinucprof)))

## Plot ####
color = colorRamp2(
  c(0, 0.1),
  c("white", "blue")
)

pdf(
  "figures/figure1/figure1_heatmap.pdf",
  width = 12,
  height = 8
)

ComplexHeatmap::Heatmap(trinucprof[74:1,]/rowSums(trinucprof[74:1,]), 
                        cluster_columns = F, 
                        cluster_rows = F, 
                        col = color, 
                        show_column_names = FALSE,
                        name = "Weight in spectrum",
                        heatmap_legend_param = list(
                          title_gp = grid::gpar(fontface = "plain")
                        )) +
  rowAnnotation('Similarity to cohort mean' = anno_barplot(as.double(similarity['sum',][74:1]), 
                                                           bar_width = 0.95, 
                                                           gp = grid::gpar(fill = "#ccc9c8", col = NA), 
                                                           border = F))

dev.off()
