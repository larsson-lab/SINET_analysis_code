# Read libraries ####
library(readr)
library(vcfR)
library(dplyr)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ggpmisc)
library(gtools)
library(ComplexHeatmap)
library(circlize)
library(lsa)
library(ggnewscale)
library(ggplot2)
library(patchwork)

# Load needed input data ####
SNVs = read.delim('input_data/WGS/SNVs.tsv', sep = '\t')
metadata_WGS = read.delim('input_data/metadata.tsv', sep = '\t')
metadata = read_delim('input_data/metadata_libraries_SINET.txt')
purity = read_delim('input_data/WGS/purity.tsv') %>% mutate(Sample = as.character(Sample))
colnames(purity)[1] = 'Library'
purity = inner_join(purity, metadata_WGS, 'Library')
source('scripts/analysis/phylogenies_trunk.R') #add number of trunk mutations to metadata_WGS

## --- SUMMARIZE STATS --- ####
samples = metadata$sample_name
genome_scaled = local({
  genome = 6300000000
  scale = 0.86 # mappable genome size retrieved from in silico genome
  genome_scaled = genome*scale
  genome_scaled
})

stats_all = data.frame()
for (sample in samples){
  path = paste('input_data/NanoSeq/', sample, '_stats.txt', sep = '')
  stats = as.data.frame(t(read_delim(path, col_names = F)))
  colnames(stats) = stats[1,]
  stats = stats[2,]
  stats = stats %>% mutate(sample_name = sample)
  stats_all = rbind(stats_all, stats)
}
rm(stats, path)
stats_all = inner_join(stats_all, metadata, by = 'sample_name')

# --- SUMMARIZE BURDEN --- ####
burden_all = data.frame()
for (sample in samples){
  path = paste('input_data/NanoSeq/', sample, '_sbs_burden.txt', sep = '')
  burden = as.data.frame(t(read_delim(path, col_names = F)))
  colnames(burden) = burden[1,]
  burden = burden[2,]
  burden = burden %>% mutate(sample_name = sample)
  burden_all = rbind(burden_all, burden)  
}
rm(burden, path)
burden_all = inner_join(burden_all, metadata, by = 'sample_name')
burden_all = burden_all %>% mutate(burden = (genome_scaled*as.numeric(`Uncorrected burden`)))
burden_all = burden_all %>% mutate(burden_lci= (genome_scaled*as.numeric(`Uncorrected burden 95% lower`)))
burden_all = burden_all %>% mutate(burden_uci = (genome_scaled*as.numeric(`Uncorrected burden 95% upper`)))
burden_all = burden_all %>% mutate(mutation_rate = (genome_scaled*as.numeric(`Uncorrected burden`))/Age)
burden_all = burden_all %>% mutate(mutation_rate_lci = (genome_scaled*as.numeric(`Uncorrected burden 95% lower`))/Age)
burden_all = burden_all %>% mutate(mutation_rate_uci = (genome_scaled*as.numeric(`Uncorrected burden 95% upper`))/Age)

burden_all = burden_all %>% mutate(B_normal = case_when(Clinical_type == 'primary' ~ mutation_rate,
                                                        Clinical_type == 'liver_met' ~ 60,
                                                        Clinical_type == 'lymph_met' ~ 52.69,
                                                        Clinical_type == 'carcinosis' ~ mutation_rate,
                                                        Clinical_type == 'blood' ~ mutation_rate,
                                                        Clinical_type == 'normal' ~ mutation_rate))

burden_all = left_join(burden_all, purity %>% select(Sample, Purity), by = 'Sample')
burden_all = burden_all %>% mutate(Purity = if_else(is.na(Purity), true = 1, false = Purity))
burden_all = burden_all %>% 
  mutate(burden_scaled = ((burden - (1-Purity)*B_normal*Age)/Purity)) %>%
  mutate(burden_uci_scaled = ((burden_uci - (1-Purity)*B_normal*Age)/Purity)) %>%
  mutate(burden_lci_scaled = ((burden_lci - (1-Purity)*B_normal*Age)/Purity)) 
burden_all = burden_all %>% 
  mutate(mutation_rate_scaled = burden_scaled/Age) %>%
  mutate(mutation_rate_uci_scaled = burden_uci_scaled/Age) %>%
  mutate(mutation_rate_lci_scaled = burden_lci_scaled/Age) 
burden_all = burden_all %>% mutate(q = 1-((1-Purity)*((B_normal*Age)/(burden))))

## --- SUMMARIZE MUTATIONS --- ####
mutations = data.frame()

for (sample in samples){
  path = paste('input_data/NanoSeq/', sample,  '_snv.vcf', sep = '')
  muts = vcfR::read.vcfR(path)
  muts = as.data.frame(muts@fix)
  muts = muts %>% mutate(sample_name = sample) %>% filter(FILTER == 'PASS')
  mutations = rbind(mutations, muts)
}
rm(muts, path)
mutations = inner_join(mutations, metadata, by = 'sample_name')

## --- ADD NUMBER OF MUTATIONS PER SAMPLE TO STATS --- ####
summary_mutations = as.data.frame(table(mutations$sample_name))
colnames(summary_mutations) = c('sample_name', 'count')
stats_all = inner_join(stats_all, summary_mutations, by = 'sample_name')

# --- MUTATION OVERLAP WITH WGS --- ####
mutations_WGS = SNVs %>% dplyr::select(chr, pos, ref, alt, sample_id) %>% mutate(sample_id = as.character(sample_id))
colnames(mutations_WGS) = c('CHROM', 'POS', 'REF', 'ALT', 'Library')
mutations_tmp = mutations %>% mutate(POS = as.double(POS)) %>% filter(!Clinical_type %in% c('blood'))

clonal = inner_join(mutations_tmp, mutations_WGS, by = c('CHROM', 'POS', 'REF', 'ALT', 'Library'))
fractions = full_join(as.data.frame(table(mutations_tmp$Sample)), as.data.frame(table(clonal$Sample)), by = 'Var1')
colnames(fractions) = c('Sample', 'mutations_tot', 'mutations_clonal')
fractions = fractions %>% mutate(mutations_clonal = if_else(is.na(mutations_clonal), true = 0, false = mutations_clonal))
fractions = fractions %>%  mutate(fraction_clonal = as.numeric(mutations_clonal)/as.numeric(mutations_tot)) %>% mutate(fraction_nonclonal = 1 - fraction_clonal)
fractions = inner_join(fractions, metadata, by = 'Sample')
fractions_scaled = inner_join(fractions %>% select(Sample, mutations_tot, mutations_clonal, fraction_clonal, fraction_nonclonal), burden_all, by = 'Sample')

fractions_scaled = fractions_scaled %>% 
  mutate(fraction_clonal_scaled = fraction_clonal/q) 

# --- PREPARE DATA FOR TRINUCPROFILE CLONAL VS NON-CLONAL --- ####
clonal = clonal %>% mutate(group = 'clonal') 
nonclonal = anti_join(mutations %>% mutate(POS = as.double(POS)) %>% filter(Clinical_type != 'blood'), clonal, by = c('CHROM', 'POS', 'REF', 'ALT', 'Library')) %>% 
  mutate(group = 'nonclonal') %>% dplyr::select('group', 'CHROM', 'POS', 'REF', 'ALT')
clonal = clonal %>% dplyr::select('group', 'CHROM', 'POS', 'REF', 'ALT')
clonal_nonclonal = rbind(clonal, nonclonal)
colnames(clonal_nonclonal) = c('group', 'chr', 'pos', 'ref', 'alt')
write_delim(clonal_nonclonal, 'input_data/clonal_nonclonal.tsv', delim = '\t') #write to memory to use in signature analysis
clonal_nonclonal_profile = mut.to.sigs.input(clonal_nonclonal, sample.id = 'group', bsg= BSgenome.Hsapiens.UCSC.hg38) 
clonal_nonclonal_profile = as.data.frame(t(clonal_nonclonal_profile))
clonal_nonclonal_profile = clonal_nonclonal_profile %>% rownames_to_column(var = 'context')
rm(clonal, nonclonal, mutations_tmp, mutations_WGS)

## --- MUTATION RATE MEAN AND SD --- ####
x = (burden_all %>% filter(!Clinical_type %in% c('blood', 'normal', 'liver_met')) %>%  dplyr::select(mutation_rate_scaled))[[1]] 
mean = mean(x)
uci = mean + 1.64*sd(x)
lci = mean - 1.64*sd(x)

mean_primary = mean((burden_all %>% filter(Clinical_type %in% c('primary')) %>% dplyr::select(mutation_rate_scaled))[[1]])
mean_lymph = mean((burden_all %>% filter(Clinical_type %in% c('lymph_met')) %>% dplyr::select(mutation_rate_scaled))[[1]])
mean_peritoneal = mean((burden_all %>% filter(Clinical_type %in% c('carcinosis')) %>% dplyr::select(mutation_rate_scaled))[[1]])
mean_liver = mean((burden_all %>% filter(Clinical_type %in% c('liver_met')) %>% dplyr::select(mutation_rate_scaled))[[1]])

sd_primary = sd((burden_all %>% filter(Clinical_type %in% c('primary')) %>% dplyr::select(mutation_rate_scaled))[[1]])
sd_lymph = sd((burden_all %>% filter(Clinical_type %in% c('lymph_met')) %>% dplyr::select(mutation_rate_scaled))[[1]])
sd_peritoneal = sd((burden_all %>% filter(Clinical_type %in% c('carcinosis')) %>% dplyr::select(mutation_rate_scaled))[[1]])
sd_liver = sd((burden_all %>% filter(Clinical_type %in% c('liver_met')) %>% dplyr::select(mutation_rate_scaled))[[1]])
sd_mean = sd(x)

mean_df = data.frame(c(mean_primary, mean_lymph, mean_peritoneal, mean_liver, mean), row.names = c('primary', 'lymph_met', 'carcinosis', 'liver_met', 'combined'))
mean_df = rownames_to_column(mean_df, var = 'Clinical_type')
mean_df = cbind(mean_df, c(sd_primary, sd_lymph, sd_peritoneal, sd_liver, sd_mean))
colnames(mean_df) = c('Clinical_type', 'mean', 'SD')
mean_df = mean_df %>% 
  mutate(uci_90 = mean + 1.64*SD) %>% 
  mutate(lci_90 = mean - 1.64*SD) 

rm(sample, samples)

## --- WGS data --- ####
mutation_count = as.data.frame(table(SNVs$sample_id))
colnames(mutation_count) = c('Library', 'count')

## --- Age---####
mutation_count = inner_join(mutation_count, purity, by = 'Library')
mutation_count = mutation_count %>% filter(!Clinical_type %in% c('blood', 'normal')) 
mutation_count = mutation_count %>% mutate(mean = if_else(Clinical_type == 'liver_met', true = mean_liver, false = mean))
mutation_count = mutation_count %>% 
  mutate(uci_90 = mean + 1.64*sd_mean) %>% 
  mutate(lci_90 = mean - 1.64*sd_mean) #use sd_mean also for liver

## --- Adjust mutation count for WGD --- ####
WGD = read_delim('input_data/WGS/WGD.tsv', delim = '\t')
mutation_count = mutation_count %>% mutate(count_WGD = count)

mutation_count[mutation_count$Sample == '3B', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '3B', 'count'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]
mutation_count[mutation_count$Sample == '3D', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '3D', 'count'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]

mutation_count[mutation_count$Sample == '5E', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '5E', 'count'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]
mutation_count[mutation_count$Sample == '5F', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '5F', 'count'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]

mutation_count[mutation_count$Sample == '6A', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '6A', 'count'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]
mutation_count[mutation_count$Sample == '6C', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '6C', 'count'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]

mutation_count[mutation_count$Sample == '6B', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '6B', 'count'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]
mutation_count[mutation_count$Sample == '6D', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '6D', 'count'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]

mutation_count[mutation_count$Sample == '16A', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '16A', 'count'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]
mutation_count[mutation_count$Sample == '16B', 'count_WGD'] = (((mutation_count[mutation_count$Sample == '16B', 'count'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]

## -- Adjust counts for purity --- ####
get_r = function(a, b, p) {
  r = (1-exp(-a*p^b))/(1-exp(-a))
  return(r)
}

get_correction_factor = function(p){
  r = get_r(a,b,p)
  cf = 1/r
  return(cf)
}

# Read parameters a and b
source('scripts/analysis/curvefit.R')

# Calculate r
mutation_count = mutation_count %>% 
  mutate(recall = get_r(a,b,Purity)) %>% 
  mutate(correction = get_correction_factor(Purity)) %>% 
  mutate(count_corrected = count_WGD*correction)

mutation_count = mutation_count %>% mutate(ce_mean = if_else(Clinical_type == 'liver_met', true = count_corrected/mean_liver, false = count_corrected/mean)) %>%
  mutate(ce_uci_90 = count_corrected/uci_90) %>% 
  mutate(ce_lci_90 = count_corrected/lci_90) 

# --- PREPARE DATA FOR HEATMAP --- ####
samples = metadata$sample_name
profile_all = data.frame(row.names = seq_len(96))

for (sample in samples){
  print(sample)
  path = paste('input_data/NanoSeq/', sample, '_sbs_96_corrected.txt', sep = '')
  profile = read_delim(path, col_names = T)
  profile = profile[,c(1,2)]
  colnames(profile) = c('context', sample)
  if (sample == samples[1]){
    profile_all = cbind(profile_all, profile)
  } else{
    profile_all = cbind(profile_all, profile[,2]) 
  }
}
rm(profile, path, sample, samples)
profile_all = profile_all %>% mutate(subs = substr(context,3,5)) %>% arrange(subs, context)
profile_all = as.data.frame(t(profile_all))
colnames(profile_all) = profile_all[1,]
profile_all = profile_all[2:41,]
profile_all = profile_all %>% rownames_to_column(var = 'sample_name')
profile_all = inner_join(profile_all, metadata, by = 'sample_name')
profile_all = profile_all %>%  
  mutate(Clinical_type = factor(Clinical_type, levels = c('blood', 'normal', 'primary', 'lymph_met', 'carcinosis', 'liver_met'))) %>% 
  arrange(Clinical_type)

mat = as.matrix(profile_all[,2:97])
mat = matrix(as.numeric(mat), nrow = nrow(mat))

blood_mat = profile_all %>% filter(Clinical_type == 'blood')
blood_mat = rbind(blood_mat[,2:97], sum = colSums(sapply(blood_mat[,2:97], as.numeric)))
blood_mat = sapply(blood_mat[,1:96], as.numeric)

normal_mat = profile_all %>% filter(Clinical_type == 'normal')
normal_mat = rbind(normal_mat[,2:97], sum = colSums(sapply(normal_mat[,2:97], as.numeric)))
normal_mat = sapply(normal_mat[,1:96], as.numeric)

SINET_mat = profile_all %>% filter(Clinical_type %in% c("primary", "liver_met", "lymph_met", "carcinosis"))
SINET_mat = rbind(SINET_mat[,2:97], sum = colSums(sapply(SINET_mat[,2:97], as.numeric)))
SINET_mat = sapply(SINET_mat[,1:96], as.numeric)

blood_similarity = cosine(t(blood_mat[,1:96]/rowSums(blood_mat[,1:96])))
normal_similarity = cosine(t(normal_mat[,1:96]/rowSums(normal_mat[,1:96])))
SINET_similarity = cosine(t(SINET_mat[,1:96]/rowSums(SINET_mat[,1:96])))
similarity = c(blood_similarity[nrow(blood_similarity),1:ncol(blood_similarity)-1], normal_similarity[nrow(normal_similarity),1:ncol(normal_similarity)-1],SINET_similarity[nrow(SINET_similarity),1:ncol(SINET_similarity)-1])

rm(SINET_mat, SINET_similarity, blood_mat, blood_similarity, normal_mat, normal_similarity, metadata, clonal_nonclonal, mutations)