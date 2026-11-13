# Important! Set the code repository as working directory before running the script.
# The paths for files from Mutect2 and Varscan are examples - see README for explanation

# Read libraries ####
library(data.table) #v1.18.6.1
library(dplyr) #v1.2.1
library(readr) #v2.1.5

# Load mutations mutect ####
load_mutations_mutect = function(sample, meta_data){
  sample_code = meta_data %>%
    filter(meta_data$Sample == as.character(sample)) %>%
    select(Library)
  # tsv output from Mutect2, only mutations marked with PASS (this file should be annotated with vep)
  file = paste('vcf/mutect/', sample_code,'.filtered.onlyPASS.vep.tsv', sep="") 
  mutations = read_tsv(file, show_col_types = FALSE)
  return(mutations)
}

# Load mutations mutect extended ####
load_mutations_mutect_extended = function(sample, meta_data){
  sample_code = meta_data %>%
    filter(meta_data$Sample == as.character(sample)) %>%
    select(Library)
  # tsv output from Mutect2 with mutations obtained after FilterMutectCalls
  file = paste('vcf/mutect/', sample_code, '.filtered.tsv', sep = '') 
  mutations = read_tsv(file, show_col_types = FALSE)
  mutations$sample_id = as.character(mutations$sample_id)
  return(mutations)
}

# Load mutations varscan extended ####
load_mutations_varscan_extended= function(sample, meta_data){
  sample_code = (meta_data %>%
                   filter(meta_data$Sample == as.character(sample)) %>%
                   select(Library))[[1]]
  # tsv output from Varscan (not hc)
  file_indels = paste('vcf/varscan/', sample_code,'.indel.Somatic', sep="") #Do not forget to rewrite varscan indels!
  file_snps = paste('vcf/varscan/', sample_code, '.snp.Somatic', sep="") 
  indels = read_tsv(file_indels, show_col_types = FALSE)
  snps = read_tsv(file_snps, show_col_types = FALSE)
  
  indels = indels %>%  mutate('sample_id' = sample_code)
  snps = snps %>%  mutate('sample_id' = sample_code)
  mutations = rbind(indels, snps)
  mutations = mutations %>%  select(sample_id, chrom, position, ref, var,
                                    tumor_reads2_plus, tumor_reads2_minus, tumor_reads1_plus, 
                                    tumor_reads1_minus, normal_reads2_plus, normal_reads2_minus, 
                                    normal_reads1_plus, normal_reads1_minus)
  
  colnames(mutations) = c("sample_id", "chr", "pos", "ref", "alt", "tumor_alt_read_depth_r1",
                          "tumor_alt_read_depth_r2",  "tumor_ref_read_depth_r1", "tumor_ref_read_depth_r2",  "normal_alt_read_depth_r1",
                          "normal_alt_read_depth_r2", "normal_ref_read_depth_r1", "normal_ref_read_depth_r2")
  
  return(mutations)
}

# Load mutations varscan ####
load_mutations_varscan= function(sample, meta_data){
  sample_code = (meta_data %>%
                   filter(meta_data$Sample == as.character(sample)) %>%
                   select(Library))[[1]]
  # tsv output from Varscan (hc)
  file_indels = paste('vcf/varscan/', sample_code,'.indel.Somatic.hc', sep="") 
  file_snps = paste('vcf/varscan/', sample_code,'.snp.Somatic.hc', sep="") 
  indels = read_tsv(file_indels, show_col_types = FALSE)
  snps = read_tsv(file_snps, show_col_types = FALSE)
  
  indels = indels %>%  mutate('sample_id' = sample_code)
  snps = snps %>%  mutate('sample_id' = sample_code)
  mutations = rbind(indels, snps)
  mutations$sample_id = as.character(mutations$sample_id)
  
  return(mutations)
}

# Filter mutect ####
filter_mutations_mutect= function(coverage_n_min, strand_support_min, tumor_support_min, normal_support_max, patient){
  samples = (metadata %>% filter(Patient == patient & Tumor != '-') %>%  select(Sample))[[1]]
  chromosomes = c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY')
  mutations_filtered_tot = data.frame()
  for (sample in samples){
    mutations = load_mutations_mutect(sample, metadata)
    mutations = mutations %>%
      filter(mutations$chr %in% chromosomes)
    
    coverage_n = mutations$normal_alt_read_depth_r1 + mutations$normal_alt_read_depth_r2 +
      mutations$normal_ref_read_depth_r1 + mutations$normal_ref_read_depth_r2
    
    mutations_filtered = mutations %>% 
      filter(coverage_n >= coverage_n_min &
               tumor_alt_read_depth_r1 >= strand_support_min &
               tumor_alt_read_depth_r2 >= strand_support_min &           
               (tumor_alt_read_depth_r1 + tumor_alt_read_depth_r2) >= tumor_support_min &
               (normal_alt_read_depth_r1 + normal_alt_read_depth_r2) <= normal_support_max) 
    
    mutations_filtered_tot = rbind(mutations_filtered_tot, mutations_filtered) 
  }
  mutations_filtered_tot = mutations_filtered_tot
  mutations_filtered_tot$sample_id = as.character(mutations_filtered_tot$sample_id)
  return(mutations_filtered_tot)
}

# Filter varscan ####
filter_mutations_varscan = function(coverage_n_min, strand_support_min, tumor_support_min, normal_support_max, patient){
  
  samples = (metadata %>% filter(Patient == patient & Tumor != '-') %>%  select(Sample))[[1]]
  chromosomes = c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY')
  mutations_filtered_tot = data.frame()
  print(samples)
  for (sample in samples){
    mutations = load_mutations_varscan(sample, metadata)
    mutations = mutations %>%
      filter(mutations$chrom %in% chromosomes)
    
    coverage_n = mutations$normal_reads1 + mutations$normal_reads2
    
    mutations_filtered = mutations %>% 
      filter(coverage_n >= coverage_n_min & 
               tumor_reads2 >= tumor_support_min & 
               normal_reads2 <= normal_support_max & 
               tumor_reads2_plus >= strand_support_min &
               tumor_reads2_minus >= strand_support_min)
    
    mutations_filtered = mutations_filtered %>%  select(chrom, position, var, ref, sample_id) %>% rename(chr=chrom, pos = position, alt = var)
    mutations_filtered_tot = rbind(mutations_filtered_tot, mutations_filtered)
  }
  return(mutations_filtered_tot) 
}

# Whitelisting ####
whitelist = function(patient, sample, whitelist){
  mutations_mutect = load_mutations_mutect_extended(sample, metadata)
  mutations_varscan = load_mutations_varscan_extended(sample, metadata)
  mutations = full_join(mutations_mutect, mutations_varscan, by = c('chr', 'pos', 'ref', 'alt', 'sample_id'))
  mutations = mutations %>% mutate(tumor_alt_read_depth_r1 = coalesce(tumor_alt_read_depth_r1.x, tumor_alt_read_depth_r1.y)) %>%
    mutate(tumor_alt_read_depth_r2 = coalesce(tumor_alt_read_depth_r2.x, tumor_alt_read_depth_r2.y)) %>% 
    mutate(tumor_ref_read_depth_r1 = coalesce(tumor_ref_read_depth_r1.x, tumor_ref_read_depth_r1.y)) %>% 
    mutate(tumor_ref_read_depth_r2 = coalesce(tumor_ref_read_depth_r2.x, tumor_ref_read_depth_r2.y)) %>% 
    mutate(normal_alt_read_depth_r1 = coalesce(normal_alt_read_depth_r1.x, normal_alt_read_depth_r1.y)) %>% 
    mutate(normal_alt_read_depth_r2 = coalesce(normal_alt_read_depth_r2.x, normal_alt_read_depth_r2.y)) %>%
    mutate(normal_ref_read_depth_r1 = coalesce(normal_ref_read_depth_r1.x, normal_ref_read_depth_r1.y)) %>%
    mutate(normal_ref_read_depth_r2 = coalesce(normal_ref_read_depth_r2.x, normal_ref_read_depth_r2.y)) %>%
    select("sample_id", "chr", "pos", "ref", "alt", "tumor_alt_read_depth_r1", "tumor_alt_read_depth_r2",  
           "tumor_ref_read_depth_r1",  "tumor_ref_read_depth_r2",  "normal_alt_read_depth_r1",
           "normal_alt_read_depth_r2", "normal_ref_read_depth_r1", "normal_ref_read_depth_r2")
  extended_sample = data.frame()
  for (i in 1:ncol(shared_mutations_final[[patient]])){
    if (shared_mutations_final[[patient]][sample, i] > 10 & colnames(shared_mutations_final[[patient]])[i] != sample){
      paired = colnames(shared_mutations_final[[patient]])[i]
      paired = (metadata[metadata$Sample == paired, 'Library'])[[1]]
      sample_ID = (metadata[metadata$Sample == sample, 'Library'])[[1]]
      mismatch = anti_join(whitelist %>% filter(sample_id == paired), whitelist %>% filter(sample_id == sample_ID), by = c('chr', 'pos', 'ref', 'alt'))
      extended = semi_join(mutations, mismatch, by = c('chr', 'pos', 'ref', 'alt'))
      extended_sample = rbind(extended_sample, extended)
    }
  }
  return(extended_sample)
}

# Overlapping mutations within patient ####
overlapping_mutations_within = function(meta_data, mutations_filtered){
  #Create empty matrix 
  samples = meta_data[1:3] %>%
    filter(meta_data$Patient == patient & meta_data$Tumor != '-')
  samples = samples[order(paste(samples$Patient, samples$Tumor, sep='')),]
  overlap_matrix = matrix(ncol=nrow(samples),nrow =nrow(samples))
  
  rownames(overlap_matrix) = paste(samples$Patient, samples$Tumor, sep = '')
  colnames(overlap_matrix) = paste(samples$Patient, samples$Tumor, sep = '')
  for (i in 1:nrow(samples)){
    for (j in i:nrow(samples)) {
      if (i == j){
        shared_mutations = mutations_filtered %>% 
          filter(sample_id == samples[[i,1]])
        overlap_matrix[i,j] = nrow(shared_mutations)
      } else {
        shared_mutations = mutations_filtered %>% 
          filter(sample_id == samples[[i,1]] | sample_id == samples[[j,1]]) %>% 
          select(chr,pos,ref,alt) %>% 
          add_count(chr,pos,ref,alt) %>% 
          filter(n>1) %>% 
          distinct() 
        overlap_matrix[i,j] = nrow(shared_mutations)
        overlap_matrix[j,i] = nrow(shared_mutations)
      }
    }
  }
  return(overlap_matrix)
}

# Main run ####
#metadata = load_metadata()
metadata = read_delim('input_data/metadata.tsv')

#Create whitelist that pass filtering criteria in both Mutect and Varscan
shared_mutations_final = vector(mode = 'list', length = 16)
all_mutations_final = data.frame()

for (patient in 1:16){
  print(patient)
  mutect = filter_mutations_mutect(20,0,1,0,patient)
  varscan = filter_mutations_varscan(20,1,1,0,patient) 
  intersected = inner_join(mutect, varscan, by = c('chr', 'pos', 'ref', 'alt', 'sample_id'))
  all_mutations_final = rbind(all_mutations_final, intersected)
}

rm(mutect, varscan, intersected)

# Remove mutations found in PON
PON = fread("input_data/PON.txt",
                   col.names = c("chr", "pos", "ref", "alt")) # Note! This file is not provided and have to be generated as described in the README

PON = as.data.frame(PON) %>% filter(chr %in% c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY')) %>% mutate(pos = as.double(pos))
all_mutations_final = anti_join(all_mutations_final, PON, by = c('chr', 'pos', 'ref', 'alt'))

for (patient in 1:16){
  samples = (metadata %>%  filter(Clinical_type != 'blood') %>% filter(Patient == patient) %>% select(Library))[[1]]
  all_shared_mutations = all_mutations_final %>%  filter(sample_id %in% samples)
  all_shared_mutations = overlapping_mutations_within(metadata, all_shared_mutations)
  shared_mutations_final[[patient]] = all_shared_mutations
}

mutations_extended = data.frame()
samples = (metadata %>%  filter(Clinical_type != 'blood') %>% filter(Patient %in% c(1:16)) %>% select(Sample))[[1]]

for (sample in samples){
  patient = (metadata[metadata$Sample==sample, 'Patient'])[[1]]
  extended = whitelist(patient, sample, all_mutations_final)
  mutations_extended = rbind(mutations_extended, extended)
}

mutations_extended = unique(mutations_extended) #Add extra data that's included in the whitelist (SYMBOL, Impact etc)
for (i in 1:nrow(mutations_extended)){
  match = all_mutations_final %>% filter(chr==mutations_extended[i,2][[1]] & 
                                           pos==mutations_extended[i,3][[1]] & 
                                           ref==mutations_extended[i,4][[1]] & 
                                           alt==mutations_extended[i,5][[1]])
  mutations_extended[i, 14:21] = match[1, 14:21]
}

mutations = unique(bind_rows(all_mutations_final, mutations_extended))

shared_mutations_SNV = vector(mode = 'list', length = 16)
for (patient in 1:16){
  intersected = mutations %>% filter(ref %in% c('A', 'C', 'T', 'G') & alt %in% c('A', 'C', 'T', 'G')) %>% filter(sample_id %in% (metadata %>% filter(Patient == patient & Clinical_type != 'blood') %>%  select(Library))[[1]])
  all_shared_mutations = overlapping_mutations_within(metadata, intersected)
  shared_mutations_SNV[[patient]] = all_shared_mutations
}

# Write SNVs to tsv file for further analyses 
SNVs = mutations %>% filter(ref %in% c('A', 'C', 'T', 'G') & alt %in% c('A', 'C', 'T', 'G'))
write_delim(SNVs, 'SNVs.tsv', delim = '\t')
write_delim(mutations, 'mutations.tsv', delim = '\t')
saveRDS(shared_mutations_SNV, file = 'shared_mutations_SNV.rds')
