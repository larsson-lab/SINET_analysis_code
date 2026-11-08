library(dplyr)
library(readr)

# Read SINET WGS ####
SNV_WGS = read_delim("input_data/WGS/SNVs.tsv", delim = '\t', guess_max = 130000)
metadata_WGS = read_delim('input_data/metadata.tsv') %>% filter(!Clinical_type %in% c('normal', 'blood'))
colnames(SNV_WGS)[1]='Library'
SNV_WGS$Library = as.character(SNV_WGS$Library)
SNV_WGS = SNV_WGS %>% inner_join(metadata_WGS, by = 'Library')

SNV_WGS = SNV_WGS %>% select(Library, chr, pos, ref, alt)
colnames(SNV_WGS) = c('Sample', 'chrom', 'pos_start', 'ref', 'alt')
SNV_WGS = SNV_WGS %>%  
  mutate(Project = 'SINET_WGS', ID = '.', Genome = 'GRCh38', mut_type = 'SNP', pos_end = pos_start, Type = 'SOMATIC') %>% 
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type) 
SNV_WGS = SNV_WGS %>%  mutate(chrom = substr(chrom, 4, 6))

# Read NanoSeq ####
metadata_nano = read_delim('input_data/metadata_libraries_SINET.txt')

SNV_nano = data.frame()
samples = metadata_nano$sample_name

for (sample in samples){
  path = paste('input_data/NanoSeq/', sample,  '_snv.vcf', sep = '')
  muts = vcfR::read.vcfR(path)
  muts = as.data.frame(muts@fix)
  muts = muts %>% mutate(sample_name = sample) %>% filter(FILTER == 'PASS')
  SNV_nano = rbind(SNV_nano, muts)
}
rm(muts, path)
SNV_nano_tissue = SNV_nano

SNV_nano = SNV_nano %>% select(sample_name, CHROM, POS, REF, ALT)
colnames(SNV_nano) = c('Sample', 'chrom', 'pos_start', 'ref', 'alt')
SNV_nano = SNV_nano %>% 
  mutate(Project = 'nanoseq', ID = '.', Genome = 'GRCh38', mut_type = 'SNP', pos_end = pos_start, Type = 'SOMATIC') %>% 
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)
SNV_nano = SNV_nano %>% mutate(chrom = substr(chrom, 4, 6))

# Clonal and non clonal mutations ####
SNV_clonal_nonclonal = read_delim('input_data/clonal_nonclonal.tsv', delim = '\t')
colnames(SNV_clonal_nonclonal) = c('Sample', 'chrom', 'pos_start', 'ref', 'alt')
SNV_clonal_nonclonal = SNV_clonal_nonclonal %>% 
  mutate(Project = 'nanoseq', ID = '.', Genome = 'GRCh38', mut_type = 'SNP', pos_end = pos_start, Type = 'SOMATIC') %>% 
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)
SNV_clonal_nonclonal = SNV_clonal_nonclonal %>% mutate(chrom = substr(chrom, 4, 6))

# Blood, SINET, normal ####
SNV_nano_tissue = inner_join(SNV_nano_tissue, metadata_nano, by = 'sample_name')
SNV_nano_tissue = SNV_nano_tissue %>% select(Clinical_type, CHROM, POS, REF, ALT)
colnames(SNV_nano_tissue) = c('Sample', 'chrom', 'pos_start', 'ref', 'alt')
SNV_nano_tissue = SNV_nano_tissue %>% 
  mutate(Project = 'nanoseq', ID = '.', Genome = 'GRCh38', mut_type = 'SNP', pos_end = pos_start, Type = 'SOMATIC') %>% 
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)
SNV_nano_tissue = SNV_nano_tissue %>% mutate(chrom = substr(chrom, 4, 6))
SNV_nano_tissue = SNV_nano_tissue %>% mutate(Sample = if_else(Sample %in% c("liver_met", "lymph_met", "carcinosis", "primary"), true = 'SI-NET', false = Sample))

input_matrixgenerator = rbind(SNV_WGS, SNV_nano, SNV_clonal_nonclonal, SNV_nano_tissue)

path = 'results/matrixgenerator/input_261102.txt'
write.table(input_matrixgenerator, path, quote = FALSE, sep = '\t', row.names = FALSE)
