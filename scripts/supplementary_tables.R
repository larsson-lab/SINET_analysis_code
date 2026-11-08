library(dplyr)
library(readr)

# --- Supplementary table 1 ####
source('scripts/analysis/extended_metadata.R')
purbayes = read_delim('input_data/WGS/purity.tsv', delim = '\t', col_types = c('c', 'c'))
colnames(purbayes) = c('Library', 'Purity PurBayes')

coverage = read_delim('input_data/WGS/output_idxstats.tsv', col_names = c('Library', 'mapped_reads'))
coverage = coverage %>% mutate(Average_sequencing_depth = (mapped_reads*150)/3100000000)

supp1 = left_join(metadata, purbayes, by = 'Library')
supp1 = inner_join(supp1, coverage, by = 'Library')
supp1 = supp1 %>% select(Library, Patient, Tumor, Clinical_type, Age, SNVs, Indels, `Purity PurBayes`, SNVs_purity_corrected, Average_sequencing_depth)
supp1 = supp1 %>% mutate(`Purity PurBayes` = round(`Purity PurBayes`, digits = 2))
supp1 = supp1 %>% mutate(Average_sequencing_depth = round(Average_sequencing_depth, digits = 2))
supp1 = supp1 %>% mutate(SNVs_purity_corrected = round(SNVs_purity_corrected, digits = 2))
supp1 = supp1 %>% mutate(Included = if_else(`Purity PurBayes` > 0.4, true = 'Yes', false = 'No'))
supp1 = supp1 %>% filter(Clinical_type != 'blood')
write_delim(supp1, 'supplementary_tables/supp1.tsv', delim = '\t')

# --- Supplementary table 2 ####
# simulated samples?

# --- Supplementary table 3 ####
library(stringr)
cgc = read.table('input_data/public_data/Census_allTue Oct 22 10_40_51 2024.tsv', sep = '\t', header = TRUE)
mutations = read_delim('input_data/WGS/mutations.tsv', delim = '\t') %>% mutate(VAF = (tumor_alt_read_depth_r1+tumor_alt_read_depth_r2)/(tumor_alt_read_depth_r1+tumor_alt_read_depth_r2+tumor_ref_read_depth_r1+tumor_ref_read_depth_r2))
colnames(mutations)[1] = 'Library'
mutations$Library = as.character(mutations$Library)

supp3 = mutations %>% 
  filter(str_detect(Consequence,
                    'upstream_gene_variant|missense_variant|3_prime_UTR_variant|5_prime_UTR_variant|stop_gained|frameshift_variant|splice')) %>% 
  filter(SYMBOL %in% cgc$Gene.Symbol) 
supp3 = inner_join(supp3, metadata, by = 'Library')
supp3 = supp3 %>% select(Patient, Tumor, SYMBOL, chr, pos, ref, alt, Consequence, Existing_variation, VAF)
colnames(supp3) = c('Patient', 'Tumor', 'Gene', 'Chromosome', 'Position', 'Reference base', 'Alternative base', 'Consequence', 'Existing variantion (dbSNP154)', 'VAF')
supp3 = supp3 %>% mutate(VAF = round(VAF, digits = 2))
write_delim(supp3, 'supplementary_tables/supp3.tsv', delim = '\t')

# --- Supplementary table 4 ####
source('scripts/analysis/figure2_figure3_read_data.R')
burden_all = burden_all %>% select(Sample, `Uncorrected burden`, `Uncorrected burden 95% lower`, `Uncorrected burden 95% upper`, mutation_rate, mutation_rate_scaled, B_normal)
stats_all = stats_all %>% select(Sample, `Effective Coverage`, `Number of Effective Read Families`, `Pass-filter Duplication Rate`, count)
supp4 = inner_join(burden_all, stats_all, by = 'Sample')
supp4 = inner_join(supp4, metadata_WGS, by = 'Sample')
supp4 = supp4 %>% 
  select(Sample, Patient, Age, Clinical_type, `Effective Coverage`, `Number of Effective Read Families`, `Pass-filter Duplication Rate`, count, `Uncorrected burden`, `Uncorrected burden 95% lower`, `Uncorrected burden 95% upper`, mutation_rate, mutation_rate_scaled, B_normal) %>% 
  rename(count = 'SNVs') %>% 
  rename(mutation_rate = 'Mutation rate (SNVs per diploid genome per year)') %>% 
  rename(mutation_rate_scaled = 'Adjusted mutation rate') %>% 
  rename(B_normal = 'Assumed normal tissue mutation rate') 
supp4 = supp4 %>% 
  mutate(`Mutation rate (SNVs per diploid genome per year)` = round(`Mutation rate (SNVs per diploid genome per year)`, digits = 2)) %>% 
  mutate(`Adjusted mutation rate` = round(`Adjusted mutation rate`, digits = 2)) %>% 
  mutate(`Assumed normal tissue mutation rate` = round(`Assumed normal tissue mutation rate`, digits = 2))
supp4 = supp4 %>% arrange(Clinical_type, Patient)
write_delim(supp4, 'supplementary_tables/supp4.tsv', delim = '\t')

# --- Supplementary table 5 ####
source('scripts/analysis/figure2_figure3_read_data.R')
fractions_scaled = left_join(fractions_scaled, mutation_count %>% select(Sample, recall), by = 'Sample') %>% 
  mutate(fraction_clonal_scaled_correction = fraction_clonal_scaled/recall)
supp5 = inner_join(fractions_scaled %>% select(Sample, mutations_tot, mutations_clonal, fraction_clonal, fraction_nonclonal, fraction_clonal_scaled_correction), burden_all, by = 'Sample') %>% arrange(Patient, Tumor)
supp5 = supp5 %>% mutate(Age_clonal_expansion = fraction_clonal_scaled_correction*Age) %>% filter(Clinical_type != 'normal')

supp5 = supp5 %>% select(Sample, Clinical_type, mutations_tot, mutations_clonal, fraction_clonal, fraction_nonclonal, Purity, fraction_clonal_scaled_correction, Age, Age_clonal_expansion) 
colnames(supp5) = c('Sample', 'Clinical type',	'Total NanoSeq SNV count',	'Clonal NanoSeq SNV count (overlap with WGS)',	'Fraction clonal',	'Fraction non-clonal',	'Purity PurBayes', 'Fraction clonal (adjusted for purity)',  'Age surgery', 'Age clonal expansion')
supp5 = supp5 %>% select(Sample, `Clinical type`, `Total NanoSeq SNV count`, `Clonal NanoSeq SNV count (overlap with WGS)`, `Fraction clonal`, `Fraction non-clonal`, `Purity PurBayes`, `Fraction clonal (adjusted for purity)`, `Age surgery`, `Age clonal expansion`)
write_delim(supp5, 'supplementary_tables/supp5.tsv', delim = '\t')

# --- Supplementary table 6 ####
source('scripts/analysis/figure2_figure3_read_data.R')
source('scripts/analysis/extended_metadata.R')
supp6 = mutation_count %>% arrange(Patient, Tumor) %>% select(Sample, Clinical_type, Age, ce_mean, ce_uci_90, ce_lci_90)
supp6 = inner_join(supp6, metadata %>% select(Sample, SNVs_purity_corrected, SNVs_WGD_purity_corrected, Purity))
supp6 = supp6 %>% mutate(SNVs_WGD_purity_corrected = if_else(SNVs_WGD_purity_corrected < SNVs_purity_corrected, true = SNVs_WGD_purity_corrected, false = NA))
supp6 = supp6 %>% select(Sample, Clinical_type, Age, SNVs_purity_corrected, SNVs_WGD_purity_corrected, ce_mean, ce_lci_90, ce_uci_90, Purity)
colnames(supp6) = c('Sample',	'Clinical type',	'Age',	'SNV burden', 'SNV burden (diploid equivalent)',	'Clonal expansion age (years)',	'90% CI lower age',	'90% CI higher age', 'Purity')
supp6 = supp6 %>% mutate(Included = if_else(Purity > 0.4, true = 'Yes', false = 'No'))
supp6 = supp6 %>% mutate(NanoSeq = if_else(Sample %in% supp5$Sample, true = 'Yes', false = 'No'))
write_delim(supp6, 'supplementary_tables/supp6.tsv', delim = '\t')

# --- Supplementary table 7 ####
source('scripts/analysis/extended_metadata.R')
colnames(WGD)[1] = 'Sample'
supp7 = metadata %>% dplyr::filter(Shared>0) %>% dplyr::select(Library, Sample, Shared, Shared_purity_corrected, Shared_WGD_purity_corrected)
supp7  = left_join(supp7, WGD, by = 'Sample')
supp7 = supp7 %>% 
  mutate(before_purity_corrected = (before/(before+after))*Shared_purity_corrected) %>% 
  mutate(after_purity_corrected = Shared_purity_corrected-before_purity_corrected) %>% 
  mutate(after_purity_WGD_corrected = after_purity_corrected/2)

supp7 = supp7 %>% 
  mutate(Age_WGD = before_purity_corrected/mean) %>% 
  mutate(Age_MRCA = Shared_WGD_purity_corrected/mean) %>% 
  mutate(Age_MRCA_lower = Shared_WGD_purity_corrected/(mean+1.64*sd_mean), Age_MRCA_higer = Shared_WGD_purity_corrected/(mean-1.64*sd_mean))

supp7 = supp7 %>% select(Sample, Shared_purity_corrected, Shared_WGD_purity_corrected, Age_MRCA, Age_MRCA_lower, Age_MRCA_higer, before_purity_corrected, after_purity_corrected, after_purity_WGD_corrected, Age_WGD)
colnames(supp7) = c('Primary', 'MRCA SNV burden (purity corrected)',	'MRCA SNV burden (purity corrected diplod equivalent)',	'Age MRCA (34.3 SNVs/year)',	'90% CI lower age',	'90% CI higher age', 'SNVs before WGD (purity corrected)',	'SNVs after WGD (purity corrected)',	'SNVs after WGD (purity corrected diplod equivalent)',	'Age WGD (34.3 SNVs/year)')
supp7 = supp7 %>% 
  mutate(`Age MRCA (34.3 SNVs/year)` = round(`Age MRCA (34.3 SNVs/year)`, digits = 2)) %>%
  mutate(`90% CI lower age` = round(`90% CI lower age`, digits = 2)) %>%
  mutate(`90% CI higher age` = round(`90% CI higher age`, digits = 2)) %>%
  mutate(`Age WGD (34.3 SNVs/year)` = round(`Age WGD (34.3 SNVs/year)`, digits = 2)) %>% 
  mutate(`MRCA SNV burden (purity corrected diplod equivalent)` = round(`MRCA SNV burden (purity corrected diplod equivalent)`, digits = 2)) %>% 
  mutate(`MRCA SNV burden (purity corrected)` = round(`MRCA SNV burden (purity corrected)`, digits = 2)) %>% 
  mutate(`SNVs before WGD (purity corrected)` = round(`SNVs before WGD (purity corrected)`, digits = 2)) %>% 
  mutate(`SNVs after WGD (purity corrected)` = round(`SNVs after WGD (purity corrected)`, digits = 2)) %>%
  mutate(`SNVs after WGD (purity corrected diplod equivalent)` = round(`SNVs after WGD (purity corrected diplod equivalent)`, digits = 2)) %>%
  mutate(`MRCA SNV burden (purity corrected diplod equivalent)` = if_else(`MRCA SNV burden (purity corrected diplod equivalent)` == `MRCA SNV burden (purity corrected)`, true = NA, false = `MRCA SNV burden (purity corrected diplod equivalent)`))
write_delim(supp7, 'supplementary_tables/supp7.tsv', delim = '\t')
  