#---Load data ####
source('scripts/analysis/phylogenies_trunk.R')
metadata = metadata_WGS
mutations = read.delim('input_data/WGS/mutations.tsv', sep = '\t') 
SNVs = read.delim('input_data/WGS/SNVs.tsv', sep = '\t')

metadata = local({
  indels = anti_join(mutations, SNVs)
  metadata = metadata %>%  mutate('SNVs' = 0, 'Indels' = 0)
  for (sample in metadata$Library){
    metadata[metadata$Library==sample, 'SNVs'] = nrow(SNVs %>% filter(sample_id==sample))
    metadata[metadata$Library==sample, 'Indels'] = nrow(indels %>% filter(sample_id==sample))
  }
  metadata
})

#---Putative driver mutations---####
metadata = local({
  CDKN1B = (mutations %>% filter(SYMBOL=='CDKN1B') %>% filter(Consequence != 'intron_variant') %>% dplyr::select(sample_id))[[1]]
  TNRC6B = (mutations %>%  filter(SYMBOL == 'TNRC6B') %>% filter(Consequence != 'intron_variant') %>% dplyr::select(sample_id))[[1]]
  
  metadata = metadata %>% 
    mutate(CDKN1B = ifelse(Library %in% CDKN1B, TRUE, FALSE))
  
  metadata = metadata %>% 
    mutate(TNRC6B = ifelse(Library %in% TNRC6B, TRUE, FALSE))
  metadata
})

#---CNVs---####
chr18 = c('1A', '1C', '1F', '1G', '1H', '1I', '1J', '1K', '2G', '2H', '2K', '2L',
          '3A', '3B', '3C', '3D', '4B', '4D', '5A', '5B', '5C', '5D',
          '5E', '5F', '6A', '6B', '6C', '6D', '8A', '8C', '8D', '10A', '10D', '11A',
          '13A', '13B', '13C', '14A', '14B', '14C', '15A', '15B',
          '16A', '16B')
chr4 = c('7A', '7B', '7C', '7D', '9A', '10B', '10E', '11A', '11B', '11D')
chr5 = c('7A', '7B', '7C', '7D', '9C', '10B', '10E', '11B', '11D', '14A',
         '14B', '14C')
chr14 = c('7A', '7B', '7C', '7D', '9A', '9C', '10B', '10E', '11B', '11C', '11D', '12A') 
chr20 = c('7A', '7B', '7D', '9A', '9C', '10B', '10E', '11B', '11D', '14C')

metadata = metadata %>% mutate(chr18_loss = ifelse(metadata$Sample %in% chr18, 'LOSS', 'FALSE'))
metadata = metadata %>% mutate(chr4_gain = ifelse(metadata$Sample %in% chr4, 'GAIN', 'FALSE'))
metadata = metadata %>% mutate(chr5_gain = ifelse(metadata$Sample %in% chr5, 'GAIN', 'FALSE'))
metadata = metadata %>% mutate(chr14_gain = ifelse(metadata$Sample %in% chr14, 'GAIN', 'FALSE'))
metadata = metadata %>% mutate(chr20_gain = ifelse(metadata$Sample %in% chr20, 'GAIN', 'FALSE'))

metadata$chr18_loss[which(metadata$Sample == '7A')] = 'CN-LOH'
metadata$chr18_loss[which(metadata$Sample == '7D')] = 'CN-LOH'

{
  metadata = metadata %>% mutate(Homolog='FALSE', Homolog4 = 'FALSE', Homolog5 = 'FALSE', Homolog14 = 'FALSE', Homolog20 = 'FALSE')
  metadata$Homolog[which(metadata$Sample %in% c('1C', '1G', '1H', '1I', '1J', '1K', '2G', '2L', '3A','3C', '4B', '4D',
                                        '5B', '5E', '5F', '6A', '6B', '6C', '6D', '7A', '7D', '10A', '10D', '13A', '13B', '13C',
                                        '14A', '14B', '14C', '15A', '15B', '16A', '16B', '8A', '8C', '8D'))] = 'A'
  metadata$Homolog[which(metadata$Sample %in% c('1A', '2H', '2K', '3D', '5A', '5C', '5D'))] = 'B'
  metadata$Homolog[which(metadata$Sample %in% c('3B'))] = 'biallelic'
  
  metadata$Homolog4[which(metadata$Sample %in% c('7B', '7C', '10B', '10E', '11B', '11D'))] = 'A'
  metadata$Homolog4[which(metadata$Sample %in% c('11A'))] = 'B'
  metadata$Homolog4[which(metadata$Sample %in% c('7A', '7D'))] = 'biallelic'
  
  metadata$Homolog5[which(metadata$Sample %in% c('7B', '7C', '10B', '10E', '11B', '11D','14A', '14C'))] = 'A'
  metadata$Homolog5[which(metadata$Sample %in% c('14B'))] = 'B'
  metadata$Homolog5[which(metadata$Sample %in% c('7A', '7D'))] = 'biallelic'
  
  metadata$Homolog14[which(metadata$Sample %in% c('7B', '7C', '9A', '9C', '10B', '10E', '11B', '11D'))] = 'A'
  metadata$Homolog14[which(metadata$Sample %in% c('11C', '7A', '7D'))] = 'B'
  
  metadata$Homolog20[which(metadata$Sample %in% c('9A', '10B', '10E', '11B', '11D'))] = 'A'
  metadata$Homolog20[which(metadata$Sample %in% c('7A', '7D', '9C'))] = 'B'
}
rm(chr14, chr18, chr20, chr4, chr5, mutations, SNVs)

# -- Adjust counts for purity and WGD --- ####

## Purity ####
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

purity = read_delim('input_data/WGS/purity.tsv') %>% mutate(Sample = as.character(Sample))
colnames(purity)[1] = 'Library'

metadata = left_join(metadata, purity, by = 'Library')
metadata$SNVs = as.numeric(metadata$SNVs)

# Calculate r
metadata = metadata %>% 
  mutate(correction = get_correction_factor(Purity)) %>% 
  mutate(SNVs_purity_corrected = SNVs*correction)

## Adjust for WGD ####
WGD = read_delim('input_data/WGS/WGD.tsv', delim = '\t')

metadata = metadata %>% mutate(SNVs_WGD_corrected = SNVs)

metadata[metadata$Sample == '3B', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '3B', 'SNVs'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]
metadata[metadata$Sample == '3D', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '3D', 'SNVs'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]

metadata[metadata$Sample == '5E', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '5E', 'SNVs'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]
metadata[metadata$Sample == '5F', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '5F', 'SNVs'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]

metadata[metadata$Sample == '6A', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '6A', 'SNVs'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]
metadata[metadata$Sample == '6C', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '6C', 'SNVs'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]

metadata[metadata$Sample == '6B', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '6B', 'SNVs'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]
metadata[metadata$Sample == '6D', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '6D', 'SNVs'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]

metadata[metadata$Sample == '16A', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '16A', 'SNVs'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]
metadata[metadata$Sample == '16B', 'SNVs_WGD_corrected'] = (((metadata[metadata$Sample == '16B', 'SNVs'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]

metadata = metadata %>% mutate(SNVs_WGD_purity_corrected = SNVs_WGD_corrected*correction)

## Adjust MRCA ####
correction_1C = metadata[metadata$Sample == '1C', 'correction']
correction_1G = metadata[metadata$Sample == '1G', 'correction']
correction_1H = metadata[metadata$Sample == '1H', 'correction']
correction_1I = metadata[metadata$Sample == '1I', 'correction']
correction_1J = metadata[metadata$Sample == '1J', 'correction']
correction_1K = metadata[metadata$Sample == '1K', 'correction']
metadata = metadata %>% mutate(Shared_purity_corrected = case_when(Sample == '1C' ~ Shared*correction_1C*correction_1G*correction_1H*correction_1I*correction_1J*correction_1K,
                                                            Sample == '2G' ~ Shared*metadata[metadata$Sample == '2G', 'correction']*metadata[metadata$Sample == '2L', 'correction'],
                                                            Sample == '3B' ~ Shared*metadata[metadata$Sample == '3B', 'correction']*metadata[metadata$Sample == '3D', 'correction'],
                                                            Sample == '4B' ~ Shared*metadata[metadata$Sample == '4B', 'correction']*metadata[metadata$Sample == '4D', 'correction'],
                                                            Sample == '5E' ~ Shared*metadata[metadata$Sample == '5E', 'correction']*metadata[metadata$Sample == '5F', 'correction'],
                                                            Sample == '6A' ~ Shared*metadata[metadata$Sample == '6A', 'correction']*metadata[metadata$Sample == '6C', 'correction'],
                                                            Sample == '6B' ~ Shared*metadata[metadata$Sample == '6B', 'correction']*metadata[metadata$Sample == '6D', 'correction'],
                                                            Sample == '7A' ~ Shared*metadata[metadata$Sample == '7A', 'correction']*metadata[metadata$Sample == '7D', 'correction'],
                                                            Sample == '7B' ~ Shared*metadata[metadata$Sample == '7B', 'correction']*metadata[metadata$Sample == '7C', 'correction'],
                                                            Sample == '8A' ~ Shared*metadata[metadata$Sample == '8A', 'correction']*metadata[metadata$Sample == '8C', 'correction']*metadata[metadata$Sample == '8D', 'correction'],
                                                            Sample == '10A' ~ Shared*metadata[metadata$Sample == '10A', 'correction']*metadata[metadata$Sample == '10D', 'correction'],
                                                            Sample == '10B' ~ Shared*metadata[metadata$Sample == '10B', 'correction']*metadata[metadata$Sample == '10E', 'correction'],
                                                            Sample == '11B' ~ Shared*metadata[metadata$Sample == '11B', 'correction']*metadata[metadata$Sample == '11D', 'correction'],
                                                            Sample == '12A' ~ Shared*metadata[metadata$Sample == '12A', 'correction']*metadata[metadata$Sample == '12B', 'correction']*metadata[metadata$Sample == '12C', 'correction'],
                                                            Sample == '13A' ~ Shared*metadata[metadata$Sample == '13A', 'correction']*metadata[metadata$Sample == '13B', 'correction']*metadata[metadata$Sample == '13C', 'correction'],
                                                            Sample == '14A' ~ Shared*metadata[metadata$Sample == '14A', 'correction']*metadata[metadata$Sample == '14B', 'correction']*metadata[metadata$Sample == '14C', 'correction'],
                                                            Sample == '15A' ~ Shared*metadata[metadata$Sample == '15A', 'correction']*metadata[metadata$Sample == '15B', 'correction'],
                                                            Sample == '16A' ~ Shared*metadata[metadata$Sample == '16A', 'correction']*metadata[metadata$Sample == '16B', 'correction']
))


metadata = metadata %>% mutate(Shared_WGD_corrected = Shared) 

metadata[metadata$Sample == '3B', 'Shared_WGD_corrected'] = (((metadata[metadata$Sample == '3B', 'Shared'])[[1]] - (WGD[WGD$primary == '3B', 'before'])[[1]])/2) + (WGD[WGD$primary == '3B', 'before'])[[1]]
metadata[metadata$Sample == '5E', 'Shared_WGD_corrected'] = (((metadata[metadata$Sample == '5E', 'Shared'])[[1]] - (WGD[WGD$primary == '5E', 'before'])[[1]])/2) + (WGD[WGD$primary == '5E', 'before'])[[1]]
metadata[metadata$Sample == '6A', 'Shared_WGD_corrected'] = (((metadata[metadata$Sample == '6A', 'Shared'])[[1]] - (WGD[WGD$primary == '6A', 'before'])[[1]])/2) + (WGD[WGD$primary == '6A', 'before'])[[1]]
metadata[metadata$Sample == '6B', 'Shared_WGD_corrected'] = (((metadata[metadata$Sample == '6B', 'Shared'])[[1]] - (WGD[WGD$primary == '6B', 'before'])[[1]])/2) + (WGD[WGD$primary == '6B', 'before'])[[1]]
metadata[metadata$Sample == '16A', 'Shared_WGD_corrected'] = (((metadata[metadata$Sample == '16A', 'Shared'])[[1]] - (WGD[WGD$primary == '16A', 'before'])[[1]])/2) + (WGD[WGD$primary == '16A', 'before'])[[1]]

metadata = metadata %>% mutate(Shared_WGD_purity_corrected = case_when(Sample == '1C' ~ Shared_WGD_corrected*correction_1C*correction_1G*correction_1H*correction_1I*correction_1J*correction_1K,
                                                  Sample == '2G' ~ Shared_WGD_corrected*metadata[metadata$Sample == '2G', 'correction']*metadata[metadata$Sample == '2L', 'correction'],
                                                  Sample == '3B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '3B', 'correction']*metadata[metadata$Sample == '3D', 'correction'],
                                                  Sample == '4B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '4B', 'correction']*metadata[metadata$Sample == '4D', 'correction'],
                                                  Sample == '5E' ~ Shared_WGD_corrected*metadata[metadata$Sample == '5E', 'correction']*metadata[metadata$Sample == '5F', 'correction'],
                                                  Sample == '6A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '6A', 'correction']*metadata[metadata$Sample == '6C', 'correction'],
                                                  Sample == '6B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '6B', 'correction']*metadata[metadata$Sample == '6D', 'correction'],
                                                  Sample == '7A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '7A', 'correction']*metadata[metadata$Sample == '7D', 'correction'],
                                                  Sample == '7B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '7B', 'correction']*metadata[metadata$Sample == '7C', 'correction'],
                                                  Sample == '8A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '8A', 'correction']*metadata[metadata$Sample == '8C', 'correction']*metadata[metadata$Sample == '8D', 'correction'],
                                                  Sample == '10A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '10A', 'correction']*metadata[metadata$Sample == '10D', 'correction'],
                                                  Sample == '10B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '10B', 'correction']*metadata[metadata$Sample == '10E', 'correction'],
                                                  Sample == '11B' ~ Shared_WGD_corrected*metadata[metadata$Sample == '11B', 'correction']*metadata[metadata$Sample == '11D', 'correction'],
                                                  Sample == '12A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '12A', 'correction']*metadata[metadata$Sample == '12B', 'correction']*metadata[metadata$Sample == '12C', 'correction'],
                                                  Sample == '13A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '13A', 'correction']*metadata[metadata$Sample == '13B', 'correction']*metadata[metadata$Sample == '13C', 'correction'],
                                                  Sample == '14A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '14A', 'correction']*metadata[metadata$Sample == '14B', 'correction']*metadata[metadata$Sample == '14C', 'correction'],
                                                  Sample == '15A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '15A', 'correction']*metadata[metadata$Sample == '15B', 'correction'],
                                                  Sample == '16A' ~ Shared_WGD_corrected*metadata[metadata$Sample == '16A', 'correction']*metadata[metadata$Sample == '16B', 'correction']
))
