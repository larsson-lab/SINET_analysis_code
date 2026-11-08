metadata_WGS = read.delim('input_data/metadata.tsv', sep = '\t')
shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')
SNVs = read.delim('input_data/WGS/SNVs.tsv')

#---Shared mutations seeding primary tumors---#
metadata_WGS = metadata_WGS %>%  mutate(Shared = 0)
seeding = c('2G', '3B', '4B', '5E', '6A', '6B', '7A', '7B', '10A', 
            '10B', '11B', '15A', '16A')

seeding_special = c('1C', '8A', '12A', '13A', '14A')

count_shared = function(patient, sample){
  for (i in 1:ncol(shared_mutations_SNV[[patient]])){
    if (shared_mutations_SNV[[patient]][sample, i] > 10 & colnames(shared_mutations_SNV[[patient]])[i] != sample){
      paired = colnames(shared_mutations_SNV[[patient]])[i]
      shared = shared_mutations_SNV[[patient]][sample, i]
      return(shared)
    }
  }
}

for (sample in seeding){
  print(sample)
  patient = (metadata_WGS[metadata_WGS$Sample==sample, 'Patient'])[[1]]
  shared = count_shared(patient, sample)
  metadata_WGS[metadata_WGS$Sample==sample, 'Shared'] = shared 
}

#Seeding 12A
shared_12A_12B = inner_join(SNVs %>% filter(sample_id == '3362'), SNVs %>% filter(sample_id == '3361'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)
shared_12A_12C = inner_join(SNVs %>% filter(sample_id == '3362'), SNVs %>% filter(sample_id == '3363'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

trunk = nrow(inner_join(shared_12A_12B, shared_12A_12C))
metadata_WGS[metadata_WGS$Sample=='12A', 'Shared'] = trunk
rm(shared_12A_12B, shared_12A_12C, trunk)

#Seeding 13A
shared_13A_13B = inner_join(SNVs %>% filter(sample_id == '3529'), SNVs %>% filter(sample_id == '3527'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)
shared_13A_13C = inner_join(SNVs %>% filter(sample_id == '3529'), SNVs %>% filter(sample_id == '3528'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

trunk = nrow(inner_join(shared_13A_13B, shared_13A_13C))
metadata_WGS[metadata_WGS$Sample=='13A', 'Shared'] = trunk
rm(shared_13A_13B, shared_13A_13C, trunk)

#Seeding 14A
shared_14A_14B = inner_join(SNVs %>% filter(sample_id == '3580'), SNVs %>% filter(sample_id == '3568'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)
shared_14A_14C = inner_join(SNVs %>% filter(sample_id == '3580'), SNVs %>% filter(sample_id == '3579'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

trunk = nrow(inner_join(shared_14A_14B, shared_14A_14C))
metadata_WGS[metadata_WGS$Sample=='14A', 'Shared'] = trunk
rm(shared_14A_14B, shared_14A_14C, trunk)

#Seeding 8A
shared_8A_8C = inner_join(SNVs %>% filter(sample_id == '3564'), SNVs %>% filter(sample_id == '3562'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)
shared_8A_8D = inner_join(SNVs %>% filter(sample_id == '3564'), SNVs %>% filter(sample_id == '3567'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

trunk = nrow(inner_join(shared_8A_8C, shared_8A_8D))
metadata_WGS[metadata_WGS$Sample=='8A', 'Shared'] = trunk
rm(shared_8A_8C, shared_8A_8D, trunk)

#Seeding 1C
shared_1J_1K = inner_join(SNVs %>% filter(sample_id == '3864'), SNVs %>% filter(sample_id == '3865'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1J_1K_1I = inner_join(SNVs %>% filter(sample_id == '3863'), shared_1J_1K, by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1G_1H = inner_join(SNVs %>% filter(sample_id == '3861'), SNVs %>% filter(sample_id == '3862'), by = c('chr', 'pos', 'ref', 'alt')) %>% 
  dplyr::select(chr, pos, ref, alt)

shared_1J_1K_1I_1G_1H = inner_join(shared_1J_1K_1I, shared_1G_1H)

trunk = nrow(inner_join(shared_1J_1K_1I_1G_1H, SNVs %>% filter(sample_id == '3857')))
metadata_WGS[metadata_WGS$Sample=='1C', 'Shared'] = trunk

rm(shared_1J_1K, shared_1J_1K_1I, shared_1J_1K_1I_1G_1H, shared_1G_1H, trunk)
rm(patient, sample, seeding, seeding_special, shared, count_shared, shared_mutations_SNV)
