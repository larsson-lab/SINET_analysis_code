## Scaling mutation counts in case with WGD ####
set.seed(123)

## Load data and define functions ####
mutations = read_delim('input_data/WGS/SNVs.tsv', delim = '\t') %>%  
  mutate(VAF = (tumor_alt_read_depth_r1 + tumor_alt_read_depth_r2)/(tumor_alt_read_depth_r1 + tumor_alt_read_depth_r2+tumor_ref_read_depth_r1 + tumor_ref_read_depth_r2)) %>% 
  filter(chr != 'chr18')
SNVs = read_delim('input_data/WGS/SNVs.tsv', delim = '\t')
metadata = read_delim('input_data/metadata.tsv')
shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')
shared = function(sample1, sample2, mutations){
  
  sample1_id = (metadata %>% filter(metadata$Sample == sample1) %>% dplyr::select(Library))[[1]]
  sample2_id = (metadata %>% filter(metadata$Sample == sample2) %>% dplyr::select(Library))[[1]]
  
  sample1_VAF = mutations %>% filter(sample_id == sample1_id) %>% filter(chr %in% c(paste('chr', 1:22, sep = '')))
  sample2_VAF = mutations %>%  filter(sample_id == sample2_id) %>% filter(chr %in% c(paste('chr', 1:22, sep = '')))
  
  shared = inner_join(sample1_VAF, sample2_VAF, by = c('chr', 'pos', 'ref', 'alt')) %>% dplyr::select(VAF.x,VAF.y)
  shared = as.data.frame(shared)
}

## WGD timing patient 3 ####
shared_3B_3D = shared('3B', '3D', mutations)
clusters_3B_3D = kmeans(as.matrix(shared_3B_3D), centers = 2, iter.max = 100, nstart = 10)

shared_3B_3D = cbind(shared_3B_3D, clusters_3B_3D$cluster)
colnames(shared_3B_3D)[3] = 'cluster'

## WGD timing patient 5 ####
shared_5E_5F = shared('5E', '5F', mutations)
clusters_5E_5F = kmeans(as.matrix(shared_5E_5F), centers = 2, iter.max = 100, nstart = 10)

shared_5E_5F = cbind(shared_5E_5F, clusters_5E_5F$cluster)
colnames(shared_5E_5F)[3] = 'cluster'

## WGD timing patient 6 ####
shared_6A_6C = shared('6A', '6C', mutations)
clusters_6A_6C = kmeans(as.matrix(shared_6A_6C), centers = 2, iter.max = 1000, nstart = 10)

shared_6A_6C = cbind(shared_6A_6C, clusters_6A_6C$cluster)
colnames(shared_6A_6C)[3] = 'cluster'

shared_6B_6D = shared('6B', '6D', mutations)
clusters_6B_6D = kmeans(as.matrix(shared_6B_6D), centers = 2, iter.max = 100, nstart = 10)

shared_6B_6D = cbind(shared_6B_6D, clusters_6B_6D$cluster)
colnames(shared_6B_6D)[3] = 'cluster'

## WGD timing patient 16 ####
shared_16A_16B = shared('16A', '16B', mutations)
clusters_16A_16B = kmeans(as.matrix(shared_16A_16B), centers = 2, iter.max = 100, nstart = 10)

shared_16A_16B = cbind(shared_16A_16B, clusters_16A_16B$cluster)
colnames(shared_16A_16B)[3] = 'cluster'

## Time WGD ####
time_before_3B_3D = clusters_3B_3D$size[2]/(clusters_3B_3D$size[1] + clusters_3B_3D$size[2])
time_before_5E_5F = clusters_5E_5F$size[1]/(clusters_5E_5F$size[1] + clusters_5E_5F$size[2])
time_before_6A_6C = clusters_6A_6C$size[2]/(clusters_6A_6C$size[1] + clusters_6A_6C$size[2])
time_before_6B_6D = clusters_6B_6D$size[1]/(clusters_6B_6D$size[1] + clusters_6B_6D$size[2])
time_before_16A_16B = clusters_16A_16B$size[2]/(clusters_16A_16B$size[1] + clusters_16A_16B$size[2])

before = c(clusters_3B_3D$size[2],clusters_5E_5F$size[1], clusters_6A_6C$size[2],clusters_6B_6D$size[1],clusters_16A_16B$size[2])
after = c(clusters_3B_3D$size[1], clusters_5E_5F$size[2], clusters_6A_6C$size[1],clusters_6B_6D$size[2],clusters_16A_16B$size[1])
primary = c('3B', '5E', '6A', '6B', '16A')

WGD = data.frame(primary, before, after)

WGD$shared[WGD$primary=='3B'] = shared_mutations_SNV[[3]]['3B','3D']
WGD$shared[WGD$primary=='5E'] = shared_mutations_SNV[[5]]['5E','5F']
WGD$shared[WGD$primary=='6A'] = shared_mutations_SNV[[6]]['6A','6C']
WGD$shared[WGD$primary=='6B'] = shared_mutations_SNV[[6]]['6B','6D']
WGD$shared[WGD$primary=='16A'] = shared_mutations_SNV[[16]]['16A','16B']

WGD = WGD %>% mutate(before_scaled = round((before/(before+after))*shared)) %>% mutate(after_scaled = round((after/(before+after))*shared))
WGD = WGD %>% dplyr::select(primary, before_scaled, after_scaled)
colnames(WGD) = c('primary', 'before', 'after')

write_delim(WGD, 'input_data/WGS/WGD.tsv', delim = '\t')
