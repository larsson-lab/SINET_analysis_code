SNVs = read.delim('input_data/WGS/SNVs.tsv')
metadata = read.delim('input_data/metadata.tsv')

library(dplyr)
library(PurBayes)
library(readr)

purity_all_samples = c()
samples = metadata %>% filter(Tumor!= '-') %>% filter(Clinical_type != 'normal') 
for (row in 1:1:nrow(samples)){
  id = as.character(samples[row,1])
  print(id)
  mutations_sample = SNVs %>% filter(sample_id == id)
  alt = mutations_sample[['tumor_alt_read_depth_r1']] + mutations_sample[['tumor_alt_read_depth_r2']]
  tot =  mutations_sample[['tumor_alt_read_depth_r1']] + mutations_sample[['tumor_alt_read_depth_r2']] +
    mutations_sample[['tumor_ref_read_depth_r1']] + mutations_sample[['tumor_ref_read_depth_r2']]
  
  purity = PurBayes(tot, alt, M=NULL, Z=NULL, pop.max=5, prior=NULL, burn.in=50000, n.post = 10000, fn.jags = "PB.jags", plot = FALSE)
  summary = summary(purity)
  purity_median = summary$purity["Median"]
  purity_median = purity_median[[1]]
  
  purity_all_samples = rbind(purity_all_samples, c(id, purity_median))
  
}

colnames(purity_all_samples) = c('Sample', 'Purity')
write.table(purity_all_samples, 'input_data/WGS/purity.tsv', quote = FALSE, row.names = FALSE, sep = '\t')