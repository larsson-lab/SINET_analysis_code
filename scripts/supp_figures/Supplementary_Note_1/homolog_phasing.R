# Note! To run this script you need the Germline and LOH output files from Varscan.
# To generate them, you have to follow step 1 in the README

library(dplyr)
library(readr)

# Load sample metadata -----------------------------------------------------------
load_metadata = function(){
  meta_data = read_delim('input_data/metadata.tsv')
  return(meta_data)
}

# Read snp and filter -----------------------------------------------------
read_snp = function(sample){
  print('Reading data...')
  sample_id = (meta_data %>% filter(Sample==as.character(sample)) %>% dplyr::select(Library))[[1]]
  file_germline = paste('vcf/varscan/', sample_id, '.snp.Germline', sep='')
  file_LOH = paste('vcf/varscan/', sample_id, '.snp.LOH', sep='')
  mutations_germline = read_tsv(file_germline)
  mutations_LOH = read_tsv(file_LOH)
  mutations = rbind(mutations_germline, mutations_LOH)
  
  mutations = sample_n(mutations, 0.2*nrow(mutations))
  chromosomes = c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY')
  mutations = mutations %>%
    dplyr::select(chrom, position, ref, var, tumor_reads1, tumor_reads2, normal_reads1, normal_reads2) %>% 
    mutate(tumor_VAF = tumor_reads2/(tumor_reads2+tumor_reads1)) %>% 
    mutate(normal_VAF = normal_reads2/(normal_reads2+normal_reads1)) %>% 
    filter(between((normal_reads1+normal_reads2), 20, 60) & between(normal_VAF, 0.3, 0.7)) 
  mutations = mutations %>% 
    filter(mutations$chrom %in% chromosomes)
  return(mutations)
}

# Create scatter plot -----------------------------------------------------
phasing_scatter = function(shared, sample_1, sample_2, chr){
  plot_height = 5
  plot_width = 5
  
  path = paste('figures/supp_figures/Supplementary_Note_1/homolog_phasing/chr_', chr, '_', sample_1, '_', sample_2, '.pdf', sep = '')
  pdf(path, width = plot_width / 2.54, height = plot_height / 2.54)
  
  par(mai = c(1, 1, 0.5, 0)/2.54)
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "")
  
  axis(1, cex.axis = 0.5, pos = c(0,1), lwd.ticks=0.5, lwd=0.5, 
       tick = TRUE, at = seq(0, 1, by = 0.5), tcl = 0.1, labels = c(0,0.5,1),
       padj = -3) # Bottom axis
  axis(2, cex.axis = 0.5, pos = c(0,1), lwd.ticks=0.5, lwd=0.5, 
       tick = TRUE, at = seq(0, 1, by = 0.5), tcl = 0.1, labels = c(0,0.5,1), 
       las = 1, hadj = 0)# Left axis
  
  par(new=TRUE)
  plot(shared$tumor_VAF.x, shared$tumor_VAF.y, xlim = c(0,1), ylim = c(0,1), pch=21, 
       bg = 'black', 
       col = 'black',
       asp = 1, cex = 0.25, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
       tcl = -0.1, axes=FALSE, main= paste('chr', chr, ' germline SNPs', sep = ''), xlab = '', ylab = '', 
       cex.main = 0.5, font.main = 1)
  title(xlab=paste('VAF', sample_1, sep = ' '), cex.lab = 0.5, line = 1)
  title(ylab=paste('VAF', sample_2, sep=' '), cex.lab =0.5, line = 1)
  dev.off()
}

# Find shared snps in CNV chromosome --------------------------------------------------------
shared_snps = function(snps_sample1, snps_sample2, chromosome){
  print('Find shared SNPs...')
  shared = inner_join(x=snps_sample1, y = snps_sample2, by = c('chrom', 'position', 'ref', 'var'))
  shared = shared %>% filter(shared$chrom == paste('chr', chromosome, sep = ''))
}

# Main - homolog phasing --------------------------------------------------
meta_data = load_metadata()

main = function(sample_1, sample_2, snps_sample1, snps_sample2, chromosome){
  shared = shared_snps(snps_sample1, snps_sample2, as.character(chromosome))
  phasing_scatter(shared, sample_1, sample_2, chromosome)
}

# Patient 1 --------------------------------------------------
snps_1A = read_snp('1A')
snps_1C = read_snp('1C')
snps_1G = read_snp('1G')
snps_1H = read_snp('1H')
snps_1I = read_snp('1I')
snps_1J = read_snp('1J')
snps_1K = read_snp('1K')

#Chr18
main('1C', '1A', snps_1C, snps_1A, 18)
rm(snps_1A)
main('1C', '1G', snps_1C, snps_1G, 18)
main('1C', '1H', snps_1C, snps_1H, 18)
main('1C', '1I', snps_1C, snps_1I, 18)
main('1C', '1J', snps_1C, snps_1J, 18)
main('1C', '1K', snps_1C, snps_1K, 18)

#Chr11
main('1C', '1G', snps_1C, snps_1G, 11)
main('1C', '1H', snps_1C, snps_1H, 11)
main('1C', '1I', snps_1C, snps_1I, 11)
main('1C', '1J', snps_1C, snps_1J, 11)
main('1C', '1K', snps_1C, snps_1K, 11)

#Chr13
main('1I', '1J', snps_1I, snps_1J, 13)
main('1K', '1J', snps_1K, snps_1J, 13)

rm(snps_1C, snps_1G, snps_1H, snps_1I, snps_1J, snps_1K)
gc()

# Patient 2 --------------------------------------------------
snps_2G = read_snp('2G')
snps_2H = read_snp('2H')
snps_2K = read_snp('2K')
snps_2L = read_snp('2L')

#Chr18
main('2G', '2H', snps_2G, snps_2H, 18)
main('2G', '2K', snps_2G, snps_2K, 18)
main('2G', '2L', snps_2G, snps_2L, 18)

rm(snps_2G, snps_2H, snps_2K, snps_2L)
gc()

# Patient 3 --------------------------------------------------
snps_3A = read_snp('3A')
snps_3B = read_snp('3B')
snps_3C = read_snp('3C')
snps_3D = read_snp('3D')

#Chr18
main('3B', '3D', snps_3B, snps_3D, 18)
main('3B', '3A', snps_3B, snps_3A, 18)
main('3B', '3C', snps_3B, snps_3C, 18)
main('3A', '3C', snps_3A, snps_3C, 18)
main('3D', '3C', snps_3D, snps_3C, 18) #before: main('3D', '3C', snps_3A, snps_3C, 18)
main('3D', '3A', snps_3D, snps_3A, 18) #before: main('3D', '3A', snps_3A, snps_3C, 18)
main('3A', '3B', snps_3A, snps_3B, 18) #before: main('3A', '3B', snps_3B, snps_3A, 18)

rm(snps_3A, snps_3B, snps_3C, snps_3D)
gc()

# Patient 4 --------------------------------------------------
snps_4A = read_snp('4A')
snps_4B = read_snp('4B')
snps_4D = read_snp('4D')

#Chr18
main('4B', '4D', snps_4B, snps_4D, 18)
main('4B', '4A', snps_4B, snps_4A, 18)
main('4D', '4A', snps_4D, snps_4A, 18)

#Chr17 focal
main('4B', '4D', snps_4B, snps_4D, 17)

rm(snps_4A, snps_4B, snps_4D)
gc()
# Patient 5 --------------------------------------------------
snps_5A = read_snp('5A')
snps_5B = read_snp('5B')
snps_5C = read_snp('5C')
snps_5D = read_snp('5D')
snps_5E = read_snp('5E')
snps_5F = read_snp('5F')

#Chr18
main('5E', '5A', snps_5E, snps_5A, 18)
main('5E', '5B', snps_5E, snps_5B, 18)
main('5E', '5C', snps_5E, snps_5C, 18)
main('5E', '5D', snps_5E, snps_5D, 18)
main('5E', '5F', snps_5E, snps_5F, 18)

rm(snps_5A,snps_5B,snps_5C,snps_5D,snps_5E,snps_5F)
gc()

# Patient 6 --------------------------------------------------
snps_6A = read_snp('6A')
snps_6B = read_snp('6B')
snps_6C = read_snp('6C')
snps_6D = read_snp('6D')

#Chr18
main('6A', '6C', snps_6A, snps_6C, 18)
main('6A', '6B', snps_6A, snps_6B, 18)
main('6B', '6D', snps_6B, snps_6D, 18)

rm(snps_6A,snps_6B,snps_6C,snps_6D)
gc()

# Patient 7 --------------------------------------------------
snps_7A = read_snp('7A')
snps_7B = read_snp('7B')
snps_7C = read_snp('7C')
snps_7D = read_snp('7D')

#Chr18
main('7A', '7D', snps_7A, snps_7D, 18)

#Chr3
main('7A', '7D', snps_7A, snps_7D, 3)

#Chr7
main('7A', '7D', snps_7A, snps_7D, 7)

#Chr8
main('7A', '7D', snps_7A, snps_7D, 8)

#Chr10
main('7A', '7D', snps_7A, snps_7D, 10)

#Chr14
main('7A', '7D', snps_7A, snps_7D, 14)

#Chr20
main('7A', '7D', snps_7A, snps_7D, 20)

#Chr22
main('7A', '7D', snps_7A, snps_7D, 22)

#Chr4
main('7A', '7D', snps_7A, snps_7D, 4)
main('7A', '7B', snps_7A, snps_7B, 4)
main('7B', '7C', snps_7B, snps_7C, 4)

#Chr5
main('7A', '7D', snps_7A, snps_7D, 5)
main('7A', '7B', snps_7A, snps_7B, 5)
main('7B', '7C', snps_7B, snps_7C, 5)

#Chr13
main('7B', '7C', snps_7B, snps_7C, 13)

#Chr14
main('7B', '7C', snps_7B, snps_7C, 14)
main('7A', '7B', snps_7A, snps_7B, 14)

rm(snps_7A,snps_7B,snps_7C,snps_7D)
gc()

# Patient 8 --------------------------------------------------
snps_8A = read_snp('8A')
snps_8C = read_snp('8C')
snps_8D = read_snp('8D')

#Chr18
main('8A', '8C', snps_8A, snps_8C, 18)
main('8A', '8D', snps_8A, snps_8D, 18)

#Chr10
main('8A', '8C', snps_8A, snps_8C, 10)

rm(snps_8A,snps_8C,snps_8D)
gc()

# Patient 9 --------------------------------------------------
snps_9A = read_snp('9A')
snps_9C = read_snp('9C')

main('9A', '9C', snps_9A, snps_9C, 14)
main('9A', '9C', snps_9A, snps_9C, 20)

# Patient 10 --------------------------------------------------
snps_10A = read_snp('10A')
snps_10B = read_snp('10B')
snps_10D = read_snp('10D')
snps_10E = read_snp('10E')

#Chr18
main('10A', '10D', snps_10A, snps_10D, 18)
main('10B', '10E', snps_10B, snps_10E, 18)

#Chr9
main('10A', '10D', snps_10A, snps_10D, 9)

main('10B', '10E', snps_10B, snps_10E, 20)
main('10B', '10E', snps_10B, snps_10E, 14)
main('10B', '10E', snps_10B, snps_10E, 15)
main('10B', '10E', snps_10B, snps_10E, 12)
main('10B', '10E', snps_10B, snps_10E, 7)
main('10B', '10E', snps_10B, snps_10E, 4)
main('10B', '10E', snps_10B, snps_10E, 5)

rm(snps_10A,snps_10B,snps_10D, snps_10E)
gc()

# Patient 11 --------------------------------------------------
snps_11A = read_snp('11A')
snps_11B = read_snp('11B')
snps_11C = read_snp('11C')
snps_11D = read_snp('11D')


main('11A', '11B', snps_11A, snps_11B, 18)

main('11B', '11C', snps_11B, snps_11C, 14)
main('11B', '11D', snps_11B, snps_11D, 14)

main('11B', '11D', snps_11B, snps_11D, 4)
main('11B', '11A', snps_11B, snps_11A, 4)

main('11B', '11D', snps_11B, snps_11D, 5)

main('11B', '11D', snps_11B, snps_11D, 9)
main('11B', '11A', snps_11B, snps_11A, 9)

main('11B', '11D', snps_11B, snps_11D, 20)

rm(snps_11A, snps_11B, snps_11C, snps_11D)
gc()

# Patient 12 --------------------------------------------------
snps_12B = read_snp('12B')
snps_12C = read_snp('12C')

main('12B', '12C', snps_12B, snps_12C, 6)
main('12B', '12C', snps_12B, snps_12C, 16)

rm(snps_12B, snps_12C)
gc()

# Patient 13 --------------------------------------------------
snps_13A = read_snp('13A')
snps_13B = read_snp('13B')
snps_13C = read_snp('13C')

#Chr9
main('13A', '13B', snps_13A, snps_13B, 9)
main('13A', '13C', snps_13A, snps_13C, 9)
main('13B', '13C', snps_13B, snps_13C, 9)

#Chr11
main('13A', '13B', snps_13A, snps_13B, 11)
main('13A', '13C', snps_13A, snps_13C, 11)
main('13B', '13C', snps_13B, snps_13C, 11)

#Chr13
main('13A', '13B', snps_13A, snps_13B, 13)
main('13A', '13C', snps_13A, snps_13C, 13)
main('13B', '13C', snps_13B, snps_13C, 13)

#Chr18
main('13A', '13B', snps_13A, snps_13B, 18)
main('13A', '13C', snps_13A, snps_13C, 18)
main('13B', '13C', snps_13B, snps_13C, 18)

#Chr19
main('13B', '13C', snps_13B, snps_13C, 19)

rm(snps_13A, snps_13B, snps_13C)
gc()

# Patient 14 --------------------------------------------------
snps_14A = read_snp('14A')
snps_14B = read_snp('14B')
snps_14C = read_snp('14C')

#Chr5
main('14A', '14B', snps_14A, snps_14B, 5)
main('14A', '14C', snps_14A, snps_14C, 5)
main('14B', '14C', snps_14B, snps_14C, 5)

#Chr9
main('14B', '14C', snps_14B, snps_14C, 9)

#Chr18
main('14A', '14B', snps_14A, snps_14B, 18)
main('14A', '14C', snps_14A, snps_14C, 18)
main('14B', '14C', snps_14B, snps_14C, 18)

rm(snps_14A, snps_14B, snps_14C)
gc()

# Patient 15 --------------------------------------------------
snps_15A = read_snp('15A')
snps_15B = read_snp('15B')

#Chr18
main('15A', '15B', snps_15A, snps_15B, 18)

rm(snps_15A, snps_15B)
gc()

# Patient 16 --------------------------------------------------
snps_16A = read_snp('16A')
snps_16B = read_snp('16B')

#Chr18
main('16A', '16B', snps_16A, snps_16B, 18)

rm(snps_16A, snps_16B)
gc()
