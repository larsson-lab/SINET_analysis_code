library(dplyr)
library(readr)

source('scripts/analysis/extended_metadata.R')
mutations = read.delim('input_data/WGS/mutations.tsv')
mutations = mutations %>% mutate(VAF = (tumor_alt_read_depth_r1+tumor_alt_read_depth_r2)/(tumor_alt_read_depth_r1+tumor_alt_read_depth_r2+tumor_ref_read_depth_r1+tumor_ref_read_depth_r2))
shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')

scatterplot = function(sample1, sample2, mutations){
  
  category_colors = c("normal" = "#C39C94", "liver_met" = "#F5B6D2", "lymph_met" = "#D379B0", "carcinosis" = "#A985AB", 'primary' = '#A4D7DD')
  
  sample1_id = (metadata %>% filter(metadata$Sample == sample1) %>% select(Library))[[1]]
  sample2_id = (metadata %>% filter(metadata$Sample == sample2) %>% select(Library))[[1]]
  
  sample1_type = metadata[metadata$Sample==sample1,'Clinical_type'][[1]]
  sample2_type = metadata[metadata$Sample==sample2,'Clinical_type'][[1]]
  
  sample1_VAF = mutations %>% filter(sample_id == sample1_id) %>% filter(chr %in% c(paste('chr', 1:22, sep = '')))
  sample2_VAF = mutations %>%  filter(sample_id == sample2_id) %>% filter(chr %in% c(paste('chr', 1:22, sep = '')))
  
  shared = full_join(sample1_VAF, sample2_VAF, by = c('chr', 'pos', 'ref', 'alt')) %>% select(VAF.x,VAF.y, chr, pos, SYMBOL.x, SYMBOL.y)
  shared$VAF.x[is.na(shared$VAF.x)] = 0
  shared$VAF.y[is.na(shared$VAF.y)] = 0
  shared = as.data.frame(shared)
  
  shared = shared %>% arrange(chr == 'chr18')
  plot_height = 5
  plot_width = 5
  
  path = paste('figures/supp_figures/Supplementary_Note_1/scatter/scatter_', sample1, '_', sample2, '.pdf', sep = '')
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
  if (metadata$chr18_loss[metadata$Sample == sample1] == 'FALSE' | metadata$chr18_loss[metadata$Sample == sample2] == 'FALSE'){
    plot(shared$VAF.x, shared$VAF.y, xlim = c(0,1), ylim = c(0,1), pch=21, 
         bg =case_when(shared$VAF.x == 0 ~ category_colors[sample2_type],
                       shared$VAF.y == 0 ~ category_colors[sample1_type],
                       shared$VAF.x>0 & shared$VAF.y>0 ~ '#9CCB88'), 
         col = case_when(shared$VAF.x == 0 ~ category_colors[sample2_type],
                         shared$VAF.y == 0 ~ category_colors[sample1_type], 
                         shared$VAF.x>0 & shared$VAF.y>0 ~ '#9CCB88'),
         asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
         tcl = -0.1, axes=FALSE, main= 'Somatic mutations', xlab = '', ylab = '', 
         cex.main = 0.5, font.main = 1)
  }
  
  if (metadata$chr18_loss[metadata$Sample == sample1] %in% c('LOSS', 'CN-LOH') | metadata$chr18_loss[metadata$Sample == sample2] %in% c('LOSS', 'CN-LOH')){
    
    plot(shared$VAF.x, shared$VAF.y, xlim = c(0,1), ylim = c(0,1), pch=21, 
         bg =case_when(shared$VAF.x == 0 ~ category_colors[sample2_type],
                       shared$VAF.y == 0 ~ category_colors[sample1_type],
                       shared$VAF.x>0 & shared$VAF.y>0 ~ '#9CCB88'), 
         col = ifelse(shared$chr == 'chr18', '#2AB8CB', case_when(shared$VAF.x == 0 ~ category_colors[sample2_type],
                                                                  shared$VAF.y == 0 ~ category_colors[sample1_type], 
                                                                  shared$VAF.x>0 & shared$VAF.y>0 ~ '#9CCB88')),
         asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
         tcl = -0.1, axes=FALSE, main= 'Somatic mutations', xlab = '', ylab = '', 
         cex.main = 0.5, font.main = 1)
  }
  
  par(new=TRUE)
  plot(shared$VAF.x[shared$SYMBOL.x == 'CDKN1B'], shared$VAF.y[shared$SYMBOL.x == 'CDKN1B'], xlim = c(0,1), ylim = c(0,1), pch=21, 
       bg = 'black', 
       col = 'black',
       asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
       tcl = -0.1, axes=FALSE, xlab = '', ylab = '', 
       cex.main = 0.5, font.main = 1)
  
  par(new=TRUE)
  plot(shared$VAF.x[shared$SYMBOL.x=='TNRC6B'], shared$VAF.y[shared$SYMBOL.x=='TNRC6B'], xlim = c(0,1), ylim = c(0,1), pch=21, 
       bg = 'black', 
       col = 'black',
       asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
       tcl = -0.1, axes=FALSE, xlab = '', ylab = '', 
       cex.main = 0.5, font.main = 1)
  
  par(new=TRUE)
  plot(shared$VAF.x[shared$SYMBOL.y == 'CDKN1B'], shared$VAF.y[shared$SYMBOL.y == 'CDKN1B'], xlim = c(0,1), ylim = c(0,1), pch=21, 
       bg = 'black', 
       col = 'black',
       asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
       tcl = -0.1, axes=FALSE, xlab = '', ylab = '', 
       cex.main = 0.5, font.main = 1)
  
  par(new=TRUE)
  plot(shared$VAF.x[shared$SYMBOL.y=='TNRC6B'], shared$VAF.y[shared$SYMBOL.y=='TNRC6B'], xlim = c(0,1), ylim = c(0,1), pch=21, 
       bg = 'black', 
       col = 'black',
       asp = 1, cex = 0.5, lwd = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.5, 
       tcl = -0.1, axes=FALSE, xlab = '', ylab = '', 
       cex.main = 0.5, font.main = 1)
  
  title(xlab=paste('VAF', sample1, sep = ' '), cex.lab = 0.5, line = 1)
  title(ylab=paste('VAF', sample2, sep=' '), cex.lab =0.5, line = 1)
  
  dev.off()
}


run_scatter = function(patient, sample){
  for (i in 1:ncol(shared_mutations_SNV[[patient]])){
    print(sample)
    if (shared_mutations_SNV[[patient]][sample, i] > 10 & colnames(shared_mutations_SNV[[patient]])[i] != sample){
      paired = colnames(shared_mutations_SNV[[patient]])[i]
      print(paired)
      scatterplot(sample, paired, mutations)
    }
  }
}

samples = (metadata %>%  filter(Clinical_type != 'blood') %>% select(Sample))[[1]]

for (sample in samples){
  patient = (metadata[metadata$Sample==sample, 'Patient'])[[1]]
  run_scatter(patient, sample)
}
