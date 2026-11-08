library(readr)
library(dplyr)

input = read_delim('/home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/input_261102.txt', delim = '\t')

###BOOTSTRAP AND CREATE INPUT TO MATRIX GENERATOR###
for (i in 1:500){
  print(paste('Starting round', i, sep=' '))
  random_mutations = data.frame(matrix(ncol = 11, nrow = 0))
  colnames(random_mutations) = c("Project", "Sample", "ID", "Genome", "mut_type", "chrom", "pos_start", "pos_end", "ref", "alt", "Type")
  
  print('Sampling...')
  for (sample in unique(input$Sample)){
    mutations = input %>% filter(Sample==sample)
    colnames(mutations) = c("Project", "Sample", "ID", "Genome", "mut_type", "chrom", "pos_start", "pos_end", "ref", "alt", "Type")
    mutations = mutations[sample(nrow(mutations),nrow(mutations),replace=TRUE),]
    random_mutations = rbind(random_mutations, mutations)
  }
  print('Writing...')
  path = paste('/home/alice/Timing_publication_remake/code_repository/results/matrixgenerator/bootstrap/run', i, '/input.txt', sep = '')
  write.table(random_mutations, path, quote = FALSE, sep = '\t', row.names = FALSE)
}
