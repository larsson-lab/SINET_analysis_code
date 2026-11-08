library(grid)
library(dplyr)
library(readr)

shared_mutations_SNV = readRDS('input_data/WGS/shared_mutations_SNV.rds')
metadata = read_delim('input_data/metadata.tsv')

heatmap = function(patient, mutations){
  # Create a color palette function from white to red with 100 levels
  palette = colorRampPalette(c("white", 'orange', "red"))(100)
  
  getcolor = function(value, min_value = 0, max_value = 4660) {
    # Normalize the value to be between 0 and 1
    normalized_value = (value - min_value) / (max_value - min_value)
    color = palette[ceiling(normalized_value * (length(palette) - 1)) + 1]
    return(color)
  }
  
  
  n_samples = nrow(mutations[[patient]])
  
  # Define the size of each square in cm
  square_size = 0.5
  
  # Calculate the dimensions of the PDF based on the grid size
  plot_width = (n_samples + 1) * square_size
  plot_height = (n_samples + 1) * square_size
  
  plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), axes =FALSE, xlab = "", ylab = "")
  # Create a PDF file
  path = paste('figures/supp_figures/Supplementary_Note_1/overlap_heatmap/heatmap_', patient, '.pdf', sep = '')
  pdf(path, width = plot_width / 2.54, height = plot_height / 2.54)
  
  
  for(i in 1:(n_samples)) {
    for(j in 1:(n_samples)) {
      grid.rect(x = unit(j * square_size, "cm"),
                y = unit((n_samples - i) * square_size, "cm"),
                width = unit(square_size, "cm"),
                height = unit(square_size, "cm"),
                just = c("left", "bottom"),
                gp = gpar(col = "white", fill = getcolor(mutations[[patient]][i,j]), lwd = 0.5))
      
      grid.text(label = as.character(mutations[[patient]][i,j]),
                x = unit(j * square_size + square_size / 2, "cm"),
                y = unit((n_samples - i) * square_size + square_size / 2, "cm"),
                just = "center",
                gp = gpar(col = "black", cex = 0.5))
    }
  }
  
  category_colors = c("normal" = "#C39C94", "liver_met" = "#F5B6D2", "lymph_met" = "#D379B0", "carcinosis" = "#A985AB", 'primary' = '#A4D7DD')
  i = 0
  for (j in 1:(n_samples)){
    sample = colnames(mutations[[patient]])[[j]]
    type = metadata[metadata$Sample==sample,'Clinical_type'][[1]]
    grid.rect(x = unit(j * square_size, "cm"),
              y = unit((n_samples - i) * square_size, "cm"),
              width = unit(square_size, "cm"),
              height = unit(square_size, "cm"),
              just = c("left", "bottom"),
              gp = gpar(col = "white", fill = category_colors[type], lwd = 0.5))
    
    grid.text(label = sample,
              x = unit(j * square_size + square_size / 2, "cm"),
              y = unit((n_samples - i) * square_size + square_size / 2, "cm"),
              just = "center",
              gp = gpar(col = "black", cex = 0.5))
  }
  
  
  j = 0
  for (i in 1:(n_samples)){
    sample = colnames(mutations[[patient]])[[i]]
    type = metadata[metadata$Sample==sample,'Clinical_type'][[1]]
    grid.rect(x = unit(j * square_size, "cm"),
              y = unit((n_samples - i) * square_size, "cm"),
              width = unit(square_size, "cm"),
              height = unit(square_size, "cm"),
              just = c("left", "bottom"),
              gp = gpar(col = "white", fill = category_colors[type], lwd = 0.5))
    
    grid.text(label = sample,
              x = unit(j * square_size + square_size / 2, "cm"),
              y = unit((n_samples - i) * square_size + square_size / 2, "cm"),
              just = "center",
              gp = gpar(col = "black", cex = 0.5))
  }
  
  # Close the PDF device
  dev.off()
  
}

for (patient in 1:16){
  heatmap(patient, shared_mutations_SNV)
}
