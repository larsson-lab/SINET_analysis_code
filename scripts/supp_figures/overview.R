library(dplyr)
library(readr)

source('scripts/analysis/extended_metadata.R')
source('scripts/analysis/phylogenies_trunk.R')
metadata = inner_join(metadata_WGS, metadata, by = c("Library", "Patient", "Tumor", "Clinical_type", "Sample", "Age", "Shared"))

# Set dimensions of A4 paper in inches
width = 8.27  # Width of A4 in inches
height = 11.69/3 # Height of A4 in inches

# Set margins
left_margin = 0.5  # Left margin in inches
right_margin = 0.5  # Right margin in inches

# Calculate plot width and height after considering margins
plot_width = width 
plot_height = height #- top_margin - bottom_margin
grid_width = plot_width - 0.1 #4cm in inches
grid_margins = (plot_width - grid_width)/2

pdf('figures/supp_figures/supp_1.pdf', width = plot_width, height = plot_height)
plot(1, type = "n", xlim = c(0, plot_width), ylim = c(0, plot_height), axes =FALSE, xlab = "", ylab = "") #Create empty plot area
line_positions_patients = c(grid_margins, cumsum(grid_width*c(11/81,12/81,
                                                              4/81, 4/81, 7/81,
                                                              5/81, 5/81, 5/81,
                                                              3/81, 5/81, 4/81,
                                                              4/81, 4/81, 4/81,
                                                              2/81, 2/81)) +grid_margins, plot_width-grid_margins)

#line_positions_patients are line positions at the xaxis proportional to the number of samples
segments(line_positions_patients, 0, line_positions_patients, plot_height, lwd = 0.5)

#Add title for each patient
for (patient in 1:16){
  text(y = plot_height, x = line_positions_patients[[patient]] +
         ((line_positions_patients[[patient+1]]-line_positions_patients[[patient]])/2), labels = paste("Pat.", patient, sep = " "), cex = 5/12)
}

#Add panel with SNVs barplot
bars_top = plot_height-0.2
bars_bottom = bars_top - 0.4

for (patient in 1:16){
  bars_shared = (as.numeric((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Shared))[[1]])/4000)*(bars_top-bars_bottom)
  bars = (as.numeric((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(SNVs))[[1]])/4000)*(bars_top-bars_bottom) - bars_shared
  bar_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(bars)+1)*0.01))/(length(bars))
  x1 = c(line_positions_patients[[patient]]+0.01, cumsum(rep(bar_width, length(bars)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(bars)-1))+0.01)
  x2 = c(line_positions_patients[[patient]] + bar_width + 0.01, cumsum(rep(bar_width, length(bars)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(bars)-1)) + bar_width+0.01)
  rect(x1, bars_bottom, x2, (bars_bottom+bars_shared), col = "#9CCB89", border = NA)
  rect(x1, bars_bottom+bars_shared, x2, (bars_bottom+bars_shared+bars), col = "#9DCB8950", border=NA)
}

#Add sample names
names_top = bars_bottom - 0.025
for (patient in 1:16){
  type = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Clinical_type))[[1]])
  name =  as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Tumor))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(type)+1)*0.01))/(length(type))
  x1 = c(line_positions_patients[[patient]] +0.01, cumsum(rep(box_width, length(type)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(type)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width +0.01, cumsum(rep(box_width, length(type)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(type)-1)) + bar_width + 0.01)
  type = gsub("primary", "#A4D7DD", type)
  type = gsub("liver_met", "#F5B6D2", type)
  type = gsub("carcinosis", "#A985AB", type)
  type = gsub("normal", "#C39C94", type)
  type = gsub("lymph_met","#D379B0", type)
  rect(x1, names_top-box_width-0.01, x2, names_top, col = type, border = type, lwd = 0.5)
  text(x1+0.5*box_width, names_top -0.5*box_width, labels = name,cex = 5/12)
}
names_bottom = names_top - box_width - 0.01

#Add panel with chr18 loss 
chr18_top = names_bottom - 0.1
for (patient in 1:16){
  CNVs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(chr18_loss))[[1]])
  homologs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Homolog))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CNVs)+1)*0.01))/(length(CNVs))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + bar_width + 0.01)
  CNVs = gsub("LOSS", "#2AB8CB", CNVs)
  CNVs = gsub("FALSE", "white", CNVs)
  CNVs = gsub("CN-LOH", "#f9ba7a", CNVs)
  homologs = gsub("FALSE", "white", homologs)
  homologs = gsub("A", "#8b1812", homologs)
  homologs = gsub("B", "#eccc0c", homologs)
  homologs = gsub("biallelic", '#2ba037', homologs)
  rect(x1, chr18_top - box_width -0.01, x2, chr18_top, col = CNVs, border = homologs, lwd = 0.75)
}
chr18_bottom = chr18_top - box_width - 0.01

#Add panel with chr4 gain 
chr4_top = chr18_bottom - 0.025
for (patient in 1:16){
  CNVs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(chr4_gain))[[1]])
  homologs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Homolog4))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CNVs)+1)*0.01))/(length(CNVs))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + bar_width + 0.01)
  CNVs = gsub("GAIN", "#e98771", CNVs)
  CNVs = gsub("FALSE", "white", CNVs)
  homologs = gsub("FALSE", "white", homologs)
  homologs = gsub("A", "#8b1812", homologs)
  homologs = gsub("B", "#eccc0c", homologs)
  homologs = gsub("biallelic", '#2ba037', homologs)
  rect(x1, chr4_top - box_width - 0.01, x2, chr4_top, col = CNVs, border = homologs, lwd = 0.75)
}
chr4_bottom = chr4_top - box_width - 0.01

#Add panel with chr5 gain
chr5_top = chr4_bottom - 0.025
for (patient in 1:16){
  CNVs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(chr5_gain))[[1]])
  homologs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Homolog5))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CNVs)+1)*0.01))/(length(CNVs))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + bar_width + 0.01)
  CNVs = gsub("GAIN", "#e98771", CNVs)
  CNVs = gsub("FALSE", "white", CNVs)
  homologs = gsub("FALSE", "white", homologs)
  homologs = gsub("A", "#8b1812", homologs)
  homologs = gsub("B", "#eccc0c", homologs)
  homologs = gsub("biallelic", '#2ba037', homologs)
  rect(x1, chr5_top - box_width -0.01, x2, chr5_top, col = CNVs, border = homologs, lwd = 0.75)
}
chr5_bottom = chr5_top - box_width - 0.01

#Add panel with chr14 gain
chr14_top = chr5_bottom - 0.025
for (patient in 1:16){
  CNVs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(chr14_gain))[[1]])
  homologs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Homolog14))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CNVs)+1)*0.01))/(length(CNVs))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + bar_width + 0.01)
  CNVs = gsub("GAIN", "#e98771", CNVs)
  CNVs = gsub("FALSE", "white", CNVs)
  homologs = gsub("FALSE", "white", homologs)
  homologs = gsub("A", "#8b1812", homologs)
  homologs = gsub("B", "#eccc0c", homologs)
  rect(x1, chr14_top - box_width - 0.01, x2, chr14_top, col = CNVs, border = homologs, lwd = 0.75)
}
chr14_bottom = chr14_top - box_width - 0.01

chr20_top = chr14_bottom - 0.025
for (patient in 1:16){
  CNVs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(chr20_gain))[[1]])
  homologs = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(Homolog20))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CNVs)+1)*0.01))/(length(CNVs))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CNVs)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CNVs)-1)) + bar_width + 0.01)
  CNVs = gsub("GAIN", "#e98771", CNVs)
  CNVs = gsub("FALSE", "white", CNVs)
  homologs = gsub("FALSE", "white", homologs)
  homologs = gsub("A", "#8b1812", homologs)
  homologs = gsub("B", "#eccc0c", homologs)
  rect(x1, chr20_top - box_width - 0.01, x2, chr20_top, col = CNVs, border = homologs, lwd = 0.75)
}
chr20_bottom = chr20_top - box_width - 0.01

#Add panel with potential driver mutations
driver_top = chr20_bottom - 0.1
for (patient in 1:16){
  CDKN1B = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(CDKN1B))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(CDKN1B)+1)*0.01))/(length(CDKN1B))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(CDKN1B)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CDKN1B)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(CDKN1B)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(CDKN1B)-1)) + bar_width + 0.01)
  CDKN1B = gsub("FALSE", "white", CDKN1B)
  CDKN1B = gsub("TRUE", "#eccc0c", CDKN1B)
  rect(x1, driver_top - box_width - 0.01, x2, driver_top, col = CDKN1B, border = CDKN1B, lwd = 0.5)
}
driver_bottom = driver_top - box_width - 0.01

driver2_top = driver_bottom - 0.025
for (patient in 1:16){
  TNRC6B = as.character((metadata %>% filter(Patient == patient) %>% filter(Clinical_type != 'blood') %>% arrange(Sample) %>% select(TNRC6B))[[1]])
  box_width = ((line_positions_patients[[patient+1]] - line_positions_patients[[patient]])-((length(TNRC6B)+1)*0.01))/(length(TNRC6B))
  x1 = c(line_positions_patients[[patient]] + 0.01, cumsum(rep(box_width, length(TNRC6B)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(TNRC6B)-1)) + 0.01)
  x2 = c(line_positions_patients[[patient]] + box_width + 0.01, cumsum(rep(box_width, length(TNRC6B)-1))+line_positions_patients[[patient]] + cumsum(rep(0.01, length(TNRC6B)-1)) + bar_width + 0.01)
  TNRC6B = gsub("FALSE", "white", TNRC6B)
  TNRC6B = gsub("TRUE", "#eccc0c", TNRC6B)
  rect(x1, driver2_top - box_width - 0.01, x2, driver2_top, col = TNRC6B, border = TNRC6B, lwd = 0.5)
}

dev.off()

