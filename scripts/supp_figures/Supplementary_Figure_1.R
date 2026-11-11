source('scripts/analysis/extended_metadata.R')

library(ggplot2) #v4.0.3
library(patchwork) #1.3.2
library(ggpmisc) #v1.0.0
library(ggpubr) #v1.0.0

pdf('figures/supp_figures/Supplementary_Figure_1.pdf', width = 180*0.0393700787, height = 100*0.0393700787)

ggplot(metadata, aes(x = Purity, y = SNVs)) +
  geom_point(
    data = metadata %>% filter(Purity > 0.4),
    size = 2, col = "darkblue"
  ) +
  geom_point(
    data = metadata %>% filter(Purity < 0.4),
    size = 2, col = "darkblue", shape = 1
  ) +
  stat_poly_line(
    data = metadata %>% filter(Purity > 0.4),
    method = "lm",
    se = FALSE,
    linewidth = 0.5,
    formula = y ~ x,
    fullrange = TRUE
  ) +
  ggpubr::stat_cor(
    data = metadata %>% filter(Purity > 0.4),
    method = "pearson",
    aes(label = paste(after_stat(r.label),
                      after_stat(p.label), sep = "~~~"))
  ) +
  ylab("Burden (SNVs)") +
  xlab("Purity") +
  ggtitle("Uncorrected") +
  theme_classic() |

ggplot(metadata, aes(x = Purity, y = SNVs_purity_corrected)) +
  geom_point(
    data = metadata %>% filter(Purity > 0.4),
    size = 2, col = "darkblue"
  ) +
  geom_point(
    data = metadata %>% filter(Purity < 0.4),
    size = 2, col = "darkblue", shape = 1
  ) +
  stat_poly_line(
    data = metadata %>% filter(Purity > 0.4),
    method = "lm",
    se = FALSE,
    linewidth = 0.5,
    formula = y ~ x,
    fullrange = TRUE
  ) +
  ggpubr::stat_cor(
    data = metadata %>% filter(Purity > 0.4),
    method = "pearson",
    aes(label = paste(after_stat(r.label),
                      after_stat(p.label), sep = "~~~"))
  ) +
  ylab("Burden (SNVs)") +
  xlab("Purity") +
  ggtitle("Corrected") +
  theme_classic()

dev.off()
