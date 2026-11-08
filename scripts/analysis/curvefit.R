#library(readr)
#library(dplyr)
#library(stats)

summary = read_delim('input_data/WGS/simulated_WGS.tsv', delim = '\t')

summary = summary %>% dplyr::filter(Sample != 'tissue')

x = c(1, 0.25, 0.5, 0.75)
y = summary$`Recall (assayable genome fraction)`
df = data.frame(x, y)

fit = nls(
  y ~ (1 - exp(-a * x^b)) / (1 - exp(-a)),
  data = df,
  start = list(a = 1, b = 1)
)

a = round(summary(fit)$coefficients['a','Estimate'], 5)
b = round(summary(fit)$coefficients['b','Estimate'], 5)

rm(summary, x, y, df, fit)
