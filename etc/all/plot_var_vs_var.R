#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5 || length(args) > 9) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf> <x_var> <y_var> [<xlim_min> <xlim_max> <ylim_min> <ylim_max>]")
}

input_tsv <- args[1]
output_pdf <- args[2]
x_var <- args[3]
y_var <- args[4]
cmpf <- args[5]
zoom_limits_provided <- length(args) == 9
if (zoom_limits_provided) {
  xlim_min <- as.numeric(args[6])
  xlim_max <- as.numeric(args[7])
  ylim_min <- as.numeric(args[8])
  ylim_max <- as.numeric(args[9])
  zoomed_pdf <- sub("\\.pdf$", "_zoomed.pdf", output_pdf)
}

# Read input data
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Filter and prepare data, including coverage
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c(!!cmpf),
    !grepl('_gt50', SNPClass),
    (ILMN_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'sent+sentd') |
    (ONT_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'ont+sentdont') |
    (UG_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'ug+sentdug')
  ) %>%
  mutate(
    Pipeline = case_when(
      ILMN_cov > 0 ~ "ILMN (sent+sentd)",
      ONT_cov > 0 ~ "ONT (ont+sentdont)",
      UG_cov > 0 ~ "UG (ug+sentdug)"
    ),
    Coverage = case_when(
      ILMN_cov > 0 ~ ILMN_cov,
      ONT_cov > 0 ~ ONT_cov,
      UG_cov > 0 ~ UG_cov
    ),
    Coverage_bin = cut(
      Coverage,
      breaks=c(0, 1, 5, 10, 15, 20, 30, Inf),
      labels=c("≤1x", "1-5x", "5-10x", "10-15x", "15-20x", "20-30x", ">30x"),
      include.lowest=TRUE
    )
  )

# Explicitly defined point sizes for coverage bins
coverage_sizes <- c("≤1x"=2, "1-5x"=2.5, "5-10x"=3, "10-15x"=3.5, "15-20x"=4, "20-30x"=4.5, ">30x"=6)

# Plot using binned coverage sizes
pdf(output_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=.data[[x_var]], y=.data[[y_var]], color=SNPClass, shape=Pipeline)) +
  geom_point(size=2.1,alpha=0.7) +
   facet_wrap(~as.factor(Coverage)) +
    scale_size_manual(values=coverage_sizes) +
  labs(
    x = x_var,
    y = y_var,
    title = paste(y_var, "vs.", x_var, "by SNP Class, Pipeline, and Coverage Bin w/in ", cmpf),
    color = "SNP Class",
    shape = "Pipeline"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.5, "cm")
  )
dev.off()

pdf(zoomed_pdf, width=1200/72, height=800/72)

# Optional zoomed-in plot with same logic
if (zoom_limits_provided) {
  ggplot(d_filtered, aes(x=.data[[x_var]], y=.data[[y_var]], color=SNPClass, shape=Pipeline)) +
    geom_point(size=3.4,alpha=0.7) +
    facet_wrap(~as.factor(Coverage)) +
    scale_size_manual(values=coverage_sizes) +
    labs(
      x = x_var,
      y = y_var,
      title = paste("Zoomed", y_var, "vs.", x_var, "by SNP Class, Pipeline, and Coverage Bin, w/in ", cmpf),
      color = "SNP Class",
      shape = "Pipeline"    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.5, "cm")
    ) +
    coord_cartesian(xlim=c(xlim_min, xlim_max), ylim=c(ylim_min, ylim_max))
  }
    
dev.off()

