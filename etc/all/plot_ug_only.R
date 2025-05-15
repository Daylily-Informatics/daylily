#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Check if two arguments are provided
if (length(args) != 2) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf>")
}

# Assign arguments to variables
input_tsv <- args[1]
output_pdf <- args[2]
zoomed_pdf <- sub("\\.pdf$", "_zoomed.pdf", output_pdf)


# Read input file
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Prepare the data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    UG_cov > 0,
    paste0(Aligner, "+", SNVCaller) == 'ug+sentdug'
  ) %>%
  group_by(CmpFootprint, SNPClass, UG_cov) %>%
  summarize(mean_Fscore = mean(Fscore, na.rm=TRUE)) %>%
  ungroup()

pdf(output_pdf, width=1200/72, height=800/72)

# Create line plot
ggplot(d_filtered, aes(x = UG_cov, y = mean_Fscore, color = SNPClass)) +
  geom_line(size=0.25) +
  geom_point(size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "UG Coverage",
    y = "F-score",
    title = "F-score vs. UG Coverage by SNP Class (sent+sentd), HG006, hg38",
    color = "SNP Class"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


# Create zoomed plot
pdf(zoomed_pdf, width=1200/72, height=800/72)

ggplot(d_filtered, aes(x = UG_cov, y = mean_Fscore, color = SNPClass)) +
  geom_line(size=0.25) +
  geom_point(size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "UG Coverage",
    y = "F-score",
    title = "F-score vs. UG Coverage by SNP Class (sent+sentd), HG006, hg38",
    color = "SNP Class"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(0.980,1.0))

dev.off()
