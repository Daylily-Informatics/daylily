#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(scales)
library(viridis)


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if two arguments are provided
if (length(args) != 3) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf> <plotvar>")
}

# Assign arguments to variables
input_tsv <- args[1]
output_pdf <- args[2]
pvar <- args[3]

# Read input file
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Baseline ILMN 30x (giabHC, sent+sentd), excluding unwanted SNP classes
baseline_ilmn <- d %>%
  filter(
    CmpFootprint == 'giabHC',
    ILMN_cov == 30,
    paste0(Aligner, "+", SNVCaller) == 'sent+sentd',
    !SNPClass %in% c('Indel_gt50', 'DEL_gt50', 'INS_gt50')
  ) %>%
  group_by(SNPClass) %>%
  summarize(
    baseline_Fscore = mean(as.numeric(as.character(.data[[pvar]])), na.rm=TRUE),
    .groups = 'drop'
  )

# Hybrid UG+ONT data (ont+senthuo)
hybrid_ug_ont <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    paste0(Aligner, "+", SNVCaller) == 'ont+sentdhuo',
    UG_cov > 0, ONT_cov > 0,
    !SNPClass %in% c('Indel_gt50', 'DEL_gt50', 'INS_gt50')
  ) %>%
  group_by(CmpFootprint, SNPClass, UG_cov, ONT_cov) %>%
  summarize(
  mean_Fscore = mean(as.numeric(as.character(.data[[pvar]])), na.rm=TRUE),
    .groups = 'drop'
  ) %>%  left_join(baseline_ilmn, by='SNPClass') %>%
  mutate(
    exceeds_baseline = mean_Fscore > baseline_Fscore,
    label = sprintf("%.4f [%.4f]", mean_Fscore, baseline_Fscore)
  )

pdf(output_pdf, width=1200/72, height=800/72)

# Plot UG+ONT hybrid with correct color scale
ggplot(hybrid_ug_ont, aes(x=factor(UG_cov), y=factor(ONT_cov), fill=mean_Fscore)) +
  geom_tile(color='grey90') +
  geom_text(aes(label=label), size=2.9, color='grey10') +
  geom_text(data=hybrid_ug_ont %>% filter(exceeds_baseline),
            aes(label='*'), color='magenta', size=9, nudge_x=0.35) +
  facet_grid(SNPClass ~ CmpFootprint) +
  scale_fill_gradientn(
    colours = viridis(256),
    values = rescale(c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.997, 1)),
    limits = c(0.5, 1),
    oob = squish,
    name = "Mean F-score"
  ) +
  labs(
    x = "UG Coverage",
    y = "ONT Coverage",
    title = paste("Hybrid UG+ONT, ", pvar ," (ont+senthuo)\n(* indicates surpassing ILMN 30x baseline (giabHC sent+sentd))")
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
