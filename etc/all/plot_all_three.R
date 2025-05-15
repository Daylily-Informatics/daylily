#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf>")
}

input_tsv <- args[1]
output_pdf <- args[2]
zoomed_pdf <- sub("\\.pdf$", "_zoomed.pdf", output_pdf)
zoomed2_pdf <- sub("\\.pdf$", "_zoomed_again.pdf", output_pdf)

# Read input data
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Filter out unwanted SNP classes and prepare data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    !grepl('_gt50', SNPClass),
    (ILMN_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'sent+sentd') |
    (ONT_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'ont+sentdont') |
    (UG_cov > 0 & paste0(Aligner, "+", SNVCaller) == 'ug+sentdug')
  ) %>%
  mutate(
    Coverage = case_when(
      ILMN_cov > 0 ~ ILMN_cov,
      ONT_cov > 0 ~ ONT_cov,
      UG_cov > 0 ~ UG_cov
    ),
    Pipeline = case_when(
      ILMN_cov > 0 ~ "ILMN (sent+sentd)",
      ONT_cov > 0 ~ "ONT (ont+sentdont)",
      UG_cov > 0 ~ "UG (ug+sentdug)"
    )
  ) %>%
  group_by(CmpFootprint, Pipeline, SNPClass, Coverage) %>%
  summarize(mean_Fscore = mean(Fscore, na.rm=TRUE), .groups='drop')
d_filtered[[plot_var]] <- as.numeric(as.character(d_filtered[[plot_var]]))

# Create overview plot
pdf(output_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=Coverage, y=mean_Fscore, color=SNPClass, shape=Pipeline, linetype=Pipeline)) +
  geom_line(size=0.75) +
  geom_point(size=1.25) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "Coverage",
    y = "Mean F-score",
    title = "Mean F-score by Coverage and SNP Class across Pipelines",
    color = "SNP Class",
    shape = "Pipeline",
    linetype = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Create zoomed-in plot
pdf(zoomed_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=Coverage, y=mean_Fscore, color=SNPClass, shape=Pipeline, linetype=Pipeline)) +
  geom_line(size=0.75) +
  geom_point(size=1.25) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "Coverage",
    y = "Mean F-score",
    title = "Zoomed Mean F-score by Coverage and SNP Class across Pipelines",
    color = "SNP Class",
    shape = "Pipeline",
    linetype = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(0.990,1.0))
dev.off()


# Create zoomed-in plot
pdf(zoomed2_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=Coverage, y=mean_Fscore, color=SNPClass, shape=Pipeline, linetype=Pipeline)) +
  geom_line(size=0.75) +
  geom_point(size=1.25) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "Coverage",
    y = "Mean F-score",
    title = "Zoomed Again, Mean F-score by Coverage and SNP Class across Pipelines",
    color = "SNP Class",
    shape = "Pipeline",
    linetype = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(0.996,1.0))
dev.off()


