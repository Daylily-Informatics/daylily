#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4 || length(args) > 8) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf> <x_var> <y_var> [<xlim_min> <xlim_max> <ylim_min> <ylim_max>]")
}

input_tsv <- args[1]
output_pdf <- args[2]
x_var <- args[3]
y_var <- args[4]

# Check if zoom limits are provided
zoom_limits_provided <- length(args) == 8
if (zoom_limits_provided) {
  xlim_min <- as.numeric(args[5])
  xlim_max <- as.numeric(args[6])
  ylim_min <- as.numeric(args[7])
  ylim_max <- as.numeric(args[8])
  zoomed_pdf <- sub("\\.pdf$", "_zoomed.pdf", output_pdf)
}

# Read input data
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Filter and prepare data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima', 'giabHC_x_clinvar_genes'),
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
    )
  )

# Prepare long format for coverage-specific analysis
coverage_long <- d_filtered %>%
  pivot_longer(cols = c(ILMN_cov, ONT_cov, UG_cov), names_to = "Coverage_Type", values_to = "Coverage") %>%
  filter(Coverage > 0)

# Create overview dotplot explicitly by coverage type
pdf(output_pdf, width=1200/72, height=800/72)
ggplot(coverage_long, aes(x=.data[[x_var]], y=.data[[y_var]], color=SNPClass, shape=Pipeline)) +
  geom_point(size=2) +
  facet_grid(CmpFootprint ~ Coverage_Type, scales="free_x") +
  labs(
    x = x_var,
    y = y_var,
    title = paste(y_var, "vs.", x_var, "by SNP Class, Pipeline, and Coverage Type"),
    color = "SNP Class",
    shape = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Optionally create zoomed-in plot
if (zoom_limits_provided) {
  pdf(zoomed_pdf, width=1200/72, height=800/72)

  ggplot(coverage_long, aes(x=.data[[x_var]], y=.data[[y_var]], color=SNPClass, shape=Pipeline)) +
    geom_point(size=2) +
    facet_grid(CmpFootprint ~ Coverage_Type, scales="free_x") +
    labs(
      x = x_var,
      y = y_var,
      title = paste("Zoomed", y_var, "vs.", x_var, "by SNP Class, Pipeline, and Coverage Type"),
      color = "SNP Class",
      shape = "Pipeline"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(xlim=c(xlim_min, xlim_max), ylim=c(ylim_min, ylim_max))

  dev.off()
}

# Fit regression model and print summary
model_formula <- as.formula(paste(y_var, "~ ILMN_cov + ONT_cov + UG_cov + SNPClass + Pipeline"))
model <- lm(model_formula, data=d_filtered)

# Output regression summary to console
cat("\nRegression Model Summary:\n")
print(summary(model))
