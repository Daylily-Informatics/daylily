#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript script_name.R <input_tsv> <output_pdf> <plot_var> <ylim_min> <ylim_max>")
}

input_tsv <- args[1]
output_pdf <- args[2]
plot_var <- args[3]
ylim_min <- as.numeric(args[4])
ylim_max <- as.numeric(args[5])
zoomed_pdf <- sub("\\.pdf$", "_zoomed.pdf", output_pdf)

# Read input data
d <- read.csv(input_tsv, sep='\t', header=TRUE)

# Filter out unwanted SNP classes and prepare data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima', 'giabHC_x_clinvar_genes'),
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
  )

d_filtered[[plot_var]] <- as.numeric(as.character(d_filtered[[plot_var]]))

# Plot overview with regression lines per SNPClass and Pipeline
pdf(output_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=Coverage, y=.data[[plot_var]], color=SNPClass, shape=Pipeline, linetype=Pipeline)) +
  geom_point(size=2) +
  geom_smooth(method="lm", se=TRUE, size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "Coverage",
    y = plot_var,
    title = paste(plot_var, "by Coverage, SNP Class, and Pipeline (with regression lines)"),
    color = "SNP Class",
    shape = "Pipeline",
    linetype = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Create zoomed-in plot
pdf(zoomed_pdf, width=1200/72, height=800/72)
ggplot(d_filtered, aes(x=Coverage, y=.data[[plot_var]], color=SNPClass, shape=Pipeline, linetype=Pipeline)) +
  geom_point(size=2) +
  geom_smooth(method="lm", se=TRUE, size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "Coverage",
    y = plot_var,
    title = paste("Zoomed", plot_var, "by Coverage, SNP Class, and Pipeline (with regression lines)"),
    color = "SNP Class",
    shape = "Pipeline",
    linetype = "Pipeline"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(ylim_min, ylim_max))
dev.off()

# Regression summaries per group
cat("\nRegression summaries per SNPClass and Pipeline:\n")
d_filtered %>%
  group_by(SNPClass, Pipeline) %>%
  do(model = summary(lm(as.formula(paste(plot_var, "~ Coverage")), data=.))) %>%
  rowwise() %>%
  do({
    cat("\nGroup:", .$SNPClass, "|", .$Pipeline, "\n")
    print(.$model)
    data.frame()
  })



# Equivalence calculation per Footprint and SNPClass
footprints <- c('giabHC_x_ultima', 'giabHC')
snp_classes <- c('SNPts', 'SNPtv', 'INS_50', 'DEL_50')

for (footprint in footprints) {
  for (snp_class in snp_classes) {
    cat("\n--- Footprint:", footprint, "| SNP Class:", snp_class, "---\n")

    subset <- d_filtered %>%
      filter(CmpFootprint == footprint, SNPClass == snp_class)

    preds <- list()

    for (pipeline in unique(subset$Pipeline)) {
      pipe_data <- subset %>% filter(Pipeline == pipeline)
      
      if (nrow(pipe_data) < 2) next
      
      model <- lm(as.formula(paste(plot_var, "~ Coverage")), data=pipe_data)
      cov_range <- seq(floor(min(pipe_data$Coverage)), ceiling(max(pipe_data$Coverage)), by=1)

      pred_df <- data.frame(Coverage=cov_range)
      preds[[pipeline]] <- cbind(Coverage=cov_range, predict(model, pred_df, interval="confidence"))
    }

    if (length(preds) == 2) {
      common_cov <- intersect(preds[[1]][,"Coverage"], preds[[2]][,"Coverage"])

      equivalence <- data.frame(
        Coverage = common_cov,
        ILMN_Lower = preds[["ILMN (sent+sentd)"]][match(common_cov, preds[["ILMN (sent+sentd)"]][,"Coverage"]),"lwr"],
        ILMN_Upper = preds[["ILMN (sent+sentd)"]][match(common_cov, preds[["ILMN (sent+sentd)"]][,"Coverage"]),"upr"],
        UG_Lower = preds[["UG (ug+sentdug)"]][match(common_cov, preds[["UG (ug+sentdug)"]][,"Coverage"]),"lwr"],
        UG_Upper = preds[["UG (ug+sentdug)"]][match(common_cov, preds[["UG (ug+sentdug)"]][,"Coverage"]),"upr"]
      )

      equivalence$Equivalent <- with(equivalence, (UG_Lower <= ILMN_Upper) & (ILMN_Lower <= UG_Upper))
      print(equivalence)
    } else {
      cat("Insufficient data for", snp_class, "in", footprint, "\n")
    }
  }
}
