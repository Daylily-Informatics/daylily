args = commandArgs(trailingOnly=TRUE)

in_tsv <- args[1]

# Install packages if not already installed
packages <- c("ggplot2", "dplyr", "tidyr", "pheatmap")
installed_packages <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed_packages)) {
    install.packages(p)
  }
}

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

# Read the data
data <- read.csv(in_tsv, sep="\t", stringsAsFactors = FALSE)

# Convert columns to numeric
numeric_cols <- c("TgtRegionSize", "TN", "FN", "TP", "FP", "Fscore",
                  "Sensitivity.Recall", "Specificity", "FDR", "PPV",
                  "Precision", "AllVarMeanDP", "CovBin")
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)

# Create Aligner_Caller variable
data$Aligner_Caller <- paste(data$Aligner, data$SNVCaller, sep = "_")

# Plot 1: Bar Plots of F1 Score
plot1 <- ggplot(data, aes(x = Aligner_Caller, y = Fscore, fill = SNVCaller)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(SNPClass ~ Sample, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "F1 Score by Aligner and Caller",
       x = "Aligner and Caller",
       y = "F1 Score")

# Save Plot 1
ggsave("F1_Score_Bar_Plot.png", plot = plot1, width = 12, height = 8, dpi = 300)




# Updated Plot 2: Precision vs. Sensitivity Scatter Plots with Sample Shapes
plot2 <- ggplot(data, aes(x = Precision, y = Sensitivity.Recall, color = Aligner_Caller, shape = Sample)) +
  geom_point(size = 3) +
  facet_wrap(~ SNPClass, scales = "free") +
  theme_bw() +
  labs(title = "Precision vs. Sensitivity by Aligner and Caller",
       x = "Precision",
       y = "Sensitivity (Recall)") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))
  
# Display Plot 2
print(plot2)

# Save Plot 2
ggsave("Precision_vs_Sensitivity_Scatter_SampleShapes.png", plot = plot2, width = 12, height = 8, dpi = 300)

library(viridis)


# Plot 3: Heatmaps of Performance Metrics
heatmap_data <- data %>%
  group_by(SNPClass, Aligner, SNVCaller) %>%
  summarize(Fscore = mean(Fscore, na.rm = TRUE)) %>%
  unite("Aligner_Caller", Aligner, SNVCaller, sep = "_") %>%
  spread(key = Aligner_Caller, value = Fscore)

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$SNPClass

# Plot and save heatmap
# Plot Heatmap with Artistic Color Palette and Save
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "F1 Score Heatmap by Aligner and Caller",
         color = viridis(100),
         filename = "F1_Score_Heatmap_new.png",
         width = 10,
         height = 8)


# Plot 4: Line Charts of F1 Score Across Samples with Rotated X-axis Labels
plot4 <- ggplot(data, aes(x = Sample, y = Fscore, color = Aligner_Caller, group = Aligner_Caller)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ SNPClass, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotates x-axis labels by 45 degrees
  labs(title = "F1 Score Across Samples",
       x = "Sample",
       y = "F1 Score")

# Display Plot 4
print(plot4)

# Save Plot 4
ggsave("F1_Score_Line_Chart_Rotated_XLabels.png", plot = plot4, width = 12, height = 8, dpi = 300)