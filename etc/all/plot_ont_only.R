library(ggplot2)
library(dplyr)


d<-read.csv('./giab_all_downsampled_proc.tsv',sep='\t',header=TRUE)

# Prepare the data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    ONT_cov > 0,
    paste0(Aligner, "+", SNVCaller) == 'ont+sentdont'
  ) %>%
  group_by(CmpFootprint, SNPClass, ONT_cov) %>%
  summarize(mean_Fscore = mean(Fscore, na.rm=TRUE)) %>%
  ungroup()

# Create line plot
ggplot(d_filtered, aes(x = ONT_cov, y = mean_Fscore, color = SNPClass)) +
  geom_line(size=0.25) +
  geom_point(size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "ONT Coverage",
    y = "F-score",
    title = "F-score vs. ONT Coverage by SNP Class (sent+sentd), HG006, hg38",
    color = "SNP Class"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))