library(ggplot2)
library(dplyr)


d<-read.csv('./giab_all_downsampled_proc.tsv',sep='\t',header=TRUE)

# Prepare the data
d_filtered <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    ILMN_cov > 0,
    paste0(Aligner, "+", SNVCaller) == 'sent+sentd'
  ) %>%
  group_by(CmpFootprint, SNPClass, ILMN_cov) %>%
  summarize(mean_Fscore = mean(Fscore, na.rm=TRUE)) %>%
  ungroup()

# Create line plot
ggplot(d_filtered, aes(x = ILMN_cov, y = mean_Fscore, color = SNPClass)) +
  geom_line(size=0.25) +
  geom_point(size=0.75) +
  facet_wrap(~CmpFootprint) +
  labs(
    x = "ILMN Coverage",
    y = "Mean F-score",
    title = "Mean F-score vs. ILMN Coverage by SNP Class (sent+sentd)",
    color = "SNP Class"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
