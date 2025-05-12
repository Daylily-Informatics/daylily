library(ggplot2)
library(dplyr)
library(scales)

d <- read.csv('./giab_all_downsampled_proc.tsv', sep='\t', header=TRUE)

# Baseline ILMN 30x (giabHC, sent+sentd), excluding unwanted SNP classes
baseline_ilmn <- d %>%
  filter(
    CmpFootprint == 'giabHC',
    ILMN_cov == 30,
    paste0(Aligner, "+", SNVCaller) == 'sent+sentd',
    !SNPClass %in% c('Indel_gt50', 'DEL_gt50', 'INS_gt50')
  ) %>%
  group_by(SNPClass) %>%
  summarize(baseline_Fscore = mean(Fscore, na.rm=TRUE), .groups='drop')

# Hybrid ILMN+ONT data (ont+sentdhio)
hybrid_ilmn_ont <- d %>%
  filter(
    CmpFootprint %in% c('giabHC', 'giabHC_x_ultima'),
    paste0(Aligner, "+", SNVCaller) == 'ont+sentdhio',
    ILMN_cov > 0, ONT_cov > 0,
    !SNPClass %in% c('Indel_gt50', 'DEL_gt50', 'INS_gt50')
  ) %>%
  group_by(CmpFootprint, SNPClass, ILMN_cov, ONT_cov) %>%
  summarize(mean_Fscore = mean(Fscore, na.rm=TRUE), .groups='drop') %>%
  left_join(baseline_ilmn, by='SNPClass') %>%
  mutate(
    exceeds_baseline = mean_Fscore > baseline_Fscore,
    label = sprintf("%.3f [%.3f]", mean_Fscore, baseline_Fscore)
  )

# Plot ILMN+ONT hybrid
ggplot(hybrid_ilmn_ont, aes(x=factor(ILMN_cov), y=factor(ONT_cov), fill=mean_Fscore)) +
  geom_tile(color='grey90') +
  geom_text(aes(label=label), size=2.2, color='grey10') +
  geom_text(data=hybrid_ilmn_ont %>% filter(exceeds_baseline),
            aes(label='*'), color='magenta', size=6, nudge_x=0.3) +
  facet_grid(SNPClass ~ CmpFootprint) +
  scale_fill_viridis_c(name="Mean F-score", limits=c(0.5,1), oob=squish, option='viridis',
                       values=rescale(c(0.5, 0.997, 1))) +
  labs(
    x="ILMN Coverage",
    y="ONT Coverage",
    title="Hybrid ILMN+ONT Mean F-score (ont+sentdhio)\n(* indicates surpassing ILMN 30x baseline (giabHC sent+sentd))"
  ) +
  theme_bw(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
