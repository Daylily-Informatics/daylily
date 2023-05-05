plt <-ggplot(data=subset(d, rule_suffix=='alNsort'))+geom_point( aes(x="GIAB Samples",y=s, color=sample))+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(2000,3000) + ggtitle("7 GIAB Samples bwa-mem2-ert Runtime on 128 Cores") + xlab("") + ylab("seconds")
> ggsave("./ert_runtimes_giab.png", plt)
Saving 7.89 x 9.49 in image

plt<-ggplot(data=d) + geom_boxplot( aes(x=s,y=rule_suffix)) + ylab('Job-Rule Name') + xlab('seconds') + ggtitle('Runtimes For All Jobs Necessary To Process 7 30x GIAB Samples')

ggsave("./ert_ALL_runtimes_giab.png", plt)

`%notin%` <- Negate(`%in%`)
rm_class = c('INS_gt50','DEL_gt50','Indel_gt50')
plt<-ggplot(subset(dc, SNPClass %notin% rm_class)) + geom_point(aes(y=Fscore, x=Sample, color=SNVCaller)) + facet_wrap(~ SNPClass) + ylim(0.97,1)+theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Fscores For 30x GIAB Samples Aligned to b37 & Called By Deepvariant and Octopus")

ggsave("./giab_fscore_oct_dv_b37.png", plt)

ggplot(data=d) + geom_boxplot( aes(x=0.05417*s/60,y=rule_suffix)) + ylab('Job-Rule Name') + xlab('seconds') + ggtitle('Runtimes For All Jobs Necessary To Process 7 30x GIAB Samples') + xlim(0,3)

ggsave("./ert_ALL_runtimes_giab_cost.png", plt)
