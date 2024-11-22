args = commandArgs(trailingOnly=TRUE)

in_tsv <- args[1]
sample <- args[2]
chrm <- args[3]
alnr <- args[4]

require(data.table);

col_names <- c('CHR','START','END','COV');
col_classes <- c('factor','integer','integer','numeric');
dcov <- fread(in_tsv,header=FALSE,sep='\t', na.strings=c(""), colClasses=col_classes);
setnames(dcov,col_names);

dcov$COVnorm <- dcov$COV/median(dcov$COV);

sdRCov <- sd(dcov$COV);
sdNCov <- sd(dcov$COVnorm);
pct0 <- length(subset(dcov, COV==0))/length(dcov$COV);
pctLT5 <- length(subset(dcov, COV <5))/length(dcov$COV);
pctLT10 <- length(subset(dcov, COV <10))/length(dcov$COV);
RCcoefofvar = sdRCov/mean(dcov$COV);
NCcoefofvar = sdNCov/mean(dcov$COVnorm);


write.table(matrix(as.character(c(sample,chrm,mean(dcov$COV),median(dcov$COV),sdRCov,RCcoefofvar,mean(dcov$COVnorm),median(dcov$COVnorm),sdNCov,NCcoefofvar,pct0,pctLT5,pctLT10,alnr)),nrow=1),sep="\t",row.names=FALSE,col.names=FALSE)