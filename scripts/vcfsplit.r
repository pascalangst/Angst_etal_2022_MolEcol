
args <- commandArgs()
library(PopGenome)
sample <- args[6]
VCF_split_into_scaffolds(paste("vcf/",sample,".SNPsonly.recode.vcf", sep=""), paste("vcf/vcf_split/",sample,"/", sep=""))