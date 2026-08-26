
date = Sys.Date()

library(data.table)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)

work_dir=getwd()

dcs_sc<-read.table(paste(work_dir,"/dcs_sc/coverage_per_sample_dcs_sc.txt",sep = ""),header=T)
sscs_sc<-read.table(paste(work_dir,"/sscs_sc/coverage_per_sample_sscs_sc.txt",sep = ""),header=T)
uniq_sc<-read.table(paste(work_dir,"/unique_dcs/coverage_per_sample_uniq_dcs.txt",sep = ""),header=T)
uncollapsed<-read.table(paste(work_dir,"/uncollapsed_raw_bam/coverage_per_sample_uncollapsed.txt",sep = ""),header=T)

merged_cov<-merge(merge(merge(dcs_sc,sscs_sc),uniq_sc),uncollapsed)

write.table(merged_cov, paste(work_dir,"/merged_coverage.txt",sep = ""), quote=F, sep = "\t", row.names = F)

