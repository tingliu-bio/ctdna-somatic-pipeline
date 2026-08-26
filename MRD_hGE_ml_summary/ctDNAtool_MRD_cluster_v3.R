##Ting Liu July-2022
#ctDNA pacakge can be installed locally
#module load R/4.1.0  (R/4.0.4 works too)
#install.packages("/resources/downloads/ctDNAtools_0.4.0.tar", repos = NULL, type="source")
#install NCBI.GRCh38 if needed
##BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

library(data.table)
library(ctDNAtools)
library(purrr)
library(tidyr)
library(plyr)
library(ggplot2)
library(stringr)

#current time stamp and working dir
date=Sys.Date()
work_dir=getwd()

#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
install.packages("/resources/downloads/ctDNAtools_0.4.0.tar", repos = NULL, type="source")
install.packages("/resources/downloads/BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000.tar", repos = NULL, type="source")
#data('mutations',package = 'ctDNAtools')

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

targets<-read.table('/ressources/target_regions/hg38/picard_tools_targets_input.target',header=T)

bam_path=paste(work_dir, "reheader_tumor_bam/", sep="/")
T1_vcf_path=paste(work_dir,"T1_mutations/", sep="/")


setwd(bam_path)

files=list.files()
bam_files=files[which(str_detect(files, "new_header.bam$"))]  # "$" here is to ensure to read the bam files only

#ini the null dataframe to collect the output
out=as.data.frame(matrix(nrow=0,ncol=10))

for (bam in bam_files){
  #get the sample id from the mutation file
  sample_id=unlist(strsplit(bam, "\\."))[1]     #get the sample id
  patient_id=unlist(strsplit(sample_id, "\\_"))[1]   #get the patient id
#  cat ("sample id is:", sample_id,"\n")
#  cat ("patient id is:", patient_id,"\n")
  T1_loc=paste(T1_vcf_path,patient_id,"_T1.sorted.sscs.sc.filtered.case_only.vcf.gz",sep="")  ##get the loc of T1 vcf file
  cat (T1_loc,"\n")
 if(file.exists(T1_loc)) {  #remove the zero vcf T1 files in advance
  T1_sampleid=paste(patient_id,"_T1",sep="")
#  vcf_T1<-vcf_to_mutations_df(T1_loc,sample_name = T1_sampleid)

#read the VCF file, only keep the first 4 cols
  vcf_T1<-suppressWarnings(vcf_to_mutations_df(T1_loc))[1:4]
#run the test, then add the output line into the out frame
 test1 <- test_ctDNA(mutations = vcf_T1, bam = paste(bam_path,bam,sep=""),
                      reference = BSgenome.Hsapiens.UCSC.hg38,
                      targets = targets,
                      informative_reads_threshold = 100)  #run2 was 10000
 out<-rbind(out,test1)
 }#if
else{
 # header=c("sample,n_mutations,n_nonzero_alt,total_alt_reads,mutations_filtered,background_rate,informative_reads,multi_support_reads,pvalue,decision")
 # zero<-c(sample_id,"0","NA","NA","NA","NA","NA","NA","NA","undetermined")
  zero_table<-data.frame(sample=sample_id,n_mutations="0",n_nonzero_alt="NA",total_alt_reads="NA",mutations_filtered="NA",background_rate="NA",informative_reads="NA",multi_support_reads="NA",pvalue="NA",decision="undetermined")
 out<-rbind(out,zero_table)
 }
} #for

#return the output as a data frame

setwd(work_dir)
write.csv(out, file=paste(date, "ctDNA_output_sscs_sc.txt", sep="_"),quote=F, row.names=F)
