#!/bin/bash
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=46000M
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH -J reheader_bam

##to accoumadate the ctDNAtools in grch38 ref_genome
##we will have to reheader the bam files to get rid of the alternate contigs, alternate scaffolds
##otherwise will be errors when you run the ctDNA tools "Error: Chromosomes in mutations and/or targets don't match the specified reference"

module load samtools/1.20

work_path=$SLURM_SUBMIT_DIR
mkdir ${work_path}/header
mkdir ${work_path}/reheader_tumor_bam

files_location=${work_path}/tumor_samples
file=${files_location}/*.sorted.bam

for i in ${file}
do
#  echo $i
if [[ $i =~ ${files_location}/(.*).sorted.bam ]] ;then
   num=${BASH_REMATCH[1]}
   echo -n $num
   samtools view -H $i > ${work_path}/header/${num}.header.sam
   grep -P 'SN:chr[1-9|XY]|^@HD|^@RG|^@PG' ${work_path}/header/${num}.header.sam >${work_path}/header/${num}.new_header.sam
   samtools reheader ${work_path}/header/${num}.new_header.sam $i > ${work_path}/reheader_tumor_bam/${num}.new_header.bam
   samtools index ${work_path}/reheader_tumor_bam/${num}.new_header.bam
 fi
done
