#!/bin/bash
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=46000M
#SBATCH -t 4-10:00 # Runtime in D-HH:MM
#SBATCH -J normal_consensus
#SBATCH --array=0-99 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

#set -o pipefail
#set -e

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

#gatk_bundle="/cluster/projects/kridelgroup/resources/GATK_bundle/homo_sapiens_GRCh38"
#resources to run the gatk, no slash at the end

file_path="${work_path}/FASTQ_merged_Normals"
#path to the fq files, no slash at the end
#file=$file_path/*.fastq.gz
##use the prefix of the fq file will be used as the output prefix


#!!#note
#2026, samtools is updated from 1.10 to 1.20
module load bwa/0.7.15
module load python3/3.7.2
#module load samtools/1.10
module load samtools/1.20

##output directory
if [ ! -d "${work_path}/consensus_Normals" ]; then
        mkdir -p "${work_path}/consensus_Normals"
fi

if [ ! -d "${work_path}/consensus_output_Normals" ]; then
        mkdir -p "${work_path}/consensus_output_Normals"
fi
#s ${file_path}/*.bam  > ${work_path}/fq_files_test.txt
#create a list to save the bam files

##assign each sample into different ARRAYID
samples=${work_path}/fq1_GL_files.txt

prefixs=($(cat $samples))
#echo  "prefixs="$prefixs
#echo "prefixs="${prefixs[@]}
input_fq1=${prefixs[${SLURM_ARRAY_TASK_ID}]}

bam_path=${work_path}/consensus_Normals/bamfiles
##setup the parameters
UMI_index=IDT_dual_Index.txt
cytobandpath=hg38_cytoBand.txt

###traverse the file_path/.fq folder, then get the file prefix from the r1 reads
###the script will not check if the r2 reads exists or not, please make sure the paired reads is placed in the same ./cleanfq folder

##get the prefix and control name
if [[ "$input_fq1" =~ ${file_path}/(.*)_R1.fastq.gz ]]
 then
	echo "input_fq1=$input_fq1"
	prefix="${BASH_REMATCH[1]}"
	echo "prefix=$prefix"

# run consensus cruncher
python ConsensusCruncher.py fastq2bam \
--fastq1 $input_fq1  --fastq2 ${file_path}/${prefix}_R2.fastq.gz \
-o ${work_path}/consensus_Normals \
-r GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-b bwa \
-s samtools \
-l  ${UMI_index} \
-g picard.jar

python /ConsensusCruncher.py consensus \
-i ${bam_path}/${prefix}.sorted.bam \
-o ${work_path}/consensus_output_Normals \
-s samtools \
-g hg38 \
-b ${cytobandpath} \
--cleanup True

fi
