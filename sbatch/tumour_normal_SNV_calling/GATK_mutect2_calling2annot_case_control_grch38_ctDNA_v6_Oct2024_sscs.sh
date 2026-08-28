#!/bin/bash
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=46000M
#SBATCH -t 0-20:00 # Runtime in D-HH:MM
#SBATCH -J sscs_somatic_call
#SBATCH --array=0-99 # job array index - number of jobs = numb of unique samples in file_path folder; e.g. 6samples then --array=0-5

#set -o pipefail
#set -e

#Oct-2024 update the commd of extracting the tumor_prefix, not only LIB, but also all ids

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

gatk_bundle="/projects/GATK_bundle/homo_sapiens_GRCh38"
#resources to run the gatk, no slash at the end
ref_genome="/projects/Homo_sapiens_BWA_ref_updated_2024/seqs_for_alignment_pipelines"


file_path="${work_path}/tumor_samples"
control_path="${work_path}/germline_controls"
#path to the fq files, no slash at the end
#file=$file_path/*.fastq.gz
##use the prefix of the fq file will be used as the output prefix

#create the sample list in advance
#ls $PWD/tumor_samples/*.bam >bam_files_sscs.txt

module load bwa/0.7.15
module load python/3.4.3
module load gatk/4.1.8.1
module load samtools/1.20
module load annovar/20180416
module load tabix
module load vt/0.577
module load vep/98
##output directory

if [ ! -d "${work_path}/VCFs" ]; then
        mkdir -p "${work_path}/VCFs"
fi

if [ ! -d "${work_path}/coverage" ]; then
        mkdir -p "${work_path}/coverage"
fi

if [ ! -d "${work_path}/VEP_annot" ]; then
        mkdir -p "${work_path}/VEP_annot"
fi

if [ ! -d "${work_path}/Annovar_annot" ]; then
        mkdir -p "${work_path}/Annovar_annot"
fi
#ls ${file_path}/tumor_samples/*.bam  > ${work_path}/bam_files_test.txt
#create a list to save the bam files

##assign each sample into different ARRAYID
samples=${work_path}/bam_files_sscs.txt

prefixs=($(cat $samples))
echo  "prefixs="${prefixs[@]}
input_bam=${prefixs[${SLURM_ARRAY_TASK_ID}]}

##get the tumor samplename to input into the MuTect2
#tumor_samplename=($(samtools view -H ${input_bam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
tumor_samplename=($(samtools view -H ${input_bam} | grep '^@RG'  |head -1 | sed "s/.*SM:\([^\t]*\).*/\1/g"))
#printout the sample name to check the files

echo "tumor_samplename=$tumor_samplename"
#coordinates of probes (coding and non-coding) used as provided by IDT and Robert

amplicon_interval_list=/ressources/target_regions/hg38/picard_tools_amps_input.bed
#amplicon intervals list that liftovered by Gabrielle

#coordinates of targets (coding and non-coding) used as provided by IDT and Robert
targets_interval_list=/ressources/target_regions/hg38/picard_tools_targets_input.bed
##target intervals list that liftovered by Gabrielle

##get the prefix and control name
if [[ "$input_bam" =~ ${file_path}/(.*).sorted.bam ]]
 then
	echo "imput_bam=$input_bam"
	prefix="${BASH_REMATCH[1]}"
	echo "prefix=$prefix"

##get the file name and sample name of the control sample
tumor_prefix=($(samtools view -H ${input_bam} | grep '^@RG'  |head -1 | sed 's/.*SM:\([^\t]*\).*/\1/g' | sed 's/\(LIB-.*\)_T[0-9]/\1/g'))
#tumor_prefix=($(echo $prefix | sed 's/\(LIB-.*\)_T[0-9]/\1/g'))
#tumor_prefix=($(echo $prefix | sed 's/\(LY_RAP_000[123]\)_.*/\1/g'))
#cont_sample="${tumor_samplename##*_}"
echo "tumor_prefix=$tumor_prefix"  #get the prefix of the sample to match the controls
#echo "control_sample=$cont_sample"
#cd $control_path
#normal_file=$(ls LY_RAP_${cont_sample}_Ctl*.bam)
control_file=$(ls $control_path/${tumor_prefix}_GL*.bam)
echo "control_file=$control_file"

#control_samplename=($(samtools view -H ${control_file} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
control_samplename=($(samtools view -H ${control_file} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo $control_samplename

###get the header info from the raw_bam files
#reorg_rg=$(samtools view -H "${work_path}/raw_bam/${prefix}.bam" | grep '@RG' | sed 's#\t#\\t#g'| head -1)

#somatic mutation calling
#https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136#ref4

##create interval files
gatk BedToIntervalList \
-I "$amplicon_interval_list" \
-O "${work_path}/coverage/${prefix}_amplicon.interval_list" \
-SD "$input_bam" && echo "** bed2intervallist $prefix done **"

gatk BedToIntervalList \
-I "$targets_interval_list" \
-O "${work_path}/coverage/${prefix}_targets.interval_list" \
-SD "$input_bam" && echo "** bed2intervallist targets $prefix done **"

#save new interval lists as variables which will be used as input for final picard function
amps=${work_path}/coverage/${prefix}_amplicon.interval_list
ints=${work_path}/coverage/${prefix}_targets.interval_list

gatk CollectTargetedPcrMetrics \
-I "$input_bam" \
-O "${work_path}/coverage/${prefix}.output_pcr_metrics.txt" \
-R "${ref_genome}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta" \
--PER_TARGET_COVERAGE "${work_path}/coverage/${prefix}.per_target_coverage.txt" \
--AMPLICON_INTERVALS "$amps" \
--TARGET_INTERVALS "$ints"  && echo "** CollectTargetedPcrMetrics $prefix done **"

##calling variants
gatk Mutect2 \
-R "${ref_genome}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta" \
-I "$input_bam"  \
-I "${control_file}"  \
-tumor ${tumor_samplename} \
-normal ${control_samplename} \
-L "${work_path}/coverage/${prefix}_targets.interval_list" \
-O "${work_path}/VCFs/${prefix}.raw.vcf.gz" \
--germline-resource "/ressources/af-only-gnomad.hg38.vcf.gz" && echo "** $prefix Mutect2 variant calling done **"

#Retrieve contamination and filter SNVs called by MuTect2
##https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136
gatk GetPileupSummaries \
-I "$input_bam" \
-V "/ressources/af-only-gnomad.hg38.vcf.gz"  \
-L "${work_path}/coverage/${prefix}_targets.interval_list" \
-O "${work_path}/VCFs/${prefix}.pileups.table" && echo "** $prefix pile up summaries done **"

gatk CalculateContamination \
-I "${work_path}/VCFs/${prefix}.pileups.table" \
-O "${work_path}/VCFs/${prefix}.contamination.table" && echo "** $prefix Calculation Contamination done **"

gatk FilterMutectCalls \
-R "${ref_genome}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta" \
-V "${work_path}/VCFs/${prefix}.raw.vcf.gz" \
--contamination-table "${work_path}/VCFs/${prefix}.contamination.table" \
-O "${work_path}/VCFs/${prefix}.un_filtered.vcf.gz"  && echo "** $prefix filtering VCFs done **"

gatk SelectVariants \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${work_path}/VCFs/${prefix}.un_filtered.vcf.gz" \
--exclude-filtered \
-O "${work_path}/VCFs/${prefix}.filtered.vcf.gz" && echo "** $prefix selectVariants done  **"

gatk SelectVariants \
-R "${ref_genome}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta" \
-V "${work_path}/VCFs/${prefix}.filtered.vcf.gz" \
-sn ${tumor_samplename} \
--exclude-non-variants true \
-O "${work_path}/VCFs/${prefix}.filtered.case_only.vcf.gz" && echo "** $prefix selectVariants done  **"


#annotate by VEP/98
vep --offline --dir_cache /projects//VEP/GRCh38/98 \
-i ${work_path}/VCFs/${prefix}.filtered.case_only.vcf.gz \
--force_overwrite --variant_class --sift b --polyphen b \
--nearest symbol --gene_phenotype --regulatory --show_ref_allele \
--tab -o ${work_path}/VEP_annot/${prefix}.grch38_vep_annot.txt && echo "** $prefix vep annnotation done  **"



#Normalizing and decomposing filtered VCF files
vt normalize "${work_path}/VCFs/${prefix}.filtered.case_only.vcf.gz" -r "${ref_genome}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta" -o "${work_path}/Annovar_annot/${prefix}.filtered_norm.vcf.gz" && echo "** $prefix Normalizing done  **"

#split multiallelic variants to biallelic
vt decompose -s "${work_path}/Annovar_annot/${prefix}.filtered_norm.vcf.gz" -o "${work_path}/Annovar_annot/${prefix}.filtered_norm_decomp.vcf.gz"  && echo "** $prefix decomposed done  **"

#Annotating variants with Annovar
#updated June 2022 --onetranscript
table_annovar.pl --buildver hg38 "${work_path}/Annovar_annot/${prefix}.filtered_norm_decomp.vcf.gz" /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic70,avsnp147,clinvar_20190305 --operation g,f,f,f,f \
--outfile "${work_path}/Annovar_annot/${prefix}.annot.norm_decomp_filtered"  --onetranscript --vcfinput && echo "** $prefix annovar annotation done  **"

fi
