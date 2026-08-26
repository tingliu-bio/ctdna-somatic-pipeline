#!/bin/bash
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=46000M
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH -J merge_dcs_sc

module load vcftools
module load samtools/1.20
module load annovar/20180416
module load tabix
module load vt/0.577
module load vep/113

work_path=$SLURM_SUBMIT_DIR
##the slurm output will be written to your current directory

cd ${work_path}

if [ ! -d "${work_path}/merged_VCF" ]; then
        mkdir -p "${work_path}/merged_VCF"
fi

merged_path=${work_path}/merged_VCF
gatk_bundle="/resources/GATK_bundle/homo_sapiens_GRCh38"

outname="dcs.sc_all_samples_combined"
file="${work_path}/VCFs/*sorted.dcs.sc.filtered.case_only.vcf.gz"

bcftools merge ${work_path}/VCFs/*sorted.dcs.sc.filtered.case_only.vcf.gz -Oz -0 -o ${merged_path}/${outname}.filtered.vcf.gz  && echo "** merge vcfs done **"

#annotate by VEP/98
vep --offline --dir_cache /resources/VEP/GRCh38/98 \
-i ${merged_path}/${outname}.filtered.vcf.gz \
--force_overwrite --variant_class --sift b --polyphen b \
--nearest symbol --gene_phenotype --regulatory --show_ref_allele \
--tab -o ${merged_path}/${outname}.grch38_vep_annot.txt && echo "** $prefix vep annnotation done  **"


#Normalizing and decomposing filtered VCF files
vt normalize "${merged_path}/${outname}.filtered.vcf.gz" -r "${gatk_bundle}/Homo_sapiens_assembly38.fasta" -o "${merged_path}/${outname}.filtered_norm.vcf.gz" && echo "** $outname Normalizing done  **"

#split multiallelic variants to biallelic
vt decompose -s "${merged_path}/${outname}.filtered_norm.vcf.gz" -o "${merged_path}/${outname}.filtered_norm_decomp.vcf.gz"  && echo "** $outname decomposed done  **"

#Annotating variants with Annovar
table_annovar.pl --buildver hg38 "${merged_path}/${outname}.filtered_norm_decomp.vcf.gz" /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic70,avsnp147,clinvar_20190305 --operation g,f,f,f,f \
--outfile "${merged_path}/${outname}.annot.norm_decomp_filtered" --onetranscript --vcfinput && echo "** $outname annovar annotation done  **"

#java -jar $picard_dir/picard.jar MergeVcfs "${sample_vcfs[*]}" O=${work_path}/${outname}.filtered.vcf.gz  && echo "** merge vcfs done **"
