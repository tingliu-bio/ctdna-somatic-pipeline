current_path=$PWD
mkdir T1_mutations
cd T1_mutations
ln -s $current_path/VCFs/*_T1.sorted.sscs.sc.filtered.case_only.vcf.gz ./
ln -s ../tumour_normal_SNV_calling ./
