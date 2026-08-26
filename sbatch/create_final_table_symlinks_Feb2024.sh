current_path=$PWD
mkdir final_table
cd final_table

mkdir dcs_sc
mkdir sscs_sc
mkdir uniq_dcs
##create the symlinks for sscs_sc
cd sscs_sc
ln -s $current_path/tumour_normal_SNV_calling/sscs_sc/VCFs ./
cp /ctDNA_final_table/sscs_sc/* ./
#create the symlink
ln -s $current_path/tumour_normal_SNV_calling/sscs_sc/tumor_samples ./

##create the symlinks for dcs_sc
cd $current_path/final_table/dcs_sc
ln -s $current_path/tumour_normal_SNV_calling/dcs_sc/VCFs ./
cp /ctDNA_final_table/dcs_sc/* ./
ln -s $current_path/tumour_normal_SNV_calling/dcs_sc/tumor_samples ./

##create the symlinks for uniq_dcs
cd $current_path/final_table/uniq_dcs
cp /ctDNA_final_table/uniq_dcs/* ./
mkdir unique_dcs
cd unique_dcs
ln -s $current_path/tumour_normal_SNV_calling/unique_dcs/intersect/*intersect_uniq_dcs_updated_VAF.vcf.gz ./
