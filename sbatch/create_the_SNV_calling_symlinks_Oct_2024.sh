#Jan 2024-add the uniq_dcs folder
current_path=$PWD
mkdir tumour_normal_SNV_calling
cd tumour_normal_SNV_calling

mkdir dcs_sc
mkdir sscs_sc
mkdir unique_dcs

##create the symlinks for sscs_sc
cd sscs_sc
mkdir germline_controls
mkdir tumor_samples
cd germline_controls
ln -s $current_path/consensus_output_Normals/*/sscs_sc/*_GL.sorted.sscs.sc.sorted.bam ./
cd ../tumor_samples
ln -s $current_path/consensus_output_Tumours/*/sscs_sc/*.sorted.sscs.sc.sorted.bam ./
cd ..

##cp the GATK script from the Ting_resources folder and create the bam list files for both GL and tumor samples
cp /projects/GATK_mutect2_calling2annot_case_control_grch38_ctDNA_v6_Oct2024_sscs.sh ./
ls $PWD/tumor_samples/*sscs.sc.sorted.bam >bam_files_sscs.txt
##update the arry id based on the number of input bam files

##create the symlinks for dcs_sc
cd $current_path/tumour_normal_SNV_calling/dcs_sc
mkdir germline_controls
mkdir tumor_samples
cd germline_controls
ln -s $current_path/consensus_output_Normals/*/dcs_sc/*_GL.sorted.dcs.sc.sorted.bam ./
cd ../tumor_samples
ln -s $current_path/consensus_output_Tumours/*/dcs_sc/*.sorted.dcs.sc.sorted.bam ./
cd ..

##cp the GATK script from the Ting_resources folder and create the bam list files for both GL and tumor samples
cp /projects/GATK_mutect2_calling2annot_case_control_grch38_ctDNA_v6_Oct2024_dcs.sh ./
ls $PWD/tumor_samples/*sorted.dcs.sc.sorted.bam >bam_files_dcs.txt


##create the symlinks for uniq_dcs
#This section aims to adjust the VAF of the DCS.
cd $current_path/tumour_normal_SNV_calling/unique_dcs
mkdir germline_controls
mkdir tumor_samples
cd germline_controls
ln -s $current_path/consensus_output_Normals/*/dcs_sc/*all.unique.dcs.sorted.bam ./
cd ../tumor_samples
ln -s $current_path/consensus_output_Tumours/*/dcs_sc/*all.unique.dcs.sorted.bam ./
cd ..

##cp the GATK script from the Ting_resources folder and create the bam list files for both GL and tumor samples
cp /projects/GATK_mutect2_calling2annot_case_control_primary_grch38_ctDNA_v6_Apr2024_force_call.sh ./
cp /projects/cal_intersect_with_VCF_header.sh ./
cp /projects/cal_avg_VAF_per_sample_updated_uniq.sh ./
ls $PWD/tumor_samples/*.sorted.bam >bam_files_dcs_unique.txt
##update the arry id based on the number of input bam files
