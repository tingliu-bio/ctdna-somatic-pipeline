
#module load vcftools/0.1.15
work_path=$PWD
#cd ${work_path}

outfile=$PWD/"coverage_per_sample_sscs_sc.txt"

#file=/tumour_normal_SNV_calling/sscs_sc/VCFs/*.sorted.sscs.sc.filtered.case_only.vcf.gz
#files_location=/tumour_normal_SNV_calling/sscs_sc/coverage
files_location=./coverage
file=${files_location}/*.per_target_coverage.txt

echo -n "sample_id"  >>$outfile
echo -n -e "\t" >>$outfile
echo "sscs_mean_coverage" >>$outfile

for i in ${file}
do
#  echo $i
if [[ $i =~ ${files_location}/(.*).sorted.sscs.sc.per_target_coverage.txt ]] ;then
  num=${BASH_REMATCH[1]}
   echo -n $num  >>$outfile
   echo -n -e "\t" >>$outfile
#    if $num
    mean_cov=`tail -721 $i | awk '{total+=$7};END {print total/721}'`
    echo $mean_cov >>$outfile

    fi
done
