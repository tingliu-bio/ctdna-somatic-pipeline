
#module load vcftools/0.1.15
work_path=$PWD
#cd ${work_path}

outfile=$PWD/"col_coverage_per_sample_sscs_sc.txt"

#file=/tumour_normal_SNV_calling/sscs_sc/VCFs/*.sorted.sscs.sc.filtered.case_only.vcf.gz
#files_location=/tumour_normal_SNV_calling/sscs_sc/coverage
#files_location=./SNV_calling_per_target_coverage_consensus_crunchered
files_location=./coverage

file=${files_location}/*.per_target_coverage.txt

#echo -n "sample_id"  >>$outfile
#echo -n -e "\t" >>$outfile
#echo "mean_coverage" >>$outfile
#header=$(echo 'chrom-pos-header.txt')
header='chrom-pos-header.txt'
#header=`cat ${header}`
#echo "${header}" > matrix
cp chrom-pos-header.txt matrix

for i in ${file}
do
#  echo $i
if [[ $i =~ ${files_location}/(.*).sorted.sscs.sc.per_target_coverage.txt ]] ;then
  num=${BASH_REMATCH[1]}
#   echo -n $num  >>$outfile
#   echo -n -e "\t" >>$outfile
#    if $num
  echo $num> sample_id
#   mean_cov_col=`tail -721 $i | awk '{print $7}'`

#    cov_col=`paste -d "\n" $num $mean_cov_col`
#echo $header
#echo $cov_col
 tail -721 ${i} | awk '{print $7}' >mean_cov
  cat sample_id mean_cov >mean_cov_col
  paste matrix mean_cov_col >matrix_tmp
  cp matrix_tmp matrix
#    header=`paste $header mean_cov_col`
#    mean_cov=`tail -721 $i | awk '{total+=$7};END {print total/721}'`
#    echo $mean_cov >>$outfile
#rm mean_cov_col
#rm sample_id
    fi
done

cp matrix col_coverage_per_sample_sscs_sc.txt
rm mean_cov_col
rm sample_id
rm matrix
rm matrix_tmp
#echo "${header}" > $outfile
