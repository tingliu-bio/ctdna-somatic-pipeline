module load samtools
work_path=$PWD

outfile="avg_VAF_per_sample_sscs.txt"
files_location=${work_path}/VCFs
file=${files_location}/*.sorted.sscs.sc.filtered.case_only.vcf.gz

for i in ${file}
do
#  echo $i
if [[ $i =~ ${files_location}/(.*).sorted.sscs.sc.filtered.case_only.vcf.gz ]] ;then
  num=${BASH_REMATCH[1]}
    bgzip -d -c $i >${files_location}/${num}.sorted.sscs.sc.filtered.case_only.vcf ##unzip the gz file but  keep original files unchanged
    echo -n $num  >>$outfile #print the sample id
    echo -n -e "\t" >>$outfile
#    if $num
     num_lines=`grep '^chr' ${files_location}/${num}.sorted.sscs.sc.filtered.case_only.vcf |wc -l`
     echo -n $num_lines >>$outfile   #print the number of variants
     echo -n -e "\t" >>$outfile
#   echo $out_num
    if [ $num_lines -ne '0' ]; then  #if the value exists, means not 0 variant
      sum_AF=`grep '^chr' ${files_location}/${num}.sorted.sscs.sc.filtered.case_only.vcf |awk '{print $10}' | awk -F ':' '{total+=$3};END {print total}'`
 #     echo -n $sum_AF
#      mean_AF=`echo "s($sum_AF/$num_lines)" | bc -l`
#      echo -n $mean_AF
#      echo -n -e "\t"
      mean_AF=`grep '^chr' ${files_location}/${num}.sorted.sscs.sc.filtered.case_only.vcf |awk '{print $10}' | awk -F ':' '{total+=$3};END {print total/NR}'`
      echo $mean_AF >>$outfile
    else
      echo "0"   >>$outfile
    fi

    fi
done
