#!/usr/bin/perl
#Ting Liu 2023
#updated Nov-2023
#check the id mapping file then combined the fastq files blong to a same sample
#if only one fq file, then create the symlink, if more than one fq, then combine them
#use strict;
use warnings;

open(FH,"ctNDA.txt");
open(OUT,">lib_tmp");

my %sample_hash=();

while(<FH>){
  chomp;
  if ($_=~ /^LIB|2522/){
    my ($patient,$time_point, $fq1_path) = split('\s');
    my   $line_index={'patient_id'=>$patient,'time_point'=>$time_point,'fq1_file'=>$fq1_path};

    $sample_id=$patient.'_'.$time_point;
#      print OUT "$patient\t$time_point\t$sample_id\n";
    if(!exists $sample_hash{$sample_id}){
      $sample_hash{$sample_id}=[];
      push(@{$sample_hash{$sample_id}}, $line_index);
     }#if
    elsif(exists $sample_hash{$sample_id}){
      push(@{$sample_hash{$sample_id}}, $line_index);
     }#elsif
  }#if
}#while


foreach my $out_lib (keys %sample_hash){  ##for each $libray
#  print "$out_lib\n";
  my @fq1_arr=();
  my @fq2_arr=();
  $len= @{$sample_hash{$out_lib}};

  if($len==1) {
     my $fq_file=$sample_hash{$out_lib}->[0]->{'fq1_file'};
     my ($header) = ($fq_file =~ /(.*)_R1.fastq.gz/);
     my  $fq2_path=$header.'_R2.fastq.gz';
     my  $fq1_path=$header.'_R1.fastq.gz';  #add this line of code in case there is ^M sign at the end of the input line
#   print OUT "$fq_file\t$header\n";
     qx/ln -s $fq1_path ${out_lib}_R1.fastq.gz/;  #only one sample , create the symplink
     qx/ln -s $fq2_path ${out_lib}_R2.fastq.gz/;
     }

  else{

  for (my $i=0; $i<$len; $i++){  ##for each lib, there are one or more lanes
    my $fq_file=$sample_hash{$out_lib}->[$i]->{'fq1_file'};

    my ($header) = ($fq_file =~ /(.*)_R1.fastq.gz/);
    my  $fq2_path=$header.'_R2.fastq.gz';
    my  $fq1_path=$header.'_R1.fastq.gz';  #add this line of code in case there is ^M sign at the end of the input line
     push(@fq1_arr, $fq1_path);  #arry to store the fq1 path
     push(@fq2_arr, $fq2_path);
#  print "$out_lib\t$sample_hash{$out_lib}->[$i]->{'sampleid'}\t$sample_hash{$out_lib}->[$i]->{'Patient_ID'}\t$sample_hash{$out_lib}->[$i]->{'run_id'}\t$sample_hash{$out_lib}->[$i]->{'time_point'}\t$fq1_path\t$fq2_path\n";
  }#for
 my $fq_len=@fq1_arr;
 my $out_fq1='';
 my $out_fq2='';

 for(my $j=0; $j<$fq_len; $j++){
   $out_fq1=$out_fq1."\t".$fq1_arr[$j];
   $out_fq2=$out_fq2."\t".$fq2_arr[$j];
    }#for
#  print "$out_lib\t$sample_hash{$out_lib}->[$i]->{'sampleid'}\t$sample_hash{$out_lib}->[$i]->{'Patient_ID'}\t$sample_hash{$out_lib}->[$i]->{'run_id'}\t$sample_hash{$out_lib}->[$i]->{'time_point'}\t";
 print OUT "$out_lib\t$len\t$fq_len\t";
 print OUT "$out_fq1\t";
 print OUT "$out_fq2\n";
# $merged_fq_path= "./comined_fqs/".$out_lib ;
 qx/cat $out_fq1 > ${out_lib}_R1.fastq.gz/;
 qx/cat $out_fq2 > ${out_lib}_R2.fastq.gz/;


 }#else more than fq
}#foreach lib
