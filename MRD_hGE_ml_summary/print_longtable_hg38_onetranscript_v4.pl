#!/usr/bin/perl
#Ting 2023-06-06
#uptate the hg38 annotables
#hg38 annovar only has gene symbol, so annotable was adopted to map the ensGene_id back
#print out the longtable of mutation calls for the sake of plotting
##adding up the hGE_ml info
##Aug2023, update the cal of hGE_ml, read the avg VAF info, cal the hGE_ml per sample insdead of per mutation

use strict;
use warnings;
use Math::Complex;

open(FH,"one_transcript.annot.norm_decomp_filtered_VAF_table_sscs.txt");
open(hGE,"Yield_mod.txt");
open(avg_VAF,"avg_VAF_per_sample_sscs.txt");

##ouput files
open(OUT,">longtable_one_transcript.annot.norm_decomp_filtered_VAF_sscs_patients_with_hGE_per_sample_v4.txt");

print OUT "Sample_id\tpatient_ID\ttime_point\t";
print OUT "CHR\tPOS\tREF\tALT\tAC\tAN\tVCF_SN\tEnsembl_id\tGene.ensGene\tFunc.ensGene\tExonicFunc.ensGene\tAAChange.ensGene\tavsnp147\tcosmic70\tVAF\tinput_yield\tinput_plasma_volume\thGE_ml_mutation\tnum_mutation\tavg_VAF\thGE_ml_sample\n";

#print OUT "CHR\tPOS\tREF\tALT\tAC\tAN\tVCF_SN\tEnsembl_id\tGene.ensGene\tFunc.ensGene\tExonicFunc.ensGene\tAAChange.ensGene\tavsnp147\tcosmic70\tSample_id\tpatient_ID\ttime_point\tVAF\tinput_yield\tinput_plasma_volume\thGE_ml_mutation\tnum_mutation\tavg_VAF\thGE_ml_sample\n";

my @sample_ids;
my $sample_num=0;
my %sample_hash=();
my %sample_avg_VAF=();

while(<avg_VAF>){
  chomp;
    my ($patient_id, $num_mutation, $avg_VAF) = split('\t');

    if(!exists $sample_avg_VAF{$patient_id}){
       $sample_avg_VAF{$patient_id}={'num_mutation'=>$num_mutation,'avg_VAF'=>$avg_VAF};
       }#if
    elsif(exists $sample_avg_VAF{$patient_id}){
      print "$patient_id exists\n";
     }#elsif
}#while

while(<hGE>){
  chomp;
    my ($patient_id, $OICR_id, $timepoint,$Input_Plasma_Volume, $DNA_Initial_Yield) = split('\t');
#    $DNA_Initial_Yield = s/(^s+|s+$)//g; #get rid of the ^M
    my $tempid=$patient_id.'_'.$timepoint;

    if(!exists $sample_hash{$tempid}){
#       print "$tempid create\n";
       $sample_hash{$tempid}={'timepoint'=>$timepoint,'patient_ID'=>$patient_id,'dna_yield'=>$DNA_Initial_Yield,'input_plasma_volume'=>$Input_Plasma_Volume};
       }#if
    elsif(exists $sample_hash{$tempid}){
      print "$tempid exists\n";
     }#elsif
}#while

while (my $seq=<FH>){
     chomp($seq);

    if($seq=~/^#CHR/){ #the header line to get the sample id
#      print OUT "$seq\n";
      my @header_line= split("\t",$seq);

      for(my $i=14; $i<=$#header_line; $i++){ #sample ids
#      print OUT "$header_line[$i]\n";
        $sample_ids[$sample_num]=$header_line[$i];
        $sample_num++;
           }#for
       }#if header

    elsif($seq!~/^#/){ ##elsif vcf records
     my @vcf_record= split("\t",$seq);
     my $col_num= @vcf_record;
     my @line_info=();

     for (my $i=0; $i<= $col_num-1; $i++){  #there is a space at the end of the sample ids

      if($i<14){
       $line_info[$i]=$vcf_record[$i];
  #     print OUT "$vcf_record[$i]\t$line_info[$i]\t";
        }#if
      elsif($i>=14){
#        if(($vcf_record[$i]!~/\./)&&($vcf_record[$i]!~/\s/)){ ##if
        if(($vcf_record[$i] ne '.')&&($vcf_record[$i] ne ' ')){ ##if genotypes
         my $sample_index=$i-14;
         my ($patient_ID)=($sample_ids[$sample_index] =~ /(.*)_T\d/);
         my ($time_point)=($sample_ids[$sample_index] =~ /.*_(T\d)/);
        print OUT "$sample_ids[$sample_index]\t$patient_ID\t$time_point\t";

         for (my $j=0; $j<=$#line_info; $j++){ #print out the variant info

          print OUT "$line_info[$j]\t";   #loop print the variants
            }#for first cols;
  #        my $sample_index=$i-14;
        #  print OUT "$sample_index\t";
          my $input_yield;
          my $input_plasma_volume;
          my $hGE_ml_ctDNA;
          my $cf_DNA;
          my $num_mutation;
          my $avg_VAF;
          my $hGE_ml_per_sample;
      if(exists $sample_hash{$sample_ids[$sample_index]}){  #if hge ML exists
          $input_yield=$sample_hash{$sample_ids[$sample_index]}->{'dna_yield'};
          $input_plasma_volume=$sample_hash{$sample_ids[$sample_index]}->{'input_plasma_volume'};

          $num_mutation=$sample_avg_VAF{$sample_ids[$sample_index]}->{'num_mutation'};
          $avg_VAF=$sample_avg_VAF{$sample_ids[$sample_index]}->{'avg_VAF'};

          if(($input_yield eq 'NA')||($input_plasma_volume eq 'NA')){  ##if any of these two is NA, then NA
            $cf_DNA='NA';
            $hGE_ml_ctDNA='NA';
           }
          else{
           $cf_DNA=($input_yield*1000)/($input_plasma_volume*0.001);
           $hGE_ml_ctDNA=log10($cf_DNA*$vcf_record[$i]/3.3);  #cal per mutation hGE_ml
           $hGE_ml_per_sample=log10($cf_DNA*$avg_VAF/3.3);  #cal the per sample hge_ml
           }#else non NA
         }#if
      elsif(!exists $sample_hash{$sample_ids[$sample_index]}){
           $input_yield='NA';
           $input_plasma_volume='NA';
           $hGE_ml_ctDNA='NA';
            }#elsif
        print OUT "$sample_ids[$sample_index]\t$patient_ID\t$time_point\t$vcf_record[$i]\t$input_yield\t$input_plasma_volume\t$hGE_ml_ctDNA\t$num_mutation\t$avg_VAF\t$hGE_ml_per_sample\n";
        print OUT "$vcf_record[$i]\t$input_yield\t$input_plasma_volume\t$hGE_ml_ctDNA\t$num_mutation\t$avg_VAF\t$hGE_ml_per_sample\n";

          }#if genotypes
         }#elsif genotypes
       }#for
   }#elsif

}#while
