#!/usr/bin/perl
#uptate the hg38 annotables
#hg38 annovar only has gene symbol, so annotable was adopted to map the ensGene_id back
#this script is to re-format the vcf ouput from annovar, also adding the ensGene_id

use strict;
use warnings;


#open (FH,"one_transcript.annot.norm_decomp_filtered.hg38_multianno.vcf") or die "cannot open file";
#open (FH,"sscs.sc_all_samples_combined_onetranscript_only.hg38_multianno.vcf") or die "cannot open file";
open (FH,"./merged_VCF/sscs.sc_all_samples_combined.annot.norm_decomp_filtered.hg38_multianno.vcf") or die "cannot open file";
#open(ANNOT,"/resources/annotables/annotable_grch38_March2022.txt");
open(ANNOT,"/master_files/annotable_grch38_March2022.txt");
##annovar hg38 doesn't have ensg id

#open (DM,">Venn_diagram_table.txt");
#open (SEE, ">see");
open(OUT,">one_transcript.annot.norm_decomp_filtered_VAF_table_sscs.txt");

my %annot_hash=();

##create the hash for annotable first
while(my $seq=<ANNOT>){
  chomp($seq);

  if($seq=~/^ENSG/){ #if its not the header line
    my @line= split(' ',$seq);
#      print SEE "$line[0]\t$line[2]\n";
      if(!exists $annot_hash{$line[2]}){
        $annot_hash{$line[2]}=$line[0];
  #      print SEE "$annot_hash{$line[0]}\t$line[2]\n";
        }#if !EXISTS
       else{
#        print SEE "duplicate records $line[0]\n";
       }#else
    }#if not header
}#while

my $vcf_sn=0;

while (my $seq=<FH>){
     chomp($seq);
     next if($seq=~/^##/);#ignore the header

    if($seq=~/^#CHROM/){ #the header line to get the sample id
#      print OUT "$seq\n";
      my @header_line= split("\t",$seq);
      print OUT "#CHR\tPOS\tREF\tALT\tAC\tAN\tVCF_SN\tEnsembl_id\tGene.ensGene\tFunc.ensGene\tExonicFunc.ensGene\tAAChange.ensGene\tavsnp147\tcosmic70\t";
      for(my $i=9; $i<=$#header_line; $i++){ #sample ids
        print OUT "$header_line[$i]\t";
#        print DM "$header_line[$i]\t"; #print the sample id for venn table
         }#for
        print OUT "\n";
#        print DM "\n";
       }#if header

    elsif($seq!~/^\#/){ ##elsif vcf records
     my @vcf_record= split("\t",$seq);
     my $chrom= $vcf_record[0];
     my  @INFO= split(";",$vcf_record[7]);
         my $INFO_len= @INFO;
     my (%INFO_item) =();
       for (my $i= 0; $i<= $INFO_len-1; $i++){ #parsing the INFO item of each ID, storing each ID=desc to a pair of hash value

         my $info_id='';
         my $info_desc='';

        if(($info_id,$info_desc) = ($INFO[$i] =~ /(.*)\=(.*)/)){
        if(!exists $INFO_item{$info_id}){ #if not exists then create the hash
            $INFO_item{$info_id}=$info_desc;
         }#if
        elsif(exists $INFO_item{$info_id}){
            #$INFO_item{$info_id}='';
                next;
            }#else
          }#if info_id
        }#for each INFO col
  $vcf_sn=$chrom.'_'.$vcf_record[1].'_'.$vcf_record[3].'_'.$vcf_record[4]; #keep the running id to cal the venn diagram

print OUT "$chrom\t$vcf_record[1]\t$vcf_record[3]\t$vcf_record[4]\t" ;
#print DM "$vcf_sn\t";
#print OUT "AC=$INFO_item{'AC'}\tAN=$INFO_item{'AN'}\tFunc.ensGene=$INFO_item{'Func.ensGene'}\tGene.ensGene=$INFO_item{'Gene.ensGene'}\tExonicFunc.ensGene=$INFO_item{'ExonicFunc.ensGene'}\tAAChange.ensGene=$INFO_item{'AAChange.ensGene'}\t";
  my $Ensembl_id=$annot_hash{$INFO_item{'Gene.ensGene'}};

 if ($INFO_item{'Gene.ensGene'}=~/\\x/){  ##some annots have more than one genes, take the first one
  my ($ensid)=($INFO_item{'Gene.ensGene'} =~ /(.*)?\\x/);
   $Ensembl_id=$annot_hash{$ensid};
#   print "$ensid\t$Ensembl_id\n";
   }
  elsif(! exists $annot_hash{$INFO_item{'Gene.ensGene'}}){ ##if there is no ennsembl id found
    $Ensembl_id='-';

    }

print OUT "$INFO_item{'AC'}\t$INFO_item{'AN'}\t$vcf_sn\t$Ensembl_id\t$INFO_item{'Gene.ensGene'}\t$INFO_item{'Func.ensGene'}\t$INFO_item{'ExonicFunc.ensGene'}\t$INFO_item{'AAChange.ensGene'}\t$INFO_item{'avsnp147'}\t$INFO_item{'cosmic70'}\t";
 for(my $j=9; $j<=$#vcf_record; $j++){
#  print OUT "$vcf_record[$i]\t";
    my @GT=split(":",$vcf_record[$j]);
  #  print OUT "$GT[0]:$GT[2]\t";
    print OUT "$GT[2]\t";
    if($GT[0]=~/0\/0/){ ##print the variant table for venn diagram
#      print DM "0\t";
     }
     else{
#      print DM "$vcf_sn\t";
     }
   }#for
   print OUT "\n";
#   print DM "\n";
   }#elsif

}#while
