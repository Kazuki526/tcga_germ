#!/usr/bin/perl
use strict;
use warnings;

my $exon_interval =3;
my $pwd = "/Volumes/areca42TB2/1000genomes/GRCh38";
(-e $pwd) or die "$pwd error\n";

#top_driver_genes = ~/git/driver_genes/onlytop105/top_driver105exon.bed
#control_genes = /Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_exon.bed
my $focus_bed = $ARGV[0];
(-e $focus_bed) or die "ERROR::inputed bed is not exist!\n";
my $focus_name = $ARGV[1];
#if(($focus_name eq "") && ($focus_bed =~ /\/(.+)\.bed$/)){$focus_name=$1;}
my ($vcf,$maf)=("$focus_name/ALL_$focus_name.vcf","$focus_name/ALL_$focus_name.maf");
if((!-e$vcf)||(!-e$maf)){die "ERROR::input error? $focus_name dir does not have vcf, like.vcf, maf\n";}

#read bed
my %bed=();
my @chr=();
open(BED,"$focus_bed");
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0] =~ s/^chr//;
		if($line[0] eq "X"){next;}
		if(defined $bed{$line[0]}){
				if($bed{$line[0]} =~ /-(\d+)$/){if($1 > $line[1]){die "please sort $focus_bed\n";}}
				$bed{$line[0]}.=":$line[1]-$line[2]";
		}else{
				$bed{$line[0]}="$line[1]-$line[2]";
		}
		if(!grep($_ eq $line[0],@chr)){push(@chr,"$line[0]");}
}
close BED;

#read maf
my %focal_maf=();
open(MAF,"$maf");
<MAF>;
my $maf_header = <MAF>;chomp $maf_header;
my %col = &header2hash($maf_header);
while(<MAF>){
		chomp;
		my @line = split(/\t/,);
		$focal_maf{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}=$_;
}
close MAF;

#read sample race information
my %race=("CHB" => "EAS", "JPT" => "EAS", "CHS" => "EAS", "CDX" => "EAS", "KHV" => "EAS",
		  "CEU" => "EUR", "TSI" => "EUR", "FIN" => "EUR", "GBR" => "EUR", "IBS" => "EUR",
		  "YRI" => "AFR", "LWK" => "AFR", "GWD" => "AFR", "MSL" => "AFR", "ESN" => "AFR",
		  "ASW" => "AFR", "ACB" => "AFR", "MXL" => "AMR", "PUR" => "AMR", "CLM" => "AMR",
		  "PEL" => "AMR", "GIH" => "SAS", "PJL" => "SAS", "BEB" => "SAS", "STU" => "SAS", "ITU" => "SAS");
my %sample=();
open(INFO,"nkf -Lu /Volumes/areca42TB2/1000genomes/GRCh38/1kg_sample_info.csv|") or die "cannot open 1kg_sample_info.csv\n";
<INFO>;
while(<INFO>){
		chomp;
		my @line=split(/,/,);
		if(!defined $race{$line[2]}){die "race:$_ is not existed\n";}
		if($race{$line[2]} eq "EUR"){$sample{$line[0]}=$race{$line[2]};}
}
close INFO;


open(OUT,"|gzip -c >$focus_name/white_by_sample_$focus_name.maf.gz");
print OUT "sample_id\tallele\t$maf_header\n";
foreach my $chr(@chr){
		if(!defined $bed{$chr}){next;}
		my @focus_region=split(/:/,$bed{$chr});
		my ($focus_st,$focus_end,$focus_region_num) =(0,0,0);
		($focus_st,$focus_end) = split(/-/,$focus_region[$focus_region_num]);
		my $file = "$pwd/genotype_vcf/ALL.chr$chr" . "_GRCh38.genotypes.20170504.vcf.gz";
		open(VCF,"gunzip -c $file|")or die "cant open $file\n";
		print "read and count $file\n";
		my @focal_sample=();
		my @header=();
		while(<VCF>){
				if($_=~/^##/){next;}
				if($_ =~ /^#CHROM/){
						@header = split(/\t/,);
						for(my$i=0;$i<scalar(@header);$i++){
								if(defined $sample{$header[$i]}){push(@focal_sample,$i);}
						}
						next;
				}
				my($l1,$l3)=("","");
				if($_ =~ /^$chr\t(\d+)\t\S+\t([ATGCatgc]+)\t/){($l1,$l3)=($1,$2);}else{die "ERROR::cannot match line\n$_\n";}
				if((scalar(@focus_region)-1 == $focus_region_num)&&($l1 > $focus_end)){last;}
				while($l1 > $focus_end){
						$focus_region_num++;
						($focus_st,$focus_end) = split(/-/,$focus_region[$focus_region_num]);
				}
				if(!(($focus_st - $exon_interval <= $l1+length($l3)-1)&&
						($l1 <= $focus_end + $exon_interval))){next;}
# have multi alt
				chomp;
				my @line=split(/\t/,);
				if($line[4] =~ /,/){
						my @genotype=();
						my @alt = split(/,/,$line[4]);
						for(my $i=0;$i<scalar(@alt);$i++){
								my($posi,$ref,$alt)=&posi_correct($line[1],$line[3],$alt[$i]);
								push(@genotype,"$chr:$posi:$ref:$alt");
						}
						foreach my $vcf_posi (@focal_sample){
								if($line[$vcf_posi] eq "0|0"){next;}
								my ($allele1,$allele2)=split(/\|/,$line[$vcf_posi]);
								if($allele1 ne "0"){print OUT "$header[$vcf_posi]\tallele1\t$focal_maf{$genotype[$allele1 -1]}\n"}
								if($allele2 ne "0"){print OUT "$header[$vcf_posi]\tallele2\t$focal_maf{$genotype[$allele2 -1]}\n"}
						}
# hve only one alt
				}else{
						my $genotype = "$chr:$line[1]:$line[3]:$line[4]";
						if((length($line[3])!=1)||(length($line[4])!=1)){
								my($posi,$ref,$alt)=&posi_correct($line[1],$line[3],$line[4]);
								$genotype="$chr:$posi:$ref:$alt";
						}
						foreach my $vcf_posi (@focal_sample){
								if($line[$vcf_posi] eq "0|0"){next;}
								my ($allele1,$allele2)=split(/\|/,$line[$vcf_posi]);
								if($allele1 ne "0"){print OUT "$header[$vcf_posi]\tallele1\t$focal_maf{$genotype}\n"}
								if($allele2 ne "0"){print OUT "$header[$vcf_posi]\tallele2\t$focal_maf{$genotype}\n"}
						}
				}
		}
		close VCF;
}

close OUT;


exit;



sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
sub posi_correct($ $ $){
		my($posi,$ref,$alt)=@_;
		while($ref and $alt and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
				$ref=substr($ref,1);$alt=substr($alt,1);$posi++;
		}
		while($ref and $alt and substr($ref,-1,1) eq substr($alt,-1,1) and $ref ne $alt){
				$ref=substr($ref,0,-1);$alt=substr($alt,0,-1);
		}
		if($ref eq ""){$ref ="-";$posi--;}
		if($alt eq ""){$alt ="-";}
		if((length($ref)!=1) && (length($alt)!=1)){print "WARNNING::position correct error??\n$posi:$ref-$alt\n$_\n";}
		return($posi,$ref,$alt);
}
