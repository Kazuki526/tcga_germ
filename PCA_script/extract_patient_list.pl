#!/usr/bin/perl
use strict;
use warnings;

my $pca_vcf = "row_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz";
-e $pca_vcf or die "ERROR::$pca_vcf donot exist. doing on wrong dir?\n";

open(VCF,"gunzip -c $pca_vcf|");
open(OUT,">patient/PCA_patient_list.txt");
print OUT "patient_id\tsample_id\n";
while(<VCF>){
		if($_ !~ /^#CHROM/){next;}
		chomp;
		my @line = split(/\t/,);
		for(my $t = 9; $t < scalar(@line);$t++){
				my $pid;
				if($line[$t] =~ /^(TCGA-..-....)-/){$pid = $1;}else{die "ERROR::what patient $line[$t]\n";}
				print OUT "$pid\t$line[$t]\n";
		}
		last;
}
close VCF;
close OUT;

