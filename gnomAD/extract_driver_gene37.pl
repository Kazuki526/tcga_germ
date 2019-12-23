#!/usr/bin/perl
use warnings;
use strict;

############### this script extract only noncancer maf file ################

my $gnomAD_version = "gnomad.exomes.r2.1.1.sites";

# this script execute in gnomAD exon raw file (by chromosome) dir
my @ls = `ls file_grch37|grep $gnomAD_version`;
if(scalar(@ls)!= 24){die "ERROR::doing on wrong dir?\n";}

my @chr =qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my $maf_path="maf37/non_cancer_maf/non_cancer_chr";
if(!-e "$maf_path"."1.maf.gz"){die "ERROR::file exist error $maf_path is correct?\n";}

#read driver genes list
my %top_driver_genes =();
open(DG,"$ENV{HOME}/git/driver_genes/driver_genes.tsv") or die "ERROR::cannot open top driver genes file.\n";
<DG>;
while(<DG>){
		chomp;
		my @line = split(/\t/,);
		$top_driver_genes{$line[0]}=$line[1]; #if ref>=4 top driver genes
}
close DG;
open(CG,"/Volumes/areca42TB/GRCh38_singlefasta/control_genes.tsv") or die "ERROR::cannot open top driver genes file.\n";
<CG>;
while(<CG>){
		chomp;
		my @line = split(/\t/,);
		$top_driver_genes{$line[0]}=0;
}
close CG;

open(TDG,">maf37/non_cancer_maf/non_cancer_top_driver_gene.maf");
open(ODG,">maf37/non_cancer_maf/non_cancer_other_driver_gene.maf");
open(CTG,">maf37/non_cancer_maf/non_cancer_control_gene.maf");
foreach my $chr(@chr){
		&pull_focal_site_from_maf($chr);
}
close TDG;
close ODG;
close CTG;



sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}

sub pull_focal_site_from_maf( $ ){
		my $chr = $_[0];
		open(MAF,"gunzip -c $maf_path$chr.maf.gz|") or die "ERROR::cannot open annotation_extract of chr$chr\n";
		my $header = <MAF>;chomp$header;
		my %col=&header2hash($header);
		if($chr eq "1"){print TDG "$header\n";}
		if($chr eq "1"){print ODG "$header\n";}
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				my $gene = $line[$col{SYMBOL}];
				if(!defined $top_driver_genes{$gene}){next;}
				if($top_driver_genes{$gene} >=4){
						print TDG "$_\n";
				}elsif($top_driver_genes{$gene} >0){
						print ODG "$_\n";
				}else{
						print CTG "$_\n";
				}
		}
		close MAF;
}
