#!/usr/bin/perl
use strict;
use warnings;

my $pwd = `pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc/pancan_atlas_germ"){die "ERROR::done on wrong dir\n";}

#read driver genes list
my %top_driver_genes =();
open(TDG,"$ENV{HOME}/git/innanlab/driver_genes.tsv") or die "ERROR::cannot open top driver genes file.\n";
<TDG>;
while(<TDG>){
		chomp;
		my @line = split(/\t/,);
		$top_driver_genes{$line[0]}=$line[1]; #if ref>=4 top driver genes
}
close TDG;

#my $vcf_file = "raw_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz";
my $vcf_file = "raw_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf";
-e $vcf_file or die "ERROR::not exist $vcf_file\n";


my @ann_list = qw(Chromosome Start_Position End_Position Strand Reference_Allele Tumor_Seq_Allele2 Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons CANONICAL SIFT PolyPhen Feature);
open(TDGANN,">annotation_extract/top_driver_genes.maf");
open(ODGANN,">annotation_extract/other_driver_genes.maf");
print TDGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
print ODGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";

for(my $i=1;$i<23;$i++){my %test=&pull_focal_site_from_maf($i);}
&pull_focal_site_from_maf("X");

close TDGANN;
close ODGANN;

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
		open(EXANN,">annotation_extract/chr$chr.maf") or die "ERROR::cannot open annotation_extract\n";
		print EXANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
		open(MAF,"test/chr$chr"."_snp.maf") or die "ERROR::cannot open chr$chr maf file";
		<MAF>;
		my $header = <MAF>;chomp$header;
		if($header !~ /^Hugo_Symbol\t/){die "ERROR::maf_snp header reading miss\n";}
		my %col=&header2hash($header);
		my %maf_focal =();
		while(<MAF>){
				chomp;
				my $gene;
				if($_=~ /^(\S+)\t/){$gene=$1;}else{next;}
				if(!defined $top_driver_genes{$gene}){next;}
				if($top_driver_genes{$gene} <4){next;}
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
#				my $gene = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{gene} = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{IMPACT} = $line[$col{IMPACT}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{Consequence} = $line[$col{Consequence}];
				my $ann="$gene";
				foreach my $ann_list (@ann_list){
						$ann.="\t$line[$col{$ann_list}]";
				}
				if(defined $top_driver_genes{$gene}){
						if($top_driver_genes{$gene} >=4){
								print TDGANN "$ann\n";
						}else{
								print ODGANN "$ann\n";
						}
				}else{
						print EXANN "$ann\n";
				}
		}
		close MAF;
		open(MAF,"test/chr$chr"."_indel.maf") or die "ERROR::cannot open chr$chr maf file";
		<MAF>;
		$header = <MAF>;chomp$header;
		if($header !~ /^Hugo_Symbol\t/){die "ERROR::maf_indel header reading miss\n";}
		%col=&header2hash($header);
		while(<MAF>){
				chomp;
				my $gene;
				if($_=~ /^(\S+)\t/){$gene=$1;}else{next;}
				if(!defined $top_driver_genes{$gene}){next;}
				if($top_driver_genes{$gene} <4){next;}
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
#				my $gene = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{gene} = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{IMPACT} = $line[$col{IMPACT}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{Consequence} = $line[$col{Consequence}];
				my $ann="$gene";
				foreach my $ann_list (@ann_list){
						$ann.="\t$line[$col{$ann_list}]";
				}
				if(defined $top_driver_genes{$gene}){
						if($top_driver_genes{$gene} >=4){
								print TDGANN "$ann\n";
						}else{
								print ODGANN "$ann\n";
						}
				}else{
						print EXANN "$ann\n";
				}
		}
		close EXANN;
		return(%maf_focal);
}

sub indel_posi_edit( $ $ $ ){
		my ( $pos, $ref, $var ) = @_;
		my ( $ref_length, $var_length ) = ( length( $ref ), length( $var ));
		while( $ref and $var and substr( $ref, 0, 1 ) eq substr( $var, 0, 1 ) and $ref ne $var ) {
				( $ref, $var ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var );
				--$ref_length; --$var_length; ++$pos;
		}
		$pos =($ref eq "-" ? $pos-1 : $pos);
		return($pos,$ref,$var);
}

sub print_out ( $ $ ){
		my($out,$gene)=@_;
		if(defined $top_driver_genes{$gene}){
				if($top_driver_genes{$gene} >= 4){
						print TDGOUT "$out";
				}else{
						print ODGOUT "$out";
				}
		}else{
				print EXOUT "$out";
		}
}
sub print_vcf ( $ $ ){
		my($vcf,$gene)=@_;
		if(defined $top_driver_genes{$gene}){
				if($top_driver_genes{$gene} >= 4){
						print TDGOUT "$vcf";
				}else{
						print ODGOUT "$vcf";
				}
		}else{
				print EXOUT "$vcf";
		}
}


