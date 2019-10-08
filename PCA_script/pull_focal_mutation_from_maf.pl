#!/usr/bin/perl
use strict;
use warnings;

my $pwd = `pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc/pancan_atlas_germ"){die "ERROR::done on wrong dir\n";}

#read driver genes list
%top_driver_genes =();
open(TDG,"$ENV{HOME}/git/innanlab/driver_genes.tsv") or die "ERROR::cannot open top driver genes file.\n";
<TDG>;
while(<TDG>){
		chomp;
		my @line = split(/\t/,);
		$top_driver_genes{$line[0]}=$line[1]; #if ref>=4 top driver genes
}
close TDG;

my $vcf_file = "raw_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz";
-e $vcf_file or die "ERROR::not exist $vcf_file\n";


open(VCF,"gunzip -c $vcf_file|");
open(TDGOUT,"|gzip -c >maf_patient/top_driver_genes_patient.maf");
open(ODGOUT,"|gzip -c >maf_patient/other_driver_genes_patient.maf");
my @out_list = qw(Chromosome Start_Position End_Position Strand Reference_Allele Tumor_Seq_Allele1 Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons CANONICAL SIFT PolyPhen Feature);
open(TDGANN,">annotation_extract/top_driver_genes.maf");
open(ODGANN,">annotation_extract/other_driver_genes.maf");
print TDGANN "Hugo_Symbol\t".join("\t",@out_list)."\n";
print ODGANN "Hugo_Symbol\t".join("\t",@out_list)."\n";
my $chr="";
my %focal_site=();
my @col_info=();
while(<VCF>){
		if($_ =~ /^##/){next;}
		chomp;
		if($_ =~ /^#/){$_ =~ s/#//;@col_info=split(/\t/,);}
		my @line = split(/\t/,);
		if($line[0] ne $chr){
				if($chr ne ""){
						print "done\n";$|=1;
						close EXVCF;
						close EXANN;
				}
				$chr =$line[0];
				%focal_site = &pull_focal_site_from_maf($chr);
				open(EXVCF,">extract_vcf/chr$chr.vcf");
				open(EXANN,">annotation_extract/chr$chr.vcf");
		}
		my($DP_posi,$AD_posi)=(3,5);
		if($line[8] ne "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"){
				my@format=split(/:/,$line[8]);
				for(my $i=0;$i<@format;$i++){
						if($format[$i] eq "DP"){$DP_posi=$i;}
						if($format[$I] eq "AD"){$AD_posi=$i;}
				}
		}
		if($line[4] =~/,/){
				my($a1,$a2)=split(/,/,$line[4]);
				my($posi1,$ref1,$alt1)=&indel_posi_edit($line[1],$line[3],$a1);
				my($posi2,$reif2,$alt2)=&indel_posi_edit($line[1],$line[3],$a2);
				my($allele1,$allele2)=("$line[0]:$posi1:$ref1:$alt1","$line[0]:$posi2:$ref2:$alt2");
				if((!defined $focal{$allele1})&&(!defined $focal{$allele2})){next;}
				my($het1,$hom1,$het2,$hom2)=(0,0,0,0);
				for(my $ind =9;$ind<scalar(@col_info);$ind++){
						if($line[$ind] =~ /^\./){next;}
						if($_=~/^0\/1/){
								$het1++;
								print 

}

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
		open(OUT,">annotation_extract/chr$chr.maf") or die "ERROR::cannot open annotation_extract\n";
		print OUT "Hugo_Symbol\t".join("\t",@out_list)."\n";
		open(MAF,"test/chr$chr"."_snp.maf") or die "ERROR::cannot open chr$chr maf file";
		<MAF>;
		my %col=&header2hush(<MAF>);
		my %maf_focal =();
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
				my $gene = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{gene} = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{IMPACT} = $line[$col{IMPACT}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"}{Consequence} = $line[$col{Consequence}];
				my $out="$gene";
				foreach my $out_list (@out_list){
						$out.="\t$line[$col{$out_list}]";
				}
				if(defined $top_driver_genes{$gene}){
						if($top_driver_genes{$gene} >=4){
								print TDGANN "$out\n";
						}else{
								print ODGANN "$out\n";
						}
				}else{
						print OUT "$out\n";
				}
		}
		close MAF;
		open(MAF,"test/chr$chr"."_indel.maf") or die "ERROR::cannot open chr$chr maf file";
		<MAF>;
		%col=&header2hush(<MAF>);
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
				my $gene = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"} = $line[$col{Hugo_Symbol}];
				$maf_focal{"$line[$col{Chromosome}]:$line[$col{Start_Position}]:$line[$col{Reference_Allele}]:$line[$col{Tumor_Seq_Allele2}]"} = $line[$col{Consequence}];
				my $out="$gene";
				foreach my $out_list (@out_list){
						$out.="\t$line[$col{$out_list}]";
				}
				if(defined $top_driver_genes{$gene}){
						if($top_driver_genes{$gene} >=4){
								print TDGANN "$out\n";
						}else{
								print ODGANN "$out\n";
						}
				}else{
						print OUT "$out\n";
				}
		}
		close OUT;
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



