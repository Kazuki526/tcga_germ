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
my @ann_list = qw(Chromosome Start_Position End_Position Strand Reference_Allele Tumor_Seq_Allele2 Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons CANONICAL SIFT PolyPhen Feature);
open(TDGANN,">annotation_extract/top_driver_genes.maf");
open(ODGANN,">annotation_extract/other_driver_genes.maf");
print TDGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
print ODGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";

open(TDGVCF,"|gzip -c >extract_vcf/top_driver_genes.vcf");
open(ODGVCF,"|gzip -c >extract_vcf/other_driver_genes.vcf");
my $vcf_list = "#chr\tposi\tid\tref\talt\tqual\tfilter\tinfo\n";
print TDGVCF "$vcf_list";
print ODGVCF "$vcf_list";

open(TDGOUT,"|gzip -c >maf_patient/top_driver_genes_patient.tsv");
open(ODGOUT,"|gzip -c >maf_patient/other_driver_genes_patient.tsv");
my $out_list = "patient_id\tchr\tposi\tref\tallele1\tallele2\tfilter\tDP\tAD\tgene\tConsequence\tIMPACT\n";
print TDGOUT "$out_lsit";
print ODGOUT "$out_lsit";

my $chr="";
my %focal_site=();
my @col_info=();
while(<VCF>){
		if($_ =~ /^##/){next;}
		chomp;
		if($_ =~ /^#/){$_ =~ s/#//;@col_info=split(/\t/,);next;}
		if($_ =~ /^$chr\t/){
				if($_=~/^$chr\t(\d+)\t\S+\t(\w+)\t(\w+)\t/){
						my($posi,$ref,$alt)=($1,$2,$3);
						if(length($ref) != length($alt)){($posi,$ref,$alt)=&indel_posi_edit($posi,$ref,$alt);}
						if(!defined $focal_site{"$chr:$posi:$ref:$alt"}){next;}
				}
		}
		my @line = split(/\t/,);
		if($line[0] ne $chr){
				if($chr ne ""){
						print "done\n";$|=1;
						close EXVCF;
						close EXOUT;
				}
				$chr =$line[0];
				print "start chr$chr =>";$|=1;
				%focal_site = &pull_focal_site_from_maf($chr);
				open(EXVCF,">extract_vcf/chr$chr.vcf");
				print EXVCF "$vcf_list";
				open(EXOUT,">maf_patient/chr$chr.tsv");
				print EXOUT "$out_list";
				print "finish read MAF file =>";$|=1;
		}
		my($DP_posi,$AD_posi)=(3,5);
		if($line[8] ne "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"){
				if($line[8] !~ /^GT:/){die "ERROR::this line FORMAT is have no genotype ".join(" ",@line[0..8])."\n";}
				my@format=split(/:/,$line[8]);
				for(my $i=0;$i<@format;$i++){
						if($format[$i] eq "DP"){$DP_posi=$i;}
						if($format[$I] eq "AD"){$AD_posi=$i;}
				}
				if((!defined $format[$DP_posi])||($format[$DP_posi] ne "DP")){$DP_posi = 100;}
				if((!defined $format[$AD_posi])||($format[$AD_posi] ne "AD")){$AD_posi = 100;}
		}
		if($line[4] =~/,/){
				my @posi=($line[1]);
				my @ref=($line[3]);
				my @alt=($line[4]);
				for

				my($a1,$a2)=split(/,/,$line[4]);
				my($posi1,$ref1,$alt1)=&indel_posi_edit($line[1],$line[3],$a1);
				my($posi2,$reif2,$alt2)=&indel_posi_edit($line[1],$line[3],$a2);
				my($allele1,$allele2)=("$line[0]:$posi1:$ref1:$alt1","$line[0]:$posi2:$ref2:$alt2");
				if((!defined $focal_site{$allele1})&&(!defined $focal_site{$allele2})){next;}
				my($het1,$hom1,$het2,$hom2)=(0,0,0,0);
				for(my $ind =9;$ind<scalar(@col_info);$ind++){
						if($line[$ind] =~ /^\./){next;}
						my @format = split(/:/,$line[$ind])
						if($format[0] eq "^0/1"){
								if(!defined $focal_site{$allele1}){next;}
								$het1++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi1\t$ref1\t$ref1\t$alt1\t$line[6]\t$format[$DP_posi]\t$format[$AD_posi]\t";
								$out.= "$focal_site{$allele1}{gene}\t$focal_site{$allele1}{Consequence}\t$focal_site{$allele1}{IMPACT}\n";
								&print_out($out,$focal_site{$allele1}{gene})
						}elsif($format[0] eq "^0/2"){
								if(!defined $focal_site{$allele2}){next;}
								$het2++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi2\t$ref2\t$ref2\t$alt2\t$line[6]\t$format[$DP_posi]\t$format[$AD_posi]\t";
								$out.= "$focal_site{$allele2}{gene}\t$focal_site{$allele2}{Consequence}\t$focal_site{$allele2}{IMPACT}\n";
								&print_out($out,$focal_site{$allele2}{gene})
						}elsif($format[0] eq "^1/1"){
								if(!defined $focal_site{$allele1}){next;}
								$hom1++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi1\t$ref1\t$alt1\t$alt1\t$line[6]\t$format[$DP_posi]\t$format[$AD_posi]\t";
								$out.= "$focal_site{$allele1}{gene}\t$focal_site{$allele1}{Consequence}\t$focal_site{$allele1}{IMPACT}\n";
								&print_out($out,$focal_site{$allele1}{gene})
						}elsif($format[0] eq "^2/2"){
								if(!defined $focal_site{$allele2}){next;}
								$het2++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi2\t$ref2\t$alt2\t$alt2\t$line[6]\t$format[$DP_posi]\t$format[$AD_posi]\t";
								$out.= "$focal_site{$allele2}{gene}\t$focal_site{$allele2}{Consequence}\t$focal_site{$allele2}{IMPACT}\n";
								&print_out($out,$focal_site{$allele2}{gene})
						}elsif($format[0] eq "1/2"){
						}
				}
		}else{
				my($posi,$ref,$alt) = ($line[1],$line[3],$line[4]);
				if(length($ref) != length($alt)){($posi,$ref,$alt)=&indel_posi_edit($posi,$ref,$alt);}
				my $allele="$chr:$posi:$ref:$alt";
				if(!defined $focal_site{$allele}){next;}
				my($vcf,$het,$hom)=("",0,0);
				for(my$ind=9;$ind<scalar(@col_info);$ind++){
						if(!defined$line[$ind]){die "ERROR::line legth error\n";}
						if($line[$ind] =~ /^\.\/\./){next;}
						my @format = split(/:/,$line[$ind]);
						my $out_format = "";
						if($DP_posi == 100){$out_format="NA\t";}else{$out_format="$format[$DP_posi]\t";}
						if($AD_posi == 100){$out_format.="NA\t";}else{$out_format.="$format[$AD_posi]\t";}
						if($format[0] eq "0/1"){
								$het++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi\t$ref\t$ref\t$alt\t$line[6]\t$out_format";
								$out.= "$focal_site{$allele}{gene}\t$focal_site{$allele}{Consequence}\t$focal_site{$allele}{IMPACT}\n";
								&print_out($out,$focal_site{$allele});
						}elsif($format[0] eq "1/1"){
								$hom++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi\t$ref\t$alt\t$alt\t$line[6]\t$out_format";
								$out.= "$focal_site{$allele}{gene}\t$focal_site{$allele}{Consequence}\t$focal_site{$allele}{IMPACT}\n";
								&print_out($out,$focal_site{$allele});
						}else{print  "\nWARNING::what genotype or allele: $allele => $ind th FORMAT: $line[$ind]\n";}
				}
				$vcf = join("\t",@line[0..7]) .";het_sample=$het;hom_sample=$hom\n";
				&print_vcf($vcf,$focal_site{$allele}{gene});
		}
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
		open(EXANN,">annotation_extract/chr$chr.maf") or die "ERROR::cannot open annotation_extract\n";
		print EXANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
		open(MAF,"test/chr$chr"."_snp.maf") or die "ERROR::cannot open chr$chr maf file";
		<MAF>;
		my $header = <MAF>;chomp$header;
		my %col=&header2hash($header);
		my %maf_focal =();
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
				my $gene = $line[$col{Hugo_Symbol}];
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
		%col=&header2hash($header);
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if(($line[$col{IMPACT}] eq "MODIFIER")||($line[$col{BIOTYPE}] ne "protein_coding")){next;}
				my $gene = $line[$col{Hugo_Symbol}];
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


