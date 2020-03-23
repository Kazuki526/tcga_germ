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
my $cont_file = "/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv";
open(CG,"$cont_file") or die "ERROR::cannot open $cont_file\n";
<CG>;
while(<CG>){
		chomp;
		my @line = split(/\t/,);
		$top_driver_genes{$line[0]}=0;
}
close CG;

my $vcf_file = "raw_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz";
-e $vcf_file or die "ERROR::not exist $vcf_file\n";


open(VCF,"gunzip -c $vcf_file|");
my @ann_list = qw(Chromosome Start_Position End_Position Strand Reference_Allele Tumor_Seq_Allele2 Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons CANONICAL SIFT PolyPhen Feature);
open(TDGANN,">annotation_extract/top_driver_genes.maf");
open(ODGANN,">annotation_extract/other_driver_genes.maf");
open(CONANN,">annotation_extract/control_genes.maf");
print TDGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
print ODGANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";
print CONANN "Hugo_Symbol\t".join("\t",@ann_list)."\n";

open(TDGVCF,"|gzip -c >extract_vcf/top_driver_genes.vcf.gz");
open(ODGVCF,"|gzip -c >extract_vcf/other_driver_genes.vcf.gz");
open(CONVCF,"|gzip -c >extract_vcf/control_genes.vcf.gz");
my $vcf_list = "#chr\tposi\tid\tref\talt\tqual\tfilter\tinfo\n";
print TDGVCF "$vcf_list";
print ODGVCF "$vcf_list";
print CONVCF "$vcf_list";

open(TDGOUT,"|gzip -c >maf_patient/top_driver_genes_patient.tsv.gz");
open(ODGOUT,"|gzip -c >maf_patient/other_driver_genes_patient.tsv.gz");
open(CONOUT,"|gzip -c >maf_patient/control_genes_patient.tsv.gz");
my $out_list = "patient_id\tchr\tposi\tref\tallele1\tallele2\tfilter\tDP\tAD\tgene\tConsequence\tIMPACT\n";
print TDGOUT "$out_list";
print ODGOUT "$out_list";
print CONOUT "$out_list";

my $chr="";
my %focal_site=();
my @col_info=();
while(<VCF>){
		if($_ =~ /^##/){next;}
		if($_ =~ /^#CHROM/){$_ =~ s/#//;chomp;@col_info=split(/\t/,);next;}
		if($_ =~ /^$chr\t/){
				if($_=~/^$chr\t(\d+)\t\S+\t(\w+)\t(\w+)\t/){
						my($posi,$ref,$alt)=($1,$2,$3);
						if(length($ref) != length($alt)){($posi,$ref,$alt)=&indel_posi_edit($posi,$ref,$alt);}
						if(!defined $focal_site{"$chr:$posi:$ref:$alt"}){next;}
				}
		}
		chomp;
		my @line = split(/\t/,);
		if($line[0] ne $chr){
				if($chr ne ""){
						print "done chr$chr\n";$|=1;
						close EXVCF;
						close EXOUT;
				}
				$chr =$line[0];
				print "start chr$chr => ";$|=1;
				%focal_site = &pull_focal_site_from_maf($chr);
				open(EXVCF,"|gzip -c >extract_vcf/chr$chr.vcf.gz");
				print EXVCF "$vcf_list";
				open(EXOUT,"|gzip -c >maf_patient/chr$chr.tsv.gz");
				print EXOUT "$out_list";
				print "finish read MAF file =>";$|=1;
		}
		my($DP_posi,$AD_posi)=(3,5);
		if($line[8] ne "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR"){
				if($line[8] !~ /^GT:/){die "ERROR::this line FORMAT is have no genotype ".join(" ",@line[0..8])."\n";}
				my@format=split(/:/,$line[8]);
				for(my $i=0;$i<@format;$i++){
						if($format[$i] eq "DP"){$DP_posi=$i;}
						if($format[$i] eq "AD"){$AD_posi=$i;}
				}
				if((!defined $format[$DP_posi])||($format[$DP_posi] ne "DP")){$DP_posi = 100;}
				if((!defined $format[$AD_posi])||($format[$AD_posi] ne "AD")){$AD_posi = 100;}
		}
		if($line[4] =~/,/){
				my $focal_exist=0;
				my @focal =(0);
				my @posi=($line[1]);
				my @ref=($line[3]);
				my @alt=($line[3]);
				my @variant=split(/,/,$line[4]);
				my @het=(0);my @hom=(0);
				for(my$i=1;$i<=scalar(@variant);$i++){
						($posi[$i],$ref[$i],$alt[$i]) = &indel_posi_edit($line[1],$line[3],$variant[$i-1]);
						if(defined$focal_site{"$line[0]:$posi[$i]:$ref[$i]:$alt[$i]"}){$focal_exist++;$focal[$i]=1;}else{$focal[$i]=0;}
						$het[$i]=0;$hom[$i]=0;
				}
				if($focal_exist==0){next;}
				for(my $ind =9;$ind<scalar(@col_info);$ind++){
						if($line[$ind] =~ /^\./){next;}
						my @format = split(/:/,$line[$ind]);
						if($format[0] !~ /^\d+\/\d+$/){die "ERROR::FORMAT error $line[0] $line[1] patient $ind:$col_info[$ind] FORMAT:$line[8] => $line[$ind]\n";}
						my @allele=split(/\//,$format[0]);
						my $out_format = "";
						if($DP_posi == 100){$out_format="NA\t";}else{$out_format="$format[$DP_posi]\t";}
						if($AD_posi == 100){$out_format.="NA\t";}else{$out_format.="$format[$AD_posi]\t";}
						if(($allele[0] != $allele[1])&&($allele[0] !=0)){
								my($a0,$a1)=@allele;
								if(($focal[$a0]==0)&&($focal[$a1]==0)){next;}
								my($allele0,$allele1)=("$line[0]:$posi[$a0]:$ref[$a0]:$alt[$a0]","$line[0]:$posi[$a1]:$ref[$a1]:$alt[$a1]");
								$het[$a0]++;$het[$a1]++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi[$a0]\t$ref[$a0]\t$alt[$a0]\tNA\t$line[6]\t$out_format\t";
								$out.= "$focal_site{$allele0}{gene}\t$focal_site{$allele0}{Consequence}\t$focal_site{$allele0}{IMPACT}\n";
								&print_out($out,$focal_site{$allele0}{gene});
								$out = "$col_info[$ind]\t$line[0]\t$posi[$a1]\t$ref[$a1]\tNA\t$alt[$a1]\t$line[6]\t$out_format\t";
								$out.= "$focal_site{$allele1}{gene}\t$focal_site{$allele1}{Consequence}\t$focal_site{$allele1}{IMPACT}\n";
								&print_out($out,$focal_site{$allele1}{gene});
						}else{
								my $a=$allele[1];
								if($focal[$a]==0){next;}
								my $allele = "$line[0]:$posi[$a]:$ref[$a]:$alt[$a]";
								if($allele[0]==0){
										$het[$a]++;
										my $out = "$col_info[$ind]\t$line[0]\t$posi[$a]\t$ref[$a]\t$ref[$a]\t$alt[$a]\t$line[6]\t$out_format\t";
										$out.= "$focal_site{$allele}{gene}\t$focal_site{$allele}{Consequence}\t$focal_site{$allele}{IMPACT}\n";
										&print_out($out,$focal_site{$allele}{gene});
								}else{
										$hom[$a]++;
										my $out = "$col_info[$ind]\t$line[0]\t$posi[$a]\t$ref[$a]\t$alt[$a]\t$alt[$a]\t$line[6]\t$out_format\t";
										$out.= "$focal_site{$allele}{gene}\t$focal_site{$allele}{Consequence}\t$focal_site{$allele}{IMPACT}\n";
										&print_out($out,$focal_site{$allele}{gene});
								}
						}
				}
				for(my$i=1;$i<=scalar(@variant);$i++){
						if($focal[$i] ==0){next;}
						my $vcf = "$line[0]\t$posi[$i]\t$line[2]\t$ref[$i]\t$alt[$i]\t$line[5]\t$line[6]\t$line[7];het_sample=$het[$i];hom_sample=$hom[$i]\n";
						&print_vcf($vcf,$focal_site{"$line[0]:$posi[$i]:$ref[$i]:$alt[$i]"}{gene});
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
								&print_out($out,$focal_site{$allele}{gene});
						}elsif($format[0] eq "1/1"){
								$hom++;
								my $out = "$col_info[$ind]\t$line[0]\t$posi\t$ref\t$alt\t$alt\t$line[6]\t$out_format";
								$out.= "$focal_site{$allele}{gene}\t$focal_site{$allele}{Consequence}\t$focal_site{$allele}{IMPACT}\n";
								&print_out($out,$focal_site{$allele}{gene});
						}#else{print  "\nWARNING::what genotype or allele: $allele => $ind th FORMAT: $line[$ind]\n";}
				}
				$vcf = join("\t",@line[0..7]) .";het_sample=$het;hom_sample=$hom\n";
				&print_vcf($vcf,$focal_site{$allele}{gene});
		}
}
print "done chr$chr\n";$|=1;
close VCF;
close TDGANN;
close ODGANN;
close TDGVCF;
close ODGVCF;
close TDGOUT;
close ODGOUT;
close EXVCF;
close EXOUT;


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
		open(MAF,"annotation/chr$chr"."_snp.maf") or die "ERROR::cannot open chr$chr maf file";
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
		open(MAF,"annotation/chr$chr"."_indel.maf") or die "ERROR::cannot open chr$chr maf file";
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
						}elsif($top_driver_genes{$gene}==0){
								print CONANN "$ann\n";
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
		if(length($ref) == length($var)){
				return($pos,$ref,$var);
		}else{
				my ( $ref_length, $var_length ) = ( length( $ref ), length( $var ));
				while( $ref and $var and substr( $ref, 0, 1 ) eq substr( $var, 0, 1 ) and $ref ne $var ) {
						( $ref, $var ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var );
						--$ref_length; --$var_length; ++$pos;
				}
				$pos =($ref eq "-" ? $pos-1 : $pos);
				return($pos,$ref,$var);
		}
}

sub print_out ( $ $ ){
		my($out,$gene)=@_;
		if(defined $top_driver_genes{$gene}){
				if($top_driver_genes{$gene} >= 4){
						print TDGOUT "$out";
				}elsif($top_driver_genes{$gene}==0){
						print CONOUT "$out";
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
						print TDGVCF "$vcf";
				}elsif($top_driver_genes{$gene}==0){
						print CONVCF "$vcf";
				}else{
						print ODGVCF "$vcf";
				}
		}else{
				print EXVCF "$vcf";
		}
}


