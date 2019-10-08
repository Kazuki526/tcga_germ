#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my $pwd = `pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc/pancan_atlas_germ"){die "ERROR::done on wrong dir\n";}

my $vcf_file = "raw_file/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz";
-e $vcf_file or die "ERROR::not exist $vcf_file\n";

open(VCF,"gunzip -c $vcf_file|");
my $chr=1;my @chr=();
print "start $chr => ";
open(OUTSNP,">annotation/chr$chr"."_snp_like.vcf");
open(OUTIND,">annotation/chr$chr"."_indel_like.vcf");
while(<VCF>){
		if($_ =~ /^#/){next;}
		chomp;
		my @line = split(/\t/,);
		if(($chr ne $line[0])&&($chr ne "other")){
				close OUTSNP;close OUTIND;
				push(@chr,$chr);
#				`bash script/vcf2maf.sh chr$chr`;
				print "done $chr\n";
				if($line[0] =~/[\dXY]{1,2}/){$chr = $line[0];
				}else{$chr ="other";}
				print "start $chr => ";
				open(OUTSNP,">annotation/chr$chr"."_snp_like.vcf");
				open(OUTIND,">annotation/chr$chr"."_indel_like.vcf");
		}
		if($line[4] =~ /,/){
				my@alt=split(/,/,$line[4]);
				my$an;my@ac;
				if($line[7] =~ /;AN=(\d+);/){$an=$1;}
				if($line[7] =~ /;AC=([\d,]+)$/){
						@ac = split(/,/,$1);
				}else{die "ERROR::this line cannot extract AC value\n".join("\t",@line[0..7])."\n";}
				for(my $i=0;$i<scalar(@alt);$i++){
						if(length($line[3]) != length($alt[$i])){
								print OUTIND join("\t",@line[0..3])."\t$alt[$i]\t$line[5]\t$line[6]\t$an\t$ac[$i]\n";
						}else{
								print OUTSNP join("\t",@line[0..3])."\t$alt[$i]\t$line[5]\t$line[6]\t$an\t$ac[$i]\n";
						}
				}
		}else{
				my($an,$ac);
				if($line[7] =~ /;AN=(\d+);AC=(\d+)$/){($an,$ac)=($1,$2);
				}else{die "ERROR::this line cannot extract AN&AC value\n".join("\t",@line[0..7])."\n";}
				if(length($line[3]) != length($line[4])){
						print OUTIND join("\t",@line[0..6])."\t$an\t$ac\n";
				}else{
						print OUTSNP join("\t",@line[0..6])."\t$an\t$ac\n";
				}
		}
}
close OUTSNP;close OUTIND;

my $max_processes = 3;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $outchr (@chr){
		$pm->start and next;
		`bash $ENV{HOME}/tcga_germ/PCA_script/git//vcf2maf.sh chr$chr`;
		$pm->finish;
}
$pm->wait_all_children;
exit;
		

