#!/usr/bin/perl
use warnings;
use strict;

my $exon_interval =3;
#my $pwd=`pwd`;chomp $pwd;
#print "doing on $pwd\n";
#if($pwd ne "/Volumes/areca42TB/1000genomes/GRCh38"){die "ERROR:doing on wrong directory\n";}
my $pwd = "/Volumes/areca42TB2/1000genomes/GRCh38";
(-e $pwd) or die "$pwd error\n";

#top_driver_genes = ~/git/driver_genes/onlytop105/top_driver105exon.bed
#control_genes = /Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_exon.bed
my $focus_bed = $ARGV[0];
(-e $focus_bed) or die "ERROR::inputed bed is not exist!\n";
my $focus_name = "$ARGV[1]";
if(($focus_name eq "") && ($focus_bed =~ /\/(.+)\.bed$/)){$focus_name=$1;}

my %bed=();
my @chr=();
open(BED,"$focus_bed");
while(<BED>){
		chomp;
		my @line=split(/\t/,);
		$line[0] =~ s/^chr//;
		if(defined $bed{$line[0]}){
				if($bed{$line[0]} =~ /-(\d+)$/){if($1 > $line[1]){die "please sort $focus_bed\n";}}
				$bed{$line[0]}.=":$line[1]-$line[2]";
		}else{
				$bed{$line[0]}="$line[1]-$line[2]";
		}
		if(!grep($_ eq $line[0],@chr)){push(@chr,"$line[0]");}
}
close BED;

mkdir "$focus_name";
open(OUT,">$focus_name/ALL_$focus_name.vcf");
my $file_num=0;
foreach my $chr(@chr){
		if(!defined $bed{$chr}){next;}
		my @focus_region=split(/:/,$bed{$chr});
		my ($focus_st,$focus_end,$focus_region_num) =(0,0,0);
		($focus_st,$focus_end) = split(/-/,$focus_region[$focus_region_num]);
		my $file = "$pwd/site_vcf/ALL.chr$chr" . "_GRCh38_sites.20170504.vcf.gz";
		$file_num++;
		open(VCF,"gunzip -c $file|")or die "cant open $file\n";
		print "read and count $file\n";
		while(<VCF>){
				if($_=~/^##/){
						if($file_num==1){print OUT "$_";}
						next;
				}
				chomp;
				my @line=split(/\t/,);
				if($line[0] eq "#CHROM"){
						if($file_num==1){
								print OUT "$_\n";
						}
						next;
				}
				if((scalar(@focus_region)-1 == $focus_region_num)&&($line[1] > $focus_end)){last;}
				while($line[1] > $focus_end){
						$focus_region_num++;
						($focus_st,$focus_end) = split(/-/,$focus_region[$focus_region_num]);
				}
				if(!(($focus_st - $exon_interval <= $line[1]+length($line[3])-1)&&
						($line[1] <= $focus_end + $exon_interval))){next;}
# have multi alt
				if($line[4] =~ /,/){
						my @alt = split(/,/,$line[4]);
						for(my $i=0;$i<scalar(@alt);$i++){
								my($posi,$ref,$alt)=&posi_correct($line[1],$line[3],$alt[$i]);
								my @info=();
								foreach my $info (split(/;/,$line[7])){
										if($info =~ /^(.+=)(.+,.+)$/){
												my@value=split(/,/,$2);
												push(@info,"$1$value[$i]");
										}else{
												push(@info,"$info");
										}
								}
								print OUT "$line[0]\t$posi\t$line[2]\t$ref\t$alt\t$line[5]\t$line[6]\t". join(";",@info) ."\n";
						}
# hve only one alt
				}else{
						if((length($line[3])!=1) && (length($line[4])!=1)){print "WARNNING::need position correct??\n$_\n";}
						print OUT "$_\n";
				}
		}
		close VCF;
}

close OUT;


exit;



sub posi_correct($ $ $){
		my($posi,$ref,$alt)=@_;
		while(length($ref)>1 and length($alt)>1 and substr($ref,0,1) eq substr($alt,0,1) and $ref ne $alt){
				$ref=substr($ref,1);$alt=substr($alt,1);$posi++;
		}
		while(length($ref)>1 and length($alt)>1 and substr($ref,-1,1) eq substr($alt,-1,1) and $ref ne $alt){
				$ref=substr($ref,0,-1);$alt=substr($alt,0,-1);
		}
		if((length($ref)!=1) && (length($alt)!=1)){print "WARNNING::position correct error??\n$posi:$ref-$alt\n$_\n";}
		return($posi,$ref,$alt);
}

