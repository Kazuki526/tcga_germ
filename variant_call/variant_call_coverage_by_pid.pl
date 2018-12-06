#!/usr/bin/perl
use strict;
use warnings;

########## perl downloading.pl P_ID PURITY ##########

#check reference of bam are exist?
my $ref="/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";

if(!defined $ARGV[1]){die "ERROR::plese imput P_ID:patient_id and PURITY:purity from ASCAT";}
my ($pid,$purity) = @ARGV;
my $bam_json = "$pid.json";

# check token file
my $token_path=`ls $ENV{HOME}/git/innanlab/gdc|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="$ENV{HOME}/git/innanlab/gdc/$token_path";
my $token=`cat $token_path`;

#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

#check gdc-clinent exist??
my $gdc_client = "$ENV{HOME}/gdc-client";
(-e $gdc_client) or die "ERROR::gdc-client does not exits!\n";

#check existing sample_bam.json
my $sample_json="$ENV{HOME}/git/tcga_germ/variant_call/sample_bam.json";
(-e $sample_json) or die "ERROR::$sample_json exist??";

mkdir "norm_bam";
mkdir "tumor_bam";
my $response = "bam_info.tsv";
`cat $sample_json|sed s/PATIENT_ID/$pid/ >$bam_json`;
`curl --request POST --header \"Content-Type: application/json\" --data \@$bam_json 'https://api.gdc.cancer.gov/files'|nkf -Lu >$response`;


# rad response and make norm & tumor manifest
my (@norm_file,@tumor_file)=((),());
my ($norm_manifest, $tumor_manifest) = ("norm_bam/norm_manifest.tsv","tumor_bam/tumor_manifest.tsv");
open(RES,"$response");
my @header = split(/\t/,<RES>);chomp @header;
my %header_num=();
for(my $i=0;@header>$i;$i++){
		if($header[$i] eq "file_name"){$header_num{'file'}=$i;
		}elsif($header[$i] eq "cases.0.submitter_id"){$header_num{'case_id'}=$i;
		}elsif($header[$i] eq "cases.0.samples.0.sample_type"){$header_num{'sample_type'}=$i;
		}elsif($header[$i] eq "id"){$header_num{'file_id'}=$i;
		}else{die "ERROR::$response have wrong header what is $header[$i]??\n";}
}
while(<RES>){
		chomp;
		my @line=split(/\t/,);
		if($line[$header_num{sample_type}] eq "Primary Tumor"){
				`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
				my $partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
				my $t=1;
				while($partial eq "true" && $t<=10){
						`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
						$partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
						$t++;
				}
				`rm -rf $line[$header_num{file_id}]/logs`;
				`mv $line[$header_num{file_id}]/* tumor_bam/`;
				`rm -rf $line[$header_num{file_id}]`;
				if($partial eq "true"){push(@tumor_file,"tumor_bam/$line[$header_num{file}]");}
		}else{
				`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
				my $partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
				my $t=1;
				while($partial eq "true" && $t<=10){
						`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
						$partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
						$t++;
				}
				`rm -rf $line[$header_num{file_id}]/logs`;
				`mv $line[$header_num{file_id}]/* norm_bam/`;
				`rm -rf $line[$header_num{file_id}]`;
				if($partial eq "true"){push(@norm_file,"norm_bam/$line[$header_num{file}]");}
		}
}
close RES;

my %bam_file=();
# merge bam file
if(scalar(@norm_file) >1){
		`samtools merge -f norm_bam/$pid.bam @norm_file`;
		`samtools index norm_bam/$pid.bam`;
		$bam_file{norm}="norm_bam/$pid.bam";
}elsif(scalar(@norm_file)>0){$bam_file{norm}=$norm_file[0];}
if(scalar(@tumor_file) >1){
		`samtools merge -f tumor_bam/$pid.bam @tumor_file`;
		`samtools index tumor_bam/$pid.bam`;
		$bam_file{tumor}="tumor_bam/$pid.bam";
}elsif(scalar(@tumor_file)>0){$bam_file{tumor}=$tumor_file[0];}

#do varscan and coverage 
if(defined$bam_file{norm} && defined$bam_file{tumor}){
		print " $pid variant call\n";
		mkdir "vcf";
		my $mpile = "samtools mpileup -q 10 -f $ref $bam_file{norm} $bam_file{tumor}";
		`zsh -c \"varscan somatic <\($mpile\) vcf/$pid --tumor-purity $purity --p-value 0.1 --output-vcf 1 --mpileup 1\"`;
		open(DP,"samtools depth -Q 10 $bam_file{norm} $bam_file{tumor}|");
		open(OUT,"|gzip -c depth.tsv.gz");
		while(<DP>){
				chomp;
				my @line=split(/\t/,);
				if($line[2] >=6 && $line[3]>=8){print OUT "$_\n";}
		}
		close DP;close OUT;
}
exit;

sub check_download( $ $ ){
		my($file_id,$file)=@_;
		my $partial="false";
		my @ls = `ls $file_id`;chomp @ls;
		# partial check
		foreach my $f (@ls){
				if($f =~/partial$/){$partial="true";}
		}
		if($partial eq "true"){
				return($partial);
		}else{
				# bai check
				my $bai=0;
				foreach my $f (@ls){
						if($f =~/bai$/){$bai++;}
				}
				if($bai==0){`samtools index $file_id/$file`;}
		}
		return($partial);
}
		

