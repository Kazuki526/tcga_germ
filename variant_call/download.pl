#!/usr/bin/perl
use strict;
use warnings;

########## perl downloading.pl P_ID ##########

if(!defined $ARGV[0]){die "ERROR::plese imput P_ID:patient_id";}
my $pid = uc $ARGV[0];
my $bam_json = "$pid.json";

# check token file
my $token_path=`ls ~/git/innanlab/gdc|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="~/git/innanlab/gdc/$token_path";
my $token=`cat $token_path`;

#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";

#check gdc-clinent exist??
my $gdc_client = "~/gdc-client";
(-e $gdc_client) or die "ERROR::gdc-client does not exits!\n";

#check existing sample_bam.json
my $sample_json="$ENV{HOME}/git/tcga/variant_call/sample_bam.json";
(-e $sample_json) or die "ERROR::$sample_json exist??";

mkdir "norm_bam";
mkdir "tumor_bam";
my $response = "bam_info.tsv"
`cat $sample_json|sed s/\\\"PATIENT_ID\\\"/$pid/ >$bam_json`;
`curl --request POST --header \"Content-Type: application/json\" --data \@$bam_json 'https://api.gdc.cancer.gov/files'|nkf -Lu >$response`;


# rad response and make norm & tumor manifest
my (@norm_file,@tumor_file)=((),());
my (@norm_id,@tumor_id)=((),());
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
				while($partial == "true"){
						`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
						$partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
				}
				`rm -rf $line[$header_num{file_id}]/logs`;
				`mv $line[$header_num{file_id}]/* tumor_bam/`;
				push(@tumor_file,$line[$header_num{file}]);
				push(@tumor_id,$line[$header_num{file_id}]);
		}else{
				`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
				my $partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
				while($partial == "true"){
						`$gdc_client download -t $token_path $line[$header_num{file_id}]`;
						$partial =&check_download($line[$header_num{file_id}],$line[$header_num{file}]);
				}
				`rm -rf $line[$header_num{file_id}]/logs`;
				`mv $line[$header_num{file_id}]/* norm_bam/`;
				push(@norm_file,$line[$header_num{file}]);
				push(@norm_id,$line[$header_num{file_id}]);
		}
}
close RES;
close OUTN;
close OUTT;

exit;

sub check_download( $ $ ){
		my($file_id,$file)=@_;
		my $partial="false";
		my @ls = `ls $file_id`;chomp @ls;
		# partial check
		foreach my $f (@ls){
				if($f =~/partial$/){$partial="true";}
		}
		if($partial=="true"){
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
		

