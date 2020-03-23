#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my $tdgdir="/Volumes/areca42TB2/gdc/top_driver_gene";
(-e $tdgdir) or die "ERROR::changed $tdgdir??";
if(!-e "$tdgdir/blood"){mkdir "$tdgdir/blood";}
#check reference of bam are exist?
my $ref="/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";
# check token file
my $token_path=`ls $ENV{HOME}/git/innanlab/gdc|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="$ENV{HOME}/git/innanlab/gdc/$token_path";
my $token=`cat $token_path`;
#check bed.json file
my $json = "$ENV{HOME}/git/driver_genes/onlytop105/top_driver105exon_json.txt";
($json     and -e "$json"    ) or die "ERROR: there arenot json file!!\n";
#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";


my %patient_inf =();
my @patient_list =();
open(PLIST,"$ENV{HOME}/git/tcga_germ/variant_call/focal_patient_list.tsv") or die "ERROR::patient list is correct??\n";
my $header=<PLIST>; chomp $header;
my %col = &header2hash($header);
while(<PLIST>){
		chomp;
		my @line = split(/\t/,);
		push(@patient_list,$line[$col{patient_id}]);
		my($nblood,$nsolid)=($line[$col{bloodnorm}],$line[$col{solidnorm}]);
		if($nblood eq "NA"){$nblood = 0;}
		if($nsolid eq "NA"){$nsolid = 0;}
		$patient_inf{$line[$col{patient_id}]}{cancer_type}=lc $line[$col{cancer_type}];
		$patient_inf{$line[$col{patient_id}]}{gender}=$line[$col{gender}];
		$patient_inf{$line[$col{patient_id}]}{nblood}=$nblood;
		$patient_inf{$line[$col{patient_id}]}{nsolid}=$nsolid;
		$patient_inf{$line[$col{patient_id}]}{nnorm} =$nblood+$nsolid;
		$patient_inf{$line[$col{patient_id}]}{ntumor}=$line[$col{tumor}];
}
close PLIST;
#debug
if(scalar(@patient_list) != scalar(keys%patient_inf)){die "ERROR::patient_list error??". scalar(@patient_list) .":". scalar(keys%patient_inf)."\n";}


open(BLIST,"$ENV{HOME}/git/tcga_germ/variant_call/bam_file_list.tsv") or die "ERROR::bam list is correct??\n";
$header=<BLIST>; chomp $header;
%col = &header2hash($header);
while(<BLIST>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{sample_type}] eq "Primary Tumor"){
				$patient_inf{$line[$col{patient_id}]}{tumor_file}.="$line[$col{file_name}],";
				$patient_inf{$line[$col{patient_id}]}{tumor_id}.="$line[$col{id}],";
				$patient_inf{$line[$col{patient_id}]}{tumor_sample_id}.="$line[$col{sample_id}],";
		}elsif($line[$col{sample_type}] eq "Blood Derived Normal"){
				$patient_inf{$line[$col{patient_id}]}{blood_file}.="$line[$col{file_name}],";
				$patient_inf{$line[$col{patient_id}]}{blood_id}.="$line[$col{id}],";
				$patient_inf{$line[$col{patient_id}]}{norm_file}.="$line[$col{file_name}],";
				$patient_inf{$line[$col{patient_id}]}{norm_id}.="$line[$col{id}],";
				$patient_inf{$line[$col{patient_id}]}{norm_type}.="blood,";
		}elsif($line[$col{sample_type}] eq "Solid Tissue Normal"){
				$patient_inf{$line[$col{patient_id}]}{solid_file}.="$line[$col{file_name}],";
				$patient_inf{$line[$col{patient_id}]}{solid_id}.="$line[$col{id}],";
				$patient_inf{$line[$col{patient_id}]}{norm_file}.="$line[$col{file_name}],";
				$patient_inf{$line[$col{patient_id}]}{norm_id}.="$line[$col{id}],";
				$patient_inf{$line[$col{patient_id}]}{norm_type}.="solid,";
		}
}
close BLIST;

# read purity いちおう、、、
open(PURITY,"/Volumes/areca42TB2/gdc/purity/by_sample/tcga_sample_purity.tsv") or die "ERROR::cannot open purity file\n";
$header=<PURITY>; chomp $header;
%col = &header2hash($header);
while(<PURITY>){
		chomp;
		my @line = split(/\t/,);
		if(!defined $patient_inf{$line[$col{patient_id}]}{tumor_sample_id}){next;}
		if($patient_inf{$line[$col{patient_id}]}{tumor_sample_id} =~ /$line[$col{tumor_sample_id}]/){
				if($line[$col{purity}] ne "NA"){
						$patient_inf{$line[$col{patient_id}]}{purity}.="$line[$col{purity}],";
				}
		}
}
close PURITY;

my %focal_noline=();
open(IN,"/Volumes/areca42TB2/gdc/top_driver_gene/no_line_maf.tsv")or die "ERROR::cannot open no_line_maf.tsv\n";
while(<IN>){
		chomp;my@line=split(/\t/,);
		$focal_noline{$line[1]}="noline";
}
close IN;

my %pj2bp = ("brca" => "breast",
			 "coad" => "colorectal", "read" => "colorectal",
			 "gbm"  => "brain", "lgg" => "brain",
			 "kich" => "kidney","kirc" => "kidney","kirp" => "kidney",
			 "luad" => "lung", "lusc" => "lung");


open(LOG,">$ENV{HOME}/git/tcga_germ/variant_call/done_patient.txt");
my $line=0;
my $max_processes = 10;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $pid (@patient_list){
		$line++;
		if(!defined $focal_noline{$pid}){next;}
		$pm->start and next; #do fork
		my $CT = $patient_inf{$pid}{cancer_type};
		print "$line:start $CT:$pid\n";
		if(defined$pj2bp{$patient_inf{$pid}{cancer_type}}){$CT=$pj2bp{$patient_inf{$pid}{cancer_type}};}
		`rm $tdgdir/blood/$CT/*/$pid*`;
		my $nsolid = $patient_inf{$pid}{nsolid};
		my ($gender,$nnorm,$ntumor) = ($patient_inf{$pid}{gender},$patient_inf{$pid}{nblood},$patient_inf{$pid}{ntumor});
		if($nnorm==0){die "$pid have no blood norm\n";}
		my @norm_file  = split(/,/,$patient_inf{$pid}{blood_file});
		my @norm_id    = split(/,/,$patient_inf{$pid}{blood_id});
		my @tumor_file = split(/,/,$patient_inf{$pid}{tumor_file});
		my @tumor_id   = split(/,/,$patient_inf{$pid}{tumor_id});
		my $purity = 0;
		if(defined $patient_inf{$pid}{purity}){
				my @purity  = split(/,/,$patient_inf{$pid}{purity});
				$purity += $_ for @purity;
				$purity = $purity/scalar(@purity);
		}else{$purity = 1;}
#debug
		if((scalar(@norm_file)==0)||(scalar(@norm_id)==0)||
						(scalar(@tumor_file)==0)||(scalar(@tumor_id)==0)||
						(scalar(@norm_file)!=$nnorm)||(scalar(@norm_id)!=$nnorm)||
						(scalar(@tumor_file)!=$ntumor)||(scalar(@tumor_id)!=$ntumor)){
				die "ERROR::file input error?? at $pid\n";
		}

		if(!-e "$tdgdir/blood/$CT"){mkdir "$tdgdir/blood/$CT";}
		if(!-e "$tdgdir/blood/$CT/maf"){mkdir "$tdgdir/blood/$CT/maf";}
		if(!-e "$tdgdir/blood/$CT/out"){mkdir "$tdgdir/blood/$CT/out";}
		if(!-e "$tdgdir/blood/$CT/vcf"){mkdir "$tdgdir/blood/$CT/vcf";}
		if(!-e "$tdgdir/blood/$CT/depth"){mkdir "$tdgdir/blood/$CT/depth";}
		my $vcf_snp ="$tdgdir/blood/$CT/out/$pid.snp.vcf";
		my $maf ="$tdgdir/blood/$CT/maf/$pid.maf";
		if($nsolid == 0){
				if(!-e $maf){`cp $tdgdir/$CT/maf/$pid.maf $maf`;}
		}else{
#merge multi bam files
				my $norm_bam="";
				if($nnorm>1){
						for(my $i=0;$i<$nnorm;$i++){$norm_file[$i]="$tdgdir/norm_bam/$CT/$norm_file[$i]";}
						$norm_bam = "$tdgdir/norm_bam/$CT/$pid"."_blood.bam";
						if(&bam_check($norm_bam) eq "error"){
								`samtools merge -f $norm_bam @norm_file`;
								`rm $tdgdir/blood/$CT/*/$pid*`;
								if(&bam_check($norm_bam) eq "error"){
										print LOG "$CT\t$pid\t$norm_bam\tmerge_error\n";
										die "$line:$CT-$pid:$norm_bam have merge error what happened??\n";
								}
								`samtools index $norm_bam`;
						}
				}else{
						$norm_bam="$tdgdir/norm_bam/$CT/$norm_file[0]";
						if(!-e "$norm_bam.bai"){`samtools index $norm_bam`;}
				}
				my $tumor_bam="";
				if($ntumor>1){
						$tumor_bam = "$tdgdir/tumor_bam/$CT/$pid.bam";
				}else{
						$tumor_bam="$tdgdir/tumor_bam/$CT/$tumor_file[0]";
				}
#varscan
				if(!-e $vcf_snp){
						my $mpile = "samtools mpileup -q 10 -f $ref $norm_bam $tumor_bam";
						print "$line:start $CT:$pid varscan\n";
						`zsh -c \"varscan somatic <\($mpile\) $tdgdir/blood/$CT/out/$pid --tumor-purity $purity --p-value 0.1 --output-vcf 1 --mpileup 1\"`;
				}
				if(!-e $maf){
						print "$line:start $CT:$pid vcf2maf\n";
						`zsh $ENV{HOME}/git/tcga_germ/variant_call/vcf2maf.sh $pid $tdgdir/blood/$CT 2>&1`;
				}
		}
		print "$line:finish $CT:$pid\n";

		$pm->finish;
}
$pm->wait_all_children;
close LOG;
exit;





sub bam_check( $ ){
		my $bam=$_[0];
		if(!-e $bam){return "error";}
		else{
				my $head = `samtools view $bam 2>&1|head -n 1`;
				if(($head !~ /EOF\smarker\sis\sabsent/) && ($head !~ /Parse\serror/)){
						return "ok";
				}else{return "error";
				}
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
