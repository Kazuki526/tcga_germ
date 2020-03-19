#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my $tdgdir="/Volumes/areca42TB2/gdc/top_driver_gene";
(-e $tdgdir) or die "ERROR::changed $tdgdir??";
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
#check top driver genes list file
my $tdg_file = "$ENV{HOME}/git/driver_genes/onlytop105/driver_genes.tsv";
-e $tdg_file or die "ERRRO::not exist $tdg_file\n";
#check top driver genes bed file
my $bed_file = "$ENV{HOME}/git/driver_genes/onlytop105/top_driver105cds.bed";
-e $bed_file or die "ERROR:: not exist $bed_file\n";
#chech gnomAD maf fiel
my $gnomAD_file = "/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf";
-e $gnomAD_file or die "ERROR::noet exist $gnomAD_file\n";


my %patient_inf =();
my @patient_list =();
open(PLIST,"$ENV{HOME}/git/tcga_germ/variant_call/focal_patient_list.tsv") or die "ERROR::patient list is correct??\n";
my $header=<PLIST>; chomp $header;
my %col = &header2hash($header);
while(<PLIST>){
		chomp;
		my @line = split(/\t/,);
		push(@patient_list,$line[$col{patient_id}]);
		my $cancer_type="";
		if($line[$col{cancer_type}] =~ /^TCGA-(\w+)$/){$cancer_type = lc $1;}else{die "ERROR::what cancer type $line[$col{cancer_type}]\n";}
		my($nblood,$nsolid)=($line[$col{bloodnorm}],$line[$col{solidnorm}]);
		if($nblood eq "NA"){$nblood = 0;}
		if($nsolid eq "NA"){$nsolid = 0;}
		$patient_inf{$line[$col{patient_id}]}{cancer_type}=$cancer_type;
		$patient_inf{$line[$col{patient_id}]}{race}=$line[$col{race}];
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



my %pj2bp = ("brca" => "breast",
			 "coad" => "colorectal", "read" => "colorectal",
			 "gbm"  => "brain", "lgg" => "brain",
			 "kich" => "kidney","kirc" => "kidney","kirp" => "kidney",
			 "luad" => "lung", "lusc" => "lung");


my %bam_list_fordepth=();
my %focal_patient=();
my $line=0;
foreach my $pid (@patient_list){
		$line++;
		my $CT = $patient_inf{$pid}{cancer_type};
		if(defined$pj2bp{$patient_inf{$pid}{cancer_type}}){$CT=$pj2bp{$patient_inf{$pid}{cancer_type}};}
		my $nsolid = $patient_inf{$pid}{nsolid};
		my $nnorm = $patient_inf{$pid}{nblood};
		if($nnorm==0){next;}
		my @norm_file  = split(/,/,$patient_inf{$pid}{blood_file});
		if(!-e "$tdgdir/blood/$CT"){mkdir "$tdgdir/blood/$CT";}
		if(!-e "$tdgdir/blood/$CT/depth"){mkdir "$tdgdir/blood/$CT/depth";}
		my $norm_bam="";
		if($nnorm==0){next;}
		$focal_patient{$CT}{$pid}="yes";
		if($nsolid==0){next;}
		if($nnorm>1){
				$norm_bam = "$tdgdir/norm_bam/$CT/$pid"."_blood.bam";
		}else{
				$norm_bam="$tdgdir/norm_bam/$CT/$norm_file[0]";
		}
		$bam_list_fordepth{$CT}{n}++;
		$bam_list_fordepth{$CT}{pid}.="$pid\t";
		$bam_list_fordepth{$CT}{norm_bam}.="$norm_bam\t";
}

# make depth file of top_driver_genes' CDS region by 50sample
$line=0;
my $max_processes = 3;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $CT (sort keys %bam_list_fordepth){
		$line++;
		$pm->start and next; #do fork
		print "$line:make depth start $CT";
		my $depth_dir = "$tdgdir/blood/$CT/depth";
		my $num = $bam_list_fordepth{$CT}{n};print " $num\n";
		my @pid = split(/\t/,$bam_list_fordepth{$CT}{pid});
		my @norm_bam = split(/\t/,$bam_list_fordepth{$CT}{norm_bam});
		if(scalar(@norm_bam) != $num){die "ERROR::$CT focal bam num and file count is different\n";}
		my $error="";
		foreach my $bam (@norm_bam){
				my $bam_check = &bam_check($bam);
				if($bam_check ne "ok"){$error=$bam;last;}
		}
		if($error=~/./){
				die "$CT:$error have error\nstop $CT make depth file\n";
		}
		if(-e "$depth_dir/out.tsv"){die "already made dpeht $CT\n";;}
		open(ND,">$depth_dir/out.tsv");
		print ND "chr\tposition\t".join("\t",@pid)."\n";
		close ND;
		`samtools depth -a -q 10 -b $bed_file @norm_bam >>$depth_dir/out.tsv`;
		print "$line:finish make depth $CT\n";
		$pm->finish;
}
$pm->wait_all_children;

#read top driver_genes role
my %tdg_role =();
open(TDG,"$tdg_file");
<TDG>;
while(<TDG>){
		chomp;
		my @line = split(/\t/,);
		if($line[2] eq "NA"){$line[2] ="TSG";}
		$tdg_role{$line[0]}=$line[2];
}
close TDG;
#read top driver genes bed file
my %posi2role=();
open(TDGBED,"$bed_file");
while(<TDGBED>){
		chomp;
		my @line = split(/\t/,);
		my $gene="";
		if($line[3] =~ /^([^:]+):ENST/){$gene=$1;}
		for(my $posi=$line[1]+1;$posi <= $line[2];$posi++){
				$posi2role{$line[0]}{$posi}=$tdg_role{$gene};
		}
}
close TDGBED;

#read gnomAD refminor position
my %gnomad_refminor =();
open(GAD,"$gnomAD_file");
$header=<GAD>; chomp $header;
%col = &header2hash($header);
while(<GAD>){
		chomp;
		my @line = split(/\t/,);
		my $af = $line[$col{AC}]/$line[$col{AN}];
		my $af_white = 0;
		if($line[$col{AN_nfe}]+$line[$col{AN_fin}]+$line[$col{AN_asj}]>0){
				$af_white =  ($line[$col{AC_nfe}]+$line[$col{AC_fin}]+$line[$col{AC_asj}]) / ($line[$col{AN_nfe}]+$line[$col{AN_fin}]+$line[$col{AN_asj}]);
		}
		my $af_black =0;
		if($line[$col{AN_afr}]>0){
				 $af_black = $line[$col{AC_afr}] / $line[$col{AN_afr}];
		}
		if(($af >0.5)||($af_white >0.5)||($af_black >0.5)){
				$gnomad_refminor{$line[0]}{$line[1]}="refminor";
		}
}
close GAD;



$line=0;
#make coverage file
if(!-e "$tdgdir/blood/all_patient"){mkdir "$tdgdir/blood/all_patient";}
open(ARM,"|gzip -c >$tdgdir/blood/all_patient/ref_minor_coverage_blood.tsv.gz");
print ARM "patient_id\tchr\tstart\tfocal\n";
foreach my $CT (sort keys %focal_patient){
		$line++;
		my $depth_dir = "$tdgdir/blood/$CT/depth";
		print "$line:coverage start $CT\n";
		my $ref_minor_file = "$depth_dir/ref_minor_coverage_blood.tsv";
		my $num = 0;
		if(defined $bam_list_fordepth{$CT}{n}){$num=$bam_list_fordepth{$CT}{n};}
#if already done or on downloading
		if(-e $ref_minor_file){
				open(RM,"$ref_minor_file");
				<RM>;
				while(<RM>){print ARM $_;}
				close RM;
				print  "$line:coverage already done $CT\n";
				next;
		}elsif((!-e "$depth_dir/out.tsv")&&($num>0)){
				print  "$line:still downloading $CT\n";
				next;
		}

#calculate coverage and depth
		open(RM,">$ref_minor_file");
		print RM "patient_id\tchr\tstart\tfocal\n";
		if($num == 0){
				open(INRM,"$tdgdir/$CT/depth/ref_minor_coverage.tsv")or die "ERROR::cannot open $CT/depth/ref_minor_coverage.tsv";
				<INRM>;
				while(<INRM>){
						chomp;
						my @line = split(/\t/,);
						if(defined $focal_patient{$CT}{$line[0]}){
								print RM  "$_\n";
								print ARM "$_\n";
						}
				}
				close INRM;
		}else{
				my %norm_depth=();
				open(NDP,"$depth_dir/out.tsv");
				my $nheader=<NDP>;
				my @dp_col = split(/\t/,$nheader);
				while(my $line = <NDP>){
						my @line = split(/\t/,$line);
						if(!defined $gnomad_refminor{$line[0]}{$line[1]}){next;}
						for(my $i =2;$i<scalar(@line);$i++){
								if($line[$i] >=8){
										$norm_depth{$dp_col[$i]}{$line[0]}{$line[1]}="ok";
								}else{
										$norm_depth{$dp_col[$i]}{$line[0]}{$line[1]}="no";
								}
						}
				}
				close NDP;
				open(DPLIST,"$tdgdir/$CT/depth/depth_file_list.tsv")or die "ERROR::$CT/depth/depth_file_list.tsv is not exist\n";
				my $header = <DPLIST>;chomp $header;
				if($header ne "group\tcolmn_num\tpatient_id\tnorm_bam\ttumor_bam"){die "ERROR::colum changed in $CT depth_file_list.tsv??\n";}
				my %tumor_depth_posi=();
				while(<DPLIST>){
						chomp;
						my @line = split(/\t/,);
						if(defined $norm_depth{$line[2]}){
								$tumor_depth_posi{$line[0]}{$line[1]+1}=$line[2];
						}
				}
				close DPLIST;
				foreach my $group(sort keys %tumor_depth_posi){
						open(TDP,"$tdgdir/$CT/depth/tout$group.tsv");
						$header = <TDP>;chomp $header;
						@dp_col = split(/\t/,$header);
						foreach my $posi(keys %{$tumor_depth_posi{$group}}){
								if($dp_col[$posi] ne $tumor_depth_posi{$group}{$posi}){die "ERROR::patient colum position changed in $CT:tout$group.tsv file:$posi\n";}
						}
						while(<TDP>){
								chomp;
								my @line = split(/\t/,);
								if(!defined $gnomad_refminor{$line[0]}{$line[1]}){next;}
								foreach my $posi(keys %{$tumor_depth_posi{$group}}){
										if($line[$posi] >=8){
												print RM  "$dp_col[$posi]\t$line[0]\t$line[1]\tyes\n";
												print ARM "$dp_col[$posi]\t$line[0]\t$line[1]\tyes\n";
										}else{
												print RM  "$dp_col[$posi]\t$line[0]\t$line[1]\tno\n";
												print ARM "$dp_col[$posi]\t$line[0]\t$line[1]\tno\n";
										}
								}
						}
						close TDP;
				}
#print out alerady caluculate in solid&blood turn
				open(INRM,"$tdgdir/$CT/depth/ref_minor_coverage.tsv")or die "ERROR::cannot open $CT/depth/ref_minor_coverage.tsv";
				<INRM>;
				while(<INRM>){
						chomp;
						my @line = split(/\t/,);
						if((defined $focal_patient{$CT}{$line[0]})&&(!defined $bam_list_fordepth{$CT}{$line[0]})){
								print RM  "$_\n";
								print ARM "$_\n";
						}
				}
				close INRM;
				close RM;
		}
}
close ARM;
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

sub add_coverage( $ $ ){
		my ($scov,$refscov)=@_;
		open(SCOV,"$scov");
		<SCOV>;
		while(<SCOV>){
				chomp;
				my @line =split(/\t/,);
				${$refscov}{$line[0]}{$line[1]}{white}+=$line[2];
				${$refscov}{$line[0]}{$line[1]}{black}+=$line[3];
				${$refscov}{$line[0]}{$line[1]}{other}+=$line[4];
		}
		close SCOV;
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
