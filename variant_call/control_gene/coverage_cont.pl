#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my $contdir="/Volumes/areca42TB2/gdc/control_region";
(-e $contdir) or die "ERROR::changed $contdir??";
#check reference of bam are exist?
my $ref="/Volumes/areca42TB/GRCh38.d1.vd1.fa";
(-e $ref)or die "ERROR:not exist ref fasta:$ref\n";
# check token file
my $token_path=`ls $ENV{HOME}/git/innanlab/gdc|grep 'gdc-user-token'`;
if(!$token_path){die "!!ERROR!!:token file not exitst!!";}
chomp $token_path;
$token_path="$ENV{HOME}/git/innanlab/gdc/$token_path";
my $token=`cat $token_path`;
#check nkf was installed??
my $nkfpath=`which nkf`;chomp $nkfpath;
($nkfpath and -e $nkfpath) or die "ERROR:nkf was not installed. please do\nbrew install nkf\n";
#check control genes bed file
my $bed_file = "/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_cds.bed";
-e $bed_file or die "ERROR:: not exist $bed_file\n";
#chech gnomAD maf fiel
my $gnomAD_file = "/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf";
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
		$patient_inf{$line[$col{patient_id}]}{age}=$line[$col{age}];
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


my %pj2bp = ("coad" => "crc", "read" => "crc",
			 "kich" => "kcc","kirc" => "kcc","kirp" => "kcc");


my %bam_list_fordepth=();
my $line=0;
foreach my $pid (@patient_list){
		$line++;
		my $CT = $patient_inf{$pid}{cancer_type};
		if(defined$pj2bp{$patient_inf{$pid}{cancer_type}}){$CT=$pj2bp{$patient_inf{$pid}{cancer_type}};}
		my ($gender,$nnorm,$ntumor) = ($patient_inf{$pid}{gender},$patient_inf{$pid}{nnorm},$patient_inf{$pid}{ntumor});
		my @norm_file  = split(/,/,$patient_inf{$pid}{norm_file});
		my @tumor_file = split(/,/,$patient_inf{$pid}{tumor_file});
		if(!-e "$contdir/$CT"){mkdir "$contdir/$CT";}
		if(!-e "$contdir/$CT/depth"){mkdir "$contdir/$CT/depth";}
		my $norm_bam="";
		if($nnorm>1){
				$norm_bam = "$contdir/$CT/norm_bam/$pid.bam";
		}else{
				$norm_bam="$contdir/$CT/norm_bam/$norm_file[0]";
		}
		my $tumor_bam="";
		if($ntumor>1){
				$tumor_bam = "$contdir/$CT/tumor_bam/$pid.bam";
		}else{
				$tumor_bam="$contdir/$CT/tumor_bam/$tumor_file[0]";
		}
		$bam_list_fordepth{$CT}{n}++;
		$bam_list_fordepth{$CT}{pid}.="$pid\t";
		$bam_list_fordepth{$CT}{norm_bam}.="$norm_bam\t";
		$bam_list_fordepth{$CT}{tumor_bam}.="$tumor_bam\t";
}


# make depth file of top_driver_genes' CDS region by 50sample
$line=0;
my $max_processes = 0;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $CT (sort keys %bam_list_fordepth){
		$line++;
		$pm->start and next; #do fork
		print "$line:make depth start $CT\n";
		my $depth_dir = "$contdir/$CT/depth";
		my $num = $bam_list_fordepth{$CT}{n};
		my @pid = split(/\t/,$bam_list_fordepth{$CT}{pid});
		my @norm_bam = split(/\t/,$bam_list_fordepth{$CT}{norm_bam});
		my $error="";
#		foreach my $bam (@norm_bam){
#				my $bam_check = &bam_check($bam);
#				if($bam_check ne "ok"){$error=$bam;last;}
#		}
		my @tumor_bam = split(/\t/,$bam_list_fordepth{$CT}{tumor_bam});
#		foreach my $bam (@tumor_bam){
#				my $bam_check = &bam_check($bam);
#				if($bam_check ne "ok"){$error=$bam;last;}
#		}
		if($error=~/./){
				if($max_processes >1){die "$CT:$error have error\nstop $CT make depth file\n";
				}else{print "$CT:$error have error\nstop $CT make depth file\n";next;
				}
		}
		open(DLIST,">$depth_dir/depth_file_list.tsv") or die "ERROR::cannnot print out $CT dpeht file list\n";
		print DLIST "group\tcolmn_num\tpatient_id\tnorm_bam\ttumor_bam\n";
		for(my $group=1;$group-1 < $num/50;$group++){
				my $group_n = $num -($group-1)*50;
				if($group_n >50){$group_n=50;}
				for(my $i=1;$i <= $group_n;$i++){
						print DLIST "$group\t$i\t$pid[($group-1)*50+$i-1]\t$norm_bam[($group-1)*50+$i-1]\t$tumor_bam[($group-1)*50+$i-1]\n";
				}
				if((-e "$depth_dir/out$group.tsv")&&(-e "$depth_dir/tout$group.tsv")){next;}
				my @now_pid = @pid[($group-1)*50..($group-1)*50+$group_n-1];
				my @now_norm = @norm_bam[($group-1)*50..($group-1)*50+$group_n-1];
				my @now_tumor = @tumor_bam[($group-1)*50..($group-1)*50+$group_n-1];
				if(!-e "$depth_dir/out$group.tsv"){
						open(ND,">$depth_dir/out$group.tsv");
						print ND "chr\tposition\t".join("\t",@now_pid)."\n";
						close ND;
						`samtools depth -a -q 10 -b $bed_file @now_norm >>$depth_dir/out$group.tsv`;
				}
				if(!-e "$depth_dir/tout$group.tsv"){
						open(TD,">$depth_dir/tout$group.tsv");
						print TD "chr\tposition\t".join("\t",@now_pid)."\n";
						close TD;
						`samtools depth -a -q 10 -b $bed_file @now_tumor >>$depth_dir/tout$group.tsv`;
				}
				print "$line:finish $CT of ". (($group-1)*50+$group_n) ."/ $num\n";
		}
		close DLIST;
		print "$line:finish make depth $CT\n";
		$pm->finish;
}
$pm->wait_all_children;

#read top driver genes bed file
my %posi2role=();
open(CGBED,"$bed_file");
while(<CGBED>){
		chomp;
		my @line = split(/\t/,);
		my $gene="";
		if($line[3] =~ /^([^:]+):ENST/){$gene=$1;}
		for(my $posi=$line[1]+1;$posi <= $line[2];$posi++){
				$posi2role{$line[0]}{$posi}="control_gene";
		}
}
close CGBED;

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
open(ARM,"|gzip -c >$contdir/all_patient/ref_minor_coverage.tsv.gz");
print ARM "patient_id\tchr\tstart\tfocal\n";
open(APCOV,">$contdir/all_patient/patient_coverage.tsv");
print APCOV "patient_id\trole\tcoverage\n";
open(APDEP,">$contdir/all_patient/patient_depth.tsv");
print APDEP "patient_id\trole\tnorm_depth\ttumor_depth\tlength_bp\n";
my %all_site_coverage = ();
foreach my $CT (sort keys %bam_list_fordepth){
		$line++;
		my $depth_dir = "$contdir/$CT/depth";
		print "$line:coverage start $CT\n";
		my $patient_cov_file = "$depth_dir/coverage_by_patient.tsv";
		my $patient_dep_file = "$depth_dir/depth_by_patient.tsv";
		my $site_cov_file = "$depth_dir/coverage_by_site.tsv";
		my $ref_minor_file = "$depth_dir/ref_minor_coverage.tsv";
#if already done or on downloading
		if((-e $patient_cov_file)&&(-e $patient_dep_file)&&(-e $site_cov_file)&&(-e $ref_minor_file)){
				&add_coverage($site_cov_file,\%all_site_coverage);
				open(PCOV,"$patient_cov_file");
				<PCOV>;
				while(<PCOV>){print APCOV $_;}
				close PCOV;
				open(PDEP,"$patient_dep_file");
				<PDEP>;
				while(<PDEP>){print APCOV $_;}
				close PDEP;
				open(RM,"$ref_minor_file");
				<RM>;
				while(<RM>){print ARM $_;}
				close RM;
				print  "$line:coverage already done $CT\n";
				next;
		}elsif(!-e "$depth_dir/out1.tsv"){
				print  "$line:still downloading $CT\n";
				next;
		}

#calculate coverage and depth
		my $num = $bam_list_fordepth{$CT}{n};
		my %patient_coverage=();
		my %patient_depth=();
		my %site_coverage=();
		open(RM,">$ref_minor_file");
		print RM "patient_id\tchr\tstart\tfocal\n";
		for(my $group=1;$group-1 < $num/50;$group++){
				open(NDP,"$depth_dir/out$group.tsv");
				open(TDP,"$depth_dir/tout$group.tsv");
				my $nheader=<NDP>;my $theader=<TDP>;
				if($nheader ne $theader){
						die "ERROR::$CT:header different between depth/out$group.tsv and dpeht/tout$group.tsv\n";
				}
				chomp $nheader;
				my @dp_col = split(/\t/,$nheader);
				while(my $nline = <NDP>){
						my $tline = <TDP>;
						my @nline = split(/\t/,$nline);
						my @tline = split(/\t/,$tline);
						if("$nline[0]:$nline[1]" ne "$tline[0]:$tline[1]"){
								die "position is different at $CT\nout$group.tsv $nline[0]:$nline[1]\ntout$group.tsv $tline[0]:$tline[1]\n";
						}
						for(my $i =2;$i<scalar(@nline);$i++){
								if($patient_inf{$dp_col[$i]}{age} eq "NA"){next;}
								$patient_depth{$dp_col[$i]}{norm}+=$nline[$i];
								$patient_depth{$dp_col[$i]}{tumor}+=$tline[$i];
								$patient_depth{$dp_col[$i]}{bp}++;
								if(($nline[$i] >=8) && ($tline[$i] >=8)){
										$site_coverage{$nline[0]}{$nline[1]}{$patient_inf{$dp_col[$i]}{race}}+=2;
										$all_site_coverage{$nline[0]}{$nline[1]}{$patient_inf{$dp_col[$i]}{race}}+=2;
										$patient_coverage{$dp_col[$i]}++;
										if(defined $gnomad_refminor{$nline[0]}{$nline[1]}){
												print RM  "$dp_col[$i]\t$nline[0]\t$nline[1]\tyes\n";
												print ARM "$dp_col[$i]\t$nline[0]\t$nline[1]\tyes\n";
										}
								}else{
										if(defined $gnomad_refminor{$nline[0]}{$nline[1]}){
												print RM  "$dp_col[$i]\t$nline[0]\t$nline[1]\tno\n";
												print ARM "$dp_col[$i]\t$nline[0]\t$nline[1]\tno\n";
										}
								}
						}

				}
				close NDP;close TDP;

		}
		close RM;
		open(PCOV,">$patient_cov_file");
		print PCOV "patient_id\trole\tcoverage\n";
		foreach my $pid(sort keys %patient_coverage){
				print PCOV  "$pid\tcontrol\t$patient_coverage{$pid}\n";
				print APCOV "$pid\tcontrol\t$patient_coverage{$pid}\n";
		}
		close PCOV;
		open(PDEP,">$patient_dep_file");
		print PDEP "patient_id\trole\tnorm_depth\ttumor_depth\tlength_bp\n";
		foreach my $pid(sort keys %patient_depth){
				print PDEP  "$pid\tcontrol\t$patient_depth{$pid}{norm}\t$patient_depth{$pid}{tumor}\t$patient_depth{$pid}{bp}\n";
				print APDEP "$pid\tcontrol\t$patient_depth{$pid}{norm}\t$patient_depth{$pid}{tumor}\t$patient_depth{$pid}{bp}\n";
		}
		close PDEP;
		open(SCOV,">$site_cov_file");
		print SCOV "chr\tstart\tan_white\tan_black\tan_other\n";
		foreach my $chr (sort keys %site_coverage){
				foreach my $posi (sort{$a <=>$b}keys %{$site_coverage{$chr}}){
						if(!defined $site_coverage{$chr}{$posi}{white}){$site_coverage{$chr}{$posi}{white}=0;}
						if(!defined $site_coverage{$chr}{$posi}{black}){$site_coverage{$chr}{$posi}{black}=0;}
						if(!defined $site_coverage{$chr}{$posi}{other}){$site_coverage{$chr}{$posi}{other}=0;}
						print SCOV "$chr\t$posi\t$site_coverage{$chr}{$posi}{white}\t$site_coverage{$chr}{$posi}{black}\t$site_coverage{$chr}{$posi}{other}\n";
				}
		}
		close SCOV;
}
close ARM;
close APCOV;
close APDEP;
#print out all cancer type coverage
open(ASCOV,">$contdir/all_patient/site_coverage_all.tsv");
print ASCOV "chr\tstart\tan_white\tan_black\tan_other\n";
foreach my $chr (sort keys %all_site_coverage){
		foreach my $posi (sort{$a <=>$b}keys %{$all_site_coverage{$chr}}){
						if(!defined $all_site_coverage{$chr}{$posi}{white}){$all_site_coverage{$chr}{$posi}{white}=0;}
						if(!defined $all_site_coverage{$chr}{$posi}{black}){$all_site_coverage{$chr}{$posi}{black}=0;}
						if(!defined $all_site_coverage{$chr}{$posi}{other}){$all_site_coverage{$chr}{$posi}{other}=0;}
				print ASCOV "$chr\t$posi\t$all_site_coverage{$chr}{$posi}{white}\t$all_site_coverage{$chr}{$posi}{black}\t$all_site_coverage{$chr}{$posi}{other}\n";
		}
}
close ASCOV;
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
