#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my ($start_time) =(time())[0];

my %patient_inf =();
my @patient_list =();
open(PLIST,"$ENV{HOME}/git/tcga_germ/variant_call/white_focal_patient.tsv") or die "ERROR::patient list is correct??\n";
<PLIST>; #header
while(<PLIST>){
		chomp;
		my @line = split(/\t/,);
		push(@patient_list,$line[1]);
		$patient_inf{$line[1]}=$line[3];
}
close PLIST;

open(LOG,">done_patient.txt");
my $max_processes = 5;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $pid (@patient_list){
		my ($now_time) =(time())[0];
		my $spent = $now_time - $start_time;
		if($spent/24/3600 > 25){last;}
#		$pm->start and next; #do fork

		mkdir $pid;
		chdir $pid;
		`perl $ENV{HOME}/git/tcga_germ/variant_call/variant_call_coverage_by_pid.pl $pid $patient_inf{$pid}`;
		`zsh $ENV{HOME}/git/tcga_germ/variant_call/vcf2maf.sh >vcf2maf.out 2>&1`;
		`rm -rf norm_bam tumor_bam`;

		print "done $pid\n";print LOG "done $pid\n";
		last;
#		$pm->finesh;
}
#$pm->wait_all_children;
exit;
