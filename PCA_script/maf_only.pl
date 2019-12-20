#!/usr/bin/perl
use warnings;
use strict;
use Parallel::ForkManager;

my $pwd = `pwd`;chomp $pwd;
if($pwd ne "/Volumes/areca42TB2/gdc/pancan_atlas_germ"){die "ERROR::done on wrong dir\n";}

my @chr=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);

my $max_processes = 3;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $outchr (@chr){
		$pm->start and next;
		print "start chr$outchr\n";
		`bash $ENV{HOME}/git/tcga_germ/PCA_script/vcf2maf.sh chr$outchr`;
		$pm->finish;
}
$pm->wait_all_children;
exit;
		

