#!/usr/bin/perl
use strict;
use warnings;

my $input_maf = "$ARGV[0]";
(-e $input_maf) or die "ERROR::please input maf file\n";
my $input="";
if($input_maf =~ /^(.+)\.maf.gz$/){$input=$1;
}elsif($input_maf =~ /^(.+)\.maf$/){$input=$1;
}else{die "ERROR:: input file is maf file?? $input_maf\n";}

my $liftover=$ENV{'HOME'}."/liftover/liftover";
(-e $liftover) or die "ERROR:liftover not exist at $liftover\n";
my $chainfile=$ENV{"HOME"}."/liftover/hg38ToHg19.over.chain";
(-e $chainfile) or die "ERROR:chain file not exist at $chainfile\n";
my $hg19="/Volumes/areca42TB/GRCh37-lite/GRCh37-lite.fa";
(-e $hg19 ) or die "ERROR:hg19 fasta not exist at $hg19\n";
my $hg38=$ENV{"HOME"}."/liftover/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
(-e $hg38 ) or die "ERROR:hg38 fasta not exist at $hg38\n";

my %variants=();
# make bed file for liftover
if($input_maf =~ /.gz$/){
		open(IN,"gunzip -c $input_maf|") or die "ERROR::cannot gunzip open $input_maf\n";
}else{
		open(IN,"$input_maf") or die "ERROR::cannot open $input_maf\n";
}
my $header=<IN>;chomp $header;
my %col=&header2hash($header);
open(BED,">$input.bed");
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my($chr,$posi,$ref) = ($line[$col{chr}],$line[$col{posi}],$line[$col{ref}]);
		$chr =~ s/^chr//;
		my $vinf = "chr$chr:$posi:$ref";
		if(defined $variants{$vinf}){next;}
		my ($vclass,$bed_region,$bed_out) = ("","","");
		if($ref eq "-"){
				$vclass="ins";
				$bed_region = "$chr:". ($posi-1) ."-". ($posi+1);
				$bed_out = "chr$chr\t". ($posi-1) ."\t". ($posi+1)."\t$vinf\n";
		}elsif(length($ref)>1){
				$vclass="del";
				my $vleng=length($ref);
				$bed_region = "$chr:". ($posi-2) ."-". ($posi+$vleng);
				$bed_out = "chr$chr\t". ($posi-2) ."\t". ($posi+$vleng) ."\t$vinf\n";
		}elsif($ref =~ /^[ATGCatgc]$/){
				$vclass="snv";
				$bed_region = "$chr:$posi-$posi";
				$bed_out = "chr$chr\t". ($posi-1) ."\t$posi\t$vinf\n";
		}else{die "what variant class??\n$_\n";}
		print BED "$bed_out";
		$variants{$vinf}{vclass}=$vclass;
		$variants{$vinf}{bed_region}=$bed_region;
}
close IN;
close BED;

my %remap = map{chomp;my @c=split("\t");$c[0]=~s/^chr//;($c[3],"$c[0]:$c[1]-$c[2]")}`$liftover $input.bed $chainfile /dev/stdout /dev/null 2>/dev/null`;

open(OUT,">$input"."_liftover_result.tsv");
print OUT "chr\tposi\tref\tliftovrer_result\tchange\n";
#check nucl difference between 
foreach my $vinf (keys%variants){
		my($chr,$posi,$ref) = split(/:/,$vinf);
		if(!defined $remap{$vinf}){print OUT "$chr\t$posi\t$ref\tliftover_error\n";next;}
		my($chr37,$start37,$end37);
		if($remap{$vinf} =~ /^(.+):(\d+)-(\d+)$/){($chr37,$start37,$end37)=($1,$2,$3);}else{print "lifover error?? $vinf => $remap{$vinf}\n";}
		if($chr ne "chr$chr37"){print OUT "$chr\t$posi\t$ref\tchr_change\n";next;}
		if($variants{$vinf}{vclass} eq "snv"){
				my $fa38 =&faidx($variants{$vinf}{bed_region},$hg38);
				if($ref ne $fa38){die "ERROR:SNV position error?? $vinf\n";}
				my $fa37 =&faidx("$chr37:$end37-$end37",$hg19);
				if($fa37 eq  $fa38){
						print OUT "$chr\t$posi\t$ref\tok\t$fa38>$fa37\n";
				}else{
						print OUT "$chr\t$posi\t$ref\tref_change\t$fa38>$fa37\n";
				}
		}else{
				my $fa38 =&faidx($variants{$vinf}{bed_region},$hg38);
				my $fa37 =&faidx("$chr37:$start37-$end37",$hg19);
				if($fa37 eq  $fa38){
						print OUT "$chr\t$posi\t$ref\tok\t$fa38>$fa37\n";
				}else{
						print OUT "$chr\t$posi\t$ref\tref_change\t$fa38>$fa37\n";
				}
		}
}
close OUT;



sub faidx ( $ $ ){
		my ($region, $fasta_file) = @_;
		my $fasta=`samtools faidx $fasta_file $region`;
		my $seq ="";
		foreach my$l (split(/\n/,$fasta)){
				if($l =~ /^>/){next;}
				$seq.=$l;
		}
		return($seq);
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
