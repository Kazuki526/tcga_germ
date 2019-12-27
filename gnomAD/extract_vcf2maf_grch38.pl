#!/usr/bin/perl
use warnings;
use strict;

my $gnomAD_version = "gnomad.exomes.r2.1.1.sites";

# this script execute in gnomAD exon raw file (by chromosome) dir
my @ls = `ls file_grch38|grep $gnomAD_version`;
if(scalar(@ls)!= 24){die "ERROR::doing on wrong dir?\n";}

#make extract dir
mkdir "maf38";
mkdir "maf38/all_maf";
mkdir "maf38/non_cancer_maf";
mkdir "maf38/control_maf";

#vep extract colum
my @vep_focal = qw(SYMBOL Consequence IMPACT Gene HGVSc HGVSp cDNA_position CDS_position Protein_position Amino_acids Codons
				   STRAND CANONICAL SIFT PolyPhen LoF LoF_filter Feature);
my @race = qw(afr sas amr eas nfe fin asj oth);
my $ac_col="AC\tAN\tnhomalt";
foreach my $race (@race){
		$ac_col .= "\tAC_$race\tAN_$race\tnhomalt_$race";
}
my @ac_col = split(/\t/,$ac_col);
my @chr = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);


#read transcript_ID to cDNA_length
open(ENST,"/Volumes/areca42TB/GRCh37-lite/ensembl_info/ensembl_transcript_length.tsv");
my %transcript=();
<ENST>;#header
while(<ENST>){
		chomp;
		my @line = split(/\t/,);
		$transcript{$line[2]}=$line[4];
}
close ENST;

# Prioritize Sequence Ontology terms in order of severity, as estimated by Ensembl:
# http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
sub GetEffectPriority {
		my ( $effect ) = @_;
		$effect = '' unless( defined $effect );
		my %effectPriority = (
						'transcript_ablation' => 1, # A feature ablation whereby the deleted region includes a transcript feature
						'exon_loss_variant' => 1, # A sequence variant whereby an exon is lost from the transcript
						'splice_donor_variant' => 2, # A splice variant that changes the 2 base region at the 5' end of an intron
						'splice_acceptor_variant' => 2, # A splice variant that changes the 2 base region at the 3' end of an intron
						'stop_gained' => 3, # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
						'frameshift_variant' => 3, # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
						'stop_lost' => 3, # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
						'start_lost' => 4, # A codon variant that changes at least one base of the canonical start codon
						'initiator_codon_variant' => 4, # A codon variant that changes at least one base of the first codon of a transcript
						'disruptive_inframe_insertion' => 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
						'disruptive_inframe_deletion' => 5, # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
						'inframe_insertion' => 5, # An inframe non synonymous variant that inserts bases into the coding sequence
						'inframe_deletion' => 5, # An inframe non synonymous variant that deletes bases from the coding sequence
						'protein_altering_variant' => 5, # A sequence variant which is predicted to change the protein encoded in the coding sequence
						'missense_variant' => 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
						'conservative_missense_variant' => 6, # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
						'rare_amino_acid_variant' => 6, # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
						'transcript_amplification' => 7, # A feature amplification of a region containing a transcript
						'splice_region_variant' => 8, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
						'stop_retained_variant' => 9, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
						'synonymous_variant' => 9, # A sequence variant where there is no resulting change to the encoded amino acid
						'incomplete_terminal_codon_variant' => 10, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
						'coding_sequence_variant' => 11, # A sequence variant that changes the coding sequence
						'mature_miRNA_variant' => 11, # A transcript variant located with the sequence of the mature miRNA
						'exon_variant' => 11, # A sequence variant that changes exon sequence
						'5_prime_UTR_variant' => 12, # A UTR variant of the 5' UTR
						'5_prime_UTR_premature_start_codon_gain_variant' => 12, # snpEff-specific effect, creating a start codon in 5' UTR
						'3_prime_UTR_variant' => 12, # A UTR variant of the 3' UTR
						'non_coding_exon_variant' => 13, # A sequence variant that changes non-coding exon sequence
						'non_coding_transcript_exon_variant' => 13, # snpEff-specific synonym for non_coding_exon_variant
						'non_coding_transcript_variant' => 14, # A transcript variant of a non coding RNA gene
						'nc_transcript_variant' => 14, # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
						'intron_variant' => 14, # A transcript variant occurring within an intron
						'intragenic_variant' => 14, # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
						'INTRAGENIC' => 14, # snpEff-specific synonym of intragenic_variant
						'NMD_transcript_variant' => 15, # A variant in a transcript that is the target of NMD
						'upstream_gene_variant' => 16, # A sequence variant located 5' of a gene
						'downstream_gene_variant' => 16, # A sequence variant located 3' of a gene
						'TFBS_ablation' => 17, # A feature ablation whereby the deleted region includes a transcription factor binding site
						'TFBS_amplification' => 17, # A feature amplification of a region containing a transcription factor binding site
						'TF_binding_site_variant' => 17, # A sequence variant located within a transcription factor binding site
						'regulatory_region_ablation' => 17, # A feature ablation whereby the deleted region includes a regulatory region
						'regulatory_region_amplification' => 17, # A feature amplification of a region containing a regulatory region
						'regulatory_region_variant' => 17, # A sequence variant located within a regulatory region
						'regulatory_region' =>17, # snpEff-specific effect that should really be regulatory_region_variant
						'feature_elongation' => 18, # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
						'feature_truncation' => 18, # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
						'intergenic_variant' => 19, # A sequence variant located in the intergenic region, between genes
						'intergenic_region' => 19, # snpEff-specific effect that should really be intergenic_variant
						'' => 20
						);
		unless( defined $effectPriority{$effect} ) {
				warn "WARNING: Unrecognized effect \"$effect\". Assigning lowest priority!\n";
				return 20;
		}
		return $effectPriority{$effect};
}

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(5);
foreach my $chr (@chr){
#fork start
		print "start chr$chr\n";
		$pm->start and next;
		&print_extracted_file("$chr");
#fork end 
		$pm->finish;
}
$pm->wait_all_children;
exit;

sub print_extracted_file( $ ){
		my $chr = $_[0];
		my $infile = "file_grch38/$gnomAD_version.$chr.liftover_grch38.vcf.bgz";
		open(VCF,"gunzip -c $infile|");
		open(ALL,"|gzip -c >maf38/all_maf/gnomAD_chr$chr.maf.gz");
		open(CAN,"|gzip -c >maf38/non_cancer_maf/non_cancer_chr$chr.maf.gz");
		open(CONT,"|gzip -c >maf38/control_maf/control_chr$chr.maf.gz");
		print ALL "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."$ac_col\n";
		print CAN "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."$ac_col\n";
		print CONT "chr\tposi\tref\talt\tfilter\t". join("\t",@vep_focal) ."$ac_col\n";
		my %vepcol=();
		my $info_check=3;
		while(<VCF>){
				if($_ =~ /^#/){
						if($_ =~ /^##INFO=<ID=vep/){
								chomp;
								%vepcol = &vep2colum($_);
						}
						next;
				}
				my @line = split(/\t/,);
				if($line[6] =~/AC0/){next;}
				my $vepout = &vepout($line[7],\%vepcol);
				if($vepout eq ""){next;}
				if(length($line[3]) != length($line[4])){
						my ( $ref_length, $var_length ) = ( length( $line[3] ), length( $line[4] ));
# Backup the VCF-style position and REF/ALT alleles, so we can use it later
						my ( $pos, $ref, $var ) = ( $line[1], $line[3], $line[4] );
# Remove any prefixed reference bps from all alleles, using "-" for simple indels
						while( $ref and $var and substr( $ref, 0, 1 ) eq substr( $var, 0, 1 ) and $ref ne $var ) {
								( $ref, $var ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var );
								--$ref_length; --$var_length; ++$pos;
						}
						$pos =($ref eq "-" ? $pos-1 : $pos);
						($line[1],$line[3],$line[4])=($pos,$ref,$var);
				}
				my $outbase = "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[6]$vepout";
				my %info = &info2hash($line[7]);
				my $allout ="";
				foreach my $out_info (@ac_col){
						$allout .= "\t$info{$out_info}";
				}
				print ALL "$outbase\t$allout\n";
				if($info{non_cancer_AC}!=0){
						my $noncanout ="";
						foreach my $out_info(@ac_col){
								my $out_info_nc=$info{"non_cancer_$out_info"};
								$noncanout .= "\t$out_info_nc";
						}
						print CAN "$outbase\t$noncanout\n";
				}
				if($info{controls_AC}!=0){
						my $contout ="";
						foreach my $out_info(@ac_col){
								my $out_info_co=$info{"controls_$out_info"};
								$contout .= "\t$out_info_co";
						}
						print CONT "$outbase\t$contout\n";
				}
		}
		close VCF;
		close ALL;
		close CAN;
		close CONT;
}



sub vep2colum( $ ){
		my $info = $_[0];
		my $vepinfo="";
		if($info =~ /^##INFO=.*Format: (Allele.*LoF_info)\">$/){$vepinfo=$1;}else{die "ERROR::INFO vep is not matched\n";}
		my @vepinfo = split(/\|/,$vepinfo);
		my %vepcol=();
		for(my $i=0;$i < scalar(@vepinfo);$i++){
				$vepcol{$vepinfo[$i]}=$i;
		}
		foreach my $vep_focal(@vep_focal){
				if(!defined $vepcol{$vep_focal}){die "ERROR::not exist $vep_focal on INFO of vep\n";}
		}
		return(%vepcol);
}

sub most_effect( $ $ ){
		my @vep_text = split(/,/,$_[0]);
		my %vepcol = %{$_[1]};
		my %all_vep =();
		for(my $i=0;$i<scalar(@vep_text);$i++){
				$vep_text[$i] =~ s/\&/,/g;
				my @vep = split(/\|/,$vep_text[$i]);
				if(($vep[$vepcol{IMPACT}] eq "MODIFIER") || ($vep[$vepcol{BIOTYPE}] ne "protein_coding")){next;}
				my ($Consequence)=split(",",$vep[$vepcol{Consequence}]);
				$all_vep{$i}{effect_prio}=&GetEffectPriority($Consequence);
				$all_vep{$i}{length}=$transcript{$vep[$vepcol{Feature}]};
		}
		if(scalar(keys %all_vep) ==0){return("");
		}else{
				my @arranged_index = sort{
						$all_vep{$a}{effect_prio} <=> $all_vep{$b}{effect_prio} ||
						$all_vep{$a}{length} <=> $all_vep{$b}{length}
				}keys %all_vep;
				return($vep_text[$arranged_index[0]]);
		}
}

sub vepout( $ $ ){
		my $info = $_[0];
		my %vepcol = %{$_[1]};
		my $vep_text;
		if($info =~ /;vep=([^;]+)$/){$vep_text=$1;
		}elsif($info =~ /;vep=([^;]+);was_mixed$/){$vep_text=$1;
		}else{die  "ERROR::this line have no vep result\n$info\n";}
		my $most_effect_vep = &most_effect($vep_text,\%vepcol);
		my $out = "";
		if($most_effect_vep eq ""){
				return($out);
		}else{
				my @vep = split(/\|/,$most_effect_vep);
				foreach my $vep_focal(@vep_focal){
						if(!defined $vep[$vepcol{$vep_focal}]){$out .="\t";
						}else{$out .= "\t$vep[$vepcol{$vep_focal}]";}
				}
				return($out);
		}
}
sub info2hash( $ ){
		my $info = $_[0];
		my @info = split(/;/,$info);
		my %out =();
		for(my $i=0;$i<scalar(@info);$i++){
				if($info[$i] =~ /^([^=]+)=(\d+)$/){$out{$1}=$2;}
		}
		return(%out);
}
		




