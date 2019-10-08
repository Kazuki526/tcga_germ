#!/usr/bin/bash

pwd=`pwd`


export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this
export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$HOME/vep/samtools/bin:$PATH 
cd ~/vep/vcf2maf

perl vcf2maf.pl --input-vcf $pwd/annotation/"$1"_snp_like.vcf --output-maf $pwd/annotation/"$1"_snp.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --ncbi-build GRCh37
perl vcf2maf.pl --input-vcf $pwd/annotation/"$1"_indel_like.vcf --output-maf $pwd/annotation/"$1"_indel.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --ncbi-build GRCh37
