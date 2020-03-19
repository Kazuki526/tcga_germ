#!/usr/bin/zsh

# $1 = patient_id
# $2 = each cancer type dir

cat $2/vcf/$1.indel.vcf|sed s/^chr// >$2/out/$1.vcf
cat $2/vcf/$1.snp.vcf  |sed -e '/^#/d' -e s/^chr// >>$2/out/$1.vcf


export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this
export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$HOME/vep/samtools/bin:$PATH 
cd ~/vep/vcf2maf

perl vcf2maf.pl --input-vcf $2/out/$1.vcf --output-maf $2/maf/$1.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38
