#!/usr/bin/zsh

patient_dir=`pwd`
mkdir maf

grep -v "^chr.._" vcf/$1.snp.vcf|grep -v "^chr._" |sed s/chr// >vcf/$1.snp2maf.vcf

export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB} #vcftools neeed this
export VEP_PATH=$HOME/vep
export VEP_DATA=$HOME/.vep
export PERL5LIB=$VEP_PATH:$PERL5LIB
export PATH=$VEP_PATH/htslib:$HOME/vep/samtools/bin:$PATH 
cd ~/vep/vcf2maf

perl vcf2maf.pl --input-vcf $patient_dir/vcf/$1.snp2maf.vcf --output-maf $patient_dir/maf/$1.snp.maf --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta ~/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --ncbi-build GRCh38
