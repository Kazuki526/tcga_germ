library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd('/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene')
write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}

######################################## driver gene infomation ##########################################
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%dplyr::rename(gene_symbol=gene)%>>%
  mutate(role=ifelse(is.na(role),"oncogene/TSG",role))%>>%
  inner_join(read_tsv("~/git/driver_genes/onlytop105/top_driver105.bed",
                      col_names = c("chr","start","end","id","info","strand"))%>>%
               tidyr::separate(id,into=c("gene_symbol","ENSG"),sep=";")%>>%
               filter(chr != "chrX")%>>%dplyr::select(gene_symbol))%>>%
  dplyr::select(gene_symbol,role)
write_df(driver_genes,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/driver_genes.tsv")

############################################# patient infomation ############################################
#patient_list = read_tsv("/Volumes/areca42TB/tcga/all_patient/patient_list.tsv")
patient_list = read_tsv("~/git/tcga_germ/variant_call/focal_patient_list.tsv")%>>%
  mutate(cancer_type=str_extract(cancer_type,"[A-Z]+$"))%>>%filter(!is.na(age))
#write_df(patient_list,"/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")

############################################# TCGA maf files #############################################
if(0){
  classify_consequence = function(.data) {
    mutate(.data,
           mutype= dplyr::recode(Consequence,
                                 downstream_gene_variant = 'flank',
                                 `3_prime_UTR_variant` = 'flank',
                                 `3_prime_UTR_variant,NMD_transcript_variant`='NMD',
                                 upstream_gene_variant = 'flank',
                                 `5_prime_UTR_variant` = 'flank',
                                 frameshift_variant = 'truncating',
                                 frameshift_variant = 'truncating',
                                 inframe_deletion = 'inframe_indel',
                                 inframe_insertion = 'inframe_indel',
                                 intron_variant = 'intron',
                                 splice_region_variant = 'splice_region',
                                 coding_sequence_variant = 'missense',
                                 missense_variant = 'missense',
                                 stop_gained = 'truncating',
                                 stop_lost = 'truncating',
                                 stop_retained_variant = 'silent',
                                 synonymous_variant = 'silent',
                                 splice_acceptor_variant = 'splice',
                                 splice_donor_variant = 'splice',
                                 protein_altering_variant = 'inframe_indel',
                                 start_lost = 'truncating',
                                 `coding_sequence_variant,3_prime_UTR_variant`='flank',
                                 `coding_sequence_variant,5_prime_UTR_variant`='flank',
                                 `stop_gained,frameshift_variant` = 'truncating',
                                 `stop_gained,frameshift_variant,start_lost` = 'truncating',
                                 `missense_variant,splice_region_variant`='missense',
                                 `intron_variant,non_coding_transcript_variant`='intron',
                                 `intron_variant,NMD_transcript_variant`='NMD',
                                 `non_coding_transcript_exon_variant,non_coding_transcript_variant`='non_coding',
                                 `frameshift_variant,splice_region_variant`='truncating',
                                 `frameshift_variant,start_lost`='truncating',
                                 `frameshift_variant,stop_lost`='truncating',
                                 `inframe_deletion,splice_region_variant`='splice_region',
                                 `inframe_insertion,splice_region_variant`='splice_region',
                                 `frameshift_variant,stop_lost,splice_region_variant`='truncating',
                                 `frameshift_variant,stop_retained_variant`='inframe_indel',
                                 `protein_altering_variant,splice_region_variant`='inframe_indel',
                                 `splice_acceptor_variant,intron_variant`='splice',
                                 `splice_acceptor_variant,coding_sequence_variant`='splice',
                                 `splice_acceptor_variant,coding_sequence_variant,intron_variant`='splice',
                                 `splice_acceptor_variant,5_prime_UTR_variant`='splice',
                                 `splice_acceptor_variant,5_prime_UTR_variant,intron_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,intron_variant`='splice',
                                 `splice_donor_variant,intron_variant`='splice',
                                 `splice_donor_variant,3_prime_UTR_variant,intron_variant`='splice',
                                 `splice_donor_variant,5_prime_UTR_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant`='splice',
                                 `splice_donor_variant,splice_acceptor_variant,intron_variant`='splice',
                                 `splice_region_variant,5_prime_UTR_variant`='flank',
                                 `splice_region_variant,3_prime_UTR_variant`='flank',
                                 `splice_region_variant,intron_variant` = 'splice_region',
                                 `splice_region_variant,synonymous_variant`='silent',
                                 `splice_region_variant,stop_retained_variant`='silent',
                                 `start_lost,inframe_deletion`='truncating',
                                 `start_lost,splice_region_variant`='truncating',
                                 `stop_gained,start_lost`='truncating',
                                 `stop_lost,inframe_deletion`='truncating',
                                 `stop_lost,splice_region_variant`='truncating',
                                 `stop_gained,protein_altering_variant`='truncating',
                                 `stop_gained,frameshift_variant,splice_region_variant`='truncating',
                                 `stop_gained,inframe_deletion`='truncating',
                                 `stop_gained,inframe_insertion`='truncating',
                                 `stop_gained,splice_region_variant`='truncating'))
  }
  ct_exchange = tibble(cancer_type=c("BRCA","COAD","READ","GBM","LGG","KICH","KIRC","KIRP","LUAD","LUSC"),
                       bp=c("breast","colorectal","colorectal","brain","brain","kidney","kidney","kidney","lung","lung"))
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                'Tumor_Seq_Allele2',"Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2",
                'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position","FILTER",
                't_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','Transcript_ID')
  .cols = .colnames %>>%
  {setNames(c('c','c','d','d','c','c', 'c','c','c', 'c','c','c','c','c','c', 'i','i','i','i','i','i','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types = .cols) %>>%
      classify_consequence()
  }
  #blood driverd + solid derived normal
  patient_list%>>%left_join(ct_exchange)%>>%mutate(bp=ifelse(is.na(bp),tolower(cancer_type),bp))%>>%
    mutate(filename=paste('/Volumes/areca42Tb2/gdc/top_driver_gene/',bp,"/maf/",patient_id,".maf",sep="")) %>>%
    dplyr::select(patient_id,cancer_type,filename)%>>%
    mutate(maf=purrr::map(filename,~strip_maf(.))) %>>%
    dplyr::select(-filename)%>>%
    unnest() %>>%#(?.)%>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,t_allele1=Tumor_Seq_Allele1,t_allele2=Tumor_Seq_Allele2,
                  n_allele1=Match_Norm_Seq_Allele1,n_allele2=Match_Norm_Seq_Allele2) %>>%
    inner_join(driver_genes%>>%dplyr::select(gene_symbol,role)) %>>%
    mutate(soma_or_germ=ifelse((t_allele1==n_allele1)&(t_allele2==n_allele2),"germline","somatic"),
           t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
           n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
    mutate(LOH=ifelse((t_genotype=="homo")&(n_genotype=="hetero"),"LOH","no"))%>>%
    mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH))%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_topdriver_maf.tsv.gz")
  #blood derived only
  patient_list%>>%left_join(ct_exchange)%>>%mutate(bp=ifelse(is.na(bp),tolower(cancer_type),bp))%>>%
    filter(!is.na(bloodnorm))%>>%
    mutate(filename=paste('/Volumes/areca42Tb2/gdc/top_driver_gene/blood/',bp,"/maf/",patient_id,".maf",sep="")) %>>%
    dplyr::select(patient_id,cancer_type,filename)%>>%
    mutate(maf=purrr::map(filename,~strip_maf(.))) %>>%
    dplyr::select(-filename)%>>%
    unnest() %>>%#(?.)%>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,t_allele1=Tumor_Seq_Allele1,t_allele2=Tumor_Seq_Allele2,
                  n_allele1=Match_Norm_Seq_Allele1,n_allele2=Match_Norm_Seq_Allele2) %>>%
    inner_join(driver_genes%>>%dplyr::select(gene_symbol,role)) %>>%
    mutate(soma_or_germ=ifelse((t_allele1==n_allele1)&(t_allele2==n_allele2),"germline","somatic"),
           t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
           n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
    mutate(LOH=ifelse((t_genotype=="homo")&(n_genotype=="hetero"),"LOH","no"))%>>%
    mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH))%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_topdriver_maf.tsv.gz")
}


norm_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_topdriver_maf.tsv.gz")%>>%
  filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
  mutate(chr = ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
  filter(FILTER=="PASS",chr!="chrX")
write_df(norm_maf_all,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/norm_maf_all.tsv.gz")
blood_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_topdriver_maf.tsv.gz")%>>%
  filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
  mutate(chr = ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
  filter(FILTER=="PASS",chr!="chrX")
write_df(blood_maf_all,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/bloodnorm_maf_all.tsv.gz")
######################################## tally each site ####################################
tally_norm_maf = norm_maf_all%>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt,chr!="chrX")%>>%
  mutate(homo=ifelse(n_genotype=="homo",1,0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
            cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
  ungroup() %>>%
  left_join(read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/all_patient/site_coverage_all.tsv"))%>>%
  mutate(an_cancer = an_white+an_black+an_other)
write_df(tally_norm_maf,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/tally_norm_maf.tsv.gz")
####################################### ref minor #############################################
ref_minor=read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/all_patient/ref_minor_coverage.tsv.gz")
write_df(ref_minor,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/ref_minor_coverage.tsv.gz")
ref_minor_blood = read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/blood/all_patient/ref_minor_coverage_blood.tsv.gz")
write_df(ref_minor_blood,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/ref_minor_coverage_blood.tsv.gz")
coverage=read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/all_patient/patient_coverage.tsv") %>>%
  tidyr::spread(role,coverage)%>>%mutate(all=TSG+oncogene+`oncogene/TSG`)
write_df(coverage,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/coverage.tsv.gz")
####################################################################
#pathogenic variantをもつ患者のlist
patient_with_ps_ = read_tsv("/Volumes/areca42TB/tcga/burden_plot/cell2018/charger/tcga/tcga_all_variants_grch38_charger.tsv") %>>%
  dplyr::rename(chr=Chromosome,start=Start,ref=Reference,alt=Alternate,gene_symbol=HUGO_Symbol,
                clinvar=ClinVar_Pathogenicity, acmg=ACMG_Classification, charger=CharGer_Classification,
                trait=ClinVar_Traits) %>>%
  filter(clinvar=="Pathogenic"|charger=="Pathogenic"|charger=="Likely Pathogenic")%>>%
  left_join(norm_maf_all) %>>%
  dplyr::select(gene_symbol,cancer_type,patient_id)
#write_df(patient_with_ps,"all_patient/pathogenic_site_list.tsv")
#この中でsignificantなのは？
patient_with_ps = tibble(cancer_type=c("BRCA", "OV",   "BRCA", "OV",   "KIRP","BRCA","LUAD","PRAD","KIRC","GBM",  "UCEC"),
                         gene_symbol=c("BRCA1","BRCA1","BRCA2","BRCA2","MET", "ATM", "ATM", "ATM", "BAP1","AXIN2","PTEN")) %>>%
  mutate(significance="significant") %>>%
  right_join(patient_with_ps_) %>>%
  count(cancer_type,significance,patient_id)%>>%
  dplyr::select(-n)
  
write_df(patient_with_ps,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")
