library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd('/Volumes/DR8TB2/tcga_rare_germ/control_gene/')
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
control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")

####################################### TCGA data ########################################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_race = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_race.tsv")
patien_all_info = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/all_patient_info.tsv")

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
                                 intron_variant = 'silent',
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
                                 `intron_variant,non_coding_transcript_variant`='silent',
                                 `intron_variant,NMD_transcript_variant`='NMD',
                                 `non_coding_transcript_exon_variant,non_coding_transcript_variant`='silent',
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
                                 `splice_donor_variant,coding_sequence_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,intron_variant`='splice',
                                 `splice_donor_variant,intron_variant`='splice',
                                 `splice_donor_variant,3_prime_UTR_variant,intron_variant`='splice',
                                 `splice_donor_variant,5_prime_UTR_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant`='splice',
                                 `splice_region_variant,5_prime_UTR_variant`='flank',
                                 `splice_region_variant,3_prime_UTR_variant`='flank',
                                 `splice_region_variant,intron_variant` = 'splice_region',
                                 `splice_region_variant,synonymous_variant`='silent',
                                 `splice_region_variant,stop_retained_variant`='silent',
                                 `start_lost,inframe_deletion`='truncating',
                                 `start_lost,splice_region_variant`='truncating',
                                 `stop_lost,inframe_deletion`='truncating',
                                 `stop_gained,protein_altering_variant`='truncating',
                                 `stop_gained,frameshift_variant,splice_region_variant`='truncating',
                                 `stop_gained,inframe_deletion`='truncating',
                                 `stop_gained,inframe_insertion`='truncating',
                                 `stop_gained,splice_region_variant`='truncating'))
  }
  options(scipen = 10) #startなどを指数表記のまま保存しないため
  extract_norm_maf=function(.bp){
    print(paste0("strat ",.bp))
    .colnames = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1',
                  'Tumor_Seq_Allele2',"Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2",
                  'Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position","FILTER",
                  't_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','Transcript_ID')
    .cols = .colnames %>>%
    {setNames(c('c','c','d','d','c','c', 'c','c','c', 'c','c','c','c','c','c', 'i','i','i','i','i','i','c'), .)} %>>%
    {do.call(readr::cols_only, as.list(.))}
    strip_maf = function(infile) {
      read_tsv(infile, comment='#', col_types=.cols) %>>%
        classify_consequence()
    }
    #.bp is bodypart of it cancer
    maf=tibble::tibble(file=list.files(paste0("/Volumes/areca42TB2/gdc/control_region/",.bp,"/maf"))) %>>%
      mutate(filename = paste0('/Volumes/areca42TB2/gdc/control_region/',.bp,"/maf/",file),
             patient_id = str_replace(file,".maf",""))%>>%
      mutate(tbl=purrr::map(filename,~strip_maf(.))) %>>%
      unnest() %>>%#(?.)%>>%
      dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                    ref=Reference_Allele,t_allele1=Tumor_Seq_Allele1,t_allele2=Tumor_Seq_Allele2,
                    n_allele1=Match_Norm_Seq_Allele1,n_allele2=Match_Norm_Seq_Allele2) %>>%
      mutate(soma_or_germ=ifelse((t_allele1==n_allele1)&(t_allele2==n_allele2),"germline","somatic"),
             t_genotype=ifelse(t_allele1==t_allele2,"homo","hetero"),
             n_genotype=ifelse(n_allele1==n_allele2,"homo","hetero")) %>>%
      mutate(LOH=ifelse((t_genotype=="homo")&(n_genotype=="hetero"),"LOH","no"))
  }
  tibble::tibble(cancer_type=c("brca","crc","gbm","hnsc","kcc","lgg",
                               "luad","lusc","ov","prad","thca","ucec")) %>>%
    mutate(purrr::map(cancer_type,~extract_norm_maf(.))) %>>%
    unnest() %>>%
    dplyr::select(-file,-filename) %>>%
    mutate(cancer_type = toupper(cancer_type)) %>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_region.maf.gz")
}

norm_maf_all_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_region.maf.gz") %>>%
  inner_join(control_genes%>>%dplyr::select(-role))%>>%
  filter(mutype!="flank",mutype!="splice_region") %>>%
  filter(Consequence!="intron_variant",Consequence!="intron_variant,non_coding_transcript_variant") %>>%
  left_join(patient_list) %>>%
  #filter(gene_symbol != "GDF2", gene_symbol != "GPRIN2") %>>%
  #filter(gene_symbol != "MRC1", gene_symbol != "CCDC168")%>>% #ExACでMAF==0なのにtcgaではいっぱいなsiteが多数
  mutate(LOH=ifelse((soma_or_germ == "somatic" & LOH =="no" & ref != n_allele2),"back_mutation",LOH))%>>%
  dplyr::select(-cancer_type,-age,-gender) %>>%left_join(patient_list)%>>%
  filter(!is.na(cancer_type))
write_df(norm_maf_all_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/norm_maf_all_cont.tsv.gz")
######  coverage file
coverage_all_by_cancer_type_cont =
  read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/coverage_all_cont.tsv.gz")
coverage_all_cont = coverage_all_by_cancer_type_cont %>>%
  mutate(an_cancer = an_white + an_black + an_other) %>>%
  group_by(chr,start) %>>%
  summarise(an_cancer = sum(an_cancer))


tally_norm_maf_cont = norm_maf_all_cont %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt)%>>%
  filter(!(chr=="chrX" & gender=="male" & allele=="n_allele1")) %>>%
  mutate(homo=ifelse(n_genotype=="homo",ifelse(chr=="X" & gender=="male",0,1),0))%>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
            Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
            cDNA_position=first(cDNA_position),CDS_position=first(CDS_position),
            Protein_position=first(Protein_position)) %>>%
  ungroup() %>>%
  left_join(coverage_all_cont)
write_df(tally_norm_maf_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/tally_nomr_maf_cont.tsv.gz")

######################################## gnomAD (control cases data) ##############################################
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))
if(0){
  ###ref_minor_list
  control_gnomad %>>%
    filter(AF > 0.5| AF_white > 0.5| (AF_black > 0.5)) %>>%
    write_df("/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_list_gnomAD.tsv")
  ### do make_coverage.pl
}
ref_minor_focal_cont = 
  read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/ref_minor_coverage_by_patient.tsv.gz")%>>%
  left_join(patient_list)


