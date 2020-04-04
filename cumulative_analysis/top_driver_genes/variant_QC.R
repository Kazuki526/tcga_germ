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
driver_genes=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/driver_genes.tsv")

####################################### TCGA data ########################################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
norm_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/norm_maf_all.tsv.gz")
blood_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/bloodnorm_maf_all.tsv.gz")
tally_norm_maf = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/tally_norm_maf.tsv.gz")
ref_minor = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/ref_minor_coverage.tsv.gz")
ref_minor_blood = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/ref_minor_coverage_blood.tsv.gz")
patient_with_ps = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")
######################################## gnomAD (control cases data) ##############################################
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",col_types = cols(LoF_filter="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
#exac = read_tsv("/Volumes/areca42TB/exac/file/extracted/exac_nontcga.maf")
#vcf_exac = read_tsv("/Volumes/areca42TB/exac/file/extracted/exac_nontcga.vcf")
#######################################################################################################
################################filtering duplicate? position #########################################
#######################################################################################################
if(0){
HWE_test_heterom = function(ac,an,hom){
  AF=ac/an
  .matrix = matrix(c(round(an*AF*(1-AF)),round(an/2*AF^2),ac-hom*2,hom), nrow = 2)
  fisher.test(.matrix,alternative = "less")$p.value
}

duplicate_site = tdg_gnomad %>>%
  mutate(hom_gad_hwe=(AC^2/AN)/2) %>>%
  filter(chr!="chrX",AC!=0,!is.na(nhomalt)) %>>%
  nest(-chr,-start,-ref,-alt) %>>%
  mutate(HWE_gnomAD = purrr::map(data,~HWE_test_heterom(.$AC,.$AN,.$nhomalt))) %>>%
  unnest() %>>%
  mutate(FDR_gnomAD=p.adjust(HWE_gnomAD)) %>>%filter(FDR_gnomAD<0.01)%>>%
  full_join(tally_norm_maf%>>%
              dplyr::select(chr,start,ref,alt,ac_cancer,hom_cancer,gene_symbol,mutype,an_cancer) %>>%
              mutate(hom_cancer_hwe = ac_cancer^2/an_cancer/2) %>>%
              filter(!is.na(an_cancer)) %>>%#filter(ac_cancer<an_cancer)%>>%
              nest(-chr,-start,-ref,-alt) %>>%
              mutate(HWE_cancer=purrr::map(data,~HWE_test_heterom(.$ac_cancer,.$an_cancer,.$hom_cancer))) %>>%
              unnest() %>>%
              mutate(FDR_cancer=p.adjust(HWE_cancer)) %>>%filter(FDR_cancer<0.01))
write_df(duplicate_site,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/duplicate_site.tsv")
}
duplicate_site = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/duplicate_site.tsv")

somatic_recurrent = norm_maf_all%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)

#normalでaltalt, tumorでrefaltとなってる際にnormalでrefのdepth=0のものだけ採用！
varscan_error = norm_maf_all %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(an_error = n()*2) %>>%
  ungroup() 
varscan_error_site = norm_maf_all %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  dplyr::select(patient_id,chr,start,ref,alt) %>>%
  mutate(varscan_error_focal="yes")


quality_filter= function(.data,.data_type="vcf",.fdr=0.01,.database="cancer",
                         .duplicate=T,.somatic=T,.varscan=F){
  .site = duplicate_site
  if(.database!="all"){
    if(.database=="cancer"){
      .site = .site %>>% filter(FDR_cancer < .fdr)
    }else if(.database=="gnomAD"){
      .site = .site %>>% filter(FDR_gnomAD < .fdr)
    }else{stop(paste0("database variable is wrong .database=",.database,"\ncancer, gnomAD is correct."))}
  }
  if(.varscan){
    if(.data_type=="vcf"){
      .data = .data %>>%
        left_join(varscan_error) %>>%
        mutate(ac_cancer = ifelse(!is.na(an_error),ac_cancer-an_error,ac_cancer),
               an_cancer = ifelse(!is.na(an_error),an_cancer-an_error,an_cancer)) %>>%
        dplyr::select(-an_error)
    }else if(.data_type=="maf"){
      .data = .data %>>%
        left_join(varscan_error_site) %>>%
        filter(is.na(varscan_error_focal)) %>>%
        group_by(patient_id,gene_symbol,Protein_position) %>>%
        mutate(same_codon=n()) %>>%
        filter(same_codon==1) %>>%
        ungroup() %>>%
        dplyr::select(-varscan_error_focal, -same_codon)
    }else{stop(paste0(".datatype variabel is wrong .datatype=",.type,"\nvcf, maf is correct"))}
  }
  if(.data_type=="maf"){
    .data = .data %>>%
      filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
      mutate(alt=n_allele2)
  }
  kb_around = function(.start){
    .head=.start-1001
    .tail=.start+1000
    data.frame(start=.head:.tail)
  }
  .remove_site = .site %>>%
    dplyr::select(chr,start) %>>%
    mutate(start = purrr::map(start,~kb_around(.)),duplicate_focal = "yes") %>>%
    unnest() %>>%distinct()
  .data %>>%
    left_join(.remove_site) %>>%
    left_join(somatic_recurrent) %>>%
    filter(if(.duplicate==T){is.na(duplicate_focal)}else{chr==chr}) %>>%
    filter(if(.duplicate==T){gene_symbol!="KMT2C"}else{chr==chr}) %>>%
    filter(if(.somatic==T){is.na(recurrent_focal)}else{chr==chr})%>>%
    dplyr::select(-recurrent_focal,-duplicate_focal)
}
QC_maf_norm =quality_filter(norm_maf_all,.data_type="maf",.varscan=T)
write_df(QC_maf_norm,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/QC_norm_maf_all.tsv.gz")
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#####################################                  完成形                   #########################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
maf_trim_for_cumulative = function(.maf=norm_maf_all,.vcf=tdg_gnomad,.race="all",.fdr=0.01,.blood=F,
                                   .database="cancer",.duplicate=T,.somatic=T,.varscan=T){
  if(.race=="white"){
    .maf = .maf %>>%left_join(patient_list%>>%dplyr::select(patient_id,race),by=c("patient_id")) %>>%filter(race==.race)
    .vcf = .vcf %>>%mutate(AF=AF_white) %>>%dplyr::select(chr,start,ref,alt,AF)
  }else if(.race=="black"){
    .maf = .maf %>>%left_join(patient_list%>>%dplyr::select(patient_id,race),by=c("patient_id")) %>>%filter(race==.race)
    .vcf = .vcf %>>%mutate(AF=AF_black) %>>%dplyr::select(chr,start,ref,alt,AF) 
  }else if(.race=="all"){
    if(substitute(.vcf)=="tdg_gnomad"){
      .vcf = .vcf %>>% dplyr::select(chr,start,ref,alt,AF)
    }}else {stop(paste0(".race is wrong .race=",.race,"\nuse all, white or black!"))}
  .ref_minor=ref_minor
  if(.blood){
    .ref_minor = ref_minor_blood
    }
  ref_minor_maf = .ref_minor %>>%filter(focal=="yes")%>>%
    left_join(patient_list%>>%dplyr::select(patient_id,race),by=c("patient_id")) %>>%
    filter(if(.race !="all"){race==.race}else{chr==chr}) %>>%
    inner_join(.vcf%>>%filter(AF>0.5)%>>%mutate(MAF=1-AF),by=c("chr","start"))%>>%
    left_join(tally_norm_maf%>>%dplyr::select(-ac_cancer,-an_cancer),by=c("chr", "start", "ref", "alt"))%>>%
    left_join(.maf) %>>%filter(alt==n_allele2)%>>%
    mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
           n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
           soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
           LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
                  mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,Protein_position,FILTER) %>>%
    quality_filter(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
    filter(MAC > 0) %>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,Protein_position,FILTER)
  
  .maf %>>%
    quality_filter(.data_type = "maf",.fdr = .fdr,.database = .database,
                   .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
    left_join(.vcf) %>>%
    mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
    filter(AF < 0.5) %>>% mutate(MAF = AF) %>>%
    mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
    dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,Protein_position,FILTER) %>>%
    rbind(ref_minor_maf)%>%
    inner_join(patient_list%>>%filter(!is.na(age))%>>%dplyr::select(patient_id))
}

all_maf_for_cumulative = maf_trim_for_cumulative()
write_df(all_maf_for_cumulative,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")
white_maf_for_cumulative = maf_trim_for_cumulative(.race = "white")
write_df(white_maf_for_cumulative,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")
black_maf_for_cumulative = maf_trim_for_cumulative(.race = "black")
write_df(black_maf_for_cumulative,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/black_maf_for_cumulative.tsv.gz")

blood_maf_for_cumulative = maf_trim_for_cumulative(.maf = blood_maf_all,.blood = T)
write_df(blood_maf_for_cumulative,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_maf_for_cumulative.tsv.gz")
blood_white_maf_for_cumulative = maf_trim_for_cumulative(.maf = blood_maf_all,.blood = T,.race = "white")
write_df(blood_white_maf_for_cumulative,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_white_maf_for_cumulative.tsv.gz")
################################################################################################################################
# filter out low coverage patients
coverage=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/coverage.tsv.gz") 
patient_tdg = patient_list %>>%
  left_join(coverage%>>%dplyr::select(patient_id,all))%>>%
  filter(all > (max(all,na.rm = T)/2))%>>%
  mutate(border=mean(all)-sd(all)*2)%>>%
  filter(all>border)%>>%dplyr::select(-all,-border)
write_df(patient_tdg,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")
