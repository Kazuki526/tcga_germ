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
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
patien_all_info = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/all_patient_info.tsv",col_types = "ccccciidci")
norm_maf_all_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/norm_maf_all_cont.tsv.gz")
tally_norm_maf_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/tally_nomr_maf_cont.tsv.gz")
ref_minor_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/ref_minor_coverage.tsv.gz")


######################################## gnomAD (control cases data) ##############################################
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))

################################################################################################################
################################################ quality filter ################################################
################################################################################################################
### HWE test ###
if(0){
  options(scipen = 10) #startなどを指数表記のまま保存しないため
  HWE_test_heterom = function(ac,an,hom){
    AF=ac/an
    .matrix = matrix(c(round(an*AF*(1-AF)),round(an/2*AF^2),ac-hom*2,hom), nrow = 2)
    fisher.test(.matrix,alternative = "less")$p.value
  }
  duplicate_site_cont =control_gnomad %>>%
    filter(alt!="-")%>>%
    dplyr::select(gene_symbol,chr,start,ref,alt,AC,AN,nhomalt) %>>%
    dplyr::rename(ac_gnomad=AC,an_gnomad=AN,hom_gnomad=nhomalt) %>>%
    mutate(hom_gnomad_hwe=(ac_gnomad^2/an_gnomad)/2) %>>%
    filter(ac_gnomad!=0,!is.na(hom_gnomad)) %>>%
    nest(-chr,-start,-ref,-alt) %>>%
    mutate(HWE_gnomad = purrr::map(data,~HWE_test_heterom(.$ac_gnomad,.$an_gnomad,.$hom_gnomad))) %>>%
    unnest() %>>%
    mutate(FDR_gnomad=p.adjust(HWE_gnomad)) %>>%filter(FDR_gnomad<0.01) %>>%
    full_join(tally_norm_maf_cont%>>%
                filter(!is.na(an_cancer),alt!="-")%>>%
                dplyr::select(chr,start,ref,alt,ac_cancer,hom_cancer,gene_symbol,mutype,an_cancer) %>>%
                mutate(hom_cancer_hwe = ac_cancer^2/an_cancer/2) %>>%
                nest(-chr,-start,-ref,-alt) %>>%
                mutate(HWE_cancer=purrr::map(data,~HWE_test_heterom(.$ac_cancer,.$an_cancer,.$hom_cancer))) %>>%
                unnest() %>>%
                mutate(FDR_cancer=p.adjust(HWE_cancer)) )%>>%
    filter(FDR_cancer<0.01|FDR_gnomad<0.01)
  write_df(duplicate_site_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/duplicate_site.tsv")
}
duplicate_site_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/duplicate_site.tsv") 
##### somaticでrecurrentなmutationはgermで起こっているとは考えにくい（様々なエラーが考えられる）
somatic_recurrent_cont = norm_maf_all_cont%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)

#normalでaltalt, tumorでrefaltとなってる際にnormalでrefのdepth=0のものだけ採用！
#また同じサイトでこのエラーが有意に多い(100patient以上の)siteは解析に使用しないことにした。(3site)
#しかし control_regionでは存在しない、、
varscan_error_cont = norm_maf_all_cont %>>%
  filter(LOH == "back_mutation")%>>%#, n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  group_by(chr,start,end,ref,alt) %>>%
  summarise(an_error = n()*2) %>>%
  ungroup() 
varscan_error_site_cont = norm_maf_all_cont %>>%
  filter(LOH == "back_mutation", n_ref_count != 0) %>>%
  dplyr::rename(alt = n_allele2) %>>%
  dplyr::select(patient_id,chr,start,ref,alt) %>>%
  mutate(varscan_error_focal="yes")


quality_filter_cont =
  function(.data,.data_type="vcf",.fdr=0.01,.database="cancer",.duplicate=T,.somatic=T,.varscan=T){
    .site = duplicate_site_cont
    if(.database!="all"){
      if(.database=="cancer"){
        .site = .site %>>% filter(FDR_cancer < .fdr)
      }else if(.database=="gnomAD"){
        .site = .site %>>% filter(FDR_gnomad < .fdr)
      }else{stop(paste0("database variable is wrong .database=",.database,
                        "\ncancer, gnomAD is correct."))}
    }
    if(.varscan){
      if(.data_type=="vcf"){
        .data = .data %>>%
          left_join(varscan_error_cont) %>>%
          mutate(ac_cancer = ifelse(!is.na(an_cancer),ac_cancer-an_error,ac_cancer),
                 an_cancer = ifelse(!is.na(an_cancer),an_cancer-an_error,an_cancer) - an_male_cancer) %>>%
          dplyr::select(-an_male_cancer)
      }else if(.data_type=="maf"){
        .data = .data %>>%
          left_join(varscan_error_site_cont) %>>%
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
    kb_sagittal = function(.start){
      .head=.start-1001
      .tail=.start+1000
      data.frame(start=.head:.tail)
    }
    .remove_site = .site %>>%
      dplyr::select(chr,start) %>>%
      mutate(start = purrr::map(start,~kb_sagittal(.)),duplicate_focal = "yes") %>>%
      unnest() %>>%distinct()
    .data = .data %>>%
      left_join(.remove_site) %>>%
      left_join(somatic_recurrent_cont) %>>%
      filter(if(.duplicate==T){is.na(duplicate_focal)}else{chr==chr}) %>>%
      filter(if(.somatic==T){is.na(recurrent_focal)}else{chr==chr})%>>%
      dplyr::select(-recurrent_focal,-duplicate_focal)
  }
QC_maf_norm_cont=quality_filter_cont(norm_maf_all_cont,.data_type = "maf")
write_df(QC_maf_norm_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/QC_norm_maf_all_cont.tsv.gz")
######################################################################################################
################################### germline mutation list ###########################################
######################################################################################################
if(0){
  maf_trim_for_cumulative_cont = function(.maf=norm_maf_all_cont,.vcf=control_gnomad,.race="all",.fdr=0.01,
                                          .database="cancer",.duplicate=T,.somatic=T,.varscan=T){
    if(.race=="white"){
      .maf = .maf %>>%left_join(patient_list%>>%dplyr::select(patient_id,race)) %>>%filter(race==.race)
      .vcf = .vcf %>>%mutate(AF=AF_white) %>>%
        dplyr::select(chr,start,ref,alt,AF)
    }else if(.race=="black"){
      .maf = .maf %>>%left_join(patient_list%>>%dplyr::select(patient_id,race)) %>>%filter(race==.race)
      .vcf = .vcf %>>%mutate(AF=AF_black) %>>%
        dplyr::select(chr,start,ref,alt,AF) 
    }else if(.race=="all"){
      .vcf = .vcf %>>%
        dplyr::select(chr,start,ref,alt,AF)
      }else {stop(paste0(".race is wrong .race=",.race,"\nuse all, white or black!"))}
    ref_minor = ref_minor_cont %>>%
      filter(focal == "yes") %>>%dplyr::select(-focal)%>>%
      left_join(patient_list%>>%dplyr::select(patient_id,race)) %>>%
      filter(if(.race !="all"){race==.race}else{chr==chr}) %>>%
      left_join(left_join(tally_norm_maf_cont,.vcf)%>>%
                  dplyr::select(-ac_cancer,-hom_cancer)) %>>%
      filter(AF > 0.5) %>>% mutate(MAF=1-AF) %>>%
      left_join(.maf) %>>%
      mutate(n_allele1 = ifelse(is.na(n_allele1),ref,n_allele1),
             n_allele2 = ifelse(is.na(n_allele2),ref,n_allele2),
             soma_or_germ =ifelse(is.na(soma_or_germ),"ref",soma_or_germ),
             LOH = ifelse(is.na(LOH),"ref",LOH))%>>%
      dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,
                    mutype,AF,MAF,n_allele1,n_allele2,soma_or_germ,LOH,Protein_position,FILTER) %>>%
      quality_filter_cont(.data_type = "maf",.fdr = .fdr,.database = .database,
                          .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
      mutate(MAC=ifelse(n_allele1 == alt,0,ifelse(n_allele2==ref,2,1))) %>>%
      filter(MAC > 0)%>>%
      dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,Protein_position,FILTER)
    .maf %>>%
      quality_filter_cont(.data_type = "maf",.fdr = .fdr,.database = .database,
                          .duplicate = .duplicate,.somatic = .somatic,.varscan = .varscan) %>>%
      left_join(.vcf) %>>%
      mutate(AF=ifelse(is.na(AF),0,AF)) %>>%
      filter(AF <= 0.5) %>>% mutate(MAF = AF) %>>%
      mutate(MAC = ifelse(n_allele1==alt,2,ifelse(n_allele2==alt,1,0))) %>>%
      dplyr::select(patient_id,cancer_type,gene_symbol,chr,start,end,ref,alt,mutype,AF,MAF,MAC,Protein_position,FILTER) %>>%
      rbind(ref_minor)
  }
all_maf_for_cumulative_cont = maf_trim_for_cumulative_cont()
write_df(all_maf_for_cumulative_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")
white_maf_for_cumulative_cont = maf_trim_for_cumulative_cont(.race = "white")
write_df(white_maf_for_cumulative_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
black_maf_for_cumulative_cont = maf_trim_for_cumulative_cont(.race = "black")
write_df(black_maf_for_cumulative_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/black_maf_for_cumulative_control.tsv.gz")
}
#all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")
#white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
#black_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/black_maf_for_cumulative_control.tsv.gz")
################################################################################################################################
# filter out low coverage patients
coverage_cont=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/coverage_cont.tsv.gz") 
patient_cont = patient_list %>>%
  left_join(coverage_cont%>>%dplyr::select(patient_id,coverage))%>>%
  filter(coverage>max(coverage)/2)%>>%
  mutate(border=mean(coverage)-sd(coverage)*2)%>>%
  filter(coverage>border)%>>%dplyr::select(-coverage,-border)
write_df(patient_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv")
