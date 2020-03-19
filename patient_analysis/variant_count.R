#load classify_consequence() from prepare_tbls.R
library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd('/Volumes/DR8TB2/tcga_rare_germ')
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
control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")
#site number
read_tsv("~/git/driver_genes/onlytop105/top_driver105cds.bed",col_names = c("chr","start","end","info","score","strand"))%>>%
  separate(info,c("gene_symbol","ensg","cds_n"),sep = ":")%>>%filter(chr!="chrX")%>>%
  left_join(driver_genes)%>>%mutate(leng=end-start+5)%>>%group_by(role)%>>%
  summarise(leng=sum(leng))%>>%(?.)%>>%{sum(.$leng)}
read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_cds.bed",col_names = c("chr","start","end","info","score","strand"))%>>%
  separate(info,c("gene_symbol","ensg","cds_n"),sep = ":")%>>%filter(chr!="chrX")%>>%
  mutate(leng=end-start+5)%>>%summarise(leng=sum(leng)) #-638*4=594,291

####################################### TCGA data ########################################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
tally_norm_maf = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/tally_norm_maf.tsv.gz")
tally_norm_maf_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/tally_nomr_maf_cont.tsv.gz")

######################################## gnomAD (control cases data) ##############################################
#load classify_consequence() from prepare_tbls.R
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",col_types = cols(LoF_filter="c"))%>>%
  classify_consequence()%>>%
  filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
  mutate(AF=AC/AN,AC_white=AC_fin+AC_nfe+AC_asj,AN_white=AN_fin+AN_nfe+AN_asj,
         AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,mutype,filter,SYMBOL,AC,AN,nhomalt,AF,AC_white,AN_white,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  classify_consequence()%>>%
  filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
  mutate(AF=AC/AN,AC_white=AC_fin+AC_nfe+AC_asj,AN_white=AN_fin+AN_nfe+AN_asj,
         AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,mutype,filter,SYMBOL,AC,AN,nhomalt,AF,AC_white,AN_white,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))
#########################################################################################################
############################################### QC process ##############################################
#########################################################################################################
duplicate_site = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/duplicate_site.tsv")
duplicate_site_cont = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/duplicate_site_control_gnomad.tsv") 
kb_around = function(.start){
  .head=.start-1001
  .tail=.start+1000
  data.frame(start=.head:.tail)
}
duplicate_site_all= full_join(duplicate_site,duplicate_site_cont)%>>%
  filter(FDR_cancer<0.01)%>>%dplyr::select(chr,start)%>>%
  mutate(start = purrr::map(start,~kb_around(.))) %>>%
  unnest() %>>%distinct()
somatic_recurrent = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/norm_maf_all.tsv.gz")%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)
somatic_recurrent_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/norm_maf_all_cont.tsv.gz")%>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%
  filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%
  dplyr::select(-n)
somatic_recurrent_all = full_join(somatic_recurrent,somatic_recurrent_cont)%>>%
  dplyr::select(chr,start,ref,alt)
count_variant = function(.data,.db="tcga",.gene="tdg",.duplicate=F,.somatic=F){
  #check
  if(.db !="tcga" & .db != "gnomAD" &.db != "mix"){stop(".db must be tcga, gnomAD or mix")}
  if(.gene != "tdg" & .gene != "cont"){stop(".gene must be tdg or cont")}

  .data %>>%
    filter(get(.db)=="ok")%>>%
    {if(.duplicate){.%>>%anti_join(duplicate_site_all)}else{.}}%>>%
    {if(.somatic){.%>>%anti_join(somatic_recurrent_all)}else{.}}%>>%
    mutate(mutation_class = ifelse(ref=="-"|alt=="-",
                                   paste(.db,.gene,"indel",sep = "_"),
                                   paste(.db,.gene,"SNV",sep = "_"))) %>>%
    count(mutation_class)
}
mixtdg_tcga_gnomad = tally_norm_maf%>>%dplyr::select(chr,start,ref,alt)%>>%mutate(tcga="ok")%>>%
  full_join(tdg_gnomad%>>%dplyr::select(chr,start,ref,alt,AF,AF_white)%>>%mutate(gnomAD="ok"))%>>%
  mutate_at(vars(AF,AF_white),.funs=~ifelse(is.na(.),0,.))%>>%mutate(mix="ok")
mixcont_tcga_gnomad = tally_norm_maf_cont%>>%dplyr::select(chr,start,ref,alt)%>>%mutate(tcga="ok")%>>%
  full_join(control_gnomad%>>%dplyr::select(chr,start,ref,alt,AF,AF_white)%>>%mutate(gnomAD="ok"))%>>%
  mutate_at(vars(AF,AF_white),.funs=~ifelse(is.na(.),0,.))%>>%mutate(mix="ok")

#non QC
mixtdg_tcga_gnomad%>>%count_variant()%>>%
  full_join(mixcont_tcga_gnomad%>>%count_variant(.gene="cont"))%>>%
  tidyr::spread(mutation_class,n)%>>%mutate(step="non_QC")%>>%
#after duplicate filter
  full_join(mixtdg_tcga_gnomad%>>%count_variant(.duplicate=T)%>>%
              full_join(mixcont_tcga_gnomad%>>%count_variant(.gene="cont",.duplicate=T))%>>%
              tidyr::spread(mutation_class,n)%>>%mutate(step="duplicate_filter"))%>>%
#after duplicate filter & somatic recurrent filter
  full_join(mixtdg_tcga_gnomad%>>%count_variant(.duplicate=T,.somatic=T)%>>%
              full_join(mixcont_tcga_gnomad%>>%count_variant(.gene="cont",.duplicate=T,.somatic=T))%>>%
              tidyr::spread(mutation_class,n)%>>%mutate(step="somatic_recurrent_filter"))%>>%
#same in gnomAD
  left_join(mixtdg_tcga_gnomad%>>%count_variant(.db="gnomAD")%>>%
              full_join(mixcont_tcga_gnomad%>>%count_variant(.db="gnomAD",.gene="cont"))%>>%
              tidyr::spread(mutation_class,n)%>>%mutate(step="non_QC")%>>%
              full_join(mixtdg_tcga_gnomad%>>%count_variant(.db="gnomAD",.duplicate=T)%>>%
                          full_join(mixcont_tcga_gnomad%>>%count_variant(.db="gnomAD",.gene="cont",.duplicate=T))%>>%
                          tidyr::spread(mutation_class,n)%>>%mutate(step="duplicate_filter"))%>>%
              full_join(mixtdg_tcga_gnomad%>>%count_variant(.db="gnomAD",.duplicate=T,.somatic=T)%>>%
                          full_join(mixcont_tcga_gnomad%>>%count_variant(.db="gnomAD",.gene="cont",.duplicate=T,.somatic=T))%>>%
                          tidyr::spread(mutation_class,n)%>>%mutate(step="somatic_recurrent_filter")))%>>%
  left_join(mixtdg_tcga_gnomad%>>%count_variant(.db="mix")%>>%
              full_join(mixcont_tcga_gnomad%>>%count_variant(.db="mix",.gene="cont"))%>>%
              tidyr::spread(mutation_class,n)%>>%mutate(step="non_QC")%>>%
              full_join(mixtdg_tcga_gnomad%>>%count_variant(.db="mix",.duplicate=T)%>>%
                          full_join(mixcont_tcga_gnomad%>>%count_variant(.db="mix",.gene="cont",.duplicate=T))%>>%
                          tidyr::spread(mutation_class,n)%>>%mutate(step="duplicate_filter"))%>>%
              full_join(mixtdg_tcga_gnomad%>>%count_variant(.db="mix",.duplicate=T,.somatic=T)%>>%
                          full_join(mixcont_tcga_gnomad%>>%count_variant(.db="mix",.gene="cont",.duplicate=T,.somatic=T))%>>%
                          tidyr::spread(mutation_class,n)%>>%mutate(step="somatic_recurrent_filter")))%>>%
  dplyr::select(step,tcga_tdg_SNV,tcga_tdg_indel,tcga_cont_SNV,tcga_cont_indel,
                gnomAD_tdg_SNV,gnomAD_tdg_indel,gnomAD_cont_SNV,gnomAD_cont_indel,
                mix_tdg_SNV,mix_tdg_indel,mix_cont_SNV,mix_cont_indel)%>>%
  write_df("~/Dropbox/work/rare_germ/screening_process.tsv")


############################################### each annotation ##############################################
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv")
all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")
all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")
make_count_tbl=function(.tbl){
  .tbl %>>%
    group_by(gene_symbol,chr,start,ref,alt,mutype)%>>%summarise(MAC=sum(MAC))%>>%ungroup()%>>%
    left_join(driver_genes)%>>%
    mutate(role=ifelse(is.na(role),"control",role),
           mutype = ifelse(ref=="-"|alt=="-",paste0("indel_",mutype),paste0("SNV_",mutype)))%>>%
    group_by(role,mutype)%>>%
    summarise(lnum=n(),vnum=sum(MAC)) %>>% #location_num & variant_num
    tidyr::pivot_wider(names_from=role,values_from=c(lnum,vnum))%>>%
    mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
    mutate(lnum_TDG=lnum_oncogene+lnum_TSG+`lnum_oncogene/TSG`,vnum_TDG=vnum_oncogene+vnum_TSG+`vnum_oncogene/TSG`)%>>%
    left_join(tibble(mutype=c("SNV_truncating","SNV_splice","SNV_missense","SNV_silent","indel_truncating","indel_splice","indel_inframe_indel"),
                     order =c(1,2,3,4,5,6,7)))%>>%arrange(order)%>>%
    dplyr::select(mutype,lnum_TSG,vnum_TSG,lnum_oncogene,vnum_oncogene,`lnum_oncogene/TSG`,`vnum_oncogene/TSG`,
                  lnum_TDG,vnum_TDG,lnum_control,vnum_control)
}


full_join(all_maf_for_cumulative,all_maf_for_cumulative_cont)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/tcga_all_variant.tsv")
full_join(all_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id)),
          all_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id)))%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/tcga_age_all_variant.tsv")
full_join(tdg_gnomad,control_gnomad)%>>%mutate(MAC=ifelse(AF<0.5,AC,AN-AC))%>>%
  filter(!is.na(AC))%>>% make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/gnomad_all_variant.tsv")
#################only MAF<0/05%###############
full_join(all_maf_for_cumulative,all_maf_for_cumulative_cont)%>>%
  filter(MAF<0.0005)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/tcga_rare_variant.tsv")
full_join(all_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id)),
          all_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id)))%>>%
  filter(MAF<0.0005)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/tcga_age_rare_variant.tsv")
full_join(tdg_gnomad,control_gnomad)%>>%
  mutate(MAC=ifelse(AF<0.5,AC,AN-AC),MAF=ifelse(AF<0.5,AF,1-AF))%>>%filter(!is.na(AC))%>>%
  filter(MAF<0.0005)%>>% make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/gnomad_rare_variant.tsv")
###################################################  only white ########################################################
white_maf_for_cumulative=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")
white_maf_for_cumulative_cont=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
full_join(white_maf_for_cumulative,white_maf_for_cumulative_cont)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_tcga_all_variant.tsv")
full_join(white_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id)),
          white_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id)))%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_tcga_age_all_variant.tsv")
full_join(tdg_gnomad,control_gnomad)%>>%mutate(MAC=ifelse(AF_white<0.5,AC_white,AN_white-AC_white))%>>%
  filter(!is.na(AC_white))%>>% make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_gnomad_all_variant.tsv")
#################only MAF<0/05%###############
full_join(white_maf_for_cumulative,white_maf_for_cumulative_cont)%>>%
  filter(MAF<0.0005)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_tcga_rare_variant.tsv")
full_join(white_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id)),
          white_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id)))%>>%
  filter(MAF<0.0005)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_tcga_age_rare_variant.tsv")
full_join(tdg_gnomad,control_gnomad)%>>%
  mutate(MAC=ifelse(AF_white<0.5,AC_white,AN_white-AC_white),MAF=ifelse(AF_white<0.5,AF_white,1-AF_white))%>>%
  filter(!is.na(AC_white),MAF<0.0005)%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_gnomad_rare_variant.tsv")

