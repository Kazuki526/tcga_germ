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
if(0){
#count site number
read_tsv("~/git/driver_genes/onlytop105/top_driver105cds.bed",col_names = c("chr","start","end","info","score","strand"))%>>%
  separate(info,c("gene_symbol","ensg","cds_n"),sep = ":")%>>%filter(chr!="chrX")%>>%
  left_join(driver_genes)%>>%mutate(leng=end-start+5)%>>%group_by(role)%>>%
  summarise(leng=sum(leng))%>>%(?.)%>>%{sum(.$leng)}
read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_cds.bed",col_names = c("chr","start","end","info","score","strand"))%>>%
  separate(info,c("gene_symbol","ensg","cds_n"),sep = ":")%>>%filter(chr!="chrX")%>>%
  mutate(leng=end-start+5)%>>%summarise(leng=sum(leng)) #-638*4=594,291
}
####################################### TCGA data ########################################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv",col_types = "cciciiiic")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv",col_types = "cciciiiic")
norm_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/norm_maf_all.tsv.gz")
norm_maf_all_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/norm_maf_all_cont.tsv.gz")


##################### tally norm in white #####################
if(0){
  tally_norm_maf_white = norm_maf_all%>>%
    filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
    tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
    filter(ref != alt,chr!="chrX")%>>%
    mutate(homo=ifelse(n_genotype=="homo",1,0))%>>%
    inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
    group_by(chr,start,end,ref,alt) %>>%
    summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
              Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
              cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
    ungroup() %>>%
    left_join(read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/all_patient/site_coverage_all.tsv"))%>>%
    mutate(an_cancer = an_white+an_black+an_other)
  write_df(tally_norm_maf_white,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/tally_norm_maf_white.tsv.gz")
  tally_norm_maf_white_cont = norm_maf_all_cont %>>%
    filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
    tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
    filter(ref != alt)%>>%
    mutate(homo=ifelse(n_genotype=="homo",1,0))%>>%
    inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
    group_by(chr,start,end,ref,alt) %>>%
    summarise(ac_cancer=n(),hom_cancer=sum(homo)/2,gene_symbol=first(gene_symbol),
              Consequence=first(Consequence),PolyPhen=first(PolyPhen),mutype=first(mutype),
              cDNA_position=first(cDNA_position),CDS_position=first(CDS_position)) %>>%
    ungroup() %>>%
    left_join(read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/site_coverage_all.tsv"))%>>%
    mutate(an_cancer = an_white+an_black+an_other)
  write_df(tally_norm_maf_white_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/tally_nomr_maf_white_cont.tsv.gz")
}
tally_norm_maf_white = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/tally_norm_maf_white.tsv.gz")
tally_norm_maf_white_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/tally_nomr_maf_white_cont.tsv.gz")
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

somatic_recurrent_all = 
  norm_maf_all %>>%
  filter(soma_or_germ=="somatic",LOH=="no")%>>%
  count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
  dplyr::rename(alt=t_allele2)%>>%filter(n>10) %>>%
  mutate(recurrent_focal="yes") %>>%dplyr::select(-n)%>>%
  full_join(norm_maf_all_cont%>>%
              filter(soma_or_germ=="somatic",LOH=="no")%>>%
              count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
              dplyr::rename(alt=t_allele2)%>>%filter(n>10) %>>%
              mutate(recurrent_focal="yes") %>>% dplyr::select(-n))%>>%
  dplyr::select(chr,start,ref,alt)

  

count_variant = function(.data,.duplicate=F,.somatic=F,.back=F,.lowcov=F,.count_by="variant_num"){
  #check
  if(.count_by!="variant_num" & .count_by!="site_num"){stop(".count_by must be variant_num or site_num")}
  .data %>>%
    mutate(mutation_class=paste(role,vclass,.count_by,sep="_"))%>>%
    {if(.duplicate){.%>>%anti_join(duplicate_site_all)%>>%filter(gene_symbol!="KMT2C")}else{.}}%>>%
    {if(.somatic){.%>>%anti_join(somatic_recurrent_all)}else{.}}%>>%
    {if(.back){.%>>%filter(!(LOH == "back_mutation" & n_ref_count != 0))}else{.}}%>>%
    {if(.lowcov){.%>>%filter(!is.na(coverage))}else{.}}%>>%
    {if(.count_by=="site_num"){.%>>%count(chr,start,ref,alt,mutation_class)}else{.}}%>>%
    count(mutation_class)%>>%
    tidyr::spread(key=mutation_class,value=n)
}

norm_variant =norm_maf_all %>>%
  filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
  tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
  filter(ref != alt,chr!="chrX")%>>%
  mutate(vclass=ifelse(ref=="-"|alt=="-","indel","SNV"),role="TDG")%>>%
  dplyr::select(patient_id,gene_symbol,chr,start,vclass,allele,ref,alt,n_ref_count,LOH,role)%>>%
  inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  left_join(patient_tdg%>>%mutate(coverage="enough")%>>%dplyr::select(patient_id,coverage))%>>%
  bind_rows(norm_maf_all_cont%>>%
              filter(!(soma_or_germ=="somatic" & LOH=="no")) %>>%
              tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
              filter(ref != alt,chr!="chrX")%>>%
              mutate(vclass=ifelse(ref=="-"|alt=="-","indel","SNV"),role="control")%>>%
              dplyr::select(patient_id,gene_symbol,chr,start,vclass,allele,ref,alt,n_ref_count,LOH,role)%>>%
              inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
              left_join(patient_cont%>>%mutate(coverage="enough")%>>%dplyr::select(patient_id,coverage)))

count_variant(norm_variant,.count_by = "site_num")%>>%bind_cols(count_variant(norm_variant))%>>%
  dplyr::select(TDG_SNV_site_num,TDG_SNV_variant_num,TDG_indel_site_num,TDG_indel_variant_num,
                control_SNV_site_num,control_SNV_variant_num,control_indel_site_num,control_indel_variant_num)%>>%
  bind_rows(count_variant(norm_variant,.count_by = "site_num",.duplicate = T)%>>%
              bind_cols(count_variant(norm_variant,.duplicate = T)))%>>%
  bind_rows(count_variant(norm_variant,.count_by = "site_num",.duplicate = T,.somatic = T)%>>%
              bind_cols(count_variant(norm_variant,.duplicate = T,.somatic = T)))%>>%
  bind_rows(count_variant(norm_variant,.count_by = "site_num",.duplicate = T,.somatic = T,.back = T)%>>%
              bind_cols(count_variant(norm_variant,.duplicate = T,.somatic = T,.back = T)))%>>%
  bind_rows(count_variant(norm_variant,.count_by = "site_num",.duplicate = T,.somatic = T,.back = T,.lowcov = T)%>>%
              bind_cols(count_variant(norm_variant,.duplicate = T,.somatic = T,.back = T,.lowcov = T)))%>>%
  write_df("~/Dropbox/work/rare_germ/screening_process.tsv")

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
############################################### each annotation ##############################################
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

if(0){
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
}
###################################################  only white ########################################################
white_maf_for_cumulative=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")
white_maf_for_cumulative_cont=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
#full_join(white_maf_for_cumulative,white_maf_for_cumulative_cont)%>>%
#  make_count_tbl%>>%
#  write_df("~/Dropbox/work/rare_germ/white_tcga_all_variant.tsv")
full_join(white_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id)),
          white_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id)))%>>%
  make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_tcga_age_all_variant.tsv")
full_join(tdg_gnomad,control_gnomad)%>>%mutate(MAC=ifelse(AF_white<0.5,AC_white,AN_white-AC_white))%>>%
  filter(!is.na(AC_white))%>>% make_count_tbl%>>%
  write_df("~/Dropbox/work/rare_germ/white_gnomad_all_variant.tsv")
#################only MAF<0/05%###############
#full_join(white_maf_for_cumulative,white_maf_for_cumulative_cont)%>>%
#  filter(MAF<0.0005)%>>%
#  make_count_tbl%>>%
#  write_df("~/Dropbox/work/rare_germ/white_tcga_rare_variant.tsv")
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

