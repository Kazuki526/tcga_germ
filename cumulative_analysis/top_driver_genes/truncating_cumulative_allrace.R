library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
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
options(scipen=1)

driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%dplyr::rename(gene_symbol=gene)%>>%
  mutate(role=ifelse(is.na(role),"TSG/oncogene",role))
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_race = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_race.tsv")
################ read MAF extracted ################
all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")

####################################################################################################################
####################################################### TSG ########################################################
####################################################################################################################

#患者ごとのtruncating な遺伝子の数
truncating_count = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%#filter(truncating_count_n<5)%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all.tsv")
ggsave("age_plot/cumulative/tsg/all_race/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT.tsv")
ggsave("age_plot/cumulative/tsg/all_race/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.0307
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/all_race/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
truncating_count_rare = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot_all = truncating_count_rare%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_rare.tsv") 
ggsave("age_plot/cumulative/tsg/all_race/rare_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_rare%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/tsg/all_race/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.001986
t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/all_race/TSG_rare.pdf",.plot,width = 14,height = 8)

if(0){
  ##### rare BRCA1,2を別にしてみたら？####
  brca_truncate_num = patient_list %>>%
    left_join(all_maf_for_cumulative %>>%
                filter((gene_symbol=="BRCA1"|gene_symbol=="BRCA2")&MAF<0.0005) %>>%  ####ここを変えてBRCA1,2のみか以外か変える
                filter((mutype=="truncating"|mutype=="splice")) %>>%
                count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
                group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
           age=round(age/365.25*100)/100)
  
  .plot_all = brca_truncate_num %>>%truncate_plot_allcantype()
  ggsave("age_plot/cumulative/all_race/brca_truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = brca_truncate_num %>>%truncate_plot_bycantype()
  ggsave("age_plot/cumulative/all_race/brca_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
  
  .plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                             .plot_by + theme(axis.title.y = element_blank()),
                             labels = "auto",rel_widths = c(0.5,1))
  ggsave("age_plot/fig/truncate/all_race/rare_brca12.pdf",.plot,width = 10,height = 10)
  
  ##filterをいじってBRCA1,2以外では？
  brca_truncate_num_ = patient_list %>>%
    left_join(all_maf_for_cumulative %>>%
                filter(!((gene_symbol=="BRCA1"|gene_symbol=="BRCA2")&MAF<0.0005)) %>>%  ####ここを変えてBRCA1,2のみか以外か変える
                left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
                filter(mutype=="truncating"|mutype=="splice") %>>%
                filter(role=="TSG"|role=="oncogene/TSG")%>>%
                count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
                group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
           age=round(age/365.25*100)/100)
  .plot = brca_truncate_num_ %>>%truncate_plot_allcantype()
  ggsave("age_plot/cumulative/brca_non_truncating.pdf",.plot,height = 5,width = 5)
  .plot = brca_truncate_num_ %>>%truncate_plot_bycantype()
  ggsave("age_plot/cumulative/brca_non_truncating_by_cancerype.pdf",.plot,height = 10,width = 10)
  
  
  ###burden test 同様singleton, doubleton除去して見たら？
  truncating_focal_site = quality_filter(norm_maf_all,.data_type="maf") %>>%
    filter(mutype == "splice"| mutype =="truncating") %>>%
    dplyr::select(-alt) %>>%
    tidyr::gather(allele,alt,n_allele1,n_allele2) %>>%
    filter(ref != alt)%>>%
    group_by(chr,start,end,ref,alt) %>>%
    summarise(ac_cancer=n(),gene_symbol=first(gene_symbol),
              Consequence=first(Consequence),mutype=first(mutype)) %>>%
    ungroup() %>>%
    left_join(coverage_all) %>>%
    full_join(quality_filter(exac)) %>>% filter(mutype == "splice"| mutype =="truncating") %>>%
    filter((is.na(AC_Adj)|AC_Adj/AN_Adj *100 <0.05), sum(AC_Adj,ac_cancer,na.rm = T) >2, chr!="chrX") %>>%
    mutate(focal = "ok") %>>%
    dplyr::select(gene_symbol,chr,start,end,ref,alt,focal)
  truncating_count_rare_cut_single = all_maf_for_cumulative %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="TSG"|role=="oncogene/TSG")%>>%
    filter(MAF < 0.0005)%>>%
    left_join(truncating_focal_site) %>>%
    filter(!is.na(focal))%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    #0個の患者も入れる
    {left_join(patient_list,.)}%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
           age=round(age/365.25*100)/100)
  truncate_plot_allcantype(truncating_count_rare_cut_single)
}
#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
###### truncating ########
truncating_count_onco_rare = all_maf_for_cumulative %>>%
  filter(MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare.tsv")
ggsave("age_plot/cumulative/oncogene/all_race/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/oncogene/all_race/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.2996
t.test(truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/all_race/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(0){
  truncating_count_onco = all_maf_for_cumulative %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(patient_list,.)}%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate.tsv")
  ggsave("age_plot/cumulative/oncogene/all_race/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT.tsv")
  ggsave("age_plot/cumulative/oncogene/all_race/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.2996
  t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
         truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                             .plot_by + theme(axis.title.y = element_blank()),
                             labels = "auto",rel_widths = c(0.5,1))
  ggsave("age_plot/fig/truncate/all_race/alloncogene.pdf",.plot,width = 10,height = 10)
}






















