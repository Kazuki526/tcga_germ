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
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
black_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/black_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")

######################################################################################################################
####################################################### white ########################################################
######################################################################################################################
############ TSG ##############
#患者ごとのtruncating な遺伝子の数
truncating_count = white_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  left_join(patient_race)%>>%filter(race=="white")%>>%dplyr::select(-race)%>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all_white.tsv")
ggsave("age_plot/cumulative/tsg/white/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT_white.tsv")
ggsave("age_plot/cumulative/tsg/white/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.01767
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/white/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
if(1){
truncating_count_rare = all_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  left_join(patient_race)%>>%filter(race=="white")%>>%dplyr::select(-race)%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot_all = truncating_count_rare%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_rare_white.tsv") 
ggsave("age_plot/cumulative/tsg/white/rare_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_rare%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT_white.tsv")
ggsave("age_plot/cumulative/tsg/white/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.0007569
t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/white/TSG_rare.pdf",.plot,width = 14,height = 8)
}

#################### oncogene ###########################
truncating_count_onco_rare = all_maf_for_cumulative %>>%
  filter(MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  left_join(patient_race)%>>%filter(race=="white")%>>%dplyr::select(-race)%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_white.tsv")
ggsave("age_plot/cumulative/oncogene/white/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT_white.tsv")
ggsave("age_plot/cumulative/oncogene/white/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value= 0.263
t.test(truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/white/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(1){
  truncating_count_onco = all_maf_for_cumulative %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(patient_list,.)}%>>%
    left_join(patient_race)%>>%filter(race=="white")%>>%dplyr::select(-race)%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate_white.tsv")
  ggsave("age_plot/cumulative/oncogene/white/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT_white.tsv")
  ggsave("age_plot/cumulative/oncogene/white/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.263
  t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
         truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                             .plot_by + theme(axis.title.y = element_blank()),
                             labels = "auto",rel_widths = c(0.5,1))
  ggsave("age_plot/fig/truncate/white/alloncogene.pdf",.plot,width = 10,height = 10)
}

######################################################################################################################
####################################################### black ########################################################患者が少なすぎる
######################################################################################################################
############ TSG ##############
#患者ごとのtruncating な遺伝子の数
truncating_count = black_maf_for_cumulative %>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  #0個の患者も入れる
  {left_join(patient_list,.)}%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  left_join(patient_race)%>>%filter(race=="black")%>>%dplyr::select(-race)%>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all_black.tsv")
ggsave("age_plot/cumulative/tsg/black/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT_black.tsv")
ggsave("age_plot/cumulative/tsg/black/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.4451
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all,
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
ggsave("age_plot/fig/truncate/black/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
if(1){
  truncating_count_rare = all_maf_for_cumulative %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="TSG"|role=="oncogene/TSG")%>>%
    filter(MAF < 0.0005)%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    #0個の患者も入れる
    {left_join(patient_list,.)}%>>%
    left_join(patient_race)%>>%filter(race=="black")%>>%dplyr::select(-race)%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
           age=round(age/365.25*100)/100)
  
  .plot_all = truncating_count_rare%>>%
    truncate_plot_allcantype(.permu_file = "TSG/truncate_rare_black.tsv") 
  ggsave("age_plot/cumulative/tsg/black/rare_truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_rare%>>%
    truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT_black.tsv")
  ggsave("age_plot/cumulative/tsg/black/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
  #あるなしのt.test p-value=0.1803
  t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
         truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::plot_grid(.plot_all,
                             .plot_by + theme(axis.title.y = element_blank()),
                             labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                             rel_widths = c(1,1.8))
  ggsave("age_plot/fig/truncate/black/TSG_rare.pdf",.plot,width = 14,height = 8)
}

###################### oncogene ##############################
truncating_count_onco_rare = all_maf_for_cumulative %>>%
  filter(MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
  {left_join(patient_list,.)}%>>%
  left_join(patient_race)%>>%filter(race=="black")%>>%dplyr::select(-race)%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_black.tsv")
ggsave("age_plot/cumulative/oncogene/black/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT_black.tsv")
ggsave("age_plot/cumulative/oncogene/black/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.0359
t.test(truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                           .plot_by + theme(axis.title.y = element_blank()),
                           labels = "auto",rel_widths = c(0.5,1))
ggsave("age_plot/fig/truncate/black/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(1){
  truncating_count_onco = all_maf_for_cumulative %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%
    {left_join(patient_list,.)}%>>%
    left_join(patient_race)%>>%filter(race=="black")%>>%dplyr::select(-race)%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate_black.tsv")
  ggsave("age_plot/cumulative/oncogene/black/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT_black.tsv")
  ggsave("age_plot/cumulative/oncogene/black/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.0359
  t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
         truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::plot_grid(.plot_all+ggtitle("all cancer type"),
                             .plot_by + theme(axis.title.y = element_blank()),
                             labels = "auto",rel_widths = c(0.5,1))
  ggsave("age_plot/fig/truncate/black/alloncogene.pdf",.plot,width = 10,height = 10)
}