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
options(scipen=100)

driver_genes=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/driver_genes.tsv")
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
patient_hicov = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list_exclude_low_coverage.tsv",col_types = "cciciiiic")
################ read MAF extracted ################
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
black_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/black_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",
                      col_types=cols(chr="c",LoF_filter="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
######################################################################################################################
####################################################### white ########################################################
######################################################################################################################
############ TSG ##############
#患者ごとのtruncating な遺伝子の数
truncating_count = white_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  #0個の患者も入れる
  right_join(patient_hicov)%>>%filter(race=="white",!is.na(age))%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all_white.tsv")
ggsave("age_plot/cumulative/tsg/white/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT_white.tsv")
ggsave("age_plot/cumulative/tsg/white/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.0006985
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
  cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                      strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                     x=0.37,y=0.07,width=0.63,height=0.93)+
  cowplot::draw_text("Number of truncating and splice variants in TSG",x=0.5,y=0.05,size=20)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
ggsave("age_plot/fig/truncate/white/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
if(1){
truncating_count_rare = white_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  #0個の患者も入れる
  right_join(patient_hicov)%>>%filter(race=="white",!is.na(age))%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot_all = truncating_count_rare%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_rare_white.tsv") 
ggsave("age_plot/cumulative/tsg/white/rare_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_rare%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT_white.tsv")
ggsave("age_plot/cumulative/tsg/white/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.0003727
t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")
if(0){
  regression_table = make_regression_tabel_truncate(.maf=white_maf_for_cumulative,.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_trunc_white.tsv")
}
regression_table = read_tsv("age_plot/cumulative/regression/TSG_trunc_white.tsv")
.plot_reglog = regression_plot_log(regression_table,.black=0.05)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.38,width=0.36,height=0.62)+
  cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                         strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/truncate/white/TSG_rare.pdf",.plot,width = 14,height = 8)
}

#################### oncogene ###########################
truncating_count_onco_rare = white_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS",MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_hicov)%>>%filter(race=="white",!is.na(age))%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_white.tsv")
ggsave("age_plot/cumulative/oncogene/white/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT_white.tsv")
ggsave("age_plot/cumulative/oncogene/white/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value= 0.1976
t.test(truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")
if(0){
  regression_table = make_regression_tabel_truncate(.maf=white_maf_for_cumulative,.max_maf = 50,.race="white",.role = "oncogene")
  write_df(regression_table,"age_plot/cumulative/regression/oncogene_trunc_white.tsv")
}
regression_table = read_tsv("age_plot/cumulative/regression/oncogene_trunc_white.tsv")
.plot_reglog = regression_plot_log(regression_table,.black=0.05)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.38,width=0.36,height=0.62)+
  cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                      strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/truncate/white/oncogene_rare.pdf",.plot,width = 14,height = 8)


if(0){
  truncating_count_onco = white_maf_for_cumulative %>>%
    filter(chr!="chrX",FILTER=="PASS")%>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
    right_join(patient_hicov)%>>%filter(race=="white",!is.na(age))%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate_white.tsv")
  ggsave("age_plot/cumulative/oncogene/white/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT_white.tsv")
  ggsave("age_plot/cumulative/oncogene/white/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.1976
  t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
         truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::ggdraw()+
    cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
    cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                        strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                       x=0.37,y=0.07,width=0.63,height=0.93)+
    cowplot::draw_text("Number of truncating and splice variants in oncogene",x=0.5,y=0.05,size=20)+
    cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
    cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
  ggsave("age_plot/fig/truncate/white/alloncogene.pdf",.plot,width = 10,height = 10)
}
if(0){
######################################################################################################################
####################################################### black ########################################################患者が少なすぎる
######################################################################################################################
############ TSG ##############
#患者ごとのtruncating な遺伝子の数
truncating_count = black_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  #0個の患者も入れる
  right_join(patient_hicov)%>>%filter(race=="black",!is.na(age))%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all_black.tsv")
ggsave("age_plot/cumulative/tsg/black/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT_black.tsv")
ggsave("age_plot/cumulative/tsg/black/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.1544
t.test(truncating_count[truncating_count$truncating_count_n>0,]$age/365.25,
       truncating_count[truncating_count$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
  cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                      strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                     x=0.37,y=0.07,width=0.63,height=0.93)+
  cowplot::draw_text("Number of truncating and splice variants in TSG",x=0.5,y=0.05,size=20)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
ggsave("age_plot/fig/truncate/black/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
if(1){
  truncating_count_rare = black_maf_for_cumulative %>>%
    filter(chr!="chrX",FILTER=="PASS")%>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="TSG"|role=="oncogene/TSG")%>>%
    filter(MAF < 0.0005)%>>%
    #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
    #0個の患者も入れる
    right_join(patient_hicov)%>>%filter(race=="black",!is.na(age))%>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
           age=round(age/365.25*100)/100)
  
  .plot_all = truncating_count_rare%>>%
    truncate_plot_allcantype(.permu_file = "TSG/truncate_rare_black.tsv") 
  ggsave("age_plot/cumulative/tsg/black/rare_truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_rare%>>%
    truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT_black.tsv")
  ggsave("age_plot/cumulative/tsg/black/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
  #あるなしのt.test p-value=0.0677
  t.test(truncating_count_rare[truncating_count_rare$truncating_count_n>0,]$age/365.25,
         truncating_count_rare[truncating_count_rare$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::ggdraw()+
    cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
    cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                        strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                       x=0.37,y=0.07,width=0.63,height=0.93)+
    cowplot::draw_text("Number of rare truncating and splice variants in TSG",x=0.5,y=0.05,size=20)+
    cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
    cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
  ggsave("age_plot/fig/truncate/black/TSG_rare.pdf",.plot,width = 14,height = 8)
}

###################### oncogene ##############################
truncating_count_onco_rare = black_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS",MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_hicov)%>>%filter(race=="black",!is.na(age))%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_black.tsv")
ggsave("age_plot/cumulative/oncogene/black/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT_black.tsv")
ggsave("age_plot/cumulative/oncogene/black/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.04298
t.test(truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n>0,]$age/365.25,
       truncating_count_onco_rare[truncating_count_onco_rare$truncating_count_n==0,]$age/365.25,alternative="less")

.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
  cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                      strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                     x=0.37,y=0.07,width=0.63,height=0.93)+
  cowplot::draw_text("Number of rare truncating and splice variants in oncogene",x=0.5,y=0.05,size=20)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
ggsave("age_plot/fig/truncate/black/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(1){
  truncating_count_onco = black_maf_for_cumulative %>>%
    filter(chr!="chrX",FILTER=="PASS")%>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
    right_join(patient_hicov)%>>%filter(race=="black",!is.na(age))%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate_black.tsv")
  ggsave("age_plot/cumulative/oncogene/black/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT_black.tsv")
  ggsave("age_plot/cumulative/oncogene/black/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.04298
  t.test(truncating_count_onco[truncating_count_onco$truncating_count_n>0,]$age/365.25,
         truncating_count_onco[truncating_count_onco$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::ggdraw()+
    cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
    cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                        strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                       x=0.37,y=0.07,width=0.63,height=0.93)+
    cowplot::draw_text("Number of truncating and splice variants in oncogene",x=0.5,y=0.05,size=20)+
    cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
    cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
  ggsave("age_plot/fig/truncate/black/alloncogene.pdf",.plot,width = 10,height = 10)
}}

