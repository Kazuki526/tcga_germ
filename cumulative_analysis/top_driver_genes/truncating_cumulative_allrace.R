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
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")
################ read MAF extracted ################
all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")

####################################################################################################################
####################################################### TSG ########################################################
####################################################################################################################

#患者ごとのtruncating な遺伝子の数
truncating_count = all_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  #0個の患者も入れる
  right_join(patient_tdg)%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n)) %>>%
  mutate(age=round(age/365.25*100)/100)


.plot_all = truncating_count%>>%#filter(truncating_count_n<5)%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_all.tsv")
ggsave("age_plot/cumulative/tsg/all_race/truncating_all.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_all_byCT.tsv")
ggsave("age_plot/cumulative/tsg/all_race/truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#TSGのtruncateあるなしのt_testでは？？p_value=0.01997
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
ggsave("age_plot/fig/truncate/all_race/allTSG.pdf",.plot,width = 14,height = 8)

############### MAF<0.05%のtruncating mutationのみでやってみたら？
truncating_count_rare = all_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="TSG"|role=="oncogene/TSG")%>>%
  filter(MAF < 0.0005)%>>%
  count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  #0個の患者も入れる
  right_join(patient_tdg)%>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n),
         age=round(age/365.25*100)/100)

.plot_all = truncating_count_rare%>>%
  truncate_plot_allcantype(.permu_file = "TSG/truncate_rare.tsv") 
ggsave("age_plot/cumulative/tsg/all_race/rare_truncating.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_rare%>>%
  truncate_plot_bycantype(.permu_file = "TSG/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/tsg/all_race/raer_truncating_by_cancerype.pdf",.plot_by,height = 10,width = 10)
#あるなしのt.test p-value=0.003491
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
ggsave("age_plot/fig/truncate/all_race/TSG_rare.pdf",.plot,width = 14,height = 8)

if(0){
#BRCA1,2を除外したら？(BRCA1:76(12.8%),BRCA2:260(43.8%), all TSG:594)
#BRCA1,2 truncating or spliceを持つ患者は
  brca_paid = all_maf_for_cumulative %>>%
    filter(chr!="chrX",FILTER=="PASS")%>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="TSG"|role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>%
    filter(gene_symbol=="BRCA1"|gene_symbol=="BRCA2")%>>%dplyr::select(cancer_type,patient_id)
  .plot_all = truncating_count%>>%anti_join(brca_paid)%>%
    truncate_plot_allcantype(.permu_file = "TSG/truncate_notbrca.tsv")
  ggsave("age_plot/cumulative/tsg/all_race/truncating_notbrca.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count%>>%anti_join(brca_paid)%>%
    truncate_plot_bycantype(.permu_file = "TSG/truncate_notbrca_byCT.tsv")
  ggsave("age_plot/cumulative/tsg/all_race/truncating_notbrca_by_cancerype.pdf",.plot_by,height = 10,width = 10)
  #TSGのtruncateあるなしのt_testでは？？p_value=0.2227
  truncating_count_nbr=truncating_count%>>%anti_join(brca_paid)
  t.test(truncating_count_nbr[truncating_count_nbr$truncating_count_n>0,]$age/365.25,
         truncating_count_nbr[truncating_count_nbr$truncating_count_n==0,]$age/365.25,alternative="less")
  
  .plot = cowplot::ggdraw()+
    cowplot::draw_plot(.plot_all+theme(axis.title.x = element_blank()),x=0,y=0.07,width=0.36,height=0.93)+
    cowplot::draw_plot(.plot_by + theme(axis.title = element_blank(),axis.text.y = element_text(size=10),
                                        strip.text = element_text(size=8,margin=margin(1,0,1,0))),
                       x=0.37,y=0.07,width=0.63,height=0.93)+
    cowplot::draw_text("Number of truncating and splice variants",x=0.5,y=0.05,size=20)+
    cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
    cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)
  ggsave("age_plot/fig/truncate/all_race/allTSG_notbrca.pdf",.plot,width = 14,height = 8)
}


#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
###### truncating ########
truncating_count_onco_rare = all_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS")%>>%
  filter(MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_tdg)%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

.plot_all = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare.tsv")
ggsave("age_plot/cumulative/oncogene/all_race/truncating_rare.pdf",.plot_all,height = 5,width = 5)
.plot_by = truncating_count_onco_rare %>>%
  truncate_plot_bycantype(.permu_file = "oncogene/truncate_rare_byCT.tsv")
ggsave("age_plot/cumulative/oncogene/all_race/truncating_rare_byCT.pdf",.plot_by,height = 10,width = 10)
#oncogeneのtruncateあるなしのt_testでは？？p-value=0.2958
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
ggsave("age_plot/fig/truncate/all_race/oncogene_rare.pdf",.plot,width = 10,height = 10)

if(1){
  truncating_count_onco = all_maf_for_cumulative %>>%
    filter(chr!="chrX",FILTER=="PASS")%>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
    filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(role=="oncogene"|role=="oncogene/TSG")%>>%
    count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
    group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
    right_join(patient_tdg)%>>%
    mutate(age=round(age/365.25*100)/100) %>>%
    mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
  
  .plot_all = truncating_count_onco %>>%
    truncate_plot_allcantype(.permu_file = "oncogene/truncate.tsv")
  ggsave("age_plot/cumulative/oncogene/all_race/truncating.pdf",.plot_all,height = 5,width = 5)
  .plot_by = truncating_count_onco %>>%
    truncate_plot_bycantype(.permu_file = "oncogene/truncate_byCT.tsv")
  ggsave("age_plot/cumulative/oncogene/all_race/truncating_byCT.pdf",.plot_by,height = 10,width = 10)
  #oncogeneのtruncateあるなしのt_testでは？？p-value=0.2958
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
  ggsave("age_plot/fig/truncate/all_race/alloncogene.pdf",.plot,width = 10,height = 10)
}






















