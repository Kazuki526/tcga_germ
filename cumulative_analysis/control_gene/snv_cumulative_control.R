library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
loadNamespace('cowplot')
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
options(scipen=100)

control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")

####################################### TCGA data ############################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv")
all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
##################################### gnomAD data #############################################
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))
######################################################################################################
truncating_count_cont_rare = all_maf_for_cumulative_cont %>>%
  filter(MAF < 0.0005,FILTER=="PASS")%>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  inner_join(control_genes%>>%dplyr::select(gene_symbol))%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>%dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_cont)%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))
#missense
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.permu_file = "all_race/missense005.tsv")
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.permu_file = "all_race/missense1.tsv")
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.permu_file = "all_race/missense10.tsv")

#silent
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "all_race/silent005.tsv")
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent1.tsv")
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent10.tsv")

######## by cancer type #######
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.facet_by_cancer_type=T,
                     .permu_file = "all_race/missense005_byCT.tsv",.title = T)
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.facet_by_cancer_type=T,
                     .permu_file = "all_race/silent005_byCT.tsv",.title = T)
#####################################################################################################
# MAF 0.01%ごとに
if(1){
 regression_table_cont = make_regression_tabel_cont() 
 write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn.tsv")
 regression_table_silent_cont = make_regression_tabel_cont(.mutype = "silent")
 write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn.tsv")
}
######################################################################################################
############################################ white ###################################################
######################################################################################################
truncating_count_cont_rare_white = white_maf_for_cumulative_cont %>>%
  filter(FILTER=="PASS",MAF < 0.0005)%>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  inner_join(control_genes%>>%dplyr::select(gene_symbol))%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_cont)%>>%filter(race=="white",!is.na(age))%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

#missense
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = T,.permu_file = "white/missense005_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 1, .more_1par = T,.permu_file = "white/missense1_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 10, .more_1par = T,.permu_file = "white/missense10_white.tsv")

#silent
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "white/silent005_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "white/silent1_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "white/silent10_white.tsv")


######## by cancer type #######
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = T,.facet_by_cancer_type=T,
                     .permu_file = "white/missense005_byCT_white.tsv",.title = T)
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = F,.mutype = "silent",.facet_by_cancer_type=T,
                     .permu_file = "white/silent005_byCT_white.tsv",.title = T)
#####################################################################################################
# MAF 0.01%ごとに
if(1){
  regression_table_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,.race="white") 
  write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn_white.tsv")
  regression_table_silent_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.mutype = "silent")
  write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn_white.tsv")
}
############################### figure 用に調整 ##############################################
if(0){
  .plot = cowplot::plot_grid(.plot005 + theme(axis.title = element_text(size =15),title = element_text(size = 20)),
                             .plot005_by + theme(axis.title.y = element_blank(),axis.title.x = element_text(size = 15))+
                               ggtitle(label = NULL),
                             labels = "auto",label_size = 25,ncol = 2,scale = 0.95,rel_widths = c(1,1.8))
  ggsave("age_plot/fig/control/maf005_nonsyn.pdf",.plot,width = 14,height = 8)
  .plot = cowplot::plot_grid(.plot05 + theme(axis.title = element_text(size =15),title = element_text(size = 20)),
                             .plot05_by + theme(axis.title.y = element_blank(),axis.title.x = element_text(size = 15),
                                                axis.text.x = element_text(size=8))+ ggtitle(label = NULL),
                             labels = "auto",label_size = 25,ncol = 2,scale = 0.95,rel_widths = c(1,1.8))
  ggsave("age_plot/fig/control/maf05_nonsyn.pdf",.plot,width = 14,height = 8)
  #####################################
  .plot_reg = cowplot::ggdraw()+
    cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+annotate("rect",xmin=0,xmax=1,ymin=-0.49,ymax=0.02,alpha=0.2)+
                         scale_y_continuous(limits = c(-0.49,0.02),expand = c(0,0)),
                       x=0.05, y=0.53, width = 0.9, height = 0.45)+
    cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
                       x=0.05, y=0   , width = 0.9, height = 0.50)+
    cowplot::draw_text("regression coefficient",size = 30, x=0.025, y=0.5, angle=90)
  .plots_reg = cowplot::ggdraw()+
    cowplot::draw_plot(.plots_reg10 + theme(axis.title = element_blank())+annotate("rect",xmin=0,xmax=1,ymin=-0.72,ymax=0.02,alpha=0.2)+
                         scale_y_continuous(limits = c(-0.72,0.02),expand = c(0,0)),
                       x=0, y=0.53, width = 0.9, height = 0.45)+
    cowplot::draw_plot(.plots_reg1  + theme(axis.title.y = element_blank()),
                       x=0, y=0   , width = 0.9, height = 0.50)
  .plot = cowplot::plot_grid(.plot_reg,.plots_reg)
  ggsave("age_plot/fig/presentation/control_reg.pdf",.plot,width = 15,height = 8)
}
#####################################################################################
############################### figure 用に調整 ##############################################
.plot_trunc = truncating_count_cont_rare %>>%
  truncate_plot_allcantype(.permu_file = "all_race/truncate_rare.tsv")
lm1=read_tsv("age_plot/cumulative/all_race/lm_missense0-1regression_all_race.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/all_race/lm_missense0-10regression_all_race.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot_cont(.MAF_end = 0.05,.more_1par = T,
                         .permu_file = "all_race/missense005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
lm1s=read_tsv("age_plot/cumulative/all_race/lm_silent0-1regression_all_race.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10s=read_tsv("age_plot/cumulative/all_race/lm_silent0-10regression_all_race.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005=cumulative_plot_cont(.MAF_end = 0.05,.more_1par = F,.mutype = "silent",
                          .permu_file = "all_race/silent005.tsv",.all_color = "darkred",.save = F,
                          .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1s$X.Intercept.,slope=lm1s$missense_num),colour="blue",linetype=lm1s$LT)+
  geom_abline(aes(intercept=lm10s$X.Intercept.,slope=lm10s$missense_num),colour="green",linetype=lm10s$LT)

regression_table=read_tsv("age_plot/cumulative/regression/nonsyn.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/syn.tsv")
.plot_reg = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_reg  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)


.plot_truncw = truncating_count_cont_rare_white %>>%
  truncate_plot_allcantype(.permu_file = "all_race/truncate_rare_white.tsv")
lm1w=read_tsv("age_plot/cumulative/white/lm_missense0-1regression_white.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10w=read_tsv("age_plot/cumulative/white/lm_missense0-10regression_white.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005w=cumulative_plot_cont(.maf= white_maf_for_cumulative_cont,.MAF_end = 0.05,.race = "white",.more_1par = T,
                          .permu_file = "white/missense005_white.tsv",.all_color = "darkred",.save = F,
                          .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1w$X.Intercept.,slope=lm1w$missense_num),colour="blue",linetype=lm1w$LT)+
  geom_abline(aes(intercept=lm10w$X.Intercept.,slope=lm10w$missense_num),colour="green",linetype=lm10w$LT)
lm1sw=read_tsv("age_plot/cumulative/white/lm_silent0-1regression_white.tsv") %>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10sw=read_tsv("age_plot/cumulative/white/lm_silent0-10regression_white.tsv") %>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005w=cumulative_plot_cont(.maf= white_maf_for_cumulative_cont,.MAF_end = 0.05,.race = "white",.mutype = "silent",
                           .permu_file = "white/silent005_white.tsv",.all_color = "darkred",.save = F,
                           .more_1par = F,.regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1sw$X.Intercept.,slope=lm1sw$missense_num),colour="blue",linetype=lm1sw$LT)+
  geom_abline(aes(intercept=lm10sw$X.Intercept.,slope=lm10sw$missense_num),colour="green",linetype=lm10sw$LT)

regression_table=read_tsv("age_plot/cumulative/regression/nonsyn_white.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/syn_white.tsv")
.plot_regw = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_regw  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_trunc+ggtitle(label = NULL),x=0.05,y=0.5,width=0.25,height=0.5)+
  cowplot::draw_plot(.plot005+ggtitle(label = NULL),x=0.3,y=0.7,width=0.35,height=0.3)+
  cowplot::draw_plot(.plot_reg,x=0.3,y=0.5,width=0.35,height=0.2)+  
  cowplot::draw_plot(.plots005+ggtitle(label = NULL),x=0.65,y=0.7,width=0.35,height=0.3)+
  cowplot::draw_plot(.plots_reg,x=0.65,y=0.5,width=0.35,height=0.2)+
  
  cowplot::draw_plot(.plot_truncw+ggtitle(label = NULL),x=0.05,y=0,width=0.25,height=0.5)+
  cowplot::draw_plot(.plot005w+ggtitle(label = NULL),x=0.3,y=0.2,width=0.35,height=0.3)+
  cowplot::draw_plot(.plot_regw,x=0.3,y=0,width=0.35,height=0.2)+  
  cowplot::draw_plot(.plots005w+ggtitle(label = NULL),x=0.65,y=0.2,width=0.35,height=0.3)+
  cowplot::draw_plot(.plots_regw,x=0.65,y=0,width=0.35,height=0.2)+
  
  cowplot::draw_plot_label("a",x=0.06,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.31,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.31,y=0.69,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("d",x=0.66,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("e",x=0.66,y=0.69,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_text("In all races",x=0.025,y=0.75,size=20,angle=90)+
  cowplot::draw_plot_label("f",x=0.06,y=0.49,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("g",x=0.31,y=0.49,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("h",x=0.31,y=0.19,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("i",x=0.66,y=0.49,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("j",x=0.66,y=0.19,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_text("In Caucasian",x=0.025,y=0.25,size=20,angle=90)


.plot
ggsave("age_plot/fig/control_ns_reg_and_violin_allwhite.pdf",.plot,width = 12,height = 8)


