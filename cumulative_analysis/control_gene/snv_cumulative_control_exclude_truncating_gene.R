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
#missense
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.permu_file = "all_race/missense005_notrunc.tsv",.no_trunc=T)
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.permu_file = "all_race/missense1_notrunc.tsv",.no_trunc=T)
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.permu_file = "all_race/missense10_notrunc.tsv",.no_trunc=T)

#silent
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "all_race/silent005_notrunc.tsv",.no_trunc=T)
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent1_notrunc.tsv",.no_trunc=T)
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent10_notrunc.tsv",.no_trunc=T)


#####################################################################################################
# MAF 0.01%ごとに
if(1){
  regression_table_cont = make_regression_tabel_cont(.no_trunc=T) 
  write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn_notrunc.tsv")
  regression_table_silent_cont = make_regression_tabel_cont(.mutype = "silent",.no_trunc=T)
  write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn_notrunc.tsv")
}
######################################################################################################
############################################ white ###################################################
######################################################################################################
#missense
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.no_trunc=T,
                     .MAF_end = 0.05, .more_1par = T,.permu_file = "white/missense005_white_notrunc.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.no_trunc=T,
                     .MAF_end = 1, .more_1par = T,.permu_file = "white/missense1_white_notrunc.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 10, .more_1par = T,.permu_file = "white/missense10_white_notrunc.tsv")

#silent
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.no_trunc=T,
                     .MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "white/silent005_white_notrunc.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.no_trunc=T,
                     .MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "white/silent1_white_notrunc.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.no_trunc=T,
                     .MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "white/silent10_white_notrunc.tsv")



#####################################################################################################
# MAF 0.01%ごとに
if(1){
  regression_table_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,.race="white",.no_trunc=T) 
  write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn_white_notrunc.tsv")
  regression_table_silent_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,.race = "white",.mutype = "silent",.no_trunc=T)
  write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn_white_notrunc.tsv")
}

#####################################################################################
############################### figure 用に調整 ##############################################
lm1=read_tsv("age_plot/cumulative/all_race/lm_missense0-1regression_notrunc.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/all_race/lm_missense0-10regression_notrunc.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot_cont(.MAF_end = 0.05,.more_1par = T,.no_trunc=T,
                              .permu_file = "all_race/missense005_notrunc.tsv",.all_color = "darkred",.save = F,
                              .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
lm1s=read_tsv("age_plot/cumulative/all_race/lm_silent0-1regression_notrunc.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10s=read_tsv("age_plot/cumulative/all_race/lm_silent0-10regression_notrunc.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005=cumulative_plot_cont(.MAF_end = 0.05,.more_1par = F,.mutype = "silent",.no_trunc=T,
                               .permu_file = "all_race/silent005_notrunc.tsv",.all_color = "darkred",.save = F,
                               .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1s$X.Intercept.,slope=lm1s$missense_num),colour="blue",linetype=lm1s$LT)+
  geom_abline(aes(intercept=lm10s$X.Intercept.,slope=lm10s$missense_num),colour="green",linetype=lm10s$LT)

regression_table=read_tsv("age_plot/cumulative/regression/nonsyn_notrunc.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/syn_notrunc.tsv")
.plot_reg = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_reg  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)


lm1w=read_tsv("age_plot/cumulative/white/lm_missense0-1regression_white_notrunc.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10w=read_tsv("age_plot/cumulative/white/lm_missense0-10regression_white_notrunc.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005w=cumulative_plot_cont(.maf= white_maf_for_cumulative_cont,.MAF_end = 0.05,.race = "white",.more_1par = T,
                               .permu_file = "white/missense005_white_notrunc.tsv",.all_color = "darkred",.save = F,
                               .regression_size = 5,.pnum_size = 3,.no_trunc=T)+
  geom_abline(aes(intercept=lm1w$X.Intercept.,slope=lm1w$missense_num),colour="blue",linetype=lm1w$LT)+
  geom_abline(aes(intercept=lm10w$X.Intercept.,slope=lm10w$missense_num),colour="green",linetype=lm10w$LT)
lm1sw=read_tsv("age_plot/cumulative/white/lm_silent0-1regression_white_notrunc.tsv") %>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10sw=read_tsv("age_plot/cumulative/white/lm_silent0-10regression_white_notrunc.tsv") %>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005w=cumulative_plot_cont(.maf= white_maf_for_cumulative_cont,.MAF_end = 0.05,.race = "white",.mutype = "silent",
                                .permu_file = "white/silent005_white_notrunc.tsv",.all_color = "darkred",.save = F,
                                .more_1par = F,.regression_size = 5,.pnum_size = 3,.no_trunc=T)+
  geom_abline(aes(intercept=lm1sw$X.Intercept.,slope=lm1sw$missense_num),colour="blue",linetype=lm1sw$LT)+
  geom_abline(aes(intercept=lm10sw$X.Intercept.,slope=lm10sw$missense_num),colour="green",linetype=lm10sw$LT)

regression_table=read_tsv("age_plot/cumulative/regression/nonsyn_white_notrunc.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/syn_white_notrunc.tsv")
.plot_regw = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_regw  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005+ggtitle(label = NULL),x=0.03,y=0.7,width=0.48,height=0.3)+
  cowplot::draw_plot(.plot_reg,x=0.03,y=0.5,width=0.48,height=0.2)+  
  cowplot::draw_plot(.plots005+ggtitle(label = NULL),x=0.52,y=0.7,width=0.48,height=0.3)+
  cowplot::draw_plot(.plots_reg,x=0.52,y=0.5,width=0.48,height=0.2)+
  
  cowplot::draw_plot(.plot005w+ggtitle(label = NULL),x=0.03,y=0.2,width=0.48,height=0.3)+
  cowplot::draw_plot(.plot_regw,x=0.03,y=0,width=0.48,height=0.2)+  
  cowplot::draw_plot(.plots005w+ggtitle(label = NULL),x=0.52,y=0.2,width=0.48,height=0.3)+
  cowplot::draw_plot(.plots_regw,x=0.52,y=0,width=0.48,height=0.2)+
  
  cowplot::draw_plot_label("a",x=0.04,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.04,y=0.69,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.53,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("d",x=0.53,y=0.69,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_text("In all races",x=0.015,y=0.75,size=20,angle=90)+
  cowplot::draw_plot_label("e",x=0.04,y=0.49,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("f",x=0.04,y=0.19,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("g",x=0.53,y=0.49,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("h",x=0.53,y=0.19,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_text("In Caucasian",x=0.015,y=0.25,size=20,angle=90)


.plot
ggsave("age_plot/fig/control_ns_reg_and_violin_allwhite_notrunc.pdf",.plot,width = 12,height = 8)


