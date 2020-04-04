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
patient_with_ps = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")

################ read MAF extracted ################
all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",
                      col_types=cols(chr="c",LoF_filter="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
###############################################################################################
######################################### All race ############################################
###############################################################################################
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

cumulative_plot(.MAF_end = 0.05,.permu_file = "oncogene/missense_005.tsv",.role = "oncogene")

cumulative_plot(.MAF_end = 1,.permu_file = "oncogene/missense_1.tsv", .role = "oncogene")

cumulative_plot(.MAF_end = 10, .permu_file = "oncogene/missense_10.tsv", .role = "oncogene")

##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.max_maf = 50,.role = "oncogene")
  write_df(regression_table,"age_plot/cumulative/regression/oncogene_nonsyn.tsv")
}
#silent
cumulative_plot(.MAF_end = 0.05,.role="oncogene",.mutype = "silent",
                .permu_file = "oncogene/silent_005.tsv",)

cumulative_plot(.MAF_end = 1,.role="oncogene",.mutype = "silent",
                .permu_file = "oncogene/silent_1.tsv")

cumulative_plot(.MAF_end = 10,.role="oncogene",.mutype = "silent",
                .permu_file = "oncogene/silent_10.tsv")
if(0){
  regression_table = make_regression_tabel(.role="oncogene",.mutype = "silent",.max_maf = 50)
  write_df(regression_table,"age_plot/cumulative/regression/oncogene_syn.tsv")
}
############################################################################################
######################################### white ############################################
############################################################################################
truncating_count_onco_rare_white = white_maf_for_cumulative %>>%
  filter(chr!="chrX",FILTER=="PASS",MAF < 0.0005)%>>%
  left_join(driver_genes %>>%dplyr::select(gene_symbol,role),by="gene_symbol") %>>%
  filter(mutype=="truncating"|mutype=="splice") %>>%
  filter(role=="oncogene"|role=="oncogene/TSG")%>>%
  #count(cancer_type,patient_id,gene_symbol) %>>% dplyr::select(-n)%>>%
  group_by(cancer_type,patient_id) %>>%summarise(truncating_count_n=n())%>>%ungroup()%>>%
  right_join(patient_tdg)%>>%filter(race=="white",!is.na(age))%>>%
  mutate(age=round(age/365.25*100)/100) %>>%
  mutate(truncating_count_n=ifelse(is.na(truncating_count_n),0,truncating_count_n))

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.role="oncogene",
                .permu_file = "oncogene/missense_005_white.tsv",.regression_size = 5,.pnum_size = 3)

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 1,.race="white",.role="oncogene",
                .permu_file = "oncogene/missense_1_white.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",.role="oncogene",
                .permu_file = "oncogene/missense_10_white.tsv")

##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.maf= white_maf_for_cumulative,.max_maf = 50,.race="white",.role="oncogene")
  write_df(regression_table,"age_plot/cumulative/regression/oncogene_nonsyn_white.tsv")
}
#silent
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.mutype = "silent",.role="oncogene",
                .permu_file = "oncogene/silent_005_white.tsv",.regression_size = 5,.pnum_size = 3)

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 1,.race="white",.mutype = "silent",.role="oncogene",
                .permu_file = "oncogene/silent_1_white.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",.mutype = "silent",.role="oncogene",
                .permu_file = "oncogene/silent_10_white.tsv")
if(0){
  regression_table = make_regression_tabel(.maf=white_maf_for_cumulative,.mutype = "silent",.max_maf = 50,.race="white",.role="oncogene")
  write_df(regression_table,"age_plot/cumulative/regression/oncogene_syn_white.tsv")
}

############################### figure 用に調整 ##############################################
if(0){
.plot_trunc = truncating_count_onco_rare %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare.tsv")
lm1=read_tsv("age_plot/cumulative/oncogene/all_race/lm_missense0-1_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/oncogene/all_race/lm_missense0-10_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.MAF_end = 0.05,.role = "oncogene",
                         .permu_file = "oncogene/missense_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
lm1s=read_tsv("age_plot/cumulative/oncogene/all_race/lm_silent0-1_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10s=read_tsv("age_plot/cumulative/oncogene/all_race/lm_silent0-10_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005=cumulative_plot(.MAF_end = 0.05,.role = "oncogene",.mutype = "silent",
                          .permu_file = "oncogene/silent_005.tsv",.all_color = "darkred",.save = F,
                          .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1s$X.Intercept.,slope=lm1s$missense_num),colour="blue",linetype=lm1s$LT)+
  geom_abline(aes(intercept=lm10s$X.Intercept.,slope=lm10s$missense_num),colour="green",linetype=lm10s$LT)

regression_table=read_tsv("age_plot/cumulative/regression/oncogene_nonsyn.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/oncogene_syn.tsv")
.plot_reg = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_reg  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)


.plot_truncw = truncating_count_onco_rare_white %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_white.tsv")
lm1w=read_tsv("age_plot/cumulative/oncogene/white/lm_missense0-1_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10w=read_tsv("age_plot/cumulative/oncogene/white/lm_missense0-10_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005w=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",.role = "oncogene",
                         .permu_file = "oncogene/missense_005_white.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1w$X.Intercept.,slope=lm1w$missense_num),colour="blue",linetype=lm1w$LT)+
  geom_abline(aes(intercept=lm10w$X.Intercept.,slope=lm10w$missense_num),colour="green",linetype=lm10w$LT)
lm1sw=read_tsv("age_plot/cumulative/oncogene/white/lm_silent0-1_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10sw=read_tsv("age_plot/cumulative/oncogene/white/lm_silent0-10_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005w=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",.mutype = "silent",
                          .permu_file = "oncogene/silent_005_white.tsv",.all_color = "darkred",.save = F,
                          .role = "oncogene",.regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1sw$X.Intercept.,slope=lm1sw$missense_num),colour="blue",linetype=lm1sw$LT)+
  geom_abline(aes(intercept=lm10sw$X.Intercept.,slope=lm10sw$missense_num),colour="green",linetype=lm10sw$LT)

regression_table=read_tsv("age_plot/cumulative/regression/oncogene_nonsyn_white.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/oncogene_syn_white.tsv")
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
ggsave("age_plot/fig/poster_fig/oncogene_ns_reg_and_violin_allwhite.pdf",.plot,width = 16,height = 12)
}
.plot_truncw = truncating_count_onco_rare_white %>>%
  truncate_plot_allcantype(.permu_file = "oncogene/truncate_rare_white.tsv")
lm1w=read_tsv("age_plot/cumulative/oncogene/white/lm_missense0-1_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10w=read_tsv("age_plot/cumulative/oncogene/white/lm_missense0-10_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005w=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",.role = "oncogene",
                          .permu_file = "oncogene/missense_005_white.tsv",.all_color = "darkred",.save = F,
                          .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1w$X.Intercept.,slope=lm1w$missense_num),colour="blue",linetype=lm1w$LT)+
  geom_abline(aes(intercept=lm10w$X.Intercept.,slope=lm10w$missense_num),colour="green",linetype=lm10w$LT)
lm1sw=read_tsv("age_plot/cumulative/oncogene/white/lm_silent0-1_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10sw=read_tsv("age_plot/cumulative/oncogene/white/lm_silent0-10_regression.tsv") %>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plots005w=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",.mutype = "silent",
                           .permu_file = "oncogene/silent_005_white.tsv",.all_color = "darkred",.save = F,
                           .role = "oncogene",.regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1sw$X.Intercept.,slope=lm1sw$missense_num),colour="blue",linetype=lm1sw$LT)+
  geom_abline(aes(intercept=lm10sw$X.Intercept.,slope=lm10sw$missense_num),colour="green",linetype=lm10sw$LT)

regression_table_trunc = read_tsv("age_plot/cumulative/regression/oncogene_trunc_white.tsv")
regression_table=read_tsv("age_plot/cumulative/regression/oncogene_nonsyn_white.tsv")
regression_table_silent=read_tsv("age_plot/cumulative/regression/oncogene_syn_white.tsv")
.plott_regw = regression_plot_log(regression_table_trunc,.black=0.05)
.plot_regw = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plots_regw  = regression_plot_log(regression_table_silent,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_truncw+ggtitle(label = NULL),x=0,y=0.4,width=0.25,height=0.6)+
  cowplot::draw_plot(.plott_regw,x=0,y=0,width=0.25,height=0.4)+
  cowplot::draw_plot(.plot005w+ggtitle(label = NULL),x=0.3,y=0.4,width=0.35,height=0.6)+
  cowplot::draw_plot(.plot_regw,x=0.3,y=0,width=0.35,height=0.4)+  
  cowplot::draw_plot(.plots005w+ggtitle(label = NULL),x=0.65,y=0.4,width=0.35,height=0.6)+
  cowplot::draw_plot(.plots_regw,x=0.65,y=0,width=0.35,height=0.4)+
  
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 30)+
  cowplot::draw_plot_label("b",x=0.01,y=0.39,hjust = 0,vjust = 1,size = 30)+
  cowplot::draw_plot_label("c",x=0.31,y=0.99,hjust = 0,vjust = 1,size = 30)+
  cowplot::draw_plot_label("d",x=0.31,y=0.39,hjust = 0,vjust = 1,size = 30)+
  cowplot::draw_plot_label("e",x=0.66,y=0.99,hjust = 0,vjust = 1,size = 30)+
  cowplot::draw_plot_label("f",x=0.66,y=0.39,hjust = 0,vjust = 1,size = 30)
.plot
ggsave("age_plot/fig/poster_fig/oncogene_ns_reg_and_violin_white.pdf",.plot,width = 16,height = 12)
