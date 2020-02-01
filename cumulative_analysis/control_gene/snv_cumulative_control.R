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
options(scipen=1)

control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")

####################################### TCGA data ############################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_race = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_race.tsv")
patien_all_info = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/all_patient_info.tsv")
all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")%>>%
  filter(FILTER=="PASS")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")%>>%
  filter(FILTER=="PASS")
##################################### gnomAD data #############################################
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))
######################################################################################################
#missense
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.permu_file = "all_race/missense005.tsv")
#cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,.permu_file = "missense05.tsv")
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.permu_file = "all_race/missense1.tsv")
##^X.intercept =60.35836,R=-0.03584669,P=0.0041
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.permu_file = "all_race/missense10.tsv")
##^X.intercept =60.55437,R=-0.01745769,P=0.0219 

#silent
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "all_race/silent005.tsv")
#cumulative_plot_cont(.MAF_end = 0.5, .more_1par = T,.mutype = "silent",.title = T,.permu_file = "all_race/silent05.tsv")
cumulative_plot_cont(.MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent1.tsv")
##^X.intercept =60.27423,R=-0.03990302,P=0.0031
cumulative_plot_cont(.MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "all_race/silent10.tsv")
##^X.intercept =60.67661,R=-0.0276724,P=0.0017

######## by cancer type #######
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,.facet_by_cancer_type=T,
                     .permu_file = "all_race/missense005_byCT.tsv",.title = T)
cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",.facet_by_cancer_type=T,
                     .permu_file = "all_race/silent005_byCT.tsv",.title = T)
#####################################################################################################
# MAF 0.01%ごとに
if(0){
 regression_table_cont = make_regression_tabel_cont() 
 write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn.tsv")
 regression_table_silent_cont = make_regression_tabel_cont(.mutype = "silent")
 write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn.tsv")
}
regression_table_cont=read_tsv("age_plot/cumulative/regression/nonsyn.tsv")
regression_plot_byside(regression_table_cont,.blue = 1,.green = 10)
regression_table_silent_cont=read_tsv("age_plot/cumulative/regression/syn.tsv")
regression_plot_byside(regression_table_silent_cont,.blue = 1,.green = 10)
######################################################################################################
############################################ white ###################################################
######################################################################################################
#missense
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = T,.permu_file = "white/missense005_white.tsv")
#cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
#                     .MAF_end = 0.5, .more_1par = T,.permu_file = "white/missense05_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 1, .more_1par = T,.permu_file = "white/missense1_white.tsv")
##^X.intercept =60.20037,R=-0.04025278,P=0.1951
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 10, .more_1par = T,.permu_file = "white/missense10_white.tsv")
##^X.intercept =59.7421,R=0.005655658,P=0.5944

#silent
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = F,.mutype = "silent",.title = T,.permu_file = "white/silent005_white.tsv")
#cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
#                     .MAF_end = 0.5, .more_1par = T,.mutype = "silent",.title = T,.permu_file = "white/silent05_white.tsv")
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 1, .more_1par = T,.mutype = "silent",.permu_file = "white/silent1_white.tsv")
##^X.intercept =60.05776,R=-0.03487504,P=0.2788
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 10, .more_1par = T,.mutype = "silent",.permu_file = "white/silent10_white.tsv")
##^X.intercept =60.10723,R=-0.009556983,P=0.367


######## by cancer type #######
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = T,.facet_by_cancer_type=T,
                     .permu_file = "white/missense005_byCT_white.tsv",.title = T)
cumulative_plot_cont(.maf=white_maf_for_cumulative_cont,.race = "white",
                     .MAF_end = 0.05, .more_1par = F,.mutype = "silent",.facet_by_cancer_type=T,
                     .permu_file = "white/silent005_byCT_white.tsv",.title = T)
#####################################################################################################
# MAF 0.01%ごとに
if(0){
  regression_table_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,.race="white") 
  write_df(regression_table_cont,"age_plot/cumulative/regression/nonsyn_white.tsv")
  regression_table_silent_cont = make_regression_tabel_cont(.maf=white_maf_for_cumulative_cont,
                                                            .race = "white",.mutype = "silent")
  write_df(regression_table_silent_cont,"age_plot/cumulative/regression/syn_white.tsv")
}
regression_table_cont=read_tsv("age_plot/cumulative/regression/nonsyn_white.tsv")
regression_plot_byside(regression_table_cont,.blue = 1,.green = 10)
regression_table_silent_cont=read_tsv("age_plot/cumulative/regression/syn_white.tsv")
regression_plot_byside(regression_table_silent_cont,.blue = 1,.green = 10)
############################### figure 用に調整 ##############################################
.plot = cowplot::plot_grid(.plot005 + theme(axis.title = element_text(size =15),
                                            title = element_text(size = 20)),
                           .plot005_by + theme(axis.title.y = element_blank(),
                                               axis.title.x = element_text(size = 15))+
                             ggtitle(label = NULL),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
.plot
ggsave("age_plot/fig/control/maf005_nonsyn.pdf",.plot,width = 14,height = 8)

.plot = cowplot::plot_grid(.plot05 + theme(axis.title = element_text(size =15),
                                           title = element_text(size = 20)),
                           .plot05_by + theme(axis.title.y = element_blank(),
                                              axis.title.x = element_text(size = 15),
                                              axis.text.x = element_text(size=8))+
                             ggtitle(label = NULL),
                           labels = "auto",label_size = 25,ncol = 2,scale = 0.95,
                           rel_widths = c(1,1.8))
.plot
ggsave("age_plot/fig/control/maf05_nonsyn.pdf",.plot,width = 14,height = 8)




##################################################################
.plot_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.49,ymax=0.02,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.49,0.02),expand = c(0,0)),
                     x=0.05, y=0.53, width = 0.9, height = 0.45)+
  cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
                     x=0.05, y=0   , width = 0.9, height = 0.50)+
  cowplot::draw_text("regression coefficient",size = 30, x=0.025, y=0.5, angle=90)
.plot_reg
.plots_reg = cowplot::ggdraw()+
  cowplot::draw_plot(.plots_reg10 + theme(axis.title = element_blank())+
                       annotate("rect",xmin=0,xmax=1,ymin=-0.72,ymax=0.02,alpha=0.2)+
                       scale_y_continuous(limits = c(-0.72,0.02),expand = c(0,0)),
                     x=0, y=0.53, width = 0.9, height = 0.45)+
  cowplot::draw_plot(.plots_reg1  + theme(axis.title.y = element_blank()),
                     x=0, y=0   , width = 0.9, height = 0.50)

.plot = cowplot::plot_grid(.plot_reg,.plots_reg)
ggsave("age_plot/fig/presentation/control_reg.pdf",.plot,width = 15,height = 8)

.plots005 = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = F,.mutype = "silent",
                                 .save = F,.regression_size = 5,.pnum_size = 3)+
  theme(title = element_text(size = 15),axis.text = element_text(size = 15),
        axis.title = element_text(size =20))
.plot005 = cumulative_plot_cont(.MAF_end = 0.05, .more_1par = T,
                                .save = F,.regression_size = 5,.pnum_size = 3)+
  theme(title = element_text(size = 15),axis.text = element_text(size = 15),
        axis.title = element_text(size =20))
.plot = cowplot::plot_grid(.plot_reg,.plots_reg,.plot005,.plots005,
                           labels = "auto",label_size = 30,ncol = 2)
.plot
ggsave("age_plot/fig/control/ns_reg_and_violin.pdf",.plot,width = 15,height = 15)
###############################################
#stage1,2でのみ
regression_table_stage12_cont = make_regression_tabel_cont(.patient_list = patient_list %>>%
                                                             left_join(all_patient_info %>>%
                                                                         dplyr::select(patient_id,stage)) %>>%
                                                             filter(stage==1)) 
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot()
ggsave("age_plot/control_region/regression_plot-1_missense_stage1_control.pdf",.plot,height = 8,width = 20)
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_missense_stage_1_control.pdf",.plot,height = 8,width = 20)

regression_table_stage12_silent_cont = make_regression_tabel_cont(.patient_list = patient_list %>>%
                                                                    left_join(all_patient_info %>>%
                                                                                dplyr::select(patient_id,stage)) %>>%
                                                                    filter(stage==1), .mutype = "silent") 
.plot=regression_table_stage12_silent_cont　%>>%
  regression_tbl_plot()
ggsave("age_plot/control_region/regression_plot-1_silent_stage1_control.pdf",.plot,height = 8,width = 20)
.plot=regression_table_stage12_cont　%>>%
  regression_tbl_plot(.maf_max = 10,.expand = 0.15)
ggsave("age_plot/control_region/regression_plot-10_silent_stage_1_control.pdf",.plot,height = 8,width = 20)


