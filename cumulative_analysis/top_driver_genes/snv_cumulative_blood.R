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
blood_patient_list = readtsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")%>>%
  filter(!is.na(bloodnorm),!is.na(age))
blood_patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")%>>%
  filter(!is.na(bloodnorm))

patient_with_ps = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")

################ read MAF extracted ################
blood_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
blood_white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/blood_white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",
                      col_types=cols(chr="c",LoF_filter="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
####全cancer_typeまとめて####
cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.blood=T,
                .patient_list = blood_patient_tdg,
                .permu_file = "TSG/missense_005_white_blood.tsv",.regression_size = 5,.pnum_size = 3)

#cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.5,.race="white",.blood=T,
#.patient_list = blood_patient_tdg,
#                .permu_file = "TSG/missense_05_white_blood.tsv")

cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 1,.race="white",.blood=T,
                .patient_list = blood_patient_tdg,
                .permu_file = "TSG/missense_1_white_blood.tsv")

#cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 5,.race="white",.blood=T,
#.patient_list = blood_patient_tdg,
#                .permu_file = "TSG/missense_5_white_blood.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",.blood=T,
                .patient_list = blood_patient_tdg,
                .permu_file = "TSG/missense_10_white_blood.tsv")


#### cancer_type ごとに #####
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,
                .patient_list = blood_patient_tdg,.blood = T,
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/missense_005_byCT_white_blood.tsv")


##### 0.01%ごとのplot
if(1){
  regression_table = make_regression_tabel(.maf=blood_white_maf_for_cumulative,
                                           .patient_list = blood_patient_tdg,
                                           .max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_white_blood.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_white_blood.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_white_blood.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-1_blood_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-10_blood_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .patient_list = blood_patient_tdg,.blood = T,
                         .permu_file = "TSG/missense_005_white.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,
                             .patient_list = blood_patient_tdg,.blood = T,
                             .facet_by_cancer_type = T,.race = "white",
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/missense_005_byCT_white.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.y = element_text(size =15),axis.title.x = element_blank())+
                       ggtitle(label = NULL),
                     x=0,y=0.38,width=0.36,height=0.62)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         strip.text = element_text(size=8,margin=margin(1,0,1,0)))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg_white_blood.pdf",.plot,width = 14,height = 8)
##########################################    pathogenc siteをtruncating 同様に除くと？ ##############################################
cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.pathogenic = T,
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/sup/missense_005_white_pathogenic_blood.tsv",.regression_size = 5,.pnum_size = 3)

cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 1,.race="white",.pathogenic = T,
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/sup/missense_1_white_pathogenic_blood.tsv")

cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 10,.race="white",.pathogenic = T,
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/sup/missense_10_white_pathogenic_blood.tsv")


#### cancer_type ごとに #####
cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.pathogenic = T,
                .patient_list = blood_patient_tdg,.blood = T,
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/sup/missense_005_byCT_white_pathogenic_blood.tsv")


##### 0.01%ごとのplot
if(1){
  regression_table = make_regression_tabel(.maf=white_maf_for_cumulative,
                                           .patient_list = blood_patient_tdg,
                                           .pathogenic = T,.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_white_pathogenic_blood.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_white_pathogenic_blood.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_white_pathogenic_blood.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-1_pathogenic_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-10_pathogenic_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .patient_list = blood_patient_tdg,.blood = T,
                         .permu_file = "TSG/sup/missense_005_white_pathogenic_blood.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.pathogenic = T)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,
                             .patient_list = blood_patient_tdg,.blood = T,
                             .facet_by_cancer_type = T,.race = "white",.pathogenic = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/sup/missense_005_byCT_white_pathogenic_blood.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.y = element_text(size =15),axis.title.x = element_blank())+
                       ggtitle(label = NULL),
                     x=0,y=0.38,width=0.36,height=0.62)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         strip.text = element_text(size=8,margin=margin(1,0,1,0)))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg_white_pathogenic_blood.pdf",.plot,width = 14,height = 8)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
########################################################### synonymou on TSG #############################################################
#silent
cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.mutype = "silent",
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/silent_005_white_blood.tsv",.regression_size = 5,.pnum_size = 3)
##X.intercept =60.15189,R=-0.264923,P=0.1202

cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 1,.race="white",.mutype = "silent",
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/silent_1_white_blood.tsv")
##X.intercept =60.25276,R=-0.1578433,P=0.126

cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 10,.race="white",.mutype = "silent",
                .patient_list = blood_patient_tdg,.blood = T,
                .permu_file = "TSG/silent_10_white_blood.tsv")
##X.intercept =60.47167,R=-0.09961537,P=0.1229


#### cancer_type ごとに #####
cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.mutype = "silent",
                .patient_list = blood_patient_tdg,.blood = T,
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/silent_005_byCT_white_blood.tsv")


##### 0.01%ごとのplot
if(1){
  regression_table = make_regression_tabel(.maf=blood_white_maf_for_cumulative,
                                           .patient_list = blood_patient_tdg,
                                           .mutype = "silent",.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_syn_white_blood.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_syn_white_blood.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
#.plot_reg=regression_plot_byside(regression_table_,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_syn_byside_white_blood.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-1_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-10_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .patient_list = blood_patient_tdg,.blood = T,
                         .permu_file = "TSG/silent_005_white_blood.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.mutype = "silent")+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.mutype="silent",
                             .patient_list = blood_patient_tdg,.blood = T,
                             .facet_by_cancer_type = T,.race = "white",
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/silent_005_byCT_white_blood.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.y = element_text(size =15),axis.title.x = element_blank())+
                       ggtitle(label = NULL),
                     x=0,y=0.38,width=0.36,height=0.62)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         strip.text = element_text(size=8,margin=margin(1,0,1,0)))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/poster_fig/maf005_syn_with_logreg_white_blood.pdf",.plot,width = 14,height = 8)



##########################################################################################################
##########################################################################################################
##########################################################################################################
########################## nonsyn とsynを並べて (nonsynはpathogenicも除いて)##############################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-1_pathogenic_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-10_pathogenic_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .patient_list = blood_patient_tdg,.blood = T,
                         .permu_file = "TSG/sup/missense_005_white_pathogenic_blood.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.pathogenic = T)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_white_pathogenic_blood.tsv")
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)


lm1=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-1_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-10_blood_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005_syn=cumulative_plot(.maf= blood_white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .patient_list = blood_patient_tdg,.blood = T,
                         .permu_file = "TSG/silent_005_white_blood.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.mutype = "silent")+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
regression_table=read_tsv("age_plot/cumulative/regression/TSG_syn_white_blood.tsv")
.plot_reglog_syn = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)

.plot = cowplot::plot_grid(.plot005,.plot005_syn+theme(axis.title.y = element_blank()),
                           .plot_reglog,.plot_reglog_syn+theme(axis.title.y = element_blank()),
                           labels = "auto",
                           rel_heights = c(1.3,0.7),scale = 0.95)
.plot
ggsave("age_plot/fig/poster_fig/maf005_syn_nonsyn_with_logreg_white_blood.pdf",.plot,width = 12,height = 8)
