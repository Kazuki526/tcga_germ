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
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",
                      col_types=cols(chr="c",LoF_filter="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
####全cancer_typeまとめて####
cumulative_plot(.MAF_end = 0.05,.permu_file = "TSG/missense_005.tsv",.regression_size = 5,.pnum_size = 3)
##X.intercept =60.49784-0.6645733,R=-0.6645733,P<0.0001
cumulative_plot(.MAF_end = 0.5,.permu_file = "TSG/missense_05.tsv")
cumulative_plot(.MAF_end = 1,.permu_file = "TSG/missense_1.tsv")
##X.intercept =60.54521,R=-0.3465629,P=0.0001
cumulative_plot(.MAF_end = 5,.permu_file = "TSG/missense_5.tsv")
cumulative_plot(.MAF_end = 10,.permu_file = "TSG/missense_10.tsv")
##X.intercept =60.79035,R=-0.1773188,P=0.0003


#### cancer_type ごとに #####
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                             .permu_file = "TSG/missense_005_byCT.tsv")

#########################################################################
# mutationのあるgene数でcount
cumulative_plot(.MAF_end = 0.05,.by_gene = T,.permu_file = "TSG/sup/missense_005_genenum.tsv")
############################################################
##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.max_maf = 50)
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/all_race/lm_missense0-1_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/all_race/lm_missense0-10_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/missense_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/missense_005_byCT.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg.pdf",.plot,width = 14,height = 8)
##########################################    pathogenc siteをtruncating 同様に除くと？ ##############################################
cumulative_plot(.MAF_end = 0.05,.pathogenic = T,
                .permu_file = "TSG/sup/missense_005_pathogenic.tsv",.regression_size = 5,.pnum_size = 3)
cumulative_plot(.MAF_end = 0.5,.pathogenic = T,.permu_file = "TSG/missense_05_pathogenic.tsv")
cumulative_plot(.MAF_end = 1,.pathogenic = T,.permu_file = "TSG/sup/missense_1_pathogenic.tsv")
cumulative_plot(.MAF_end = 5,.pathogenic = T,.permu_file = "TSG/sup/missense_5_pathogenic.tsv")
cumulative_plot(.MAF_end = 10,.pathogenic = T,.permu_file = "TSG/sup/missense_10_pathogenic.tsv")

#### cancer_type ごとに #####
.plot005_by =cumulative_plot(.MAF_end = 0.05,.pathogenic = T,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                             .permu_file = "TSG/sup/missense_005_byCT_pathogenic.tsv")


##### 0.01%ごとのplot
if(1){
  regression_table = make_regression_tabel(.pathogenic = T,.max_maf = 50)
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_pathogenic.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_pathogenic.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_pathogenic.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/all_race/lm_missense0-1_pathogenic_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/all_race/lm_missense0-10_pathogenic_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/sup/missense_005_pathogenic.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.pathogenic = T)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,.pathogenic = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/sup/missense_005_byCT_pathogenic.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg_pathogenic.pdf",.plot,width = 14,height = 8)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
########################################################### synonymou on TSG #############################################################
#silent
cumulative_plot(.MAF_end = 0.05,.mutype = "silent",.permu_file = "TSG/silent_005tsv",.regression_size = 5,.pnum_size = 3)
cumulative_plot(.MAF_end = 1,.mutype = "silent",.permu_file = "TSG/silent_1.tsv")
cumulative_plot(.MAF_end = 10,.mutype = "silent",.permu_file = "TSG/silent_10.tsv")

#### cancer_type ごとに #####
cumulative_plot(.MAF_end = 0.05,.mutype = "silent",
                .facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/silent_005_byCT.tsv")


##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.mutype = "silent",.max_maf = 50)
  write_df(regression_table,"age_plot/cumulative/regression/TSG_syn.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_syn.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_syn_byside.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1 =read_tsv("age_plot/cumulative/TSG/all_race/lm_silent0-1_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/all_race/lm_silent0-10_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/silent_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.mutype = "silent")+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.MAF_end = 0.05,.mutype="silent",
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/silent_005_byCT.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_syn_with_logreg.pdf",.plot,width = 14,height = 8)


