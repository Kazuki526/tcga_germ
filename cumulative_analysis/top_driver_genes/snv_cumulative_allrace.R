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
  mutate(role=ifelse(is.na(role),"TSG",role))%>>%
  dplyr::select(gene_symbol,role)
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")
patient_race = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_race.tsv")
patient_with_ps = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")

################ read MAF extracted ################
all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",col_types=cols(chr="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
####全cancer_typeまとめて####
cumulative_plot(.MAF_end = 0.05,
                .permu_file = "TSG/missense_005.tsv",.regression_size = 5,.pnum_size = 3)
##X.intercept =60.32619,R=-0.3803452,P=0.0361

cumulative_plot(.MAF_end = 0.5,
                .permu_file = "TSG/missense_05.tsv")
##X.intercept =60.45501,R=-0.3311621,P=0.0291

cumulative_plot(.MAF_end = 1,
                .permu_file = "TSG/missense_1.tsv")
##X.intercept =60.73801,R=-0.3175257,P=0.0001

cumulative_plot(.MAF_end = 5,
                .permu_file = "TSG/missense_5.tsv")
##X.intercept =60.15509,R=-0.04198409,P=0.3403

cumulative_plot(.MAF_end = 10,
                .permu_file = "TSG/missense_10.tsv")
##X.intercept =60.78216,R=-0.1493859,P=0.0036


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
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/missense_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  ##MAF<1  :X.intercept =60.73801,R=-0.3175257,P=0.0001
  geom_abline(aes(intercept=60.73801,slope=--0.3175257),colour="blue")+
  ##MAF<10 :X.intercept =60.78216,R=-0.1493859,P=0.0036
  geom_abline(aes(intercept=60.78216,slope=-0.1493859),colour="green")
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/missense_005_byCT.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.x = element_text(size =15),
                                      axis.title.y = element_text(size =15))+
                       ggtitle(label = NULL),
                     x=0,y=0.35,width=0.36,height=0.65)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_text(size = 15))+
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
                .permu_file = "TSG/sup/missense_005_pathogenic.tsv",.regression_size = 5,.pnum_size = 3,.save = F)
##X.intercept =60.33791,R=-0.3592016,P=0.0473

cumulative_plot(.MAF_end = 0.5,.pathogenic = T,
                .permu_file = "TSG/missense_05_pathogenic.tsv",.save = F)
##X.intercept =60.46035,R=-0.3130085,P=0.0364

cumulative_plot(.MAF_end = 1,.pathogenic = T,
                .permu_file = "TSG/sup/missense_1_pathogenic.tsv",.save = F)
##X.intercept =60.7622,R=-0.3187213,P=0.0002

cumulative_plot(.MAF_end = 5,.pathogenic = T,
                .permu_file = "TSG/sup/missense_5_pathogenic.tsv",.save = F)
##X.intercept =60.15763,R=-0.03358689,P=0.3694

cumulative_plot(.MAF_end = 10,.pathogenic = T,
                .permu_file = "TSG/sup/missense_10_pathogenic.tsv",.save = F)
##X.intercept =60.78392,R=-0.1459481,P=0.0041


#### cancer_type ごとに #####
.plot005_by =cumulative_plot(.MAF_end = 0.05,.pathogenic = T,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/sup/missense_005_byCT_pathogenic.tsv")


##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.pathogenic = T,.max_maf = 50)
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_pathogenic.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_pathogenic.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_pathogenic.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/sup/missense_005_pathogenic.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.pathogenic = T)+
  ##MAF<1  :X.intercept =60.7622,R=-0.3187213,P=0.0002
  geom_abline(aes(intercept=60.7622,slope=-0.3187213),colour="blue")+
  ##MAF<10 :X.intercept =60.78392,R=-0.1459481,P=0.0041
  geom_abline(aes(intercept=60.78392,slope=-0.1459481),colour="green")
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,.pathogenic = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/sup/missense_005_byCT_pathogenic.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.x = element_text(size =15),
                                      axis.title.y = element_text(size =15))+
                       ggtitle(label = NULL),
                     x=0,y=0.35,width=0.36,height=0.65)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_text(size = 15))+
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
cumulative_plot(.MAF_end = 0.05,.mutype = "silent",
                .permu_file = "TSG/silent_005tsv",.regression_size = 5,.pnum_size = 3)
##X.intercept =60.15526,R=-0.2331806,P=0.1329

cumulative_plot(.MAF_end = 1,.mutype = "silent",
                .permu_file = "TSG/silent_1.tsv")
##X.intercept =60.41023,R=-0.1507266,P=0.001

cumulative_plot(.MAF_end = 10,.mutype = "silent",
                .permu_file = "TSG/silent_10.tsv")
##X.intercept =60.73014,R=-0.1080266,P=0.0003


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
.plot005=cumulative_plot(.MAF_end = 0.05,
                         .permu_file = "TSG/silent_005.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.mutype = "silent")+
  ##MAF<1  :X.intercept =60.41023,R=-0.1507266,P=0.001
  geom_abline(aes(intercept=60.41023,slope=-0.1507266),colour="blue")+
  ##MAF<10 :X.intercept =60.73014,R=-0.1080266,P=0.0003
  geom_abline(aes(intercept=60.73014,slope=-0.1080266),colour="green")
.plot005_by =cumulative_plot(.MAF_end = 0.05,
                             .facet_by_cancer_type = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/silent_005_byCT.tsv",.all_color = "darkred")
#reg_plot logで
.plot_reglog = regression_plot_log(regression_table,.dred=0.05,.blue=1,.green = 10)
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot005 + theme(axis.title.x = element_text(size =15),
                                      axis.title.y = element_text(size =15))+
                       ggtitle(label = NULL),
                     x=0,y=0.35,width=0.36,height=0.65)+
  cowplot::draw_plot(.plot005_by + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_text(size = 15))+
                       ggtitle(label = NULL),
                     x=0.36,y=0,width=0.64,height=1)+
  cowplot::draw_plot(.plot_reglog,x=0,y=0,width=0.36,height=0.35)+
  cowplot::draw_plot_label("a",x=0.01,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("b",x=0.36,y=0.99,hjust = 0,vjust = 1,size = 20)+
  cowplot::draw_plot_label("c",x=0.01,y=0.35,hjust = 0,vjust = 0,size = 20)
.plot
ggsave("age_plot/fig/poster_fig/maf005_syn_with_logreg.pdf",.plot,width = 14,height = 8)


#