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
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",col_types=cols(chr="c"))%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
####################################################################################################################
#####################################################  TSG  ########################################################
####################################################################################################################
####全cancer_typeまとめて####
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",
                .permu_file = "TSG/missense_005_white.tsv",.regression_size = 5,.pnum_size = 3)

#cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.5,.race="white",
#                .permu_file = "TSG/missense_05_white.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 1,.race="white",
                .permu_file = "TSG/missense_1_white.tsv")

#cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 5,.race="white",
#                .permu_file = "TSG/missense_5_white.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",
                .permu_file = "TSG/missense_10_white.tsv")


#### cancer_type ごとに #####
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/missense_005_byCT_white.tsv")

#########################################################################
# mutationのあるgene数でcount
.plot=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.by_gene = T,.save = F,
                        .permu_file = "TSG/sup/missense_005_genenum_white.tsv",.regression_size = 5,.pnum_size = 3)
ggsave("age_plot/cumulative/TSG/white/missense0-0.05_genenum.pdf")
############################################################
##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.maf=white_maf_for_cumulative,.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_white.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_white.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_white.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-1_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-10_regression.tsv")%>>%mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .permu_file = "TSG/missense_005_white.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,
                             .facet_by_cancer_type = T,.race = "white",
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/missense_005_byCT_white.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg_white.pdf",.plot,width = 14,height = 8)
 ##########################################    pathogenc siteをtruncating 同様に除くと？ ##############################################
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.pathogenic = T,
                .permu_file = "TSG/sup/missense_005_white_pathogenic.tsv",.regression_size = 5,.pnum_size = 3)

#cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.5,.race="white",.pathogenic = T,
#                .permu_file = "TSG/sup/missense_05_white_pathogenic.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 1,.race="white",.pathogenic = T,
                .permu_file = "TSG/sup/missense_1_white_pathogenic.tsv")

#cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 5,.race="white",.pathogenic = T,
#                .permu_file = "TSG/sup/missense_5_white_pathogenic.tsv")

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",.pathogenic = T,
                .permu_file = "TSG/sup/missense_10_white_pathogenic.tsv")


#### cancer_type ごとに #####
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.pathogenic = T,
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/sup/missense_005_byCT_white_pathogenic.tsv")


##### 0.01%ごとのplot
if(1){
  regression_table = make_regression_tabel(.maf=white_maf_for_cumulative,.pathogenic = T,.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_nonsyn_white_pathogenic.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_nonsyn_white_pathogenic.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_nonsyn_byside_white_pathogenic.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-1_pathogenic_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_missense0-10_pathogenic_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .permu_file = "TSG/sup/missense_005_white_pathogenic.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.pathogenic = T)+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,
                             .facet_by_cancer_type = T,.race = "white",.pathogenic = T,
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/sup/missense_005_byCT_white_pathogenic.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_nonsyn_with_logreg_white_pathogenic.pdf",.plot,width = 14,height = 8)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
########################################################### synonymou on TSG #############################################################
#silent
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race="white",.mutype = "silent",
                .permu_file = "TSG/silent_005_white.tsv",.regression_size = 5,.pnum_size = 3)
##X.intercept =60.15189,R=-0.264923,P=0.1202

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 1,.race="white",.mutype = "silent",
                .permu_file = "TSG/silent_1_white.tsv")
##X.intercept =60.25276,R=-0.1578433,P=0.126

cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 10,.race="white",.mutype = "silent",
                .permu_file = "TSG/silent_10_white.tsv")
##X.intercept =60.47167,R=-0.09961537,P=0.1229


#### cancer_type ごとに #####
cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.mutype = "silent",
                .race = "white",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = "TSG/silent_005_byCT_white.tsv")


##### 0.01%ごとのplot
if(0){
  regression_table = make_regression_tabel(.maf=white_maf_for_cumulative,.mutype = "silent",.max_maf = 50,.race="white")
  regression_table_ = make_regression_tabel(.mutype = "silent",.max_maf = 50,.race="white")
  write_df(regression_table,"age_plot/cumulative/regression/TSG_syn_white.tsv")
}
regression_table=read_tsv("age_plot/cumulative/regression/TSG_syn_white.tsv")
.plot_reg=regression_plot_byside(regression_table,.blue = 1,.green = 10)
.plot_reg=regression_plot_byside(regression_table_,.blue = 1,.green = 10)
ggsave("age_plot/fig/regression/TSG_syn_byside_white.pdf",.plot_reg,height = 6,width = 12)

##########################################################################################################
###################################### figrure用に調整 #######################################
lm1=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-1_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
lm10=read_tsv("age_plot/cumulative/TSG/white/lm_silent0-10_regression.tsv")%>>%
  mutate(LT=ifelse(p_value<0.05,"solid","dashed"))
.plot005=cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,.race = "white",
                         .permu_file = "TSG/silent_005_white.tsv",.all_color = "darkred",.save = F,
                         .regression_size = 5,.pnum_size = 3,.mutype = "silent")+
  geom_abline(aes(intercept=lm1$X.Intercept.,slope=lm1$missense_num),colour="blue",linetype=lm1$LT)+
  geom_abline(aes(intercept=lm10$X.Intercept.,slope=lm10$missense_num),colour="green",linetype=lm10$LT)
.plot005_by =cumulative_plot(.maf= white_maf_for_cumulative,.MAF_end = 0.05,
                             .facet_by_cancer_type = T,.race = "white",
                             .pnum_size = 2.5,.regression_size = 3.5,.width = 12,.save = F,
                             .permu_file = "TSG/silent_005_byCT_white.tsv",.all_color = "darkred")
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
ggsave("age_plot/fig/poster_fig/maf005_syn_with_logreg_white.pdf",.plot,width = 14,height = 8)


##############################################################################################
#どの遺伝子が1番効いてる？
get_lm_coef = function(.tbl){
  if(first(.tbl$focal =="yes")){
    .patient_count = left_join(patient_list,.tbl, by = "patient_id") %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num),
             age=signif(age/365.25,4))
    model=lm(age ~ missense_num, data = .patient_count)
    as.data.frame(as.list(coef(model))) %>>%
      mutate(lm_p_value = 1 - pf(summary(model)$fstatistic["value"],summary(model)$fstatistic["numdf"],
                                 summary(model)$fstatistic["dendf"]),
             t_test_p_value = t.test(.patient_count[.patient_count$missense_num > 0,]$age/365.25,
                                     .patient_count[.patient_count$missense_num ==0,]$age/365.25,
                                     alternative= "less")$p.value,
             all_missense_num = sum(.patient_count$missense_num))
  }else{
    data.frame(X.Intercept.=NA,missense_num=NA,lm_p_value=NA,t_test_p_value=NA,
               all_missense_num=sum(.tbl$missense_num,na.rm = T))
  }
}
driver_genes %>>%dplyr::select(gene,role) %>>%dplyr::rename(gene_symbol=gene) %>>%
  filter((role=="TSG"|role=="oncogene/TSG"),gene_symbol!="KMT2C") %>>%
  left_join(all_maf_for_cumlative %>>%filter(MAF<0.005)%>>%
              group_by(gene_symbol,patient_id) %>>% summarise(missense_num=sum(MAC))) %>>%
  mutate(focal = ifelse(is.na(patient_id),"no","yes")) %>>%
  nest(-gene_symbol,-role) %>>%
  mutate(regression =purrr::map(data, ~get_lm_coef(.))) %>>%
  dplyr::select(-role,-data) %>>%unnest() %>>%
  arrange(missense_num) %>>%
  write_df("~/Dropbox/install/tvz/TSG_by_gene_regression_coef.tsv")


#################################################################################################################
################################################## oncogene #####################################################
#################################################################################################################
####### missense ########
cumulative_plot(.MAF_end = 0.05,.path = "onco",.role = "oncogene",
                .permu_file = "oncogene/missense005.tsv",.test_tail = "two")

####### missense by cancer type ########
cumulative_plot(.MAF_end = 0.05,.path = "onco/by_cancer_type",
                .role = "oncogene",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = )
cumulative_plot(.MAF_end = 3,.path = "onco/by_cancer_type",
                .role = "oncogene",.facet_by_cancer_type = T,
                .pnum_size = 2.5,.regression_size = 3.5,.width = 12,
                .permu_file = )

##################################################################
regression_table_onco = make_regression_tabel(.maf=white_maf_for_cumulative,
                                              .role = "oncogene",.max_maf = 50)
write_df(regression_table_onco,"age_plot/cumulative/regression/onco_nonsyn.tsv")
regression_table_onco = read_tsv("age_plot/cumulative/regression/onco_nonsyn.tsv")
###################################### figrure用に調整 #######################################
# .plot_reg = cowplot::ggdraw()+
#   cowplot::draw_plot(.plot_reg10 + theme(axis.title = element_blank())+
#                        annotate("rect",xmin=0,xmax=1,ymin=-0.098,ymax=0.28,alpha=0.2)+
#                        scale_y_continuous(limits = c(-0.098,0.28),expand = c(0,0)),
#                      x=0.05, y=0.53, width = 0.9, height = 0.47)+
#   cowplot::draw_plot(.plot_reg1  + theme(axis.title.y = element_blank()),
#                      x=0.05, y=0   , width = 0.9, height = 0.53)+
#   cowplot::draw_text("regression coefficient",size = 15, x=0.025, y=0.5, angle=90)
# 
# .plot25  = cumulative_plot(.MAF_end = 2.5, .role = "oncogene",
#                             .save = F,.regression_size = 5,.pnum_size = 3)
# .plot005 = cumulative_plot(.MAF_end = 0.05, .role = "oncogene",
#                             .save = F,.regression_size = 5,.pnum_size = 3)
# .plot = cowplot::plot_grid(.plot_reg,
#                             cowplot::plot_grid(.plot005 + theme(axis.title.x = element_blank(),
#                                                                  axis.title.y = element_text(size = 15),
#                                                                  axis.text = element_text(size = 12),
#                                                                  title = element_text(size = 20)),
#                                                .plot25 + theme(title = element_text(size = 20),
#                                                                 axis.title = element_text(size = 15),
#                                                                 axis.text = element_text(size = 12)),
#                                                labels = c("b","c"),label_size = 30,
#                                                ncol = 1,rel_heights = c(1,1.1)),
#                             labels = c("a",""),ncol=2,rel_widths = c(1.5,1),label_size = 30)
# .plot
# ggsave("age_plot/fig/presentation/oncogene_nonsyn_reg_and_violin.pdf",.plot,width = 15,height =7 )


###############################################################################
#silent
cumulative_plot(.MAF_end = 0.05,.path = "oncogene/all_race",.role = "oncogene",.mutype = "silent",
                .permu_file = "oncogene/silent005.tsv")
regression_table_silent_onco = make_regression_tabel(.maf=white_maf_for_cumulative,
                                                     .mutype = "silent",.role = "oncogene",.max_maf = 50)
write_df(regression_table_silent_onco,"age_plot/cumulative/regression/onco_syn.tsv")
regression_table_silent_onco = read_tsv("age_plot/cumulative/regression/onco_syn.tsv")
###################################### figrure用に調整 #######################################
# .plots_reg = cowplot::ggdraw()+
#   cowplot::draw_plot(.plots_reg10 + theme(axis.title = element_blank())+
#                        annotate("rect",xmin=0,xmax=1,ymin=-0.37,ymax=0.02,alpha=0.2)+
#                        scale_y_continuous(limits = c(-0.37,0.02),expand = c(0,0)),
#                      x=0.05, y=0.53, width = 0.9, height = 0.47)+
#   cowplot::draw_plot(.plots_reg1  + theme(axis.title.y = element_blank())+
#                        scale_y_continuous(breaks = c(-0.3,-0.2,-0.1),labels = c(-0.3,-0.2,-0.1)),
#                      x=0.05, y=0   , width = 0.9, height = 0.53)+
#   cowplot::draw_text("regression coefficient",size = 15, x=0.025, y=0.5, angle=90)
# 
# .plots002  = cumulative_plot(.MAF_end = 0.02, .mutype = "silent",.role = "oncogene",
#                             .save = F,.regression_size = 5,.pnum_size = 3)
# .plots005 = cumulative_plot(.MAF_end = 0.05, .mutype = "silent",.role = "oncogene",
#                             .save = F,.regression_size = 5,.pnum_size = 3)
# .plots = cowplot::plot_grid(.plots_reg,
#                             cowplot::plot_grid(.plots002 + theme(axis.title.x = element_blank(),
#                                                                  axis.title.y = element_text(size = 15),
#                                                                  axis.text = element_text(size = 12),
#                                                                  title = element_text(size = 20)),
#                                                .plots005 + theme(title = element_text(size = 20),
#                                                                 axis.title = element_text(size = 15),
#                                                                 axis.text = element_text(size = 12)),
#                                                labels = c("b","c"),label_size = 30,
#                                                ncol = 1,rel_heights = c(1,1.1)),
#                             labels = c("a",""),ncol=2,rel_widths = c(1.5,1),label_size = 30)
# .plots
# ggsave("age_plot/fig/presentation/oncogene_syn_reg_and_violin.pdf",.plots,width = 15,height =7 )
