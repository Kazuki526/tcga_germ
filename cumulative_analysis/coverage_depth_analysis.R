library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
loadNamespace('cowplot')
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

setwd('/Volumes/DR8TB2/tcga_rare_germ/')
######################################## gene infomation ##########################################
driver_genes=read_tsv("~/git/driver_genes/driver_genes.tsv")%>>%
  filter(refs>3) %>>%dplyr::rename(gene_symbol=gene)%>>%
  mutate(role=ifelse(is.na(role),"TSG/oncogene",role))%>>%
  dplyr::select(gene_symbol,role)
control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")

####################################### TCGA data ########################################################
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")#%>>%filter(!is.na(age))%>>%mutate(age = age/365.25)
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv")

all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%filter(chr!="chrX",FILTER=="PASS")
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%filter(chr!="chrX",FILTER=="PASS")
all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")%>>%filter(FILTER=="PASS")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")%>>%filter(FILTER=="PASS")

coverage = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/coverage.tsv.gz") 
coverage_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/coverage_cont.tsv.gz")
depth = read_tsv("/Volumes/areca42TB2/gdc/top_driver_gene/all_patient/patient_depth.tsv")
depth_cont = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/patient_depth.tsv")

#################################################################################################################
#coverage and onset age 
coverage %>>%mutate(TSG=TSG+`oncogene/TSG`,oncogene=oncogene+`oncogene/TSG`)%>>%
  dplyr::select(patient_id,TSG,oncogene)%>>%
  tidyr::gather(role,coverage,TSG,oncogene) %>>%
  right_join(patient_tdg)%>>%mutate(age=age/365.25)%>>%
  ggplot(aes(x=coverage,y=age))+
  geom_point()+
  facet_grid(.~role,scales = "free")+
  stat_smooth(method=lm, formula = y ~  +x,se=F,colour="black")
#negatively correlated with onset age,,,

###################################################################################################################
#depth and onset age
.plot = depth %>>%tidyr::pivot_longer(cols=c(norm_depth,tumor_depth),names_to="tn",values_to="depth")%>>%
  tidyr::pivot_wider(names_from=role,values_from=c(depth,length_bp))%>>%
  mutate(TSG=(depth_TSG+`depth_oncogene/TSG`)/(length_bp_TSG+`length_bp_oncogene/TSG`),
         oncogene=(depth_oncogene+`depth_oncogene/TSG`)/(length_bp_oncogene+`length_bp_oncogene/TSG`))%>>%
  dplyr::select(patient_id,tn,TSG,oncogene)%>>%
  tidyr::pivot_longer(cols=c(TSG,oncogene),names_to = "role",values_to = "mean_depth")%>>%
  right_join(patient_tdg)%>>%
  bind_rows(depth_cont%>>%
              tidyr::pivot_longer(cols=c(norm_depth,tumor_depth),names_to="tn",values_to="depth")%>>%
              mutate(mean_depth=depth/length_bp)%>>%dplyr::select(-depth,-length_bp)%>>%
              right_join(patient_cont))%>>%
  mutate(age=age/365.25,tn=ifelse(tn=="norm_depth","Normal sample","Tumor sample"))%>>%
  ggplot(aes(x=mean_depth,y=age))+
  geom_point(size=0.1)+
  facet_grid(tn ~ role)+theme_bw()+
  stat_smooth(method = lm,formula = y ~ +x,se=F)
.plot  
ggsave("by_race_anaysis/depth_age_correlation.pdf",.plot,width = 8,height = 4)
####################################################################################################################
# variantnum/coverage and age correlation
make_coverage_regression_tabel = function(.role = "TSG",.race="white",.mutype="missense"){
  regression_out = function(.class,.maf,.role,.patient_list){
    if((.class*10000) %% 1000 == 0){print(paste0("doing MAF=",.class*100))}
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .class)%>>%
      group_by(patient_id) %>>%
      summarise(MAC=sum(MAC)) %>>% ungroup() %>>%
      right_join(.patient_list,by = c("patient_id")) %>>%
      mutate(missense_num = ifelse(is.na(MAC),0,MAC)/coverage)
    #相関直線を
    lm=lm(age/365.25 ~ missense_num, data=missense_count)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"]))
  }
  if(.race!="white"){stop(paste0(".race is wrong .race=",.race,"\nwe can use white only!"))}
  .patient_list=tibble()
  .maf=tibble()
  .coverage=tibble()
  if(.role=="TSG"|.role=="oncogene"){
    truncating_patients = white_maf_for_cumulative %>>%
      filter(mutype=="truncating"|mutype=="splice") %>>%
      left_join(driver_genes %>>%dplyr::select(gene_symbol,role), by = "gene_symbol") %>>%
      filter(role==.role | role=="oncogene/TSG")%>>%
      count(patient_id) %>>% dplyr::select(-n)
    .patient_list = patient_tdg %>>%filter(race==.race)%>>%
      anti_join(truncating_patients,by="patient_id")%>>%
      left_join(coverage,by="patient_id")%>>%
      mutate(age=age/365.25,coverage=(get(.role)+`oncogene/TSG`)/1000000)%>>%
      dplyr::select(patient_id,age,coverage)
    .maf = white_maf_for_cumulative %>>%
      filter(mutype==.mutype) %>>%
      left_join(driver_genes %>>%dplyr::select(gene_symbol,role), by = "gene_symbol") %>>%
      filter(role==.role | role=="oncogene/TSG")%>>%
      dplyr::select(patient_id,MAF,MAC)
  }else if(.role=="control"){
    .patient_list = patient_cont %>>%filter(race == .race) %>>%
      left_join(coverage_cont,by="patient_id")%>>%
      mutate(age=age/365.25,coverage=coverage/1000000)%>>%
      dplyr::select(patient_id,age,coverage)
    .maf = white_maf_for_cumulative_cont %>>%
      filter(mutype==.mutype) %>>%
      inner_join(control_genes%>>%dplyr::select(gene_symbol),by="gene_symbol")%>>%
      dplyr::select(patient_id,MAF,MAC)
  }else{stop(paste0(".role error:: .role==",.role,"\n .role must be TSG, oncogene or control"))}
  tibble::tibble(MAF=1:(50*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.role,.patient_list)))%>>%
    unnest()
}
if(0){
  TSG_mis_cove_regtbl=make_coverage_regression_tabel()
  write_df(TSG_mis_cove_regtbl,"by_race_anaysis/coverage_regression/TSG_nonsyn.tsv")
  cont_mis_cove_regtbl=make_coverage_regression_tabel(.role="control")
  write_df(cont_mis_cove_regtbl,"by_race_anaysis/coverage_regression/control_nonsyn.tsv")
  cont_sil_cove_regtbl=make_coverage_regression_tabel(.role = "control",.mutype = "silent")
  write_df(cont_sil_cove_regtbl,"by_race_anaysis/coverage_regression/control_syn.tsv")
}
TSG_mis_cove_regtbl = read_tsv("by_race_anaysis/coverage_regression/TSG_nonsyn.tsv")
cont_mis_cove_regtbl = read_tsv("by_race_anaysis/coverage_regression/control_nonsyn.tsv")
cont_sil_cove_regtbl = read_tsv("by_race_anaysis/coverage_regression/control_syn.tsv")

.plot_reg_tn=regression_plot_log(TSG_mis_cove_regtbl)+theme(axis.title = element_blank())
.plot_reg_cn=regression_plot_log(cont_mis_cove_regtbl)+theme(axis.title = element_blank())
.plot_reg_cs=regression_plot_log(cont_sil_cove_regtbl)+theme(axis.title = element_blank())
.plot = cowplot::ggdraw()+
  cowplot::draw_plot(.plot_reg_tn+ggtitle("TSG nonsynonymous"),
                     x=0.03,y=0.68,width = 0.97,height =0.32 )+
  cowplot::draw_plot(.plot_reg_cn+ggtitle("Control genes nonsynonymous"),
                     x=0.03,y=0.36,width = 0.97,height =0.32 )+
  cowplot::draw_plot(.plot_reg_cs+ggtitle("Control genes nonsynonymous"),
                     x=0.03,y=0.04,width = 0.97,height =0.32 )+
  cowplot::draw_text("MAF (%)",x=0.5,y=0.02)+
  cowplot::draw_text("Regression coefficient (onset age VS variants num/covered length (Mb)",x=0.015,y=0.5,angle=90)
.plot
ggsave("by_race_anaysis/coverage_correct_regression.pdf",.plot,height = 8,width = 6)
