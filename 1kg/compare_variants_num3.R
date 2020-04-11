library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
loadNamespace('cowplot')
setwd('/Volumes/DR8TB2/tcga_rare_germ/')
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
.MAF=0.0005
######################################## gene infomation ##########################################
driver_genes=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/driver_genes.tsv")
tdg_bed = read_tsv("~/git/driver_genes/onlytop105/top_driver105exon.bed",col_names = c("chr","bed_start","bed_end", "info","none","nonen"))%>>%
  mutate(gene_symbol=str_extract(info,"^[^:]+"))%>>%dplyr::select(chr,bed_start,bed_end,gene_symbol)
control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")
cont_bed = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_exon.bed",
                    col_names = c("chr","bed_start","bed_end", "info","none","nonen"))%>>%
  mutate(gene_symbol=str_extract(info,"^[^:]+"))%>>%dplyr::select(chr,bed_start,bed_end,gene_symbol)
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
######################################## cancer patient data #######################################
patient_hicov = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list_exclude_low_coverage.tsv",col_types = "cciciiiic")
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
###################################### 1000 genomes and gnomAD #######################################
inconsistent_site=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/inconsistent_site.tsv")
sample_num_1kg=503
tdg_1kg=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/1kg_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
control_1kg=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene//1kg_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
tdg_gnomad_white=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/gnomAD_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
control_gnomad_white=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/gnomAD_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)

################################### compare TDG variants nums corrected by CG  ####################################
driver_genes_2type = driver_genes %>>%
  mutate(role=ifelse(role=="oncogene/TSG","TSG",role))%>>%
  bind_rows(driver_genes%>>%filter(role=="oncogene/TSG")%>>%mutate(role="oncogene"))
sample_tbl=
  white_maf_for_cumulative%>>%inner_join(patient_hicov%>>%dplyr::select(patient_id))%>>%
  filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
  mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
  inner_join(driver_genes_2type)%>>%
  group_by(patient_id,mutype,role)%>>%
  summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC")%>>%
  right_join(patient_hicov%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
  tidyr::pivot_longer(col=-patient_id,names_to = c("role","mutype"),values_to = "MAC",names_sep = "_")%>>%
  mutate(Database = "TCGA")%>>%
  bind_rows(white_maf_for_cumulative_cont%>>%inner_join(patient_hicov%>>%dplyr::select(patient_id))%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              inner_join(control_genes)%>>%
              group_by(patient_id,mutype)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
              right_join(patient_hicov%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
              mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
              tidyr::pivot_longer(col=-patient_id,names_to = "mutype",values_to = "MAC")%>>%
              mutate(Database = "TCGA",role="control"))%>>%
  bind_rows(tdg_1kg%>>%filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              inner_join(driver_genes_2type)%>>%
              group_by(sample_id,mutype,role)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC")%>>%
              right_join(tdg_1kg%>>%count(sample_id)%>>%dplyr::select(sample_id)) %>>%
              mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
              tidyr::pivot_longer(col=-sample_id,names_to = c("role","mutype"),values_to = "MAC",names_sep = "_")%>>%
              dplyr::rename(patient_id=sample_id)%>>%mutate(Database="1000 genomes"))%>>%
  bind_rows(control_1kg%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              inner_join(control_genes)%>>%
              group_by(sample_id,mutype)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
              right_join(tdg_1kg%>>%count(sample_id)%>>%dplyr::select(sample_id)) %>>%
              mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
              tidyr::pivot_longer(col=-sample_id,names_to = "mutype",values_to = "MAC")%>>%
              dplyr::rename(patient_id=sample_id)%>>%mutate(Database="1000 genomes",role="control"))%>>%
  ungroup()
# before correct
.plot_tsg_bef=sample_tbl %>>%filter(role!="oncogene")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(ifelse(role=="TSG",role,"Control genes"),levels=c("TSG","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  facet_wrap(.~mutype,scales = "free")+scale_fill_manual(values=c(TSG="#beaed4",`Control genes`="#7fc97f"))+
  theme_bw()
.plot_oncg_bef=sample_tbl %>>%filter(role!="TSG")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes")),
         role=factor(ifelse(role=="oncogene","Oncogene","Control genes"),levels=c("Oncogene","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  facet_wrap(.~mutype,scales = "free")+scale_fill_manual(values=c(Oncogene="#fdc086",`Control genes`="#7fc97f"))+
  theme_bw()
.plot_bef=cowplot::ggdraw()+
  cowplot::draw_plot(.plot_tsg_bef+theme(axis.title = element_blank(),axis.text.x = element_text(angle = -15,vjust = 0.5)),
                     x=0.04,y=0.5,width=0.96,height=0.5)+
  cowplot::draw_plot(.plot_oncg_bef+theme(axis.title = element_blank(),axis.text.x = element_text(angle = -15,vjust = 0.5)),
                     x=0.04,y=0,width=0.96,height=0.5)+
  cowplot::draw_text("Average number of variants",x=0.02,angle=90)
# corrected
.plot_tsg=sample_tbl %>>%filter(role!="oncogene")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes")),
         role=factor(ifelse(role=="TSG",role,"Control genes"),levels=c("TSG","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  group_by(Database,mutype)%>>%
  mutate(sd= sd/max(mean),mean=mean/max(mean))%>>%ungroup()%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  facet_grid(.~mutype)+scale_fill_manual(values=c(TSG="#beaed4",`Control genes`="#7fc97f"))+
  theme_bw()+geom_hline(yintercept = 1,color="red")
.plot_oncg=sample_tbl %>>%filter(role!="TSG")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes")),
         role=factor(ifelse(role=="oncogene","Oncogene","Control genes"),levels=c("Oncogene","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  group_by(Database,mutype)%>>%
  mutate(sd= sd/max(mean),mean=mean/max(mean))%>>%ungroup()%>>%                                filter(Database!="gnomAD non-cancer")%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  facet_grid(.~mutype)+scale_fill_manual(values=c(Oncogene="#fdc086",`Control genes`="#7fc97f"))+
  theme_bw()+geom_hline(yintercept = 1,color="red")
.plot_compare=cowplot::ggdraw()+
  cowplot::draw_plot(.plot_tsg+theme(axis.title = element_blank(),axis.text.x = element_text(angle = -15,vjust = 0.5)),
                     x=0.04,y=0.5,width=0.96,height=0.5)+
  cowplot::draw_plot(.plot_oncg+theme(axis.title = element_blank(),axis.text.x = element_text(angle = -15,vjust = 0.5)),
                     x=0.04,y=0,width=0.96,height=0.5)+
  cowplot::draw_text("Correction value of average number of variants",x=0.02,angle=90)
.plot=cowplot::plot_grid(.plot_bef,NULL,.plot_compare,ncol = 1,rel_heights = c(1,0.2,1),labels = c("a",NA,"b"))
.plot
ggsave("~/Dropbox/work/rare_germ/compare_analysis/compare_variants3.pdf",.plot,height = 12,width = 7)

######### distribution ##########################
.control_ave=sample_tbl%>>%filter(role=="control")%>>%
  group_by(Database,mutype)%>>%summarise(cont_ave=mean(MAC)) %>>%ungroup
sample_tbl %>>%filter(role=="TSG")%>>%
  left_join(.control_ave)%>>%mutate(corrected_MAC=MAC/cont_ave) %>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%
  ggplot(aes(x=corrected_MAC))+
  geom_bar(aes(y=..prop..))+
  facet_grid(Database ~ mutype,scales = "free")+
  xlab("Correction value number of variants")+ylab("Proportion")+
  theme_bw()
ggsave("~/Dropbox/work/rare_germ/compare_analysis/MAC_correct_distribution.pdf",height = 4,width = 8)




#for table
sample_tbl %>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  bind_rows(tdg_gnomad_white%>>%inner_join(driver_genes_2type)%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              group_by(role,mutype)%>>%summarise(mean=sum(AF_white))%>>%
              full_join(control_gnomad_white%>>%inner_join(control_genes)%>>%
                          filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
                          mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
                          group_by(role,mutype)%>>%summarise(mean=sum(AF_white)))%>>%
              mutate(Database="gnomAD non-cancer"))%>>%(?.)%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(ifelse(role=="control","Control genes",role),levels=c("TSG","oncogene","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%
  arrange(Database,role,mutype)%>>%
  write_df("~/Dropbox/work/rare_germ/compare_analysis/compare_variants.tsv")
sample_tbl %>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  bind_rows(tdg_gnomad_white%>>%inner_join(driver_genes_2type)%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              group_by(role,mutype)%>>%summarise(mean=sum(AF_white))%>>%
              full_join(control_gnomad_white%>>%inner_join(control_genes)%>>%
                          filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
                          mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
                          group_by(role,mutype)%>>%summarise(mean=sum(AF_white)))%>>%
              mutate(Database="gnomAD non-cancer"))%>>%(?.)%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(ifelse(role=="control","Control genes",role),levels=c("TSG","oncogene","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%
  arrange(Database,role,mutype)%>>%
  mutate(mean_sd=paste0(round(mean,digits=3)," (",round(sd,digits=3),")"))%>>%dplyr::select(-mean,-sd)%>>%
  tidyr::pivot_wider(names_from = "mutype",values_from = "mean_sd")%>>%
  write_df("~/Dropbox/work/rare_germ/compare_analysis/compare_variants_spread.tsv")


###### permutation test #######
control_variants = sample_tbl%>>%
  filter(role=="control")%>>%
  group_by(Database,mutype)%>>%summarise(control_MAC=mean(MAC))
permutation=function(.times,.sample_tbl){
  if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
  .sample_tbl%>>%dplyr::select(-patient_id)%>>%
    mutate(Database=sample(Database,length(Database)))%>>%
    group_by(Database)%>>%summarise_all(~mean(.))%>>%ungroup()%>>%
    dplyr::select(-Database)%>>%summarise_all(~diff(.))
}
######### permutation doin on crrected value ##########
.sample_tbl=sample_tbl %>>%filter(role!="control")%>>%left_join(control_variants)%>>%
  mutate(MAC_corrected=MAC/control_MAC)%>>%dplyr::select(-MAC,-control_MAC) %>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%ungroup()
permutation_tbl=tibble(times=1:10000)%>>%
  mutate(result=purrr::map(times,~permutation(.,.sample_tbl)))%>>%unnest()
observed_data=sample_tbl %>>%filter(role!="control")%>>%left_join(control_variants)%>>%
  mutate(MAC_corrected=MAC/control_MAC)%>>%dplyr::select(-MAC,-control_MAC) %>>%
  group_by(Database,role,mutype)%>>%summarise(MAC_corrected=mean(MAC_corrected))%>>%ungroup()%>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%(?.)%>>%
  dplyr::select(-Database)%>>%summarise_all(~diff(.))%>>%
  tidyr::pivot_longer(cols=-NULL,names_to = c("role","mutype"),values_to = "observed_dif",names_sep = "_")
permutation_tbl%>>%
  pivot_longer(cols=-times,names_to = c("role","mutype"),values_to = "perm_dif",names_sep = "_") %>>%
  left_join(observed_data)%>>%
  filter(perm_dif<observed_dif)%>>%
  count(role,mutype)
#1 TSG      missense    8691
#2 TSG      silent      2036
#3 TSG      truncating  2307
#4 oncogene missense    7866
#5 oncogene silent      3794
#6 oncogene truncating  8010


########## permutation doin on crrected value ###########
.sample_tbl_raw=sample_tbl %>>%filter(role!="control")%>>%rename(MAC_corrected=MAC)%>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%ungroup()
permutation_tbl_raw=tibble(times=1:10000)%>>%
  mutate(result=purrr::map(times,~permutation(.,.sample_tbl_raw)))%>>%unnest()
observed_data_raw=sample_tbl %>>%filter(role!="control")%>>%rename(MAC_corrected=MAC)%>>%
  group_by(Database,role,mutype)%>>%summarise(MAC_corrected=mean(MAC_corrected))%>>%ungroup()%>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%(?.)%>>%
  dplyr::select(-Database)%>>%summarise_all(~diff(.))%>>%
  tidyr::pivot_longer(cols=-NULL,names_to = c("role","mutype"),values_to = "observed_dif",names_sep = "_")
permutation_tbl_raw%>>%
  pivot_longer(cols=-times,names_to = c("role","mutype"),values_to = "perm_dif",names_sep = "_") %>>%
  left_join(observed_data_raw)%>>%
  filter(perm_dif<observed_dif)%>>%
  count(role,mutype)
#1 TSG      missense    9613
#2 TSG      silent      4109
#3 TSG      truncating  9459
#4 oncogene missense    9029
#5 oncogene silent      5811
#6 oncogene truncating  7955

###########################################################################################################
############ by gene ###########
count_by_gene = function(.tbl,.database="tcga"){
  if(.database=="tcga"){
    patient_hicov%>>%filter(race=="white")%>>%dplyr::select(patient_id)%>>%
      left_join(.tbl)%>>%mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
      summarise(`TCGA truncating`   =paste0(round(mean(truncating),digits=3),"(",round(sd(truncating),digits=3),")"),
                `TCGA nonsynonymous`=paste0(round(mean(missense),digits=3),  "(",round(sd(missense),digits=3),  ")"),
                `TCGA synonymous`   =paste0(round(mean(silent),digits=3),    "(",round(sd(silent),digits=3),    ")"))
  }else{
    tdg_1kg%>>%count(sample_id)%>>%dplyr::select(sample_id)%>>%
      left_join(.tbl)%>>%mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
      summarise(`1000 genomes truncating`   =paste0(round(mean(truncating),digits=3),"(",round(sd(truncating),digits=3),")"),
                `1000 genomes nonsynonymous`=paste0(round(mean(missense),digits=3),  "(",round(sd(missense),digits=3),  ")"),
                `1000 genomes synonymous`   =paste0(round(mean(silent),digits=3),    "(",round(sd(silent),digits=3),    ")"))
  }
}
sample_gene_tbl =
  white_maf_for_cumulative%>>%inner_join(patient_hicov%>>%dplyr::select(patient_id))%>>%
  filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
  mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
  inner_join(driver_genes)%>>%
  group_by(patient_id,mutype,gene_symbol,role)%>>%
  summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
  ungroup()%>>%nest(patient_id, silent, missense, truncating)%>>%
  mutate(tcga_count=purrr::map(data,~count_by_gene(.)))%>>%dplyr::select(-data)%>>%
  unnest()%>>%
  left_join(tdg_1kg%>>%filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              inner_join(driver_genes_2type)%>>%
              group_by(sample_id,mutype,role,gene_symbol)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
              ungroup()%>>%nest(sample_id, silent, missense, truncating)%>>%
              mutate(kg_count=purrr::map(data,~count_by_gene(.,.database="1kg")))%>>%dplyr::select(-data)%>>%
              unnest())%>>%
  left_join(driver_genes)%>>%
  mutate_all(~ifelse(is.na(.),"0(0)",.))%>>%
  arrange(role,gene_symbol)
write_df(sample_gene_tbl,"~/Dropbox/work/rare_germ/compare_analysis/compare_variants_spread_bygene.tsv")
  