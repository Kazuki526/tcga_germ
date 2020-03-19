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
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv")#%>>%filter(!is.na(age))%>>%
  mutate(age = age/365.25)
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv")


all_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/all_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")%>>%
  filter(chr!="chrX",FILTER=="PASS")
all_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/all_maf_for_cumulative_control.tsv.gz")%>>%
  filter(FILTER=="PASS")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")%>>%
  filter(FILTER=="PASS")

######################################## gnomAD (control cases data) ##############################################
tdg_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)
control_gnomad = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate(AF=AC/AN,AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj),AF_black=AC_afr/AN_afr) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC,AN,nhomalt,AF,AF_white,AF_black) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi) %>>%
  inner_join(control_genes%>>%dplyr::select(-role))

############################################################################################################################
############################################## rare variantはrace ごとに違う？ #############################################
############################################################################################################################
#################
focal_MAF=0.0005
################
mutation_num = all_maf_for_cumulative %>>%
  filter(mutype=="missense"|mutype=="silent", MAF<=focal_MAF)%>>%
  inner_join(driver_genes%>>%filter(role!="oncogene"))%>>%
  group_by(patient_id,mutype)%>>%summarise(ac=sum(MAC))%>>%ungroup()%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  right_join(patient_list%>>%dplyr::select(patient_id,race))%>>%
  mutate_all(.funs = funs(ifelse(is.na(.),0,.)))%>>%
  tidyr::gather(key = mutype,value=ac,missense,silent)%>>%
  mutate(role="TSG")%>>%
  full_join(all_maf_for_cumulative_cont%>>%
              filter(mutype=="missense"|mutype=="silent", MAF<=focal_MAF)%>>%
              group_by(patient_id,mutype)%>>%summarise(ac=sum(MAC))%>>%ungroup()%>>%
              tidyr::spread(key=mutype, value=ac)%>>%
              right_join(patient_list%>>%dplyr::select(patient_id,race))%>>%
              mutate_all(.funs = funs(ifelse(is.na(.),0,.)))%>>%
              tidyr::gather(key = mutype,value=ac,missense,silent)%>>%
              mutate(role="control"))
mutation_num%>>%filter(role=="TSG",race=="white")%>>%group_by(mutype)%>>%summarise(n=n(),ac=mean(ac))%>>%View

.plot=mutation_num %>>%
  filter(role=="TSG")%>>%inner_join(patient_tdg%>>%select(patient_id))%>>%
  ggplot(aes(x=race,y=ac))+
  geom_violin(scale="count",bw=0.5)+
  stat_summary(fun.y = mean, geom = "point",fill="black",shape=21,size=2)+
  facet_grid(. ~ mutype)+theme_bw()+ylim(c(0,20))+
  geom_signif(comparisons = list(c("black","white"),c("other","white")),
              y_position = c(19,17),test = "t.test")
.plot
ggsave(paste0("by_race_anaysis/by_race_TSGrare",focal_MAF,"_count.pdf"),.plot,width = 10,height = 8)

##じゃあonset ageは人種によって違う？？
patient_list %>>%filter((!is.na(age))) %>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  group_by(race)%>>%summarise(age=mean(age/365.25))
.plot=patient_list%>>%filter(!is.na(age))%>>%mutate(age=age/365.25)%>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  ggplot(aes(x=race,y=age))+
  geom_violin(scale="count",bw=0.5)+
  stat_summary(fun.y = mean, geom = "point",fill="black",shape=21,size=2)+theme_bw()+
  geom_signif(comparisons = list(c("black","white"),c("other","white")),
              y_position = c(105,95),test = "t.test")
.plot
ggsave("by_race_anaysis/by_race_onset_age.pdf")
if(0){
#MAF<focal_MAFのmutationが1個以上ある患者では？
patient_list%>>%
  left_join(mutation_num)%>>%filter(ac>0)%>>%
  ggplot(aes(x=race,y=age))+
  geom_violin(scale="count",bw=0.5)+
  stat_summary(fun.y = mean, geom = "point",fill="black",shape=21,size=2)+
  facet_grid(role ~ mutype)+theme_bw()
patient_list%>>%
  left_join(mutation_num)%>>%filter(ac>0)%>>%
  group_by(role,mutype,race)%>>%summarise(age=mean(age))%>>%View
}
#missense silenceの相関 (all races)
.plot=mutation_num %>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  filter(role=="TSG")%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  ggplot(aes(x=missense,y=silent))+
  geom_violin(aes(group=missense),scale = "count",bw=0.5)+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  stat_smooth(method=lm, formula = y ~  +x,se=F,colour="black")+
  #stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
  xlab("Number of Nonsynonumous Variants")+ylab("Number of Synonymous Variants")+
  #ggtitle("TSG synonymous and nonsynonymous correlation")+
  theme_bw()+#facet_grid(. ~ race)+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12))
.plot
ggsave("by_race_anaysis/TSGrarevar_correlation.pdf",.plot,width = 14,height = 6)

#人種ごとにmissense silenceの相関
.plot=mutation_num %>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  filter(role=="TSG")%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  ggplot(aes(x=missense,y=silent))+
  geom_violin(aes(group=missense),scale = "count",bw=0.5)+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  stat_smooth(method=lm, formula = y ~  +x,se=F,colour="black")+
  #stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
  xlab("Number of Nonsynonumous Variants")+ylab("Number of Synonymous Variants")+
  #ggtitle("TSG synonymous and nonsynonymous correlation")+
  theme_bw()+facet_grid(. ~ race)+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12))
.plot
ggsave("by_race_anaysis/by_race_TSGrarevar_correlation.pdf",.plot,width = 14,height = 6)

#white TSGでしっかりやると？
.plot=white_maf_for_cumulative %>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  filter(mutype=="missense"|mutype=="silent", MAF<=focal_MAF)%>>%
  inner_join(driver_genes%>>%filter(role!="oncogene"))%>>%
  group_by(patient_id,mutype)%>>%summarise(ac=sum(MAC))%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  right_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = funs(ifelse(is.na(.),0,.)))%>>%
  ggplot(aes(x=missense,y=silent))+
  geom_violin(aes(group=missense),scale = "count",bw=0.5)+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
  stat_smooth(method=lm, formula = y ~  +x,se=F,colour="black")+
  #stat_summary(fun.y=mean,geom = "point", fill="white",shape=21,size=2)
  xlab("Number of Nonsynonumous Variants")+ylab("Number of Synonymous Variants")+
  #ggtitle("TSG synonymous and nonsynonymous correlation")+
  theme_bw()+
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12))
.plot
ggsave("by_race_anaysis/white_correlation.pdf",.plot,width = 8,height = 4)
white_maf_for_cumulative %>>%
  inner_join(patient_tdg%>>%select(patient_id))%>>%
  filter(mutype=="missense"|mutype=="silent", MAF<=focal_MAF)%>>%
  inner_join(driver_genes%>>%filter(role!="oncogene"))%>>%
  group_by(patient_id,mutype)%>>%summarise(ac=sum(MAC))%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  right_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = funs(ifelse(is.na(.),0,.)))%>>%
  {lm(silent ~ missense, data=.)}

white_maf_for_cumulative %>>%
  filter(mutype=="missense"|mutype=="silent", MAF<=focal_MAF)%>>%
  inner_join(driver_genes%>>%filter(role!="oncogene"))%>>%
  group_by(patient_id,mutype)%>>%summarise(ac=sum(MAC))%>>%
  tidyr::spread(key=mutype, value=ac)%>>%
  right_join(patient_list%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = funs(ifelse(is.na(.),0,.)))%>>%
  {cor.test(.$missense,.$silent,method = "pearson")}
