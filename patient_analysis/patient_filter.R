library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd('~/git/tcga_germ/patient_analysis/')
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

bam_list = read_tsv("response.tsv")%>>%
  rename(race = cases.0.demographic.race, ethnicity = cases.0.demographic.ethnicity,
         gender = cases.0.demographic.gender, age = cases.0.diagnoses.0.age_at_diagnosis,
         cancer_type = cases.0.project.project_id, sample_type = cases.0.samples.0.sample_type,
         patient_id = cases.0.submitter_id,sample_id = cases.0.samples.0.submitter_id)%>>%
  filter(state == "released")

patient_list = bam_list %>>%
  count(patient_id,cancer_type,sample_type,age,gender,race,ethnicity)%>>%
  tidyr::spread(key=sample_type,value=n)%>>%
  filter(!(is.na(`Blood Derived Normal`) &is.na(`Solid Tissue Normal`)),
          !is.na(`Primary Tumor`), is.na(`Primary Blood Derived Cancer - Peripheral Blood`))%>>%
  rename(bloodnorm=`Blood Derived Normal`,solidnorm=`Solid Tissue Normal`,tumor=`Primary Tumor`)%>>%
  mutate(race_=ifelse(race=="white" &ethnicity!="hispanic or latino","white",
                      ifelse(race=="black or african american" &ethnicity!="hispanic or latino","black","other"))) %>>%
  mutate(race_ = ifelse(is.na(race_),"other",race_)) %>>%
  dplyr::select(-race,-ethnicity,-`Primary Blood Derived Cancer - Peripheral Blood`)%>>%
  dplyr::rename(race=race_)
write_df(patient_list,"~/git/tcga_germ/variant_call/focal_patient_list.tsv")

bam_list%>>%
  dplyr::select(patient_id,file_name,id,sample_type,sample_id)%>>%
  inner_join(focal_patient%>>%dplyr::select(patient_id))%>>%
  write_df("~/git/tcga_germ/variant_call/bam_file_list.tsv")



#######################################################################################################
##################################### make paient info table ##########################################
patient_info = bam_list %>>%
  count(cancer_type,patient_id,gender,race,ethnicity,age)%>>%dplyr::select(-n)%>>%
  right_join(patient_list%>>%dplyr::select(patient_id))%>>%
  filter(!is.na(age))
#check
patient_info %>>% count(patient_id) %>>% filter(n>1)
#gender
patient_info %>>%count(cancer_type,gender) %>>%
  mutate(gender = ifelse(is.na(gender),"not_reported",gender)) %>>%
  tidyr::spread(gender,n) %>>%
#race
  left_join(patient_info %>>%mutate(race=ifelse(is.na(race),"not reported",race))%>>%
              count(cancer_type,race)%>>%
              tidyr::spread(race,n),by="cancer_type")%>>%
  
#ethnicity
  left_join(patient_info %>>%mutate(ethnicity=ifelse(ethnicity=="not reported",NA,ethnicity))%>>%
              count(cancer_type,ethnicity)%>>%
              tidyr::spread(ethnicity,n),by="cancer_type")%>>%
#age
  left_join(patient_info %>>%filter(!is.na(age))%>>%
              group_by(cancer_type)%>>%
              summarise(`mean age of oncet`=mean(age/365.25)),by="cancer_type")%>>%
  left_join(patient_info%>>%count(cancer_type),by="cancer_type")%>>%
  dplyr::select(cancer_type,female,male,`american indian or alaska native`,asian,`black or african american`,
                `native hawaiian or other pacific islander`,white,`not reported`,`hispanic or latino`,`not hispanic or latino`,
                `<NA>`,`mean age of oncet`,n)%>>%
  write_df("~/Dropbox/work/rare_germ/patient_info_table.tsv",na="-")

#patient list 
patient_list %>>%
  mutate(cancer_type=str_extract(cancer_type,"[A-Z]+$"), age = age/365.25,
         Caucasian=ifelse(race=="white","yes","-")) %>>%
  mutate_all(~ifelse(is.na(.),0,.)) %>>%dplyr::select(-race)%>>%
  left_join(patient_info%>>%dplyr::select(patient_id,race,ethnicity))%>>%
  dplyr::select(cancer_type,patient_id,gender,race,ethnicity,Caucasian,age)%>>%
  arrange(cancer_type,patient_id)%>>%
  write_df("~/Dropbox/work/rare_germ/patient_list_S3.tsv")


#age distribution plot
pcount=patient_list%>>%filter(!is.na(age))%>>%
  group_by(cancer_type)%>>%summarise(n=n(),age=mean(age))%>>%
  mutate(arrange_order=row_number(age))
patient_list%>>%filter(!is.na(age))%>>%
  left_join(pcount%>>%dplyr::select(cancer_type,arrange_order))%>>%
  mutate(age=age/365.25)%>>%
  ggplot(aes(x=reorder(cancer_type,arrange_order),y=age))+
  geom_violin(scale = "count")+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",colour="black",shape=21,size=2)+
  geom_text(data = pcount,aes(x=cancer_type,y=10,label=n),angle=90,size=8)+
  ylab("Age of Onset")+xlab("TCGA Cancer Type")+
  theme_bw()
ggsave("~/Dropbox/work/rare_germ/figs/cancer_type_age_plot.pdf",width = 20,height = 8)

#################################################################################################################
########## by patient coverage plot
coverage=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/coverage.tsv.gz")%>>%
  inner_join(patient_list%>>%filter(!is.na(age))%>>%dplyr::select(patient_id))
coverage_cont=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/coverage_cont.tsv.gz") %>>%
  inner_join(patient_list%>>%filter(!is.na(age))%>>%dplyr::select(patient_id))

#filter out lower than half coverage
coverage = coverage%>>%filter(all > max(all)/2)%>>%mutate(border=mean(all)-sd(all)*2) #exclude 42 patient
coverage_cont = coverage_cont%>>%filter(coverage>max(coverage)/2)%>>%mutate(border=mean(coverage)-sd(coverage)*2) #exclude 65 patient

coverage%>>%summarise(mean(all),sd(all)) #mean=313193, sd=10439
.plot=coverage %>>%
  ggplot()+
  geom_histogram(aes(x=all),bins = 100)+
  geom_vline(xintercept=first(coverage$border),color="red")+
  xlab("TDG Coverage (bp)")+ylab("Number of patient")+
  theme_bw()
.plot
ggsave("~/Dropbox/work/rare_germ/figs/tdg_coverage.pdf",.plot,width = 10,height = 5)

coverage_cont%>>%summarise(mean(coverage),sd(coverage)) #mean=518595 sd=26156
.plot=coverage_cont %>>%
  ggplot()+
  geom_histogram(aes(x=coverage),bins = 100)+
  geom_vline(xintercept = first(coverage_cont$border),color="red")+
  xlab("Control genes Coverage (bp)")+ylab("Number of patient")+
  theme_bw()
.plot
ggsave("~/Dropbox/work/rare_germ/figs/cont_coverage.pdf",.plot,width = 10,height = 5)
