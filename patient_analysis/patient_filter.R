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
          !(is.na(`Primary Tumor`)&is.na(`Primary Blood Derived Cancer - Peripheral Blood`)))%>>%
  rename(bloodnorm=`Blood Derived Normal`,solidnorm=`Solid Tissue Normal`,
         tumor=`Primary Tumor`,bloodtumor=`Primary Blood Derived Cancer - Peripheral Blood`)%>>%
  mutate(race_=ifelse(race=="white" &ethnicity!="hispanic or latino","white",
                      ifelse(race=="black or african american" &ethnicity!="hispanic or latino","black","other"))) %>>%
  mutate(race_ = ifelse(is.na(race_),"other",race_)) %>>%
  dplyr::select(-race,-ethnicity)%>>%
  dplyr::rename(race=race_)%>>%
  mutate(cancer_type=str_extract(cancer_type,"[A-Z]+$"))
write_df(patient_list,"~/git/tcga_germ/variant_call/focal_patient_list.tsv")

bam_list%>>%
  dplyr::select(patient_id,file_name,id,sample_type,sample_id)%>>%
  inner_join(patient_list%>>%dplyr::select(patient_id))%>>%
  write_df("~/git/tcga_germ/variant_call/bam_file_list.tsv")



#######################################################################################################
##################################### make paient info table ##########################################
patient_info = bam_list %>>%
  count(cancer_type,patient_id,gender,race,ethnicity,age)%>>%dplyr::select(-n)%>>%
  right_join(patient_list%>>%filter(race=="white")%>%dplyr::select(patient_id))%>>%
  filter(!is.na(age))
#check
patient_info %>>% count(patient_id) %>>% filter(n>1)
#gender
patient_info %>>%count(cancer_type,gender) %>>%
  mutate(gender = ifelse(is.na(gender),"not_reported",gender)) %>>%
  bind_rows(patient_info%>>%count(gender)%>>%mutate(cancer_type="All cancer type"))%>>%
  tidyr::spread(gender,n) %>>%
#age
  left_join(patient_info %>>%filter(!is.na(age))%>>%
              group_by(cancer_type)%>>%
              summarise(`mean age of onset`=round(mean(age/365.25),digits=2),
                        `sd age of onset`=round(sd(age/365.25),digits=2))%>>%
              bind_rows(patient_info%>>%
                          summarise(`mean age of onset`=round(mean(age/365.25),digits=2),
                                    `sd age of onset`=round(sd(age/365.25),digits=2))%>>%
                          mutate(cancer_type="All cancer type")),by="cancer_type")%>>%
  left_join(patient_info%>>%count(cancer_type)%>>%
              bind_rows(patient_info%>>%summarise(n=n())%>>%
                          mutate(cancer_type="All cancer type")),by="cancer_type")%>>%
  dplyr::select(cancer_type,female,male,`mean age of onset`,`sd age of onset`,n)%>>%
  write_df("~/Dropbox/work/rare_germ/patient_info_table.tsv",na="-")

#patient list 
patient_list %>>%filter(!is.na(age))%>>%
  mutate(cancer_type=str_extract(cancer_type,"[A-Z]+$"), age = round(age/365.25,digit=2)) %>>%
  mutate_all(~ifelse(is.na(.),0,.)) %>>%filter(race=="white")%>>%
  dplyr::select(cancer_type,patient_id,gender,age,bloodnorm,solidnorm,tumor,bloodtumor)%>>%
  arrange(cancer_type,patient_id)%>>%
  write_df("~/Dropbox/work/rare_germ/patient_list_S3.tsv")


#age distribution plot
pcount=patient_list%>>%filter(!is.na(age),race=="white")%>>%
  group_by(cancer_type)%>>%summarise(n=n(),age=mean(age))%>>%
  mutate(arrange_order=row_number(age))
patient_list%>>%filter(!is.na(age),race=="white")%>>%
  left_join(pcount%>>%dplyr::select(cancer_type,arrange_order))%>>%
  mutate(age=age/365.25)%>>%
  ggplot(aes(x=reorder(cancer_type,arrange_order),y=age))+
  geom_violin(scale = "count")+
  #geom_boxplot(width=.3,fill="black")+ 
  stat_summary(fun.y=mean,geom = "point", fill="black",colour="black",shape=21,size=2)+
  geom_text(data = pcount,aes(x=cancer_type,y=10,label=n),angle=45,size=4)+
  ylab("Age of Onset")+xlab("TCGA Cancer Type")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90))
ggsave("~/Dropbox/work/rare_germ/figs/cancer_type_age_plot.pdf",width = 10,height = 4)

#################################################################################################################
########## by patient coverage plot
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv",col_types = "cciciiiic")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv",col_types = "cciciiiic")
coverage=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/coverage.tsv.gz")%>>%
  inner_join(patient_list%>>%filter(!is.na(age))%>>%dplyr::select(patient_id))
coverage_cont=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/coverage_cont.tsv.gz") %>>%
  inner_join(patient_list%>>%filter(!is.na(age))%>>%dplyr::select(patient_id))

#filter out lower than half coverage
coverage = coverage%>>%filter(all > max(all)/2)%>>%mutate(border=mean(all)-sd(all)*2) #exclude 42 patient(caucasian 9)
coverage_cont = coverage_cont%>>%filter(coverage>max(coverage)/2)%>>%mutate(border=mean(coverage)-sd(coverage)*2) #exclude 65 patient (caucasian 30)

coverage%>>%inner_join(patient_list%>>%filter(race=="white"))%>>%summarise(mean(all),sd(all)) #mean=312957, sd=10246
.plot=coverage %>>%
  inner_join(patient_list%>>%filter(race=="white"))%>>%
  ggplot()+
  geom_histogram(aes(x=all),bins = 100)+
  geom_vline(xintercept=first(coverage$border),color="red")+
  xlab("TDG Coverage (bp)")+ylab("Number of patient")+
  theme_bw()
.plot
ggsave("~/Dropbox/work/rare_germ/figs/tdg_coverage.pdf",.plot,width = 10,height = 5)

coverage_cont%>>%inner_join(patient_list%>>%filter(race=="white"))%>>%summarise(mean(coverage),sd(coverage)) #mean=517728 sd=24995
.plot=coverage_cont %>>%
  inner_join(patient_list%>>%filter(race=="white"))%>>%
  ggplot()+
  geom_histogram(aes(x=coverage),bins = 100)+
  geom_vline(xintercept = first(coverage_cont$border),color="red")+
  xlab("Control genes Coverage (bp)")+ylab("Number of patient")+
  theme_bw()
.plot
ggsave("~/Dropbox/work/rare_germ/figs/cont_coverage.pdf",.plot,width = 10,height = 5)
