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
patient_tdg = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_list_forTGD.tsv",col_types = "cciciiiic")
white_maf_for_cumulative = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/white_maf_for_cumulative.tsv.gz")
patient_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/patient_list_forcont.tsv",col_types = "cciciiiic")
white_maf_for_cumulative_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/white_maf_for_cumulative_control.tsv.gz")
if(0){
##################################### gnomAD data #############################################
tdg_gnomad_white = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_top_driver_gene.maf",
                      col_types=cols(chr="c",LoF_filter="c"))%>>%
  mutate( AC_white=AC_fin+AC_nfe+AC_asj, AN_white=AN_fin+AN_nfe+AN_asj, AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj)) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC_white,AN_white,AF_white,Consequence) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)%>>%
  classify_consequence()%>>%
  #filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
  left_join(tdg_bed) %>>%
  filter((bed_start-3<=start & start<=bed_end+2)|(bed_start-3<=start+str_length(ref)-1 & start+str_length(ref)-1 <=bed_end+2),chr!="chrX")%>>%
  dplyr::select(-bed_start,-bed_end)

control_gnomad_white = read_tsv("/Volumes/DR8TB2/gnomAD/maf38/non_cancer_maf/non_cancer_control_gene.maf")%>>%
  mutate( AC_white=AC_fin+AC_nfe+AC_asj, AN_white=AN_fin+AN_nfe+AN_asj, AF_white=(AC_fin+AC_nfe+AC_asj)/(AN_fin+AN_nfe+AN_asj)) %>>%
  dplyr::select(chr,posi,ref,alt,filter,SYMBOL,AC_white,AN_white,AF_white,Consequence) %>>%
  dplyr::rename(gene_symbol =SYMBOL,start = posi)%>>%
  inner_join(control_genes%>>%dplyr::select(-role))%>>%
  classify_consequence()%>>%
  #filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
  left_join(cont_bed) %>>%
  filter((bed_start-3<=start & start<=bed_end+2)|(bed_start-3<=start+str_length(ref)-1 & start+str_length(ref)-1 <=bed_end+2),chr!="chrX")%>>%
  dplyr::select(-bed_start,-bed_end)
###################################### 1000 genomes data ####################################
  #for QC
  norm_maf_all = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/norm_maf_all.tsv.gz")
  norm_maf_all_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/norm_maf_all_cont.tsv.gz")
  duplicate_site = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/duplicate_site.tsv")
  duplicate_site_cont = read_tsv("/Volumes/areca42TB2/gdc/control_region/all_patient/duplicate_site_control_gnomad.tsv") 
  kb_around = function(.start){
    .head=.start-1001
    .tail=.start+1000
    data.frame(start=.head:.tail)
  }
  duplicate_site_all= full_join(duplicate_site,duplicate_site_cont)%>>%
    filter(FDR_cancer<0.01)%>>%dplyr::select(chr,start)%>>%
    mutate(start = purrr::map(start,~kb_around(.))) %>>%
    unnest() %>>%distinct()
  
  somatic_recurrent_all = 
    norm_maf_all %>>%
    filter(soma_or_germ=="somatic",LOH=="no")%>>%
    count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
    dplyr::rename(alt=t_allele2)%>>%filter(n>10) %>>%
    mutate(recurrent_focal="yes") %>>%dplyr::select(-n)%>>%
    full_join(norm_maf_all_cont%>>%
                filter(soma_or_germ=="somatic",LOH=="no")%>>%
                count(gene_symbol,chr,start,end,ref,t_allele2)%>>%
                dplyr::rename(alt=t_allele2)%>>%filter(n>10) %>>%
                mutate(recurrent_focal="yes") %>>% dplyr::select(-n))%>>%
    dplyr::select(chr,start,ref,alt)
  #for 1000genomes maf data
  .colnames = c('sample_id','allele','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele',
                'Tumor_Seq_Allele2','Consequence','PolyPhen',"cDNA_position","CDS_position","Protein_position",'Transcript_ID')
  .cols = .colnames %>>%
  {setNames(c('c','c','c','c','d','d','c','c', 'c','c','c','c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types=.cols) %>>%
      classify_consequence() %>>%
      filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")
  }
  maf_1kg_all=strip_maf("/Volumes/areca42TB2/1000genomes/GRCh38/top_driver_genes/white_by_sample_top_driver_genes.maf.gz")%>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
    mutate(chr=paste0("chr",chr))
  write_df(maf_1kg_all,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/1kg_white_maf_all.tsv.gz")
  maf_1kg_cont=strip_maf("/Volumes/areca42TB2/1000genomes/GRCh38/control_genes//white_by_sample_control_genes.maf.gz")%>>%
    dplyr::rename(gene_symbol=Hugo_Symbol,chr=Chromosome,start=Start_Position,end=End_Position,
                  ref=Reference_Allele,alt=Tumor_Seq_Allele2) %>>%
    mutate(chr=paste0("chr",chr))
  write_df(maf_1kg_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/1kg_white_maf_all.tsv.gz")
  sample_1kg_list=count(maf_1kg_all,sample_id)%>>%dplyr::select(-n)
  ref_minor_convert=function(.data){
    sample_1kg_list%>>%mutate(allele1="a",allele2="a")%>>%
      tidyr::pivot_longer(cols=c("allele1","allele2"),names_to="allele",values_to="a")%>>%dplyr::select(-a)%>>%
      left_join(.data,by = c("sample_id", "allele"))%>>%
      mutate(MAC=ifelse(is.na(end),1,0))%>>%filter(MAC==1)%>>%
      dplyr::select(sample_id,allele,MAC)
  }
  tdg_1kg=tdg_gnomad_white%>>%filter(AF_white>0.5)%>>%
    dplyr::select(chr,start,ref,alt,AF_white)%>>%
    inner_join(maf_1kg_all)%>>%group_by(chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%nest()%>>%
    mutate(ref_minor=purrr::map(data,~ref_minor_convert(.)))%>>%dplyr::select(-data)%>>%unnest()%>>%
    bind_rows(maf_1kg_all%>>%anti_join(tdg_gnomad_white%>>%filter(AF_white>0.5)%>>%
                                         dplyr::select(chr,start,ref,alt,AF_white))%>>%
                left_join(tdg_gnomad_white%>>%dplyr::select(chr,start,ref,alt,AF_white))%>>%
                dplyr::select(sample_id,allele,chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%
                mutate(MAC=1))%>>%
    group_by(sample_id,chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%summarise(MAC=sum(MAC))%>>%
    mutate(MAF=ifelse(AF_white>0.5,1-AF_white,AF_white))
  control_1kg=control_gnomad_white%>>%filter(AF_white>0.5)%>>%
    dplyr::select(chr,start,ref,alt,AF_white)%>>%
    inner_join(maf_1kg_cont)%>>%group_by(chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%nest()%>>%
    mutate(ref_minor=purrr::map(data,~ref_minor_convert(.)))%>>%dplyr::select(-data)%>>%unnest()%>>%
    bind_rows(maf_1kg_cont%>>%anti_join(control_gnomad_white%>>%filter(AF_white>0.5)%>>%
                                         dplyr::select(chr,start,ref,alt,AF_white))%>>%
                left_join(control_gnomad_white%>>%dplyr::select(chr,start,ref,alt,AF_white))%>>%
                dplyr::select(sample_id,allele,chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%
                mutate(MAC=1))%>>%
    group_by(sample_id,chr,start,ref,alt,gene_symbol,mutype,AF_white)%>>%summarise(MAC=sum(MAC))%>>%
    mutate(MAF=ifelse(AF_white>0.5,1-AF_white,AF_white))
  #QC
  tdg_1kg%>>%anti_join(duplicate_site_all)%>>%filter(gene_symbol!="KMT2C")%>>%
    anti_join(somatic_recurrent_all)%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/1kg_white_maf.tsv.gz")
  control_1kg%>>%anti_join(duplicate_site_all)%>>%filter(gene_symbol!="KMT2C")%>>%
    anti_join(somatic_recurrent_all)%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/control_gene/1kg_white_maf.tsv.gz")
  # QC gnomAD
  tdg_gnomad_white%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
    anti_join(duplicate_site_all)%>>%filter(gene_symbol!="KMT2C")%>>%
    anti_join(somatic_recurrent_all)%>>%
    mutate(MAF=ifelse(AF_white>0.5,1-AF_white,AF_white))%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/gnomAD_white_maf.tsv.gz")
  control_gnomad_white%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
    anti_join(duplicate_site_all)%>>%filter(gene_symbol!="KMT2C")%>>%
    anti_join(somatic_recurrent_all)%>>%
    mutate(MAF=ifelse(AF_white>0.5,1-AF_white,AF_white))%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/control_gene/gnomAD_white_maf.tsv.gz")
  #variant count by patient 
  sample_num_cancer_tdg=count(patient_tdg%>>%filter(race=="white"))$n
  sample_num_cancer_cont=count(patient_cont%>>%filter(race=="white"))$n
  tdg_cancer = white_maf_for_cumulative %>>%inner_join(patient_tdg)%>>%inner_join(driver_genes)%>>%
    dplyr::select(patient_id,gene_symbol,chr,start,ref,alt,mutype,MAF,MAC)
  control_cancer = white_maf_for_cumulative_cont %>>%inner_join(patient_cont)%>>%inner_join(control_genes)%>>%
    dplyr::select(patient_id,gene_symbol,chr,start,ref,alt,mutype,MAF,MAC)
  
  # inconsistent site between Database
  #30 site and cancer patients have no site 
  bind_rows(tdg_1kg,control_1kg) %>>%group_by(chr,start,ref,alt,gene_symbol,mutype)%>>%
    summarise(MAC_1kg=sum(MAC),MAF_1kg=sum(MAC)/sample_num_1kg/2)%>>%
    full_join(white_maf_for_cumulative %>>%group_by(chr,start,ref,alt,gene_symbol,mutype,MAF)%>>%
                summarise(MAC=sum(MAC),MAF_cancer=sum(MAC)/sample_num_cancer_tdg/2)%>>%
                bind_rows(white_maf_for_cumulative_cont %>>%group_by(chr,start,ref,alt,gene_symbol,mutype,MAF)%>>%
                            summarise(MAC=sum(MAC),MAF_cancer=sum(MAC)/sample_num_cancer_cont/2)))%>>%
    mutate(MAF=ifelse(is.na(MAF),0,MAF),
           MAF_cancer=ifelse(is.na(MAF_cancer),0,MAF_cancer),
           MAF_1kg=ifelse(is.na(MAF_1kg),0,MAF_1kg))%>>%
    filter((MAF<.MAF & abs(MAF_1kg-MAF)>0.01)|(MAF<.MAF & abs(MAF-MAF_cancer)>0.01))%>>%(?.)%>>%
    dplyr::select(chr,start,ref,alt)%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/inconsistent_site.tsv")
}
inconsistent_site=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/inconsistent_site.tsv")
sample_num_1kg=503
tdg_1kg=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/1kg_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
control_1kg=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene//1kg_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
tdg_gnomad_white=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/gnomAD_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)
control_gnomad_white=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/gnomAD_white_maf.tsv.gz")%>>%anti_join(inconsistent_site)

######################################## variants count #############################################
#variant count by patient 
sample_num_cancer_tdg=count(patient_tdg%>>%filter(race=="white"))$n
sample_num_cancer_cont=count(patient_cont%>>%filter(race=="white"))$n
tdg_cancer = white_maf_for_cumulative %>>%inner_join(patient_tdg)%>>%inner_join(driver_genes)%>>%
  dplyr::select(patient_id,gene_symbol,chr,start,ref,alt,mutype,MAF,MAC)%>>%anti_join(inconsistent_site)
control_cancer = white_maf_for_cumulative_cont %>>%inner_join(patient_cont)%>>%inner_join(control_genes)%>>%
  dplyr::select(patient_id,gene_symbol,chr,start,ref,alt,mutype,MAF,MAC)%>>%anti_join(inconsistent_site)
count_variant=function(.MAF){
  if((.MAF*10000) %% 100 == 0){print(paste0("doing MAF=",.MAF))}
  tibble(cancer=tdg_cancer%>>%filter(MAF<=.MAF)%>>%{sum(.$MAC)/sample_num_cancer_tdg},
         kg = tdg_1kg%>>%filter(MAF<=.MAF)%>>%{sum(.$MAC)/sample_num_1kg},
         gnomad=tdg_gnomad_white%>>%filter(MAF<=.MAF)%>>%{sum(.$AF_white)})%>>%
    tidyr::pivot_longer(cols = c("cancer","kg","gnomad"),names_to="Database",values_to="Average number of variants")%>>%
    mutate(role="Top driver genes")%>>%
    bind_rows(tibble(cancer=control_cancer%>>%filter(MAF<=.MAF)%>>%{sum(.$MAC)/sample_num_cancer_cont},
                     kg = control_1kg%>>%filter(MAF<=.MAF)%>>%{sum(.$MAC)/sample_num_1kg},
                     gnomad=control_gnomad_white%>>%filter(MAF<=.MAF)%>>%{sum(.$AF_white)})%>>%
                tidyr::pivot_longer(cols = c("cancer","kg","gnomad"),names_to="Database",values_to="Average number of variants")%>>%
                mutate(role="Control genes"))
}
if(0){
  count_variantions_10 = tibble::tibble(MAF=1:1000) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(variants = purrr::map(MAF,~count_variant(.)))
  count_variantions_10_50 = tibble::tibble(MAF=100:500) %>>%
    mutate(MAF = MAF/1000) %>>%
    mutate(variants = purrr::map(MAF,~count_variant(.)))
  bind_rows(count_variantions_10,count_variantions_10_50)%>>%unnest()%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/counts_variant.tsv")
}
##### compare variants nums by database ####
variant_dif=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/counts_variant.tsv")%>>%
  mutate(Database = ifelse(Database=="kg","1000 genomes",
                           ifelse(Database=="cancer","TCGA","gnomAD non-cancer")))%>>%
  mutate(Database = factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")))
.plot_all=variant_dif%>>%
  mutate(MAF=MAF*100,facet=factor(role,levels=c("Top driver genes","Control genes")))%>>%
  ggplot()+
  geom_point(aes(x=MAF,y=`Average number of variants`,color=Database))+
  facet_wrap(~facet,scales = "free")+
  theme_bw()
.legend=cowplot::get_legend(.plot_all)
.plot5tdg=variant_dif%>>%mutate(MAF=MAF*100)%>>%
  filter(MAF<1, role=="Top driver genes")%>>%
  ggplot()+
  geom_point(aes(x=MAF,y=`Average number of variants`,color=Database))+
  theme_bw()+theme(axis.title = element_blank(),legend.position = "none")
.plot5cont=variant_dif%>>%mutate(MAF=MAF*100)%>>%
  filter(MAF<1, role=="Control genes")%>>%
  ggplot()+
  geom_point(aes(x=MAF,y=`Average number of variants`,color=Database))+
  theme_bw()+theme(axis.title = element_blank(),legend.position = "none")
.plot=cowplot::plot_grid(.plot_all+theme(axis.title = element_blank(),legend.position = "none"),
                         cowplot::plot_grid(.plot5tdg,NULL,.plot5cont,.legend,
                                            nrow=1,rel_widths = c(3,2,3,2)),
                         nrow=2)
cowplot::ggdraw()+
  cowplot::draw_plot(.plot,x=0.05,y=0.05,height = 0.95,width = 0.95)+
  cowplot::draw_text("MAF cutoff (%)",x=0.5,y=0.025)+
  cowplot::draw_text("Average number of variants",x=0.025,y=0.5,angle=90)
ggsave("~/Dropbox/work/rare_germ/compare_analysis/variants_num_dif.pdf",height = 6,width = 8)
################################### compare TDG variants nums corrected by CG  ####################################
################################### compare TDG variants nums corrected by CG  ####################################
################################### compare TDG variants nums corrected by CG  ####################################
driver_genes_2type = driver_genes %>>%
  mutate(role=ifelse(role=="oncogene/TSG","TSG",role))%>>%
  bind_rows(driver_genes%>>%filter(role=="oncogene/TSG")%>>%mutate(role="oncogene"))
sample_tbl=
  white_maf_for_cumulative%>>%inner_join(patient_tdg%>>%dplyr::select(patient_id))%>>%
  filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
  mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
  inner_join(driver_genes_2type)%>>%
  group_by(patient_id,mutype,role)%>>%
  summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC")%>>%
  right_join(patient_tdg%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
  tidyr::pivot_longer(col=-patient_id,names_to = c("role","mutype"),values_to = "MAC",names_sep = "_")%>>%
  mutate(Database = "TCGA")%>>%
  bind_rows(white_maf_for_cumulative_cont%>>%inner_join(patient_cont%>>%dplyr::select(patient_id))%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              inner_join(control_genes)%>>%
              group_by(patient_id,mutype)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
              right_join(patient_tdg%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
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
              dplyr::rename(patient_id=sample_id)%>>%mutate(Database="1000 genomes",role="control"))
#work flow 
work_flow_tbl=sample_tbl %>>%mutate(role=ifelse(role=="control","Control genes","Top driver genes"))%>>%
  group_by(Database,patient_id,role)%>>%summarise(MAC=sum(MAC))%>>%ungroup()%>>%
  group_by(Database,role)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  bind_rows(tibble(Database="gnomAD non-cancer",role=c("Top driver genes","Control genes"),
                   mean=c(tdg_gnomad_white%>>%filter(mutype!="inframe_indel",AF_white<.MAF)%>>%{sum(.$AF_white)*2},
                          control_gnomad_white%>>%filter(mutype!="inframe_indel",AF_white<.MAF)%>>%{sum(.$AF_white)*2})))%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(role,levels=c("Top driver genes","Control genes")))
.plot_befcor = ggplot(work_flow_tbl,aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  scale_fill_manual(values=c(`Top driver genes`="#d95f02",`Control genes`="#7fc97f"))+
  theme_bw()+ylab("Average of individual variants number")
.plot_corect = work_flow_tbl %>>%group_by(Database)%>>%
  mutate(sd= sd/max(mean),mean=mean/max(mean))%>>%ungroup()%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  scale_fill_manual(values=c(`Top driver genes`="#d95f02",`Control genes`="#7fc97f"))+
  theme_bw()+geom_hline(yintercept = 1,color="red")+
  ylab("Correction value")

.legend=cowplot::get_legend(.plot_corect)
.plot_workflow=cowplot::plot_grid(.plot_befcor+theme(legend.position = "none",axis.title.x = element_blank(),
                                                     axis.text.x = element_text(angle = -15,vjust = 0.5)),NULL,
                                  .plot_corect+theme(legend.position = "none",axis.title.x = element_blank(),
                                                     axis.text.x = element_text(angle = -15,vjust = 0.5)),.legend,
                                  nrow=1,rel_widths = c(2.1,1,2,1))
.plot_workflow

.plot_tsg=sample_tbl %>>%filter(role!="oncogene")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  bind_rows(tdg_gnomad_white%>>%inner_join(driver_genes_2type%>>%filter(role=="TSG"))%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              group_by(role,mutype)%>>%summarise(mean=sum(AF_white))%>>%
              full_join(control_gnomad_white%>>%inner_join(control_genes)%>>%
                          filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
                          mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
                          group_by(role,mutype)%>>%summarise(mean=sum(AF_white)))%>>%
              mutate(Database="gnomAD non-cancer"))%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(ifelse(role=="TSG",role,"Control genes"),levels=c("TSG","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  group_by(Database,mutype)%>>%
  mutate(sd= sd/max(mean),mean=mean/max(mean))%>>%ungroup()%>>%
  ggplot(aes(x=Database,y=mean,fill=role))+
  geom_bar(stat= "identity",position="dodge")+
  geom_errorbar(aes(ymin=ifelse(mean-sd<0,0,mean-sd), ymax=mean + sd, width=0.2),position=position_dodge(width = 0.9))+
  facet_grid(.~mutype)+
  facet_grid(.~mutype)+scale_fill_manual(values=c(TSG="#beaed4",`Control genes`="#7fc97f"))+
  theme_bw()+geom_hline(yintercept = 1,color="red")
.plot_oncg=sample_tbl %>>%filter(role!="TSG")%>>%
  group_by(Database,role,mutype)%>>%summarise(mean=mean(MAC),sd=sd(MAC))%>>%ungroup()%>>%
  bind_rows(tdg_gnomad_white%>>%inner_join(driver_genes_2type%>>%filter(role=="oncogene"))%>>%
              filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              group_by(role,mutype)%>>%summarise(mean=sum(AF_white))%>>%
              full_join(control_gnomad_white%>>%inner_join(control_genes)%>>%
                          filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
                          mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
                          group_by(role,mutype)%>>%summarise(mean=sum(AF_white)))%>>%
              mutate(Database="gnomAD non-cancer"))%>>%
  mutate(Database=factor(Database,levels = c("TCGA","1000 genomes","gnomAD non-cancer")),
         role=factor(ifelse(role=="oncogene","Oncogene","Control genes"),levels=c("Oncogene","Control genes")),
         mutype=factor(ifelse(mutype=="missense","nonsynonymous",ifelse(mutype=="silent","synonymous",mutype)),
                       levels=c("truncating","nonsynonymous","synonymous")))%>>%ungroup()%>>%
  group_by(Database,mutype)%>>%
  mutate(sd= sd/max(mean),mean=mean/max(mean))%>>%ungroup()%>>%
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
.plot_compare
.plot=cowplot::plot_grid(.plot_workflow,.plot_compare,nrow = 2,rel_heights = c(1,2))
.plot
ggsave("~/Dropbox/work/rare_germ/compare_analysis/compare_variants.pdf",.plot,height = 10,width = 8)


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


###### permutation test #######
control_variants = sample_tbl%>>%
  filter(role=="control")%>>%
  group_by(Database,mutype)%>>%summarise(control_MAC=mean(MAC))
observed_data=sample_tbl %>>%filter(role!="control")%>>%left_join(control_variants)%>>%
  mutate(MAC_corrected=MAC/control_MAC)%>>%dplyr::select(-MAC,-control_MAC) %>>%
  group_by(Database,role,mutype)%>>%summarise(MAC_corrected=mean(MAC_corrected))%>>%ungroup()%>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%(?.)%>>%
  dplyr::select(-Database)%>>%summarise_all(~diff(.))%>>%
  tidyr::pivot_longer(cols=-NULL,names_to = c("role","mutype"),values_to = "observed_dif",names_sep = "_")
.sample_tbl=sample_tbl %>>%filter(role!="control")%>>%left_join(control_variants)%>>%
  mutate(MAC_corrected=MAC/control_MAC)%>>%dplyr::select(-MAC,-control_MAC) %>>%
  tidyr::pivot_wider(names_from = c("role","mutype"),values_from = "MAC_corrected")%>>%ungroup()
permutation=function(.times){
  if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
  .sample_tbl%>>%dplyr::select(-patient_id)%>>%
    mutate(Database=sample(Database,length(Database)))%>>%
    group_by(Database)%>>%summarise_all(~mean(.))%>>%ungroup()%>>%
    dplyr::select(-Database)%>>%summarise_all(~diff(.))
}

permutation_tbl=tibble(times=1:10000)%>>%
  mutate(result=purrr::map(times,~permutation(.)))%>>%unnest()
permutation_tbl%>>%
  pivot_longer(cols=-times,names_to = c("role","mutype"),values_to = "perm_dif",names_sep = "_") %>>%
  left_join(observed_data)%>>%
  filter(perm_dif<observed_dif)%>>%
  count(role,mutype)
#1 TSG      missense    9218
#2 TSG      silent      2671
#3 TSG      truncating  3018
#4 oncogene missense    8253
#5 oncogene silent      4359
#6 oncogene truncating  8103
if(0){
white_maf_for_cumulative%>>%inner_join(patient_cont%>>%dplyr::select(patient_id))%>>%
  filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
  mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
  group_by(patient_id,mutype)%>>%
  summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
  right_join(patient_tdg%>>%filter(race=="white")%>>%dplyr::select(patient_id))%>>%
  mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
  tidyr::pivot_longer(col=-patient_id,names_to = "mutype",values_to = "MAC")%>>%
  mutate(Database = "TCGA")%>>%
  bind_rows(tdg_1kg%>>%filter(mutype!="inframe_indel",MAF<=.MAF)%>>%
              mutate(mutype=ifelse(mutype=="splice","truncating",mutype))%>%
              group_by(sample_id,mutype)%>>%
              summarise(MAC=sum(MAC))%>>%tidyr::pivot_wider(names_from = "mutype",values_from = "MAC")%>>%
              right_join(tdg_1kg%>>%count(sample_id)%>>%dplyr::select(sample_id)) %>>%
              mutate_all(.funs = ~ifelse(is.na(.),0,.))%>>%
              tidyr::pivot_longer(col=-sample_id,names_to = "mutype",values_to = "MAC")%>>%
              dplyr::rename(patient_id=sample_id)%>>%mutate(Database="1000 genomes"))%>>%
  left_join(control_variant_num)%>>%
  mutate(MAC_correct=MAC/average_num, facet=factor(mutype,levels = c("truncating","missnse", "silent")))%>>%
  ggplot(aes(x=Database,y=MAC_correct))+
  geom_violin(scale="count",bw=0.3)+
  stat_summary(fun.y = mean, geom = "point",fill="black",shape=21,size=2)+
  facet_wrap(~facet)+theme_bw()+
  geom_signif(comparisons = list(c("1000 genomes","TCGA")),y_position = 20,test = "t.test")
}
