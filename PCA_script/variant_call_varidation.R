library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
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

######################################## gene infomation ##########################################
driver_genes=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/driver_genes.tsv")
control_genes = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_genes.tsv")
patient_list = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list.tsv",col_types = "cciciiiic")
############################################# datas manipurate #############################################
if(0){
  classify_consequence = function(.data) {
    mutate(.data,
           mutype= dplyr::recode(Consequence,
                                 downstream_gene_variant = 'flank',
                                 `3_prime_UTR_variant` = 'flank',
                                 `3_prime_UTR_variant,NMD_transcript_variant`='NMD',
                                 upstream_gene_variant = 'flank',
                                 `5_prime_UTR_variant` = 'flank',
                                 frameshift_variant = 'truncating',
                                 frameshift_variant = 'truncating',
                                 inframe_deletion = 'inframe_indel',
                                 inframe_insertion = 'inframe_indel',
                                 intron_variant = 'intron',
                                 splice_region_variant = 'splice_region',
                                 coding_sequence_variant = 'missense',
                                 missense_variant = 'missense',
                                 stop_gained = 'truncating',
                                 stop_lost = 'truncating',
                                 stop_retained_variant = 'silent',
                                 synonymous_variant = 'silent',
                                 splice_acceptor_variant = 'splice',
                                 splice_donor_variant = 'splice',
                                 protein_altering_variant = 'inframe_indel',
                                 start_lost = 'truncating',
                                 `coding_sequence_variant,3_prime_UTR_variant`='flank',
                                 `coding_sequence_variant,5_prime_UTR_variant`='flank',
                                 `stop_gained,frameshift_variant` = 'truncating',
                                 `stop_gained,frameshift_variant,start_lost` = 'truncating',
                                 `missense_variant,splice_region_variant`='missense',
                                 `intron_variant,non_coding_transcript_variant`='intron',
                                 `intron_variant,NMD_transcript_variant`='NMD',
                                 `non_coding_transcript_exon_variant,non_coding_transcript_variant`='non_coding',
                                 `frameshift_variant,splice_region_variant`='truncating',
                                 `frameshift_variant,start_lost`='truncating',
                                 `frameshift_variant,stop_lost`='truncating',
                                 `inframe_deletion,splice_region_variant`='splice_region',
                                 `inframe_insertion,splice_region_variant`='splice_region',
                                 `frameshift_variant,stop_lost,splice_region_variant`='truncating',
                                 `frameshift_variant,stop_retained_variant`='inframe_indel',
                                 `protein_altering_variant,splice_region_variant`='inframe_indel',
                                 `splice_acceptor_variant,intron_variant`='splice',
                                 `splice_acceptor_variant,coding_sequence_variant`='splice',
                                 `splice_acceptor_variant,coding_sequence_variant,intron_variant`='splice',
                                 `splice_acceptor_variant,5_prime_UTR_variant`='splice',
                                 `splice_acceptor_variant,5_prime_UTR_variant,intron_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,intron_variant`='splice',
                                 `splice_donor_variant,intron_variant`='splice',
                                 `splice_donor_variant,3_prime_UTR_variant,intron_variant`='splice',
                                 `splice_donor_variant,5_prime_UTR_variant`='splice',
                                 `splice_donor_variant,coding_sequence_variant,3_prime_UTR_variant`='splice',
                                 `splice_donor_variant,splice_acceptor_variant,intron_variant`='splice',
                                 `splice_region_variant,5_prime_UTR_variant`='flank',
                                 `splice_region_variant,3_prime_UTR_variant`='flank',
                                 `splice_region_variant,intron_variant` = 'splice_region',
                                 `splice_region_variant,synonymous_variant`='silent',
                                 `splice_region_variant,stop_retained_variant`='silent',
                                 `start_lost,inframe_deletion`='truncating',
                                 `start_lost,splice_region_variant`='truncating',
                                 `stop_gained,start_lost`='truncating',
                                 `stop_lost,inframe_deletion`='truncating',
                                 `stop_lost,splice_region_variant`='truncating',
                                 `stop_gained,protein_altering_variant`='truncating',
                                 `stop_gained,frameshift_variant,splice_region_variant`='truncating',
                                 `stop_gained,inframe_deletion`='truncating',
                                 `stop_gained,inframe_insertion`='truncating',
                                 `stop_gained,splice_region_variant`='truncating'))
  }
  ct_exchange = tibble(cancer_type=c("BRCA","COAD","READ","GBM","LGG","KICH","KIRC","KIRP","LUAD","LUSC"),
                       bp=c("breast","colorectal","colorectal","brain","brain","kidney","kidney","kidney","lung","lung"))
  ct_exchange_cont = tibble(cancer_type=c("COAD","READ","KICH","KIRC","KIRP"),
                       bp=c("crc","crc","kcc","kcc","kcc"))
  .colnames = c('Hugo_Symbol','Chromosome','Start_Position','Reference_Allele',"Match_Norm_Seq_Allele1",
                "Match_Norm_Seq_Allele2",'Consequence',"FILTER")
  .cols = .colnames %>>%
  {setNames(c('c','c','d','c','c', 'c','c','c'), .)} %>>%
  {do.call(readr::cols_only, as.list(.))}
  strip_maf = function(infile) {
    read_tsv(infile, comment='#', col_types = .cols) %>>%
      classify_consequence()
  }
  #blood driverd + solid derived normal
  patient_list%>>%filter(!is.na(age),race=="white")%>>%
    left_join(ct_exchange)%>>%mutate(bp=ifelse(is.na(bp),tolower(cancer_type),bp))%>>%
    mutate(filename=paste('/Volumes/areca42Tb2/gdc/top_driver_gene/',bp,"/maf/",patient_id,".maf",sep="")) %>>%
    dplyr::select(patient_id,cancer_type,filename)%>>%
    mutate(maf=purrr::map(filename,~strip_maf(.))) %>>%
    dplyr::select(-filename)%>>%
    unnest() %>>%#(?.)%>>%
    dplyr::rename(gene=Hugo_Symbol,chr=Chromosome,posi=Start_Position,ref=Reference_Allele,
                  allele1=Match_Norm_Seq_Allele1,allele2=Match_Norm_Seq_Allele2) %>>%
    inner_join(driver_genes%>>%dplyr::select(gene_symbol),by=c("gene"="gene_symbol")) %>>%
    filter(!(ref==allele1 & ref==allele2))%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
    mutate(chr = ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
    filter(chr!="chrX")%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/variantcall_check/top_driver.maf.gz")
  patient_list%>>%filter(!is.na(age),race=="white")%>>%
    left_join(ct_exchange_cont)%>>%mutate(bp=ifelse(is.na(bp),tolower(cancer_type),bp))%>>%
    mutate(filename=paste('/Volumes/areca42Tb2/gdc/control_region/',bp,"/maf/",patient_id,".maf",sep="")) %>>%
    dplyr::select(patient_id,cancer_type,filename)%>>%
    mutate(maf=purrr::map(filename,~strip_maf(.))) %>>%
    dplyr::select(-filename)%>>%
    unnest() %>>%#(?.)%>>%
    dplyr::rename(gene=Hugo_Symbol,chr=Chromosome,posi=Start_Position,ref=Reference_Allele,
                  allele1=Match_Norm_Seq_Allele1,allele2=Match_Norm_Seq_Allele2) %>>%
    inner_join(control_genes%>>%dplyr::select(gene_symbol),by=c("gene"="gene_symbol")) %>>%
    filter(!(ref==allele1 & ref==allele2))%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent") %>>%
    mutate(chr = ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
    filter(,chr!="chrX")%>>%
    write_df("/Volumes/DR8TB2/tcga_rare_germ/control_gene/variantcall_check/control_region.maf.gz")
  
#extract Lu et al. cell 2018 data
  Lu_maf = read_tsv("/Volumes/areca42TB2/gdc/pancan_atlas_germ/maf_patient/lifover_maf/top_driver_genes_patient_hg38.maf.gz",
                    col_types = cols(chr="c"))%>>%filter(chr!="X")%>>%
    classify_consequence()%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
    mutate(patient_id=str_extract(patient_id,"^TCGA-..-...."),
           chr=ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
    inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))
  write_df(Lu_maf,"/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/variantcall_check/top_driver_Lu2018cell_maf.tsv.gz")
  Lu_maf_cont = read_tsv("/Volumes/areca42TB2/gdc/pancan_atlas_germ/maf_patient/lifover_maf/control_genes_patient_hg38.maf.gz",
                         col_types = cols(chr="c"))%>>%filter(chr!="X")%>>%
    classify_consequence()%>>%
    filter(mutype=="inframe_indel" | mutype=="truncating"| mutype=="splice"| mutype=="missense"| mutype=="silent")%>>%
    mutate(patient_id=str_extract(patient_id,"^TCGA-..-...."),
           chr=ifelse(str_detect(chr,"^chr"),chr,paste0("chr",chr)))%>>%
    inner_join(patient_list%>>%filter(race=="white")%>>%dplyr::select(patient_id))
  write_df(Lu_maf_cont,"/Volumes/DR8TB2/tcga_rare_germ/control_gene/variantcall_check/control_region_Lu2018cell_maf.tsv.gz")
  patient_for_vccheck=Lu_maf%>>%count(patient_id)%>>%dplyr::select(-n)%>>%
    full_join(Lu_maf_cont%>>%count(patient_id)%>>%dplyr::select(-n))
  write_df(patient_for_vccheck,"/Volumes/DR8TB2/tcga_rare_germ/patient_list_for_variant_check.tsv")
}
patient_for_vccheck=read_tsv("/Volumes/DR8TB2/tcga_rare_germ/patient_list_for_variant_check.tsv")
Lu_maf = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/variantcall_check/top_driver_Lu2018cell_maf.tsv.gz")
Lu_maf_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/variantcall_check/control_region_Lu2018cell_maf.tsv.gz")
our_maf = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/variantcall_check/top_driver.maf.gz")%>>%filter(FILTER=="PASS")
our_maf_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/variantcall_check/control_region.maf.gz")%>>%filter(FILTER=="PASS")

####################################################################################
#compare in same region
tdg_bed = read_tsv("~/git/driver_genes/onlytop105/top_driver105exon.bed",col_names = c("chr","start","end", "info","none","nonen"))%>>%
  mutate(gene=str_extract(info,"^[A-Z0-9]+"))%>>%dplyr::select(chr,start,end,gene)
cont_bed = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/control_gene_exon.bed",
                    col_names = c("chr","start","end", "info","none","nonen"))%>>%
  mutate(gene=str_extract(info,"^[A-Z0-9]+"))%>>%dplyr::select(chr,start,end,gene)
Lu_maf_exon=Lu_maf%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
  left_join(tdg_bed)%>>%
  filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end)%>>%
  bind_rows(Lu_maf_cont%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
              left_join(cont_bed)%>>%
              filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end))%>>%
  dplyr::select(patient_id,chr,posi,ref,allele1,allele2,filter,gene,mutype)%>>%
  tidyr::pivot_longer(cols = c(allele1,allele2),names_to = "allele",values_to = "alt")%>>%filter(ref!=alt)%>>%
  mutate(mutation_class=ifelse(ref=="-"|alt=="-","indel","SNV"))
#before QC 
if(0){
our_maf_exon = our_maf%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
  left_join(tdg_bed)%>>%
  filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end)%>>%
  bind_rows(our_maf_cont%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
              left_join(cont_bed)%>>%
              filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end))%>>%
  inner_join(patient_for_vccheck) %>>% dplyr::select(patient_id,chr,posi,ref,allele1,allele2,FILTER,gene,mutype)%>>%
  tidyr::pivot_longer(cols = c(allele1,allele2),names_to = "allele",values_to = "alt")%>>%filter(ref!=alt)%>>%
  mutate(mutation_class=ifelse(ref=="-"|alt=="-","indel","SNV"))

pre_QC=Lu_maf_exon %>>%full_join(our_maf_exon)#%>>%filter(allele=="allele2")
pre_QC%>>%mutate(filter=ifelse(is.na(filter),NA,"called"))%>>%
  mutate(FILTER=ifelse(is.na(FILTER),NA,"called"))%>>%
  count(mutation_class,filter,FILTER)
#1 SNV            called called 4897829
#2 SNV            called NA      473330
#3 SNV            NA     called   43738
#4 indel          called called   65862
#5 indel          called NA       32334
#6 indel          NA     called   12614
##maybe indel calling was silipped, so we check whithout position infomation
pre_QC_indel=Lu_maf_exon%>>%filter(mutation_class=="indel")%>>%dplyr::select(-posi) %>>%
  full_join(our_maf_exon%>>%filter(mutation_class=="indel")%>>%dplyr::select(-posi) )
pre_QC_indel%>>%mutate(filter=ifelse(is.na(filter),NA,"called"))%>>%
  mutate(FILTER=ifelse(is.na(FILTER),NA,"called"))%>>%#inner_join(driver_genes,by=c("gene"="gene_symbol"))%>>%
  count(mutation_class,filter,FILTER)
#1 indel          called called 77214
#2 indel          called NA     20982
#3 indel          NA     called  1272
}
### QC ######
QC_maf_norm = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/QC_norm_maf_all.tsv.gz")%>>%inner_join(patient_for_vccheck)
QC_maf_norm_cont = read_tsv("/Volumes/DR8TB2/tcga_rare_germ/control_gene/QC_norm_maf_all_cont.tsv.gz")%>>%inner_join(patient_for_vccheck)
our_maf_QC = our_maf %>>%inner_join(QC_maf_norm%>>%rename(gene=gene_symbol,posi=start)%>>%
                         dplyr::select(patient_id,gene,chr,posi,ref,alt))
our_maf_cont_QC = our_maf_cont %>>%inner_join(QC_maf_norm_cont%>>%rename(gene=gene_symbol,posi=start)%>>%
                                    dplyr::select(patient_id,gene,chr,posi,ref,alt))

our_maf_exon_QC = our_maf_QC%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
  left_join(tdg_bed)%>>%
  filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end)%>>%
  bind_rows(our_maf_cont_QC%>>%mutate(posi_end=posi+str_length(ref)-1)%>>%
              left_join(cont_bed)%>>%
              filter((start-3<=posi & posi<=end+2)|(start-3<=posi_end & posi_end<=end+2))%>>%dplyr::select(-start,-end))%>>%
  dplyr::select(patient_id,chr,posi,ref,allele1,allele2,FILTER,gene,mutype)%>>%
  tidyr::pivot_longer(cols = c(allele1,allele2),names_to = "allele",values_to = "alt")%>>%filter(ref!=alt)%>>%
  mutate(mutation_class=ifelse(ref=="-"|alt=="-","indel","SNV"))
QC_maf=Lu_maf_exon %>>%full_join(our_maf_exon_QC)
QC_maf%>>%mutate(filter=ifelse(is.na(filter),NA,"called"))%>>%
  mutate(FILTER=ifelse(is.na(FILTER),NA,"called"))%>>%
  count(mutation_class,filter,FILTER)
"1 SNV            called called 4258334
 2 SNV            called NA     1112825
 3 SNV            NA     called   21936
 4 indel          called called   29580
 5 indel          called NA       68616
 6 indel          NA     called     548"
##maybe indel calling was silipped, so we check whithout position infomation
QC_indel=Lu_maf_exon%>>%filter(mutation_class=="indel")%>>%dplyr::select(-posi) %>>%
  full_join(our_maf_exon_QC%>>%filter(mutation_class=="indel")%>>%dplyr::select(-posi) )
QC_indel%>>%#mutate(filter=ifelse(is.na(filter),NA,"called"))%>>%
  mutate(FILTER=ifelse(is.na(FILTER),NA,"called"))%>>%#inner_join(driver_genes,by=c("gene"="gene_symbol"))%>>%
  count(mutation_class,filter,FILTER)
#1 indel          called called 29710
#2 indel          called NA     68486
#3 indel          NA     called   418


#########################################################################################
# how about pathogenic variants?
Lu_maf_tdg37=read_tsv("/Volumes/areca42TB2/gdc/pancan_atlas_germ/maf_patient/top_driver_genes_patient.tsv.gz",
                      col_types = cols(chr="c"))
ps_site37=read_tsv("/Volumes/areca42TB/tcga/burden_plot/cell2018/supply_data/pathogenic_site.tsv")%>>%
  dplyr::rename(chr=Chromosome,posi=Start,ref=Reference,alt=Alternate,gene=HUGO_Symbol,cancer_type=cancer,
                clinvar=ClinVar_Pathogenicity, acmg=ACMG_Classification, charger=CharGer_Classification,
                trait=ClinVar_Traits) %>>%
  #filter(clinvar=="Pathogenic"|charger=="Pathogenic"|charger=="Likely Pathogenic")%>>%
  inner_join(driver_genes%>>%dplyr::select(gene_symbol),by=c("gene"="gene_symbol"))

Lu_maf_tdg37 %>>% mutate(patient_id=str_extract(patient_id,"TCGA-..-...."),gene=ifelse(is.na(gene),Consequence,gene)) %>>%
  #left_join(patient_list%>>%dplyr::select(patient_id,cancer_type))%>>%
  right_join(ps_site37%>>%mutate(AD=str_match(Genotype,"^./.:([0-9]+,[0-9]+):")[,2])%>>%
               dplyr::select(chr,posi,ref,alt,cancer_type,gene,AD,charger))%>>%dplyr::select(-posi,-chr)%>>%
  inner_join(patient_for_vccheck)%>>%
  inner_join(Lu_maf)%>>%(?.)%>>% #231 variants
  left_join(our_maf_exon_QC)%>>%filter(is.na(FILTER)) #8 variants

Lu_maf_tdg37 %>>% mutate(patient_id=str_extract(patient_id,"TCGA-..-...."),gene=ifelse(is.na(gene),Consequence,gene)) %>>%
  #left_join(patient_list%>>%dplyr::select(patient_id,cancer_type))%>>%
  right_join(ps_site37%>>%mutate(AD=str_match(Genotype,"^./.:([0-9]+,[0-9]+):")[,2])%>>%
               dplyr::select(chr,posi,ref,alt,cancer_type,gene,AD,charger))%>>%dplyr::select(-posi,-chr)%>>%
  inner_join(patient_for_vccheck)%>>%
  inner_join(Lu_maf)%>>%(?.)%>>% #231 variants
  left_join(our_maf_exon_QC)%>>%filter(!is.na(FILTER)) %>>%
  dplyr::select(cancer_type,patient_id)%>>%mutate(significance="significant")%>>%
  write_df("/Volumes/DR8TB2/tcga_rare_germ/top_driver_gene/patient_with_ps.tsv")


