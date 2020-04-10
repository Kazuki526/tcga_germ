#####################################################################################################################
#violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
truncating_filling = function(.truncating_count){
  .truncating_count = .truncating_count %>>%
    dplyr::select(patient_id,cancer_type,age,truncating_count_n)
  .max_count=max(.truncating_count$truncating_count_n)
  .outcome = tibble(truncating_count_n=0:.max_count) %>>%
    left_join(.truncating_count%>>%count(truncating_count_n)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$truncating_count_n) != 0){
    .outcome = .outcome %>>%
      mutate(age=-50, cancer_type="BRCA", patient_id="no_patient")
    return(bind_rows(.truncating_count,.outcome))
  }else{return(.truncating_count)}
}
#指数表記
exponent_notation = function(.num){
  .log=trunc(log10(.num))
  .log=ifelse(.log>0,.log+1,.log-1)
  paste0(.num*(10^-.log)," %*% 10^",.log)
}
#####################################################################################################################
perm_pvalue = function(.regression,.CT=NA,.tbl,.tail){
  if(!is.na(.CT)){.tbl=filter(.tbl,cancer_type==.CT)}
  if(.tail == "one"){
    .tbl = .tbl %>>%
      filter(regression < .regression )
  }else if(.tail == "two"){
    .tbl = .tbl %>>%
      filter(regression > abs(.regression))
  }
  length(.tbl$regression) /10000
}

truncate_plot_allcantype= function(.tbl,.permu=T,.test_tail="one",.permu_do=F,.permu_file=NA){
  .tbl=.tbl%>>%filter(!is.na(age))
  .max_count=max(.tbl$truncating_count_n)
  .coef_posi=ifelse(.max_count==1,2.25,ifelse(.max_count==2,3,.max_count))
  lm=lm(age ~ truncating_count_n, data=.tbl)
  regression=as.data.frame(as.list(coef(lm)))
  trunc_permu = function(.times,.tbl){
    if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
    .tbl_sample=.tbl%>>%mutate(age=sample(age,length(age)))
    as.tibble(as.list(coef(lm(age ~ truncating_count_n,data = .tbl_sample))))
  }
  if(.permu){
    if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
      regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
    }else{
      regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
        mutate(tbl=purrr::map(times,~trunc_permu(.times = .,.tbl = .tbl))) %>>%
        unnest() %>>%
        rename(regression = truncating_count_n)
      write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
    }
    .p=perm_pvalue(.tbl = regression_tbl,
                   .regression = regression$truncating_count_n,
                   .tail = .test_tail)
    regression = mutate(regression,p_value = .p)%>>%(?.)%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",p_value))))
  }else{
    regression = regression %>>%
      mutate(p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                             summary(lm)$fstatistic["dendf"]))) %>>%(?.)%>>%
      mutate(p_value=ifelse(p_value>0.5,1-p_value,p_value))%>>%
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  }
  regression=regression%>>%mutate(abintercept=X.Intercept.-truncating_count_n)
  .plot = .tbl %>>%truncating_filling()%>>%(?.)%>>%
    mutate(cancer_type ="All Cancer Types") %>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    geom_text(data =.tbl %>>%count(truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=5,label=n),size=3,position="stack")+
    geom_abline(data = regression %>>%filter(p_value > 0.05),
                aes(intercept = abintercept,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05),
                aes(intercept = abintercept,slope = truncating_count_n))+
    geom_text(data = regression, aes(x=.coef_posi,y=10,label=out_reg),size=5,hjust=1)+
    facet_wrap( ~ cancer_type)+
    xlab("Number of Truncated Gene")+ylab("Age at Diagnosis")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}

##############################################################################################
truncate_plot_bycantype = function(.tbl,.permu=T,.test_tail="one",.permu_do=F,.permu_file=NA){
  .tbl=.tbl%>>%filter(!is.na(age))
  .max_count=max(.tbl$truncating_count_n)
  .coef_posi=ifelse(.max_count == 1,2.4,.max_count+1)
  lm_p=function(.data,.permu){
    if(max(.data$truncating_count_n) == 0){
      .dtbl = tibble::tibble(X.Intercept. =0, truncating_count_n =0)
      if(.permu){return(.dtbl)}else{return(mutate(.dtbl,p_value=1))}
    }else{
      lm=lm(age ~ truncating_count_n, data=.data)
      .dtbl = as.data.frame(as.list(coef(lm)))
      if(.permu){return(.dtbl)
      }else{
        return(mutate(.dtbl,p_value=(1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                                            summary(lm)$fstatistic["dendf"]))))
      }
    }
  }
  trunc_permu_byCT = function(.times,.tbl){
    if(.times %% 1000 == 0){print(paste0("permutation ",.times," times now"))}
    .tbl%>>%mutate(age=sample(age,length(age)))%>>%
      nest(-cancer_type)%>>%dplyr::rename(data_=data)%>>%
      mutate(data_ = purrr::map(data_,~lm_p(.,.permu = T)))%>>%unnest()
  }
  regression = .tbl %>>%
    tidyr::nest(-cancer_type) %>>%
    mutate(data = purrr::map(data,~lm_p(.,.permu = .permu))) %>>%
    unnest()
  if(.permu){#if do permutation
    if(file.exists(paste0("age_plot/permute_tbl/",.permu_file))){
      regression_tbl = read_tsv(paste0("age_plot/permute_tbl/",.permu_file))
    }else{
      regression_tbl=tibble::tibble(times=seq(1,10000,by=1)) %>>%
        mutate(data = purrr::map(times,~trunc_permu_byCT(.times =.,.tbl = .tbl))) %>>%
        unnest() %>>%
        rename(regression = truncating_count_n)
      write_df(regression_tbl,paste0("age_plot/permute_tbl/",.permu_file))
    }
    regression = regression %>>%
      mutate(p_value = pmap_dbl(.,function(cancer_type,truncating_count_n,...){
        perm_pvalue(.regression =truncating_count_n,.CT=cancer_type,
                    .tbl = regression_tbl,.tail = .test_tail)})) %>>%(?.)%>>%
      mutate(p_value=ifelse(p_value>0.5,1-p_value,p_value))%>>%
      mutate(out_reg =ifelse(truncating_count_n==0,"",
               paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2))))))
  }
  regression=regression%>>%mutate(abintercept=X.Intercept.-truncating_count_n)
  .plot = .tbl %>>%truncating_filling()%>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    facet_wrap( ~ cancer_type, ncol = 5)+
    geom_text(data =.tbl %>>%count(cancer_type,truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=2.5,label=n),size=2.5,position="stack",angle=45,hjust=0)+
    geom_abline(data = regression %>>%filter(p_value > 0.05,X.Intercept. >0),
                aes(intercept = abintercept,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05,X.Intercept. >0),
                aes(intercept = abintercept,slope = truncating_count_n))+
    geom_text(data = regression, aes(x=.coef_posi,y=15,label =out_reg),size=3,hjust=1)+
    xlab("Number of truncating and splice variants")+ylab("Age at Diagnosis")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}
#################################################################################################################
#0.01%ごとにregressionしてみる
make_regression_tabel_truncate = function(.maf=all_maf_for_cumulative,.vcf=tdg_gnomad,.race="all",.role = "TSG",
                                 .fdr=0.01,.max_maf=50,.maf_filter=F,
                                 .database="all",.duplicate=T,.somatic=T,.varscan=T,.pathogenic = F,
                                 .patient_list = patient_hicov,.remove0=F){
  regression_out = function(.class,.maf,.role,.patient_list){
    if((.class*10000) %% 100 == 0){print(paste0("doing MAF=",.class*100))}
    ##missense の数
    missense_count = .maf %>>%
      filter(MAF <= .class)%>>%
      group_by(patient_id) %>>%
      summarise(missense_num=sum(MAC)) %>>% ungroup() %>>%
      {left_join(.patient_list,.,by = c("patient_id"))} %>>%
      mutate(missense_num = ifelse(is.na(missense_num),0,missense_num))
    #相関直線を
    lm=lm(age/365.25 ~ missense_num, data=missense_count)
    as.data.frame(as.list(coef(lm))) %>>%
      mutate(p_value = 1 - pf(summary(lm)$fstatistic["value"],summary(lm)$fstatistic["numdf"],
                              summary(lm)$fstatistic["dendf"]))
  }
  if((.race=="all")&(substitute(.vcf)=="tdg_gnomad")){
    .patient_list=.patient_list%>>%dplyr::select(patient_id,age)%>>%filter(!is.na(age))
    .vcf = .vcf %>>%
      dplyr::select(chr,start,ref,alt,AF)
  }
  if(.race!="all"){
    .patient_list = .patient_list %>>%
      filter(race==.race,!is.na(age))%>>%
      dplyr::select(patient_id,age)
    if(.race == "white"){
      .vcf = .vcf %>>%
        mutate(AF = AF_white) %>>%
        dplyr::select(chr,start,ref,alt,AF)
    }else if(.race=="black"){
      .vcf = .vcf %>>%
        mutate(AF = AF_black) %>>%
        dplyr::select(chr,start,ref,alt,AF)
    }else{stop(paste0(".race is wrong .race=",.race,"\nuse all, white or black!"))}
  }
  if(.maf_filter){.maf=maf_trim_for_cumulative(.vcf=.vcf,.race=.race,.fdr=.fdr,.database=.database,
                                               .duplicate=.duplicate,.somatic=.somatic,.varscan=.varscan)
  }
  .maf = .maf %>>%filter(mutype=="truncating"|mutype=="splice") %>>%
    filter(if(.remove0){MAF!=0}else{chr==chr}) %>>%
    left_join(driver_genes %>>%dplyr::select(gene_symbol,role), by = "gene_symbol") %>>%
    filter(role==.role | role=="oncogene/TSG")%>>%
    dplyr::select(patient_id,MAF,MAC)
  tibble::tibble(MAF=1:(.max_maf*100)) %>>%
    mutate(MAF = MAF/10000) %>>%
    mutate(regression = purrr::map(MAF,~regression_out(.,.maf,.role,.patient_list)))%>>%
    unnest()
}
#### regressionの図のplot 3pattern
regression_plot_log = function(.reg_tbl,.max_maf=50,.min=NA,.dred=NA,.blue=NA,.green=NA,.black=NA){
  .max=max(.reg_tbl$missense_num)
  .min=ifelse(is.na(.min),min(.reg_tbl$missense_num),.min)
  if(.max <0){.max=0}else if(.min >0){.min=0}
  .expand = (.max - .min)/10
  .plot = .reg_tbl %>>%
    mutate(MAF=100*MAF) %>>%
    ggplot()+
    geom_point(aes(x=MAF,y=missense_num))+
    #geom_vline(xintercept = 0,size =1)+
    #geom_hline(yintercept = 0,size =1)+
    labs(y="Regression Coefficient",x = "MAF (%)")+
    scale_x_log10(limits = c(0.01,.max_maf), expand = c(0,0.03),
                  breaks=c(0.01,0.1,1,10,50))+
    scale_y_continuous(limits = c(.min,.max), expand = c(0,.expand))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          panel.border = element_blank(),axis.line = element_line())
  if(!is.na(.black)){.plot=.plot+geom_vline(xintercept = .black,colour="black")}
  if(!is.na(.dred)){.plot=.plot+geom_vline(xintercept = .dred,colour="darkred")}
  if(!is.na(.blue)){.plot=.plot+geom_vline(xintercept = .blue,colour="blue")}
  if(!is.na(.green)){.plot=.plot+geom_vline(xintercept = .green,colour="green")}
  plot(.plot)
  return(.plot)
}
