#####################################################################################################################
#violin plot する際変異の最大数までに間が空いてしまう時に空の行を挿入するfunction
truncating_filling = function(.truncating_count){
  .max_count=max(.truncating_count$truncating_count_n)
  .outcome = data.frame(truncating_count_n=1:.max_count) %>>%
    left_join(.truncating_count%>>%count(truncating_count_n)) %>>%
    filter(is.na(n)) %>>%dplyr::select(-n)
  if(length(.outcome$truncating_count_n) != 0){
    .outcome = .outcome %>>%
      mutate(age=-50, cancer_type="BRCA", patient_id="no_patient",gender="male")
    return(rbind(.truncating_count,.outcome))
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
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  }
  .plot = .tbl %>>%truncating_filling()%>>%
    mutate(cancer_type ="All Cancer Types") %>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    geom_text(data =.tbl %>>%count(truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=5,label=n),size=3,position="stack")+
    geom_abline(data = regression %>>%filter(p_value > 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05),
                aes(intercept = X.Intercept.,slope = truncating_count_n))+
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
      mutate(out_reg =paste0("R=",signif(truncating_count_n,2)," P",
                             ifelse(p_value==0,"<0.0001",paste0("=",signif(p_value,2)))))
  }
  .plot = .tbl %>>%truncating_filling()%>>%
    ggplot(aes(x=as.factor(truncating_count_n), y=age))+
    geom_violin(scale = "count")+
    #geom_boxplot(width=.3,fill="black")+ 
    stat_summary(fun.y=mean,geom = "point", fill="black",shape=21,size=2)+
    scale_y_continuous(breaks = c(0,20,40,60,80),limits = c(0,90))+
    facet_wrap( ~ cancer_type, ncol = 5)+
    geom_text(data =.tbl %>>%count(cancer_type,truncating_count_n),
              aes(x=as.factor(truncating_count_n),y=2.5,label=n),size=2.5,position="stack",angle=-45)+
    geom_abline(data = regression %>>%filter(p_value > 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n),linetype="dashed")+
    geom_abline(data = regression %>>%filter(p_value <= 0.05,X.Intercept. >0),
                aes(intercept = X.Intercept.,slope = truncating_count_n))+
    geom_text(data = regression, aes(x=.coef_posi,y=15,label =out_reg),size=3.5,hjust=1)+
    xlab("Number of Truncated Gene")+ylab("Age at Diagnosis")+
    theme(panel.grid.minor.x = element_blank(),panel.background = element_rect(fill="transparent",colour="black"),
          panel.grid.major.y = element_line(colour = "gray"),panel.grid.major.x = element_line(colour = "gray95"),
          axis.line = element_line(colour = "black"),axis.ticks.y = element_blank(),
          axis.title = element_text(size=15), axis.text = element_text(size=15),strip.text = element_text(size=12),
          strip.background = element_rect(fill="transparent", colour = "black"))
  plot(.plot)
  return(.plot)
}

