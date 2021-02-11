#!/usr/bin/Rscript

library(tidyverse)
library(lme4)
library(arm)

# read in data 
snipre_data <- read_tsv("snipre_data.tsv", delim = '\t', col_names = T)

# functions for MK calculations and DOS 
mk_test <- function(dn=dn,ds=ds,pn=pn,ps=ps){
  dnds_mat <- matrix(data=c(ds,dn,ps,pn),nrow=2,byrow = F)
  pval = fisher.test(dnds_mat)$p.value
  alpha = 1-((ds*pn)/(dn*ps))
  return(list(pval,alpha))
} 

mk_tibble_calc <- function(snipre_res_obj){
  snipre_res_obj_new <- snipre_res_obj %>%
    rowwise %>%
    mutate(mk_pval = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[1]],
           alpha = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[2]]) %>%
    ungroup %>%
    dplyr::mutate(mk_pval_fdr = p.adjust(mk_pval,method="BH"))
  return(snipre_res_obj_new)
}

# run MK test
MKtest <- mk_tibble_calc(snipre_data) %>%
  mutate(dos = FR/(FR + FS) - PR/(PR + PS),
         total_poly = PR + PS,
         total_div = FR + FS) 

write.table(MKtest, "mk_output.tsv", sep = "\t", quote = F, row.names = F)

# run SnIPRE 
source("helper_scripts/SnIPRE_source.R")
source("helper_scripts/my.jags2.R")

snipre.res <- SnIPRE(MKtest)
snipre.qres <- snipre.res$new.dataset
snipre.model <- snipre.res$model
snipre.table <- table(snipre.qres$SnIPRE.class)

write.table(snipre.qres, "snipre_output.tsv", sep = "\t", quote = F, row.names = F)
