#!/usr/bin/Rscript

library(tidyverse)
library(lme4)
library(arm)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

cds <- read_tsv("data/mk_tests/onlyCDS.genes.bed", col_names = c("chr", "start", "end", "gene")) %>%
  mutate(cds.temp = end - start) %>%
  group_by(gene) %>% 
  summarise(cds.len = sum(cds.temp))

callable <- read_tsv("data/mk_tests/callable.cds.bed", col_names = c("chr", "start", "end", "gene")) %>%
  mutate(call.temp = end - start) %>%
  group_by(gene) %>% 
  summarise(call.len = sum(call.temp))

ingroup <- read_tsv(args[1], col_names = c("chr", "start", "end", "effect", "gene"))
mis.in <- ingroup %>% 
  group_by(gene) %>% 
  tally(effect ==  "missense_variant") %>%
  set_names(c("gene","PR"))
syn.in <- ingroup %>%
  group_by(gene) %>%
  tally(effect == "synonymous_variant") %>%
  set_names(c("gene", "PS"))

outgroup <- read_tsv(args[2], col_names = c("chr", "start", "end", "effect", "gene"))
mis.out <- outgroup %>% 
  group_by(gene) %>% 
  tally(effect ==  "missense_variant") %>%
  set_names(c("gene","FR"))
syn.out <- outgroup %>%
  group_by(gene) %>%
  tally(effect == "synonymous_variant") %>%
  set_names(c("gene", "FS"))

# create full table, replacing any NAs with 0s
MKtable <- full_join(mis.out, syn.out, by = "gene") %>% full_join(mis.in, by = "gene") %>% full_join(syn.in, by = "gene") %>%
  mutate(FR = replace_na(FR, 0)) %>%
  mutate(FS = replace_na(FS, 0)) %>%
  mutate(PR = replace_na(PR, 0)) %>%
  mutate(PS = replace_na(PS, 0)) 

# qc: check for complete NA conversion in each column
MKtable %>% summarise(count = sum(is.na(FR)))
MKtable %>% summarise(count = sum(is.na(FS)))
MKtable %>% summarise(count = sum(is.na(PR)))
MKtable %>% summarise(count = sum(is.na(PS)))

# qc: make sure no ratios > 1
check <- full_join(callable, cds, by = "gene") %>%
  mutate(check = call.len/cds.len) %>%
  mutate(call.len = replace_na(call.len, 0)) %>%
  mutate(cds.len = replace_na(cds.len, 0)) %>%
  mutate(check = replace_na(check, 0))
max(check$check)
check %>% filter(check > 1)

# Prep for SnIPRE
in.miss <- read_delim(args[3], delim = "\t") %>% 
  mutate(indv = N_DATA - N_MISS) %>%
  dplyr::select(indv)
out.miss <- read_delim(args[4], delim = "\t") %>% 
  mutate(indv = N_DATA - N_MISS) %>%
  dplyr::select(indv)

snipre_data <- full_join(MKtable, callable, by = "gene") %>%
  mutate(Tsil = call.len * 0.33) %>%
  mutate(Trepl = call.len * 0.66) %>%
  add_column(npop = median(in.miss$indv)/2) %>%
  add_column(nout = median(out.miss$indv)/2) %>%
  dplyr::select(-c(call.len)) %>%
  filter((Trepl/Tsil)<5) %>%
  filter((PR+FR+PS+FS)>1) %>% 
  write_delim(., "data/mk_tests/snipre_data.tsv", delim = '\t', col_names = T) 
