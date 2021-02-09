#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

ingroup <- read_delim(args[1], ) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(ingroup$missing)
remove <- ingroup %>% 
  filter(missing >= threshold) %>%
  select(INDV) %>%
  write_delim(., "ingroup.remove.indv", delim = '\t', col_names = T)

outgroup <- read_delim(args[2]) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(outgroup$missing)
remove <- outgroup %>% 
  filter(missing >= threshold) %>%
  select(INDV) %>%
  write_delim(., "outgroup.remove.indv", delim = '\t', col_names = T)
