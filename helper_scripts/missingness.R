#!/usr/bin/Rscript

library(tidyverse)

# read in data and calculate missingness
ingroup <- read.delim("ingroup_missing_data.txt", delim = "\t", col_names = T) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(ingroup$missing)
remove <- ingroup %>% 
  filter(missing >= threshold)
write.csv(remove %>% select(INDV), "ingroup.remove.indv", row.names = F, quote = F)

outgroup <- read.delim("outgroup_missing_data.txt", delim = "\t", col_names = T) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(outgroup$missing)
remove <- outgroup %>% 
  filter(missing >= threshold)
write.csv(remove %>% select(INDV), "outgroup.remove.indv", row.names = F, quote = F)

