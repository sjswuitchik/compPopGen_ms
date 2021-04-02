#this is a simple script to take a BioProject downloaded from RunSelector, and make it compatible with the clean metadata pipeline.

#ideally this should be refactored into functions shared with the SRA_parser.R script

library(tidyverse)
library(purrr)
library(stringr)

not_all_na <- function(x) {!all(is.na(x))}

read_sra_full <- function(file, path) {
  df<-read_csv(paste0(path, "/", file), cols(.default="c"), col_names = TRUE)
}


not_all_na <- function(x) {!all(is.na(x))}

write_sample_metadata <- function(df) {
  df %>%
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>% distinct() %>%
  mutate(sex = case_when(
    is.na(sex) ~ "missing",
    sex == "male and female" | sex == "mixed" | sex == "pooled male and female" ~ "pooled",
    tolower(sex) == "na" | tolower(sex) == "not applicable"  ~ "missing",
    tolower(sex) == "not collected" | tolower(sex) == "not determined" ~ "unknown",
    TRUE ~ sex
  )) %>%
  mutate(`Library Name` = case_when(
    is.na(`Library Name`) ~ str_c("library_", `Sample Name`),
    TRUE ~ `Library Name`
  )) %>%
  split(., .$BioProject) %>%
  imap(~ write_tsv(select_if(as.data.frame(.x), not_all_na), file = str_c(path_to_write, '/SRA-sample-metadata/SRA_Metadata_', .y, '.tsv')))
}

temp<-read_sra_full("SRA_Metadata_PRJEB2984.txt", "~/Downloads")
write_sample_metadata(temp, "~/Projects/popgen/compPopGen_ms/SRA")
