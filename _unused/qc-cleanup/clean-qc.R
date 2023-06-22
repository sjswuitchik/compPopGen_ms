library(tidyverse)
library(parzer)

setwd("~/Projects/popgen/compPopGen_ms/qc-cleanup/")

runs<-read_csv("all_run_data.csv", col_names = T) %>% distinct()

lat_lon<-read_csv("lat_lon_input.csv", col_names = T)
lat_sep<-read_csv("Lat_Long_separate_input.csv", col_names=T)

#clean up lat_sep

lat_sep <- lat_sep %>% filter(!is.na(Latitude), !is.na(Longitude)) %>%
  mutate(lat = parse_lat(Latitude), long = parse_lon(Longitude)) %>%
  select(BioSample, lat, long)

#clean up lat lon

lat_lon <- lat_lon %>% filter(!is.na(lat_lon)) %>%
  mutate(lat_lon_clean1 = str_replace_all(lat_lon, "[\r\n]", " ") %>% 
           str_trim()) %>%
  mutate(lat_lon_clean = ifelse(str_starts(lat_lon_clean1, "[NSEW]"),
                                str_replace_all(lat_lon_clean1, "(\\S)\\s+([NSEW])", "\\1 , \\2"),
                                str_replace_all(lat_lon_clean1, "([NSEW])\\s+([:digit:])", "\\1 , \\2"))
           ) %>%
  #need to clean up one bioproject where minutes symbol was convered to 0
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "8016.77", "8'16.77")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "23027.52", "23'27.52")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "22025.59", "22'25.59")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "20047.10", "20'47.10")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "59027.60", "59'27.60")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "26045.60", "26'45.60")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "7051.77", "7'51.77")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "22018.18", "22'18.18")) %>%
  mutate(lat_lon_clean = str_replace_all(lat_lon_clean, "36024.07", "36'24.07")) %>%
  
  mutate(lat = parse_llstr(lat_lon_clean)$lat,
         long =  parse_llstr(lat_lon_clean)$lon) %>%
  select(BioSample, lat, long)

qc %>% filter(is.na(qc)) %>% select(BioSample) %>% write_tsv("need_qc_biosamples.txt")
rs <- read_tsv("runselector_cleaned.tsv", col_names = T)

new_qc <- rs %>% mutate(lat_lon = str_trim(lat_lon)) %>% 
  mutate(lat_lon_clean = str_replace_all(lat_lon, "[\r\n\"\'\ ́]", " ") %>%
           str_replace_all("([NSEW])\\s+([:digit:])", "\\1 , \\2")) %>% 
  mutate(lat = parse_llstr(lat_lon_clean)$lat,
         long =  parse_llstr(lat_lon_clean)$lon) %>%
  select(BioSample, lat, long)


all_qc <- bind_rows(lat_sep,lat_lon,new_qc)

all_runs <- left_join(runs, all_qc) %>% filter(refGenome != "GCA_900092035.1")
table(all_runs$refGenome)

table(is.na(all_runs$lat))

all_runs %>% distinct(BioSample,refGenome,.keep_all=TRUE) %>% write_csv("final_QC_metadata.csv")

genome_key <- all_runs %>% distinct(refGenome, Organism)

#read in bam sumstats info

ss <- read_csv("all_sumstats_merged.csv", col_names=T, trim_ws = TRUE) %>%
  mutate(source = str_remove_all(File, fixed("sumstats/")) %>% str_remove_all("_bam_sumstats.txt")) %>%
  mutate(Organism = str_remove_all(source, "_GC[AF]_[0-9]+\\.[0-9]")) %>%
  select(-File, -source) %>% left_join(genome_key, by=c("Organism" = "Organism"))

final_samples <- all_runs %>% distinct(BioSample,refGenome,.keep_all=TRUE) %>% full_join(ss, by=c("BioSample" = "Sample", "refGenome" = "refGenome", "Organism" = "Organism"))

sumstats <- final_samples %>% rename(Sample = BioSample) %>% select(-Organism) %>% select(-lat) %>% select(-long)

# output separate sum stats files for each genome

# Split the tibble by the refGenome column
split_tibbles <- split(sumstats, sumstats$refGenome)

# Write each split tibble to a separate output file
for (split_value in names(split_tibbles)) {
  tibble <- split_tibbles[[split_value]] %>% select(-refGenome)
  # Use the paste0 function to construct the output filename
  filename <- paste0(split_value, "/summary_stats/final_bam_sumstats.txt")
  dir <- dirname(filename)
  if (!dir.exists(dir)) {
    dir.create(dir,recursive = TRUE)
  }
  # Use the write_tsv function from the write_tsv function from the readr package to write the tibble to a TSV file
  write_tsv(tibble, filename)
}


