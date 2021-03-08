library(tidyverse)

#read SRA searches
#kinda annoying as searches are manual but could be updated later

read_sra_clean <- function(file, path) {
  df<-read_csv(paste0(path, "/", file), cols(.default="c"), col_names = TRUE) %>%
    select(Run, BioSample, Experiment, Instrument, LibrarySelection, LibrarySource, Organism, Platform, 
           SampleName = `Sample Name`, SRAStudy = `SRA Study`, Bases, AvgSpotLen, BioProject, sex, Isolate,
           Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Ecotype, Strain,
           lat_lon)
}

files<-c("SRA-Agnatha.txt", "SRA-Amphibia.txt", "SRA-Aves.txt", "SRA-Chondrichthyes.txt", 
         "SRA-Reptilia.txt", "SRA-Sarcopterygii.txt", "SRA-fish-1.txt", "SRA-fish-2..txt")

sra_list<-lapply(files, read_sra_clean, path="~/Projects/popgen/compPopGen_ms/SRA")

#may be some parsing errors due to heading issues but shouldn't matter

#flatten to single tibble, with distinct in case of duplicates in searches

sra<-bind_rows(sra_list) %>% distinct()

#get list of all species with one BioProject with at least 10 BioSamples

species_list <- sra %>% select(Organism, BioProject, BioSample) %>% distinct() %>% 
  group_by(Organism, BioProject) %>% count() %>% filter(n >= 10) %>% ungroup() %>% select(Organism) %>% distinct()

#view

sra %>% filter(Organism %in% species_list$Organism) %>% select(Organism, BioProject, BioSample) %>%
  distinct() %>% group_by(Organism, BioProject) %>% count() %>% View()

#one species
sra %>% filter(Organism == "Anguilla anguilla") %>% select(Organism, BioProject, BioSample) %>%
  distinct() %>% group_by(Organism, BioProject) %>% count() %>% View()

#check publication, again using manual searches
#output bioproject list

sra %>% filter(Organism %in% species_list$Organism) %>% select(BioProject) %>%
  distinct() %>% write_tsv("bioprojects.txt")
