library(tidyverse)

#read SRA searches
#kinda annoying as searches are manual but could be updated later

#this parsing function may need to be updated to add more detail once we are making sample sheets

read_sra_clean <- function(file, path) {
  df<-read_csv(paste0(path, "/", file), cols(.default="c"), col_names = TRUE) %>%
    select(Run, BioSample, Experiment, Instrument, LibrarySelection, LibrarySource, Organism, Platform, 
           SampleName = `Sample Name`, SRAStudy = `SRA Study`, Bases, AvgSpotLen, BioProject, sex, Isolate,
           Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Ecotype, Strain,
           lat_lon)
}

read_assembly_clean <- function(file, path) {
  df<-read_tsv(paste0(path, "/", file, "_data_summary.tsv"), col_names = TRUE) %>%
    rename_with(~ gsub(" ", "", .x, fixed=TRUE))
}

files<-c("SRA-Agnatha.txt", "SRA-Amphibia.txt", "SRA-Aves.txt", "SRA-Chondrichthyes.txt", 
         "SRA-Reptilia.txt", "SRA-Sarcopterygii.txt", "SRA-fish-1.txt", "SRA-fish-2..txt")

sra_list<-lapply(files, read_sra_clean, path="~/Projects/popgen/compPopGen_ms/SRA")

#may be some parsing errors due to heading issues but shouldn't matter

#flatten to single tibble, with distinct in case of duplicates in searches

sra<-bind_rows(sra_list) %>% distinct()

#further clean up SRA

table(sra$LibrarySource)
table(sra$LibrarySelection)

#filter metagenomic and suspcicious library selection methods, this may miss a few things but should be cleaner

sra <- sra %>% filter(LibrarySource == "GENOMIC", 
              LibrarySelection == "RANDOM" | LibrarySelection == "unspecified" | LibrarySelection == "PCR" | LibrarySelection == "other" | LibrarySelection == "RANDOM PCR")

##loading genome info

genome_search_list <- c("amphibia", "sauropsids", "chondrichthyes", "cyclostomata", "elopocephalai", "otomorpha", "neoteleostei")

genome_list <- lapply(genome_search_list, read_assembly_clean, path="~/Projects/popgen/compPopGen_ms/SRA")
assemblies <- bind_rows(genome_list) %>% distinct() %>% 
  arrange(Taxonomyid, desc(Source), Level, desc(ContigN50)) %>%
  distinct(Taxonomyid, .keep_all = TRUE)

assemblies %>% ggplot(aes(log10(ContigN50))) + geom_histogram()

#some arbitrary filtering

assemblies <- assemblies %>% filter(ContigN50 > 1e4)


#merge

popgen <- right_join(sra, assemblies, by=c("Organism" = "OrganismScientificName"), 
                     suffix = c(".popgen", ".assembly")) %>% 
  select(-Run, -Experiment) %>% group_by(BioSample.popgen) %>%
  mutate(bases_total = sum(as.numeric(Bases))) %>% 
  select(-Bases, -AvgSpotLen) %>%
  distinct() %>% 
  mutate(coverage = as.numeric(bases_total) / as.numeric(Size))

#some quick analysis

popgen %>% mutate(covplot = ifelse(coverage < 100, coverage, 100)) %>% ggplot(aes(covplot)) + geom_histogram()

#collate with stuff we've processed

processed<-read_csv("~/Downloads/Final_MS_data - processed.csv")

#write out 

popgen %>% filter(coverage > 5) %>% group_by(Organism, BioProject.popgen) %>% count() %>% filter(n > 10) %>% 
  full_join(processed, by=c("Organism" = "species")) %>% 
  write_tsv(path="~/Projects/popgen/compPopGen_ms/SRA/bioprojects.tsv")


