library(tidyverse)
library(purrr)
library(stringr)

#read SRA searches
#kinda annoying as searches are manual but could be updated later

#this parsing function may need to be updated to add more detail once we are making sample sheets

path_to_write = "~/Projects/popgen/compPopGen_ms/SRA"

read_sra_clean <- function(file, path) {
  df<-read_csv(paste0(path, "/", file), cols(.default="c"), col_names = TRUE) %>%
    select(Run, BioSample, Experiment, Instrument, LibrarySelection, LibrarySource, Organism, Platform, 
           SampleName = `Sample Name`, SRAStudy = `SRA Study`, Bases, AvgSpotLen, BioProject, sex, Isolate,
           Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Locality = geo_loc_name, lat_lon,
           Ecotype, Strain)
}

read_assembly_clean <- function(file, path) {
  df<-read_tsv(paste0(path, "/", file, "_data_summary.tsv"), col_names = TRUE) %>%
    rename_with(~ gsub(" ", "", .x, fixed=TRUE))
}

files<-c("SRA-Agnatha.txt", "SRA-Amphibia.txt", "SRA-Aves.txt", "SRA-Chondrichthyes.txt", 
         "SRA-Reptilia.txt", "SRA-Sarcopterygii.txt", "SRA-fish-1.txt", "SRA-fish-2..txt", "SRA-stickleback.txt")

sra_list<-lapply(files, read_sra_clean, path="~/Projects/popgen/compPopGen_ms/SRA")

#may be some parsing errors due to heading issues but shouldn't matter

#flatten to single tibble, with distinct in case of duplicates in searches

sra<-bind_rows(sra_list) %>% distinct()

#further clean up SRA

table(sra$LibrarySource)
table(sra$LibrarySelection)

#filter metagenomic and suspcicious library selection methods, this may miss a few things but should be cleaner

sra_clean <- sra %>% filter(LibrarySource == "GENOMIC", 
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

popgen <- right_join(sra_clean, assemblies, by=c("Organism" = "OrganismScientificName"), 
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

# manual curation happens here

datasets <- read_csv("~/Downloads/Final_MS_data - sra.csv") %>% filter(use == "yes") %>%
  select(Organism, bioproject1 = BioProject.popgen, bioproject2 = bioproject.orig, publication)


#clean up -- ugly code

datasets_clean <- rbind(tibble(Organism = datasets$Organism, 
                               BioProject = datasets$bioproject1, 
                               Pub = datasets$publication),
                        tibble(Organism = datasets$Organism, 
                               BioProject = datasets$bioproject2, 
                               Pub = datasets$publication)) %>%
  filter(!is.na(BioProject)) %>% 
  separate_rows(BioProject) %>%
  separate_rows(Pub, sep="[,; ]") %>%
  arrange(Organism) %>%
  filter(Pub != "") %>%
  distinct()

#write this out to manually clean up publications/SRA link

datasets_clean %>% write_tsv("~/Projects/popgen/compPopGen_ms/SRA/datasets_initial.tsv")

#after cleanup

bioprojects <- read_tsv("~/Projects/popgen/compPopGen_ms/SRA/final_bioprojects.tsv")

# now, let's append our popgen data to bioprojects
# this gets all the samples for each bioproject, including species we didn't select
# for each species, mark if it is in the "use" list, for info
# then append assembly information, which might be NA

datasets_final <- bioprojects %>% 
  select(BioProject) %>% distinct() %>% 
  left_join(sra, by=c("BioProject" = "BioProject")) %>%
  mutate(use_species = ifelse(Organism %in% bioprojects$Organism, 1, 0)) %>%
  left_join(assemblies, by=c("Organism" = "OrganismScientificName"), 
            suffix = c(".popgen", ".assembly"))

# a little data exploration

datasets_final %>% group_by(BioProject.popgen) %>% summarise(used_sp = mean(use_species)) %>% filter(used_sp < 1) %>% arrange(used_sp) %>% print.data.frame()

#looks reasonable, now let's make the final "big" metadata

#from published bioprojects, here are our plausible ingroups, defined as species for which the genome assembly organism == the resequencing organism

ingroups <- datasets_final %>% group_by(Organism) %>% filter(!is.na(AssemblyAccession)) %>% count() %>% filter(n >= 5)

#so, to clean up, let's make the following
#sra download / mapping metadata, which should be a file with a column for bioproject, biosample, run, one file per species

datasets_final %>% filter(Organism %in% ingroups$Organism) %>% 
  select(Organism, BioProject = BioProject.popgen, BioSample = BioSample.popgen, Run) %>% 
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>%
  split(., .$Organism) %>%
  imap(~ write_tsv(as.data.frame(.x), path = str_c(path_to_write, '/SRA_Run_', .y, '.tsv')))

#sample metadata -- much of this will be incomplete since it needs to be linked to the publication, but it is a start

not_all_na <- function(x) {!all(is.na(x))}

datasets_final %>% filter(Organism %in% ingroups$Organism) %>% 
  select(Organism, BioProject = BioProject.popgen, BioSample = BioSample.popgen, SampleName, Sex = sex, Isolate, Country, Continent, Locality, Ecotype, Strain) %>% 
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>% distinct() %>%
  split(., .$Organism) %>%
  imap(~ write_tsv(select_if(as.data.frame(.x), not_all_na), path = str_c(path_to_write, '/SRA_Metadata_', .y, '.tsv')))

#organism metadata - genome assembly accession and annotation information

datasets_final %>% filter(Organism %in% ingroups$Organism) %>% 
  select(Organism, AssemblyName, AssemblyAccession, Annotation, Level, ContigN50, Size) %>%
  mutate(Organism = str_replace_all(Organism, " ", "_")) %>% distinct() %>% 
  write_tsv(path = str_c(path_to_write, '/Organism_Metadata.tsv'))
