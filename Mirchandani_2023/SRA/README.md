### Get SRA metadata

Testing out using Entrez searches and the SRA run selector to pull down metadata for processing in R.

Procedure is to run the searches in sra_searches.txt in sequence, saving each to a new file (for simplicity, and because the SRA Run Selector has a 20k record limit).

Then, combined and analyze in R.

### Get assembly data

To combined with genome data, pull information from https://www.ncbi.nlm.nih.gov/datasets/genomes/, searching by taxid: amphibia, sauropsids, chondrichthyes, cyclostomata, elopocephalai, otomorpha, neoteleostei

Note this strategy skips some odd fish clades: coelcanths, lungfish, Cladistia, Chondrostei, bowfins and gars, a few other genomes-poor clades.

### Pick samples

A list of possible candidates was identified by filtering all BioProjects identified from Entrez searches to require at least 10 BioSamples with coverage > 5x to be associated with a single BioProject for a species. While this will eliminate some specices that have either mostly very-low-coverage data, or have many small studies divided among multiple BioProjects, this requirement helps ensures a core of samples produced with consistent technical methods and by a single group, and greatly reduces false positive records.

Once a list of BioProjects that are of potential interest is generated in R, these were manually curated by searching for the BioProject ID in the NCBI "All Databases" search. If a PMC manuscript is linked to the BioProject ID, I skimmed the abstract and methods to determine if it was appropriate. Domesticated / captive samples, resequencing families for linkage map construction or QTL analyis, pooled sequencing projects, and ancient DNA samples were excluded.

If no published manuscript was linked to the BioProject ID, I searched the BioProject title to see if an obvious paper could be identified; if one was found, the linkage was confirmed by the Data Availablity statement referencing the correct BioProject or BioSample IDs. If no obvious manuscript could be identified by searching the BioProject title, the data was presumed to be unpublished and therefore excluded.

#### NOTE: still need to pick samples / check publications for stickleback due to an SRA search error ####

### Final cleanup

After manual curation, the final list of species to use is downloaded and imported into R, where the datasets can be cleaned to produce a species / BioProject / reference table that will be the basis of metadata creation.

However, because of historical quirks in how we tracked references, the BioProject / reference links need to be manually verified and filtered.

Any BioProject that only occurs once is assumed to have the correct publication information.

For BioProjects that are linked to more than one publication, initial check is to open the full text of the linked publication and search for the BioProject id. If it is found, keep it. If it is not found, check to see if a data availablity statement is present. If another BioProject is listed, then delete the pub/bioProject link. If no BioProject is listed, manually investigate.

### Oddball specices

A few species don't fit into the general rules above, for various reasons. These are detailed below.

The cod genus (Gadus species) are excluded, as these are the focus of a separate project using the same pipelines and any cod data will be synced from that project.

Cyanistes cyanus and Parkesia noveboracensis were sequenced for this project and do not have current BioProject IDs, and so were removed.

### Sample Metadata

Finally, we'll take the list of BioProjects and go back to our SRA data to make more detailed metadata files.
