Testing out using Entrez searches and the SRA run selector to pull down metadata for processing in R.

Procedure is to run the searches in sra_searches.txt in sequence, saving each to a new file (for simplicity, and because the SRA Run Selector has a 20k record limit).

Then, combined and analyze in R.

To combined with genome data, pull information from https://www.ncbi.nlm.nih.gov/datasets/genomes/, searching by taxid: amphibia, sauropsids, chondrichthyes, cyclostomata, elopocephalai, otomorpha, neoteleostei

Note this strategy skips some odd fish clades: coelcanths, lungfish, Cladistia, Chondrostei, bowfins and gars, a few other genomes-poor clades.
