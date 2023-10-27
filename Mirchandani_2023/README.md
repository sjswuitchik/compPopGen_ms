This repository stores scripts and metadata used to produce the [Comparative Population Genomics Data collection](https://snparcher.readthedocs.io/en/latest/datasets.html), described in the [snpArcher paper](https://www.biorxiv.org/content/10.1101/2023.06.22.546168v1).

### SRA
Scripts and documentation for searching the SRA for suitable datasets can be found in the `SRA` directory. 

### Final Sample Sheets
The `final_sample_sheets` directory stores CSV sample sheets used to run snpArcher. These are organized by reference genome accession.

### Post Process Sample Sheets
The `postprocess_sample_sheets` directory stores CSV sample sheets to run the postprocess module of snpArcher on each dataset.

### MK Tests
Scripts and data to produce the MK test figure in the manuscript are located in the `mk` directory.

### Dataset results
The Rmarkdown file `datasets.rmd` reports statisitcs, such as number of SNPs and Wattersons Theta, from each of the processed datasets.
