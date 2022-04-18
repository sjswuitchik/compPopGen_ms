Each file contains the set of individuals to remove for MK analyses. Individuals are filtered for:
 
-average coverage/missingness (outlier individuals) and coverage < 5x removed
-relatedness (any individuals with relatedness values > 0.03 - the value of a second degree relationship, one of the related individuals is selected for removal)
-population structure - if substantial structuring in the PCA and Admix plots, one population will be retained.

Note - species with no individuals to remove will have an empty file.

