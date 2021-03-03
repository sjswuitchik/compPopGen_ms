library(tidyverse)

miss <- read_delim("missing_data_per_ind.txt", delim = '\t', col_names = T) %>%
  mutate(missing = N_MISS/N_DATA) %>%
  select(INDV, missing) %>%
  write_delim(., "relative_missing_per_ind.txt", delim = '\t', col_names = T)

pdf("relative_missing.pdf")
dotchart(x = sort(miss$missing), labels = miss$INDV)
dev.off()
