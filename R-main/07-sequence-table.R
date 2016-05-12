

library(readxl)
library(dplyr)
library(readr)

x <- read_excel("data/sequence/assay_sequences.xls", skip = 1) %>%
  tbl_df()
x$`dbSNP Accession Number` = "xxx-xxx"

nc <- read_csv(file = "outputs/number-of-genotype-clusters.csv") %>%
  rename(`Locus name` = assay.name)
nc$nClusters[nc$nClusters == 0] <- "NC"


left_join(x, nc) %>%
  mutate(`No. of genotype categories` = nClusters) %>%
  select(-nClusters) %>%
  write_csv(path = "outputs/supp-tab-1-green-sturgeon-snps.csv")


# then, open that in excel (yuck!) and paste the appropriate parts into
# ./supplements/Supp-Tab-1-SNP-sequences.xls