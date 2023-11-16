
library(openxlsx)

protein_data_example <- read.xlsx("proteins.xlsx")
meta_pathways_example <- read.xlsx("meta_kegg_pathways.xlsx")

usethis::use_data(protein_data_example,overwrite = T)
usethis::use_data(meta_pathways_example,overwrite = T)

rm(list=ls())

data(protein_data_example) # a list of protein Uniprot ID
data(meta_pathways_example) # getting from http://impala.molgen.mpg.de/

