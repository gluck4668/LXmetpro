
library(openxlsx)

protein_data_example <- read.xlsx("proteins.xlsx")
meta_pathways_impala <- read.csv("meta_pathways_impala.csv")
meta_pathway_metaboAnlyst <- read.csv("meta_pathway_metaboAnlyst.csv")

usethis::use_data(protein_data_example,overwrite = T)
usethis::use_data(meta_pathways_impala,overwrite = T)
usethis::use_data(meta_pathway_metaboAnlyst,overwrite = T)

rm(list=ls())

data(protein_data_example) # a list of protein Uniprot ID
data(meta_pathways_impala) # obtained from http://impala.molgen.mpg.de/
data(meta_pathway_metaboAnlyst) # obtained from https://www.metaboanalyst.ca/

