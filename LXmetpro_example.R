
if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetpro")
library(LXmetpro)

??LXmetpro
#---------------------------------
data(protein_data_example) # a list of protein Uniprot ID
data(meta_pathways_example) # getting from http://impala.molgen.mpg.de/
#------------------------------
rm(list=ls())

devtools::load_all()

protein_data <- c("proteins.xlsx")
meta_pathways <- c("meta_kegg_pathways.xlsx")

LXmetpro(protein_data,meta_pathways)



