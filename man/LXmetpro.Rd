\name{LXmetpro}
\alias{LXmetpro}
\title{LX::LXmetpro}
\usage{
Lxmetpro(protein_data,meta_pathways)
}
\description{
It can be used to analyze the correlation between the proteins KEGG pathways and metabolites enriched pathways.
}
\examples{
if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXmetpro")
library(LXmetgen2)

??LXmetpro
#---------------------------------
data(protein_data_example) # a list of protein Uniprot ID
data(meta_pathways_example) # getting from http://impala.molgen.mpg.de/
#------------------------------
rm(list=ls())

protein_data <- c("proteins.xlsx")
meta_pathways <- c("meta_kegg_pathways.xlsx")

Lxmetpro(protein_data,meta_pathways)

}
