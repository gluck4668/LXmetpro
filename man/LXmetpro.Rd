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
library(LXmetpro)

??LXmetpro
#---------------------------------
data(protein_data_example) # a list of protein Uniprot ID
data(meta_pathways_impala) # obtained from http://impala.molgen.mpg.de/
data(meta_pathway_metaboAnlyst) # obtained from https://www.metaboanalyst.ca/
#------------------------------
rm(list=ls())

devtools::load_all()

protein_data <- c("proteins.xlsx")
meta_pathways <- c("meta_pathways_impala.csv")
# meta_pathways <- c("meta_pathway_metaboAnlyst.csv")

LXmetpro(protein_data,meta_pathways)


}
