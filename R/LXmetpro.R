
Lxmetpro <- function(protein_data,meta_pathways){

pro_pathways <- pro_pathways(protein_data)

protein_kegg_pathways_all <- pro_pathways$kegg_all_pathways

protein_metabolism_pathways <- pro_pathways$kegg_meta

dir.file <- pro_pathways$dir.file

species <- pro_pathways$species

joint_pathways <- joint_pathways(meta_pathways,protein_kegg_pathways_all,species,dir.file)

print("--------------------------------------------------------------")

print(paste("The results can be found in the folder of",dir.file))

joint_pathways$join_paht

}


