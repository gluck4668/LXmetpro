
LXmetpro <- function(protein_data,meta_pathways){

prot_path <- pro_pathways(protein_data)

species <-prot_path$species

joint_pathways <- joint_pathways(meta_pathways,protein_kegg_pathways_all,species,dir.file)

print("--------------------------------------------------------------")

print(paste("The results can be found in the folder of '",dir.file,"'"))

joint_pathways$join_paht

}


