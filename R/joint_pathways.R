
joint_pathways <- function(meta_pathways,protein_kegg_pathways_all,species,dir.file){

  #---------Jiont analysis of the gene and metabolite enriched pathways---------#

  meta_path <- read.xlsx(meta_pathways)
  meta_path <- meta_path[,c(1,6)]
  colnames(meta_path) <- c("pathways","meta_pvalue")


  desc <- grep("Descrip",colnames(protein_kegg_pathways_all),ignore.case = T)
  pvl <- grep("pvalu",colnames(protein_kegg_pathways_all),ignore.case = T)
  nn <- c(desc,pvl)
  pro_path <- protein_kegg_pathways_all[,nn]
  colnames(pro_path) <- c("pathways","pro_pvalue")

  path_venn <- inner_join(meta_path,pro_path,by="pathways")

  colnames(path_venn) <- c("Pathway","P_meta","P_pro")


  #------------------------
  pro_venn_path <- path_venn[c("Pathway","P_pro")]
  pro_venn_path$P_pro <- -log2(pro_venn_path$P_pro)
  pro_venn_path$types <- rep("proteins",nrow(pro_venn_path))
  colnames(pro_venn_path) <- c("Pathway","-log2(Pvalue)","types")

  meta_venn_path <- path_venn[c("Pathway","P_meta")]
  meta_venn_path$P_meta <- -log2(meta_venn_path$P_meta)
  meta_venn_path$types <- rep("metabolites",nrow(meta_venn_path))
  colnames(meta_venn_path) <- c("Pathway","-log2(Pvalue)","types")

  pro_meta_path <- bind_rows(pro_venn_path,meta_venn_path)
  colnames(pro_meta_path) <- c("Pathways","minus_log2_Pvalue","types")

  g_m_name <- paste0("The protein_meta_pathways_data"," (",species,")",".xlsx")
  g_m_name <-paste0(dir.file,"/", g_m_name)

  write.xlsx(pro_meta_path, g_m_name)

  y_p <- max(pro_meta_path$minus_log2_Pvalue)

  height_y <- y_p*1.2

  nrow_path <- nrow(pro_meta_path)/2

  joint_title_size <- case_when(nrow_path>=30 ~12,
                                nrow_path>=20 ~12,
                                TRUE ~14)

  joint_x_size <- case_when(nrow_path>=20 ~9,
                            nrow_path>=10 ~10,
                            TRUE ~12)

  joint_y_size <- case_when(nrow_path>=20 ~12,
                            nrow_path>=10 ~12,
                            TRUE ~12)

  joint_legend_size <- case_when(nrow_path>=30 ~12,
                                 nrow_path>=20 ~12,
                                 TRUE ~12)


  bar_width <- case_when(nrow_path>=30 ~0.9,
                         nrow_path>=20 ~0.8,
                         TRUE ~0.7)

  f1 <- ggplot(pro_meta_path, aes(x = Pathways, y = minus_log2_Pvalue,fill=types))+
    geom_bar(position = "dodge",stat = "identity",width = bar_width)+
    scale_fill_manual(values=c("#008b8b","#f08080"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          text = element_text(colour = "black", face="bold",size=12)) +
    labs(x="",y = ("-log2(Pvalue)"),title = paste('Protein-Metabolite Joint Pathways'))+
    scale_y_continuous(expand = c(0, 0),limits = c(0, height_y))

  f1

  log05 <- -log2(0.05)
  log01 <- -log2(0.01)

  line1 <- geom_hline(yintercept = c(log05),
                      linewidth = 0.6,
                      color = "blue",
                      lty = "dashed")
  line2 <- geom_hline(yintercept = c(log01),
                      linewidth = 0.6,
                      color = "red",
                      lty = "dashed")

  y1 <- geom_text(x=nrow(pro_meta_path)/2-2.5,y=log05+0.8,label = c("p<0.05"),
                  size=4,color="blue",fontface="italic")
  y2 <- geom_text(x=nrow(pro_meta_path)/2-2.5,y=log01+0.8,label = c("p<0.01"),
                  size=4,color="blue",fontface="italic")

  f2 <- f1+line1+line2+y1+y2
  f2

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =joint_title_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"),
          plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "cm")
    )+
    theme(plot.title = element_text(hjust = 0.5))


  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=joint_x_size,angle =45,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=joint_y_size))

  legend_theme <- theme(
    legend.title = element_blank(),
    legend.text = element_text(size = joint_legend_size, face = "bold"),
    legend.direction = "vertical",
    #legend.position = c(0.5,0.9),
    legend.background = element_blank()
  )

  f3 <- f2+mytheme+xytheme

  f3

  f3_name <- paste0("Protein-metabolite Joint pathways 01","(",species,")",".png")
  f3_name <-paste0(dir.file,"/", f3_name)

  ggsave(f3_name,f3,width=1200, height =1000, dpi=150,units = "px")


  f4 <- f3+
    labs(fill="")+
    theme(legend.direction = "horizontal",
          legend.position = c(0.5,0.92),
          legend.text = element_text(size=14,face = "bold") )

  f4

  f4_name <- paste0("Protein-metabolite Joint pathways 02","(",species,")",".png")
  f4_name <-paste0(dir.file,"/", f4_name)

  ggsave(f4_name,f4,width=1200, height =1000, dpi=150,units = "px")


  joint_path_result <- list(join_paht=f4)

  return(joint_path_result)

}
