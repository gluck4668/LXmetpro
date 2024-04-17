
joint_pathways <- function(meta_pathways,protein_kegg_pathways_all,species,dir.file){

#---------Jiont analysis of the gene and metabolite enriched pathways---------#
meta_data_type <- str_extract(meta_pathways,"(?<=[.]).*")%>% tolower()
if(meta_data_type=="txt")
  meta_path_0 <- read_table(meta_pathways) else
    meta_path_0 <- eval(str2expression(paste0("read.",meta_data_type,"(meta_pathways)")))

#----数据来源于IMPaLA------------------
is.impala <- grepl("overlap",colnames(meta_path_0),ignore.case = T) %>% any()
 if(is.impala){
      meta_path <- dplyr::filter(meta_path_0,pathway_source=="KEGG")
      num_all <- str_extract(meta_path$num_all_pathway_metabolites,".*(?= \\()") %>% as.numeric()
      num_meta <- meta_path$num_overlapping_metabolites
      meta_path <- mutate(meta_path,meta_Ratio=num_meta/num_all)
      meta_path <- meta_path[,c("pathway_name","P_metabolites","meta_Ratio")]
      colnames(meta_path) <- c("pathway","log2p_meta","meta_Ratio")
      meta_path$pathway <- str_extract(meta_path$pathway,".*(?= -)") %>% trimws()
      meta_path$log2p_meta = -log2(meta_path$log2p_meta)
    }

#----数据来源于MetabolAnalyst------------------
is.metabo <- grepl("impact",colnames(meta_path_0),ignore.case = T) %>% any()
if(is.metabo){
      meta_path <- meta_path_0[,c(1,5)] %>% mutate(meta_ratio=meta_path_0[,4]/meta_path_0[,2])
      colnames(meta_path) <- c("pathway","log2p_meta","meta_Ratio")
      meta_path$log2p_meta = -log2(meta_path$log2p_meta)
    }

#----protein pathways----------------------------
pro_path <- protein_kegg_pathways_all[,c("Description","pvalue","ProteinRatio")]
colnames(pro_path) <- c("pathway","log2p_pro","protein_Ratio")
pro_path$log2p_pro <- -log2(pro_path$log2p_pro)
pro_path <- data.frame(pro_path)

#---- inner_join-------------------------------
path_venn <- inner_join(meta_path,pro_path,by="pathway")
pro_meta_path <- data.frame(Pathway=rep(path_venn$pathway,2),
                            log2pvalue=c(path_venn$log2p_meta,path_venn$log2p_pro),
                            Ratio=c(path_venn$meta_Ratio,unlist(path_venn$protein_Ratio)),
                            Type=c(rep("metabolite",nrow(path_venn)),rep("protein",nrow(path_venn))))

#------准备作图------------------
g_m_name <- paste0("The protein_meta_pathways_data"," (",species,")",".xlsx")
g_m_name <-paste0(dir.file,"/", g_m_name)
write.xlsx(pro_meta_path, g_m_name)

y_p <- max(pro_meta_path$log2pvalue)
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

f1 <- ggplot(pro_meta_path, aes(x =Pathway, y = log2pvalue,fill=Type))+
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

#---------ggplot point----------------
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
  theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=10,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=10))+
  theme(legend.text=element_text(face="bold",color="black",size=10))

f5 <- ggplot(pro_meta_path)+
  geom_point(aes(x=Ratio,
                 y=Pathway,
                 shape=Type,
                 color=log2pvalue,
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'Ratio', y = 'Pathways',title=c('Protein-Metabolite Joint Pathways'),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme

f5

f5_name <-paste0("Protein-metabolite Joint pathways 03","(",species,")",".png")
f5_name <-paste0(dir.file,"/", f5_name)
ggsave(f5_name,f5,width=1450, height =1200, dpi=150,units = "px")

#------------------
line05 <- geom_vline(xintercept = c(log05),
                     linewidth = 0.5,
                     color = "black",
                     lty = "dashed")

txt05 <- geom_text(x=log05+4,y=2.5,label = c("p<0.05"),
                   size=5,color="blue",fontface="italic")

f6 <- ggplot(pro_meta_path)+
  geom_point(aes(x=log2pvalue,
                 y=Pathway,
                 shape=Type,
                 color=log2pvalue,
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = '-log2(pvalue)', y = 'Pathways',title=c('Protein-Metabolite Joint Pathways'),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme+
  line05+txt05

f6

f6_name <-paste0("Protein-metabolite Joint pathways 04","(",species,")",".png")
f6_name <-paste0(dir.file,"/", f6_name)
ggsave(f6_name,f6,width=1450, height =1200, dpi=150,units = "px")

#-------------------------------------
joint_path_result <- list(join_paht=f4)
return(joint_path_result)

}
