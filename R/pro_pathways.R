pro_pathways <- function(protein_data){

  #----------安装及引用相关普通R包-----------------------------------------
inst_packs <- function(){
  all_packages <- data.frame(installed.packages()) # 查看已安装的R包

  pack <- c("devtools","BiocManager","ggnewscale","R.utils", "ggtext",#需要安装的R包
            "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych",
            "ggplot2","ggrepel","RColorBrewer", "ggthemes","rticles",
            "grid","patchwork","Hmisc","pak","readr")

  is_pack <- pack[!pack %in% all_packages$Package] # 筛选出未安装的R包

  fun_install <- function(x) {install.packages(x,update = F,ask = F)} #安装R包函数
  sapply(is_pack,fun_install,simplify = T)  #用sapply批量安装R包

  # 批量library
  fun_library <- function(x){library(x, character.only = T)}
  sapply(pack,fun_library,simplify = T)

  #-----tidyverse比较难安装，需要pak函数来装------------------------------
  if(!"tidyverse" %in% all_packages$Package)
    pak::pak("tidyverse/tidyverse")
  library(tidyverse)

#----安装BiocManager相关R包----------------------------------------------
  Biopack <- c("DOSE","clusterProfiler","do","enrichplot",
               "pathview","BiocParallel","GO.db","KEGGREST",
               "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","purrr")
  # human: "org.Hs.eg.db"
  # mouse: "org.Mm.eg.db"
  # rat: "org.Rn.eg.db"
  # purrr： map函数
  is_Biopack <- Biopack[!Biopack %in% all_packages$Package] # 筛选出未安装的R包

  fun_Bioinstall <- function(x) {BiocManager::install(x,update = F,ask = F)} #安装R包函数
  sapply(is_Biopack,fun_Bioinstall,simplify = T)  #用sapply批量安装R包

# 批量library
  fun_library <- function(x){library(x, character.only = T)}
  sapply(Biopack,fun_library,simplify = T)
  }

  inst_packs()


#-----------读取数据，并筛选出差异表达蛋白-------------------------
pro_data_type <- str_extract(protein_data,"(?<=[.]).*")%>% tolower()
if(pro_data_type=="txt")
  pro_df <- read_table(protein_data)[,1] %>% data.frame() %>% na.omit() else
  pro_df <- eval(str2expression(paste0("read.",pro_data_type,"(protein_data)")))[,1] %>%
            data.frame() %>% na.omit()

  colnames(pro_df)[1] <- c("Protein_ID")
  pro_df <- distinct(pro_df,Protein_ID,.keep_all = T)

#----------- 蛋白enrichKEGG通路富集分析---------------------------
  kegg_pro <- enrichKEGG(pro_df[,1], organism ='hsa',
                         #universe,
                         keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                         pvalueCutoff = 0.05,
                         pAdjustMethod = 'BH',
                         qvalueCutoff = 0.2,
                         minGSSize = 3,
                         maxGSSize = 3500,
                         use_internal_data = F)
  if(!is.null(kegg_pro))
    species <- "human" #判断物种

  if(is.null(kegg_pro)){
    {kegg_pro <- enrichKEGG(pro_df[,1], organism ='rno',
                            #universe,
                            keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH',
                            qvalueCutoff = 0.2,
                            minGSSize = 3,
                            maxGSSize = 3500,
                            use_internal_data = F)
    species <- "rat"} #判断物种

    if(is.null(kegg_pro)){
      {kegg_pro <- enrichKEGG(pro_df[,1], organism ='mmu',
                              #universe,
                              keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                              pvalueCutoff = 0.05,
                              pAdjustMethod = 'BH',
                              qvalueCutoff = 0.2,
                              minGSSize = 3,
                              maxGSSize = 3500,
                              use_internal_data = F)
      species <- "mouse" } #判断物种
    }

  }


#---------------以物种为名字建立相应的文件夹--------------------
#species <<- spec # 错误: 无法改变被锁定的联编'species'的值.解决方案是新建环境变量

dir.file <<- dplyr::case_when ( species== "human" ~ "analysis results (human)",
                                 species== "mouse" ~ "analysis results (mouse)",
                                 species== "rat" ~ "analysis results (rat)",
                                 TRUE ~ "analysis results")

  if (dir.exists(dir.file)==FALSE)
    dir.create(dir.file)

#--------------蛋白KEGG通路可视化-----------------------------------
  pathways_kegg <- kegg_pro@result
  gr <- grep("GeneRatio",colnames(pathways_kegg),ignore.case = T)
  colnames(pathways_kegg)[gr] <- c("ProteinRatio")

#-----把GeneRatio分数字符串变小数-------
  num_fun <- function(i){eval(str2expression(i)) %>%
                         as.numeric() %>%
                         round(.,4)
                        }

  pathways_kegg$ProteinRatio <- map(pathways_kegg$ProteinRatio,num_fun)

#-----拆分Description-------------------
pathways_kegg$Description <- str_extract(pathways_kegg$Description,".*(?= -)")

file_path_name <- paste(species,"protein_KEGG_pathways_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(pathways_kegg,file_path_name)

row_n <- nrow(pathways_kegg)

  if(row_n>=30)
    title_text <- c("Top 30 protein KEGG pathways") else
      title_text <- c("Protein KEGG pathways")

  title_size <- case_when(row_n>30 ~12,
                          row_n>20 ~12,
                          TRUE ~11)

  xy_size <- case_when(row_n>30 ~9,
                       row_n>20 ~10,
                       TRUE ~10)

  kegg_mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

  kegg_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
    theme(legend.text=element_text(face="bold",color="black",size=xy_size))

  if(row_n<30)
    path_n <- row_n else
      path_n <- 30

  kegg_df <- pathways_kegg[1:path_n,]

#--------------kegg all pathways------------------------#
kegg_df$Count <- as.numeric(kegg_df$Count)
kegg_df$ProteinRatio <- as.numeric(kegg_df$ProteinRatio)

kegg_pathways <- ggplot(kegg_df)+
    geom_point(aes(x=ProteinRatio,
                   y=fct_reorder(Description,ProteinRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
    labs(x = 'proteinRatio', y = '',title=title_text)+
    kegg_mytheme+kegg_xytheme

kegg_pathways

path_name <- paste(species,"protein_KEGG_pathways.png")
path_name <-paste0(dir.file,"/",path_name)
ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")

#------------kegg metabolism pathways----------------#
#筛选出含有"metabol"字段的行
is.category <- grepl("category",colnames(pathways_kegg),ignore.case = T) %>% any()
if(is.category)
   kegg_meta <- dplyr::filter(pathways_kegg,grepl("metabol",pathways_kegg$category,ignore.case = T)) else
     kegg_meta <- dplyr::filter(pathways_kegg,grepl("metabol",pathways_kegg$Description,ignore.case = T))

meta_path <- grep("Metabolic pathways",pathways_kegg$Description,ignore.case = T)
if(length(meta_path)>0)
  kegg_meta <- kegg_meta[-c(meta_path),] #去掉总的Metabolic pathways

if(nrow(kegg_meta)>0)
  {kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

  kegg_meta$ProteinRatio <- as.numeric(kegg_meta$ProteinRatio)

if(nrow(kegg_meta)>=30)
   {meta_n<- 30
    title_meta_text <-c("Top 30 Protein enriched metabolism pathways")} else
    {meta_n <- nrow(kegg_meta)
    title_meta_text <-c("Protein enriched metabolism pathways")}

kegg_meta_pathways <- ggplot(kegg_meta[c(1:meta_n),])+
    geom_point(aes(x= ProteinRatio,y=reorder(Description,ProteinRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'proteinRatio', y = '',title=title_meta_text)+
    kegg_mytheme+kegg_xytheme

kegg_meta_pathways

meta_path_name <- paste(species,"metabolism pathways.png")
meta_path_name <-paste0(dir.file,"/",meta_path_name)
ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

kegg_meta_file <- dplyr::arrange(kegg_meta,desc(ProteinRatio))
meta_name <- paste(species,"metabolism pathways_data.xlsx")
meta_name <-paste0(dir.file,"/",meta_name)
write.xlsx(kegg_meta_file,meta_name)} else

print("There is no metabolism pathway!")


#------------kegg signaling pathways----------------#
#筛选出含有"signal"字段的行
kegg_signal <- dplyr::filter(pathways_kegg,grepl("signal",pathways_kegg$Description,ignore.case = T))

if(nrow(kegg_signal)>0)
  {kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

if(nrow(kegg_signal)>=30)
{signal_n <- 30
title_signal_text <- c("Top 30 Protein enriched signaling pathways")}else
{signal_n <- nrow(kegg_signal)
  title_signal_text <- c("Protein enriched signaling pathways")}

kegg_signal$ProteinRatio <- as.numeric(kegg_signal$ProteinRatio)

kegg_signal_pathways <- ggplot(kegg_signal[c(1:signal_n),])+
    geom_point(aes(x=ProteinRatio,y=reorder(Description,ProteinRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'proteinRatio', y = '',title=title_signal_text)+
    kegg_mytheme+kegg_xytheme

kegg_signal_pathways

sign_path_name <- paste(species,"Protein_signaling_pathways.png")
sign_path_name <-paste0(dir.file,"/",sign_path_name)
ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

kegg_signal_file <- dplyr::arrange(kegg_signal,desc(ProteinRatio))
sign_name <- paste(species,"Protein_signaling_pathways_data.xlsx")
sign_name <-paste0(dir.file,"/",sign_name)
write.xlsx(kegg_signal_file,sign_name)} else

  print("There is no signaling pathway!")

#-----------全局变量------------------------------------
protein_kegg_pathways_all <<- pathways_kegg

prot_result <- list(species=species)
return(prot_result)

}
