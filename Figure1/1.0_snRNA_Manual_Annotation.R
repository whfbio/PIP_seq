# ------------------------------------------------------------------------------
# Title: Manual Cell Type Annotation Based on Seurat Clusters
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script assigns manual annotations to snRNA cell clusters based on known marker
# genes and cluster relationships. Outputs annotation levels (lv1 to lv4) for
# p0 and p7 timepoints and saves annotated objects.
#
# Dependencies:
# - Seurat
# - dplyr
# - tidyverse
# - cowplot
#
# Usage:
# Load `p0_sct.rds` and `p7_sct.rds` Seurat objects before running.
# ------------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)


wd <- "/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation"
setwd(wd)

celltype.lv1.levels <- c("CM","Fib","Endo","Blood", "Peri", "SMC", "Epic", "Neur", "Mela")

celltype.lv2.levels <- c("v_CM", "v.d_CM", "v.Slit2_CM", "a_CM","p_CM",
                         "Fib", "Fib.Kirrel3",
                         "Endo.Nrxn3","Endo.Dach1","Blood",
                         "Peri", "SMC",
                         "Epic", "Neur", "Mela")                       

celltype.lv3.levels <- c("v_CM", "v.d_CM", "v.Slit2_CM", 
                         "a.right_CM","a.left_CM","av_CM",
                         "p_CM",
                         "Fib","Fib.Kirrel3",
                         "Endo.Nrxn3","Endo.Dach1",
                         "Blood",
                         "Epic",
                         "Peri", "SMC", 
                         "Neur", "Mela") 

cell_annot_order <- c("CM", "v_CM", "v.d_CM", "v.Slit2_CM","v.Ddc_CM", "v.Gm12381_CM",
                      "a_CM", "a.right_CM","a.left_CM","av_CM",
                      "p_CM", "p.G2_CM","p.M_CM",
                      "Fib", "Fib.Kirrel3","Fib.Gfpt2","Fib.Stng1", "Fib.p",
                      "Endo","Endo.Nrxn3","Endo.Dach1",
                      "Epic",
                      "Blood",
                      "Peri","Peri.p",
                      "SMC", 
                      "Neur", 
                      "Mela")


#Seurat_obj_p0 <- readRDS("./IKAP/p0_IKAP/Seurat_obj_p0.rds") 
#p0[["seurat_clusters"]] <- Idents(Seurat_obj_p0)

#this file we redo annotations, not create them
#the methodology and purpose of this file has varied overtime.
#currently it is intended to simply give clusters the celltype annotation and their sub-notations
#we have learned the significant clusters over large time and this file helps organize this complicated and messy process.
#the various levels, with level 1 being the broadest, are the classifications and sub classification.
#this file also serves to document some of the choices for certain marcation of specific clusters based off their marker gene

#p0
###############################
p0             <- readRDS("./rObjs/SCT_applied/p0_sct.rds")

#p0_celltype.df <- readRDS("./Manual_Annotation/prev_man_annot/12.14.2022/p0_celltype.df.rds")

p0 <- FindSubCluster(p0, cluster = 12, graph.name = "SCT_nn", resolution = .5,subcluster.name =  "sub_12")
p0 <- FindSubCluster(p0, cluster = 19, graph.name = "SCT_nn", resolution = .5,subcluster.name =  "sub_19")
p0 <- FindSubCluster(p0, cluster = 8, graph.name = "SCT_nn", resolution = .3,subcluster.name =  "sub_8")
                     

p0_celltype.df <- data.frame( seurat_clusters = Idents(p0)) #the idents should be the seurat clusters

p0_celltype.df <- p0_celltype.df %>%
  mutate( 
    celltype.lv1 =case_when
    (
      seurat_clusters %in% c(8,4,5,1,6,0,2,10,17,11) ~ "CM", 
      seurat_clusters %in% c(13,18,14,3)     ~ "Fib",      
      seurat_clusters %in% c(16)        ~ "Blood",    
      seurat_clusters %in% c(7,9)        ~ "Endo", 
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_12 %in% c("12_0","12_1")]   ~ "Peri",    
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_12 %in% c("12_2","12_3")]   ~ "SMC",     
      seurat_clusters %in% c(15)        ~ "Epic",       
      row.names(p0_celltype.df) %in% colnames(p0)[p0$sub_19 %in% c("19_0")]        ~ "Neur",       #FeaturePlot(p0, c("Kcnh8"), label = TRUE)
      row.names(p0_celltype.df) %in% colnames(p0)[p0$sub_19 %in% c("19_1")]        ~ "Mela"        #FeaturePlot(p0, c("Oca2"), label = TRUE)
    
    )) %>%
  mutate( 
    celltype.lv2 =case_when
    (
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_8 %in% c("8_0", "8_1","8_2")]   ~ "a_CM",  
      
      seurat_clusters %in% c(10) ~ "v.Slit2_CM",
      seurat_clusters %in% c(5)  ~ "v.d_CM",
      seurat_clusters %in% c(11,17) ~ "p_CM",
      seurat_clusters %in% c(6,4,0,2,1)  ~ "v_CM",
     

      seurat_clusters %in% c(13) ~ "Fib.Kirrel3",
      
      seurat_clusters %in% c(7)  ~ "Endo.Nrxn3",
      seurat_clusters %in% c(9)  ~ "Endo.Dach1",
      
      TRUE ~ celltype.lv1
    )) %>%
  mutate( 
    celltype.lv3 =case_when(
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_8 %in% c("8_2")]   ~ "av_CM", 
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_8 %in% c("8_0")]   ~ "a.left_CM",  
      row.names(p0_celltype.df) %in%   colnames(p0)[p0$sub_8 %in% c("8_1")]   ~ "a.right_CM",  
      
      TRUE ~ celltype.lv2
    )) %>%
  mutate( 
    celltype.lv4 =case_when(
      seurat_clusters %in% c(4)          ~ "v.Ddc_CM",
      seurat_clusters %in% c(18)         ~ "Fib.p",
      seurat_clusters %in% c(11)         ~ "p.M_CM",
      seurat_clusters %in% c(17)         ~ "p.G2_CM",
      TRUE ~ celltype.lv3
    ))


Idents(p0) <- "celltype.lv2"

p0[["celltype.lv1"]] <- p0_celltype.df$celltype.lv1
p0[["celltype.lv2"]] <- p0_celltype.df$celltype.lv2
p0[["celltype.lv3"]] <- p0_celltype.df$celltype.lv3
p0[["celltype.lv4"]] <- p0_celltype.df$celltype.lv4

DimPlot(p0, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p0, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p0, group.by = "celltype.lv3", label = TRUE) 
DimPlot(p0, group.by = "celltype.lv4", label = TRUE) 




p0[["seurat_clusters"]] <-  p0_celltype.df$seurat_clusters
p0[["celltype.lv1"]] <-  factor(p0_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p0[["celltype.lv2"]] <-  factor(p0_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p0[["celltype.lv3"]] <-  factor(p0_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p0[["celltype.lv4"]] <-  p0_celltype.df$celltype.lv4


p0 <- subset(p0, subset = celltype.lv1 != c("Mela"))
p0_celltype.df <- p0_celltype.df[p0_celltype.df$celltype.lv1 != "Mela",]

p0[["celltype.lv1"]] <-  factor(p0_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p0[["celltype.lv2"]] <-  factor(p0_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p0[["celltype.lv3"]] <-  factor(p0_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p0[["celltype.lv4"]] <-  p0_celltype.df$celltype.lv4


p0_CM <-subset(p0, subset = celltype.lv1 =="CM")

p0_CM <- subset(p0, subset = celltype.lv1 %in% c("CM"))
p0_CM <- SCTransform(p0_CM, method = "glmGamPoi")
p0_CM <- RunPCA(p0_CM)
p0_CM <- RunUMAP(p0_CM, dims = 1:15)
p0_CM <- FindNeighbors(p0_CM, dims = 1:15)
p0_CM <- FindClusters(p0_CM, resolution = .8)
DimPlot(p0_CM, group.by ="orig.ident")
DimPlot(p0_CM, label = TRUE)
DimPlot(p0_CM, label = TRUE, group.by = "celltype.lv4")


pdf("./Manual_Annotation/p0_celltype_plots.pdf", width = 10, onefile = TRUE)
DimPlot(p0,   group.by = "celltype.lv1", label = TRUE) 
DimPlot(p0,   group.by = "celltype.lv2", label = TRUE) 
DimPlot(p0,   group.by = "celltype.lv3", label = TRUE)
DimPlot(p0,   group.by = "celltype.lv4", label = TRUE)
DimPlot(p0,   group.by = "seurat_clusters"   , label = TRUE)

DimPlot(p0_CM, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p0_CM, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p0_CM, group.by = "celltype.lv3", label = TRUE)
DimPlot(p0_CM, group.by = "celltype.lv4", label = TRUE)
DimPlot(p0_CM, group.by = "seurat_clusters", label = TRUE)
dev.off()


saveRDS(p0_celltype.df,"./Manual_Annotation/p0_celltype.df.rds")
saveRDS(p0            ,"./rObjs/Man_annot/p0_man.rds")
saveRDS(p0_CM            ,"./rObjs/Man_annot/p0_CM_man.rds")




DimPlot(p0_CM, group.by = "seurat_clusters", label = TRUE)
DimPlot(p0_CM, group.by = "celltype.lv4", label = TRUE) 

#p7
####################################
p7             <- readRDS("./rObjs/SCT_applied/p7_sct.rds")

p7_celltype.df <- data.frame(seurat_clusters = p7$seurat_clusters)

p7 <- FindSubCluster(p7, cluster = 17, graph.name = "SCT_nn", resolution = .5,subcluster.name =  "sub_17")
p7 <- FindSubCluster(p7, cluster = 10, graph.name = "SCT_snn", resolution = .4,subcluster.name =  "sub_10")

p7.CM_ava <- readRDS("./Manual_Annotation/Manual_Detection/p7.CM_ava.rds")
p7.CM     <- readRDS("./Manual_Annotation/Manual_Detection/p7.CM.rds")

DimPlot(p7, label = TRUE)
p7_celltype.df <- p7_celltype.df %>%
  mutate( 
    celltype.lv1 =case_when
    (
      seurat_clusters %in% c(0,11,4,2,12,10) ~ "CM", 
      seurat_clusters %in% c(6,1,15,14)    ~ "Fib",      
      seurat_clusters %in% c(7)          ~ "Blood",    
      seurat_clusters %in% c(8,3,9)      ~ "Endo", 
      seurat_clusters %in% c(5)          ~ "Peri",     
      seurat_clusters %in% c(13)         ~ "SMC",      
      seurat_clusters %in% c(16)         ~ "Epic",  
      row.names(p7_celltype.df) %in%   colnames(p7)[p7$sub_17 %in% c("17_0")]   ~ "Neur",    
      row.names(p7_celltype.df) %in%   colnames(p7)[p7$sub_17 %in% c("17_1")]   ~ "Mela",   
        )) %>%
  mutate( 
    celltype.lv2 =case_when(

      seurat_clusters %in% c(10) ~ "a_CM",
      seurat_clusters %in% c(11,4)  ~ "p_CM",
      seurat_clusters %in% c(12) ~ "v.Slit2_CM",
      row.names(p7_celltype.df) %in% names(p7.CM)[p7.CM==1]  ~ "v.d_CM",
      celltype.lv1 %in% c("CM") ~ "v_CM",
      
      seurat_clusters %in% c(14)    ~ "Fib.Kirrel3",
      seurat_clusters %in% c(4,18)  ~ "Fib",
      
      seurat_clusters %in% c(8)  ~ "Endo.Nrxn3",
      seurat_clusters %in% c(3,9)   ~ "Endo.Dach1", #higher prolif here 
      
      TRUE ~ celltype.lv1
    )) %>%
  mutate( 
    celltype.lv3 =case_when(
      
      row.names(p7_celltype.df) %in% names(p7.CM_ava)[p7.CM_ava==1] ~ "a.left_CM" ,
      row.names(p7_celltype.df) %in% names(p7.CM_ava)[p7.CM_ava==0] ~ "a.right_CM" ,
      row.names(p7_celltype.df) %in% names(p7.CM_ava)[p7.CM_ava==2] ~ "av_CM" ,
      
      
      TRUE ~ celltype.lv2
      
      
    )) %>%
  mutate( 
    celltype.lv4 =case_when(
      seurat_clusters %in% c(15)         ~ "Fib.p",
      row.names(p7_celltype.df) %in% names(p7.CM)[p7.CM==8]  ~  "v.Gm12381_CM",
      row.names(p7_celltype.df) %in% names(p7.CM)[p7.CM==7]   ~  "v.Ddc_CM",
      TRUE ~ celltype.lv3
    ))


Idents(p7) <- "celltype.lv2"

#p7_celltype.df <- readRDS("./Manual_Annotation/prev_man_annot/12.14.2022/p7_celltype.df.rds")
p7_celltype.df <- p7_celltype.df[ colnames(p7),] 



p7[["celltype.lv1"]] <- p7_celltype.df$celltype.lv1
p7[["celltype.lv2"]] <- p7_celltype.df$celltype.lv2
p7[["celltype.lv3"]] <- p7_celltype.df$celltype.lv3
p7[["celltype.lv4"]] <- p7_celltype.df$celltype.lv4

DimPlot(p7, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p7, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p7, group.by = "celltype.lv3", label = TRUE)
DimPlot(p7, group.by = "celltype.lv4", label = TRUE)  


p7[["seurat_clusters"]] <-  p7_celltype.df$seurat_clusters
p7[["celltype.lv1"]] <-  factor(p7_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p7[["celltype.lv2"]] <-  factor(p7_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p7[["celltype.lv3"]] <-  factor(p7_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p7[["celltype.lv4"]] <-  p7_celltype.df$celltype.lv4

Idents(p7) <- "celltype.lv2"


p7 <- subset(p7, subset = celltype.lv1 != c("Mela"))

p7_celltype.df <- p7_celltype.df[p7_celltype.df$celltype.lv1 != "Mela",]
p7[["celltype.lv1"]] <-  factor(p7_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p7[["celltype.lv2"]] <-  factor(p7_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p7[["celltype.lv3"]] <-  factor(p7_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p7[["celltype.lv4"]] <-  p7_celltype.df$celltype.lv4




p7_CM <- subset(p7, subset = celltype.lv1 %in% c("CM"))
p7_CM <- SCTransform(p7_CM, method = "glmGamPoi",vars.to.regress = c("version"))
p7_CM <- RunPCA(p7_CM)
p7_CM <- RunUMAP(p7_CM, dims = 1:30)
p7_CM <- FindNeighbors(p7_CM, dims = 1:30)
p7_CM <- FindClusters(p7_CM, resolution = .8)
DimPlot(p7_CM, group.by ="version")
DimPlot(p7_CM, label = TRUE)
DimPlot(p7_CM, label = TRUE, group.by = "celltype.lv4")

Idents(p7_CM) <- "celltype.lv4"

pdf("./Manual_Annotation/p7_celltype_plots.pdf", width = 10, onefile = TRUE)
DimPlot(p7,   group.by = "celltype.lv1", label = TRUE) 
DimPlot(p7,   group.by = "celltype.lv2", label = TRUE) 
DimPlot(p7,   group.by = "celltype.lv3", label = TRUE) 
DimPlot(p7,   group.by = "celltype.lv4", label = TRUE)
DimPlot(p7,   group.by = "seurat_clusters"   , label = TRUE)

DimPlot(p7_CM, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p7_CM, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p7_CM, group.by = "celltype.lv3", label = TRUE)
DimPlot(p7_CM, group.by = "celltype.lv4", label = TRUE)
DimPlot(p7_CM, group.by = "seurat_clusters", label = TRUE)
dev.off()


saveRDS(p7_celltype.df,"./Manual_Annotation/p7_celltype.df.rds")
saveRDS(p7            ,"./rObjs/Man_annot/with_melanocytes/p7_man.rds")
saveRDS(p7_CM            ,"./rObjs/Man_annot/p7_CM_man.rds")






#p14
####################################
p14             <- readRDS("./rObjs/SCT_applied/p14_sct.rds")

p14_celltype.df <- data.frame(seurat_clusters = p14$seurat_clusters)

#p14_celltype.df <- readRDS("./Manual_Annotation/prev_man_annot/12.14.2022/p14_celltype.df.rds")
#p14_celltype.df <- p14_celltype.df[ colnames(p14),]
#p14_celltype.df$seurat_clusters <-  p14$seurat_clusters

p14_mela       <- readRDS("./Manual_Annotation/Manual_Detection/p14_mela.rds")
p14.a.left_CM  <- readRDS("./Manual_Annotation/Manual_Detection/p14.a.left_CM.rds")
p14.a.right_CM <- readRDS("./Manual_Annotation/Manual_Detection/p14.a.right_CM.rds")
p14.av_CM      <- readRDS("./Manual_Annotation/Manual_Detection/p14.av_CM.rds")
p14.v.Gm12381_CM <- readRDS("./Manual_Annotation/Manual_Detection/p14.v.Gm12381_CM.rds")



DimPlot(p14, label = TRUE) + DimPlot(p14, label = TRUE, group.by = "celltype.lv2")

p14_celltype.df <- p14_celltype.df %>%
  mutate( 
    celltype.lv1 =case_when
    (
      row.names(p14_celltype.df) %in% p14_mela  ~ "Mela",
      seurat_clusters %in% c(11,12,18,4,7,0,2,3,10) ~ "CM", #FeaturePlot(p14, c("Myh7b","Myl7"), label = TRUE)
      seurat_clusters %in% c(8,21,1,6,17)     ~ "Fib",      
      seurat_clusters %in% c(13)        ~ "Blood",    
      seurat_clusters %in% c(15,9)     ~ "Endo",      
      seurat_clusters %in% c(5,20)         ~ "Peri",    
      seurat_clusters %in% c(14)        ~ "SMC",     
      seurat_clusters %in% c(16)        ~ "Epic",       
      seurat_clusters %in% c(19)        ~ "Neur"              #FeaturePlot(p14, c("Kcnh8"), label = TRUE)
    )) %>%
  mutate( 
    celltype.lv2 =case_when(
      
      row.names(p14_celltype.df) %in% c(p14.a.left_CM, p14.a.right_CM, p14.av_CM) ~ "a_CM",
      seurat_clusters %in% c(10) ~ "v.d_CM",
      seurat_clusters %in% c(12) ~ "v.Slit2_CM",
      seurat_clusters %in% c(18) ~ "p_CM",
      
      celltype.lv1    %in% c("CM") ~ "v_CM",
      
      seurat_clusters %in% c(17) ~ "Fib.Kirrel3",
      
      row.names(p14_celltype.df) %in% p14_mela  ~ "Mela",
      seurat_clusters %in% c(15) ~ "Endo.Nrxn3",
      seurat_clusters %in% c(9)  ~ "Endo.Dach1",

      
      TRUE ~ celltype.lv1
    )) %>%
  mutate( 
    celltype.lv3 =case_when(
      
      row.names(p14_celltype.df) %in% c(p14.a.left_CM)  ~ "a.left_CM",
      row.names(p14_celltype.df) %in% c(p14.a.right_CM) ~ "a.right_CM",
      row.names(p14_celltype.df) %in% c(p14.av_CM)      ~ "av_CM",
      
      TRUE ~ celltype.lv2
    )) %>%
  mutate( 
    celltype.lv4 =case_when(
      
      row.names(p14_celltype.df) %in% p14.v.Gm12381_CM ~ "v.Gm12381_CM",
      
      seurat_clusters %in% c(21)         ~ "Fib.p",
      seurat_clusters %in% c(20)         ~ "Peri.p",
      TRUE ~ celltype.lv3
    ))




p14[["celltype.lv1"]] <- p14_celltype.df$celltype.lv1
p14[["celltype.lv2"]] <- p14_celltype.df$celltype.lv2
p14[["celltype.lv3"]] <- p14_celltype.df$celltype.lv3
p14[["celltype.lv4"]] <- p14_celltype.df$celltype.lv4

DimPlot(p14, group.by = "celltype.lv1", label = TRUE)  +DimPlot(p14, label = TRUE)
DimPlot(p14, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p14, group.by = "celltype.lv3", label = TRUE)
DimPlot(p14, group.by = "celltype.lv4", label = TRUE)  


p14[["seurat_clusters"]] <-  p14_celltype.df$seurat_clusters
p14[["celltype.lv1"]] <-  factor(p14_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p14[["celltype.lv2"]] <-  factor(p14_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p14[["celltype.lv3"]] <-  factor(p14_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p14[["celltype.lv4"]] <-  p14_celltype.df$celltype.lv4

Idents(p14) <- "celltype.lv2"

DimPlot(p14, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p14, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p14, group.by = "celltype.lv3", label = TRUE) 
DimPlot(p14, group.by = "celltype.lv4", label = TRUE) 


p14 <- subset(p14, subset = celltype.lv1 != c("Mela"))

p14_celltype.df <- p14_celltype.df[p14_celltype.df$celltype.lv1 != "Mela",]
p14[["celltype.lv1"]] <-  factor(p14_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p14[["celltype.lv2"]] <-  factor(p14_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p14[["celltype.lv3"]] <-  factor(p14_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p14[["celltype.lv4"]] <-  p14_celltype.df$celltype.lv4


p14_CM <- subset(p14, subset = celltype.lv1 %in% c("CM"))
p14_CM <- SCTransform(p14_CM, method = "glmGamPoi",vars.to.regress = c("orig.ident"))
p14_CM <- RunPCA(p14_CM)
p14_CM <- RunUMAP(p14_CM, dims = 1:15)
p14_CM <- FindNeighbors(p14_CM, dims = 1:15)
p14_CM <- FindClusters(p14_CM, resolution = .8)
DimPlot(p14_CM, group.by ="version")
DimPlot(p14_CM, label = TRUE)
DimPlot(p14_CM, label = TRUE, group.by = "orig.ident")
DimPlot(p14_CM, label = TRUE, group.by = "celltype.lv4")



saveRDS(p14_celltype.df,"./Manual_Annotation/p14_celltype.df.rds")
saveRDS(p14            ,"./rObjs/Man_annot/p14_man.rds")
saveRDS(p14_CM         ,"./rObjs/Man_annot/p14_CM_man.rds")




pdf("./Manual_Annotation/p14_celltype_plots.pdf", width = 10, onefile = TRUE)
DimPlot(p14,   group.by = "celltype.lv1", label = TRUE) 
DimPlot(p14,   group.by = "celltype.lv2", label = TRUE) 
DimPlot(p14,   group.by = "celltype.lv3", label = TRUE)
DimPlot(p14,   group.by = "celltype.lv4", label = TRUE)
DimPlot(p14,   group.by = "seurat_clusters"   , label = TRUE)

DimPlot(p14_CM, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p14_CM, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p14_CM, group.by = "celltype.lv3", label = TRUE)
DimPlot(p14_CM, group.by = "celltype.lv4", label = TRUE)
DimPlot(p14_CM, group.by = "seurat_clusters", label = TRUE)
dev.off()




#p21
##########################################################
p21             <- readRDS("./rObjs/SCT_applied/p21_sct.rds")


p21_celltype.df <- p21_celltype.df[ colnames(p21),]



p21.v.Slit2_CM <- readRDS("./Manual_Annotation/Manual_Detection/p21.v.Slit2_CM.rds")
p21.a.left_CM <- readRDS("./Manual_Annotation/Manual_Detection/p21.a.left_CM.rds")
p21.a.right_CM <- readRDS("./Manual_Annotation/Manual_Detection/p21.a.right_CM.rds")
p21.av_CM <- readRDS("./Manual_Annotation/Manual_Detection/p21.av_CM.rds")


#p21 <- FindSubCluster(p21, cluster = 11, graph.name = "SCT_snn", resolution = .4,subcluster.name =  "sub_11")
#DimPlot(p21, group.by = "sub_11", label = TRUE)


DimPlot(p21, label = TRUE) + DimPlot(p21, group.by = "celltype.lv3", label = TRUE)

p21_celltype.df <- data.frame( seurat_clusters = Idents(p21))

p21_celltype.df <- p21_celltype.df %>%
  mutate( 
    celltype.lv1 =case_when
    (
      seurat_clusters %in% c(11,14,18,4,5,1,2,3,6) ~ "CM",#FeaturePlot(p21, c("Myh7b","Myl7"), label = TRUE)
      seurat_clusters %in% c(0,10,8,17) ~ "Fib"  ,        #FeaturePlot(p21, c("Postn","Fbn1", "Col5a1"), label = TRUE)
      seurat_clusters %in% c(12)        ~ "Blood",        #FeaturePlot(p21, c("Mrc1","Fyb", "Cd86"), label = TRUE)
      seurat_clusters %in% c(13,19,9)   ~ "Endo" ,        #FeaturePlot(p21, c("Pecam1","Flt1", "Tie1"), label = TRUE)
      seurat_clusters %in% c(7)         ~ "Peri" ,        #FeaturePlot(p21, c("Pdgfrb","Notch3","Cspg4"), label = TRUE)
      seurat_clusters %in% c(16)        ~ "SMC"  ,     
      seurat_clusters %in% c(15)        ~ "Epic" ,       
      seurat_clusters %in% c(20)        ~ "Neur"          #FeaturePlot(p21, c("Kcnh8"), label = TRUE)
      #seurat_clusters %in% c(25)       ~ "Mela" ,        #FeaturePlot(p21, c("Oca2"), label = TRUE)
    )) %>%
  mutate( 
    celltype.lv2 =case_when(

      seurat_clusters %in% c(17) ~ "Fib.Kirrel3",
      
      seurat_clusters %in% c(13) ~ "Endo.Nrxn3",
      seurat_clusters %in% c(9,19) ~ "Endo.Dach1",
      
      row.names(p21_celltype.df) %in% c(p21.a.left_CM,p21.a.right_CM,p21.av_CM )  ~ "a_CM",

      
      row.names(p21_celltype.df) %in%  p21.v.Slit2_CM  ~ "v.Slit2_CM",
      seurat_clusters %in% c(4)  ~ "v.d_CM",
      celltype.lv1    %in% c("CM") ~ "v_CM",
      
      TRUE ~ celltype.lv1
    )) %>%
  mutate( 
    celltype.lv3 =case_when(
      
      row.names(p21_celltype.df) %in%  p21.a.left_CM   ~ "a.left_CM",
      row.names(p21_celltype.df) %in%  p21.a.right_CM   ~ "a.right_CM", 
      row.names(p21_celltype.df) %in%  p21.av_CM       ~ "av_CM", 
    
    
      TRUE ~ celltype.lv2
    )) %>%
  mutate( 
    celltype.lv4 =case_when(
      
      seurat_clusters %in% c(3)  ~ "v.Gm12381_CM",
    
      TRUE ~ celltype.lv3
    ))



p21[["celltype.lv1"]] <- p21_celltype.df$celltype.lv1
p21[["celltype.lv2"]] <- p21_celltype.df$celltype.lv2
p21[["celltype.lv3"]] <- p21_celltype.df$celltype.lv3
p21[["celltype.lv4"]] <- p21_celltype.df$celltype.lv4


p14_celltype.df <- readRDS("/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation/Manual_Annotation/p14_celltype.df.rds")
p14[["seurat_clusters"]] <-  p14_celltype.df$seurat_clusters
p14[["celltype.lv1"]] <-  factor(p14_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p14[["celltype.lv2"]] <-  factor(p14_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p14[["celltype.lv3"]] <-  factor(p14_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p14[["celltype.lv4"]] <-  p14_celltype.df$celltype.lv4

p21_celltype.df <- readRDS("/proj/liulab/users/mdcolon/Ref_Data_Analysis_PerturbSeq/Data_Set_Creation/Manual_Annotation/p21_celltype.df.rds")
p21[["seurat_clusters"]] <-  p21_celltype.df$seurat_clusters
p21[["celltype.lv1"]] <-  factor(p21_celltype.df$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
p21[["celltype.lv2"]] <-  factor(p21_celltype.df$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
p21[["celltype.lv3"]] <-  factor(p21_celltype.df$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
p21[["celltype.lv4"]] <-  p21_celltype.df$celltype.lv4


DimPlot(p21, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p21, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p21, group.by = "celltype.lv3", label = TRUE) 
DimPlot(p21, group.by = "celltype.lv4", label = TRUE) 

Idents(p21) <- "celltype.lv2"


p21_CM <- subset(p21, subset = celltype.lv1 %in% c("CM"))
p21_CM <- SCTransform(p21_CM, method = "glmGamPoi") %>% RunPCA()
p21_CM <- FindNeighbors(p21_CM, dims = 1:15)
p21_CM <- FindClusters(p21_CM, resolution = .5)
p21_CM <- RunUMAP(p21_CM, dims = 1:15)
DimPlot(p21_CM, label = TRUE) +DimPlot(p21_CM, label = TRUE, group.by = "celltype.lv4")

Idents(p21_CM) <- "celltype.lv3"



pdf("./Manual_Annotation/p21_celltype_plots.pdf", width = 10, onefile = TRUE)
DimPlot(p21, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p21, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p21, group.by = "celltype.lv3", label = TRUE)
DimPlot(p21, group.by = "celltype.lv4", label = TRUE)
DimPlot(p21, group.by = "seurat_clusters", label = TRUE)

DimPlot(p21_CM, group.by = "celltype.lv1", label = TRUE) 
DimPlot(p21_CM, group.by = "celltype.lv2", label = TRUE) 
DimPlot(p21_CM, group.by = "celltype.lv3", label = TRUE)
DimPlot(p21_CM, group.by = "celltype.lv4", label = TRUE)
DimPlot(p21_CM, group.by = "seurat_clusters", label = TRUE)
dev.off()


#############finally
saveRDS(p0_celltype.df,"./Manual_Annotation/p0_celltype.df.rds")
saveRDS(p7_celltype.df,"./Manual_Annotation/p7_celltype.df.rds")
saveRDS(p14_celltype.df,"./Manual_Annotation/p14_celltype.df.rds")
saveRDS(p21_celltype.df,"./Manual_Annotation/p21_celltype.df.rds")


saveRDS(p0             ,"./rObjs/Man_annot/with_melanocytes/p0_man.rds")
saveRDS(p7             ,"./rObjs/Man_annot/with_melanocytes/p7_man.rds")
saveRDS(p14            ,"./rObjs/Man_annot/with_melanocytes/p14_man.rds")
saveRDS(p21            ,"./rObjs/Man_annot/with_melanocytes/p21_man.rds")

#we remove the Melanocytes because we do not use them in following analysis
p0 <- subset(p0, subset = celltype.lv1 != c("Mela"))
p7 <- subset(p7, subset = celltype.lv1 != c("Mela"))
p14 <- subset(p14, subset = celltype.lv1 != c("Mela"))
p21 <- subset(p21, subset = celltype.lv1 != c("Mela"))

saveRDS(p0            ,"./rObjs/Man_annot/p0_man.rds")
saveRDS(p7            ,"./rObjs/Man_annot/p7_man.rds")
saveRDS(p14            ,"./rObjs/Man_annot/p14_man.rds")
saveRDS(p21            ,"./rObjs/Man_annot/p21_man.rds")



saveRDS(p0_CM            ,"./rObjs/Man_annot/p0_CM_man.rds")
saveRDS(p7_CM            ,"./rObjs/Man_annot/p7_CM_man.rds")
saveRDS(p14_CM            ,"./rObjs/Man_annot/p14_CM_man.rds")
saveRDS(p21_CM            ,"./rObjs/Man_annot/p21_CM_man.rds")

####the final step is to combine all these datasets 

#note that these 
p0  <-readRDS("./rObjs/Man_annot/p0_man.rds" )
p7  <-readRDS("./rObjs/Man_annot/p7_man.rds" )
p14 <-readRDS("./rObjs/Man_annot/p14_man.rds")
p21 <-readRDS("./rObjs/Man_annot/p21_man.rds")


p0$time_point  <- "p0"
p7$time_point  <- "p7"
p14$time_point <- "p14"
p21$time_point <- "p21"

ref <- merge(p0, y = c(p7,p14,p21), project = "PSeq_Ref_Data")
ref <- SCTransform(ref, method = "glmGamPoi")
ref <- RunPCA(ref)
ref <- RunUMAP(ref, dims = 1:20)
ref <- FindNeighbors(ref, dims = 1:15)
ref <- FindClusters(ref, resolution = .8)

ref[["celltype.lv1"]] <-  factor(ref$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
ref[["celltype.lv2"]] <-  factor(ref$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
ref[["celltype.lv3"]] <-  factor(ref$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
ref[["celltype.lv4"]] <-  ref$celltype.lv4
ref$time_point <- factor(ref$time_point,levels = c("p0","p7","p14","p21"))


DimPlot(ref)
DimPlot(ref, group.by = "orig.ident")
DimPlot(ref, group.by = "version")
DimPlot(ref, group.by = "time_point")
DimPlot(ref, group.by = "time_point") + DimPlot(ref, group.by = "celltype.lv2", label= TRUE)
DimPlot(ref, group.by = "celltype.lv3", label= TRUE)


saveRDS(ref, "./rObjs/Man_annot/ref_man.rds")


ref_CM <- subset(ref, subset = celltype.lv1 =="CM")
ref_CM <- SCTransform(ref_CM, method = "glmGamPoi")
ref_CM <- RunPCA(ref_CM)
ref_CM <- RunUMAP(ref_CM, dims = 1:15)
ref_CM <- FindNeighbors(ref_CM, dims = 1:15)
ref_CM <- FindClusters(ref_CM, resolution = .8)



DimPlot(ref_CM, label = TRUE)
DimPlot(ref_CM, group.by = "orig.ident")
DimPlot(ref_CM, group.by = "celltype.lv1") + DimPlot(ref_CM, group.by = "celltype.lv2") + DimPlot(ref_CM, group.by = "celltype.lv3")
DimPlot(ref_CM, group.by = "time_point") + DimPlot(ref_CM, group.by = "celltype.lv3") +DimPlot(ref_CM, group.by = "celltype.lv4")



ref_CM$time_point<- factor(ref_CM$time_point, levels = c("p0","p7","p14","p21"))
ref_CM[["celltype.lv1"]] <-  factor(ref_CM$celltype.lv1,  levels= celltype.lv1.levels) %>% droplevels()
ref_CM[["celltype.lv2"]] <-  factor(ref_CM$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()
ref_CM[["celltype.lv3"]] <-  factor(ref_CM$celltype.lv3,  levels= celltype.lv3.levels) %>% droplevels()
ref_CM[["celltype.lv4"]] <-  ref_CM$celltype.lv4


saveRDS(ref_CM, "./rObjs/Man_annot/ref_CM_man.rds")

ref_CM@reductions$umap <-ref_CM_man@reductions$umap


####The final annotation is based on the addtional spatial info provided from the xenium data from the corresponding timepoint
ref<-readRDS("./robjs/ref_lv5.rds")

old_annotation <- c("p_CM", "Endo.2", "v_CM", "Fib.4", "v.Ddc_CM", "Fib.1", "Fib.5", "Epic", "v.Slit2_CM", 
               "v.d_CM", "a.right_CM", "SMC", "a.left_CM", "Endo.0", "Fib.7", "Peri", "Blood", 
               "av_CM", "Endo.5", "Endo.4", "Endo.3", "Neur", "Fib.6", "Endo.1", "Fib.2", 
               "Endo.6", "Fib.0", "Fib.8","Fib.3", "v.Gm12381_CM", "Peri.p")

new_annotation <- c("p.CM", "cap.EC", "v.CM", "Postn.Fib", "Ddc.CM", "Postn.Fib", "valve.Fib", "Epic", "Slit2.CM", "Ankrd1.CM", 
                "a.CM", "SMC", "a.CM", "endocardial.EC", "p.Fib", "Peri", "Blood", "av.CM", "arteriole.EC", "p.EC", "cap.EC(High.Myo)", 
                "Neur", "Gfpt2.Fib", "cap.EC", "Col8a1.Fib", "lym.EC", "Lsamp.Fib", "Fib.8","Gfpt2.Fib","v.CM","p.Peri")


annotation_df <- data.frame(old_names = as.character(old_annotation), new_names = new_annotation)

meta<-ref@meta.data
df<-meta[,c("celltype.lv4","celltype.lv5")]
df <- df %>%
  left_join(annotation_df, by = c("celltype.lv5" = "old_names"))
# Rename the new_names column to man_annotation
df <- df %>%
  rename(final_annotation = new_names)
ref<-AddMetaData(ref, metadata = df, col.name = c("final_annotation"))
ref.subset<-subset(ref, final_annotation != "Fib.8")
ncol(ref.subset)
saveRDS(ref.subset,"./robjs/ref_man_final.rds")


#END

#####################
#
#the following scripts below are not necessary for producing the annotation
#they are maintained here for their utility for answering specific questions throughtou the generation of these annotaions

#####################

cell_annot_order <- c("CM", "v_CM", "v.d_CM", "v.Slit2_CM","v.Ddc_CM", "v.Gm12381_CM",
                      "a_CM", "a.right_CM","a.left_CM","av_CM",
                      "p_CM", "p.G2_CM","p.M_CM",
                      "Fib", "Fib.Kirrel3","Fib.Gfpt2","Fib.Stng1", "Fib.p",
                      "Endo","Endo.Nrxn3","Endo.Dach1",
                      "Epic",
                      "Blood",
                      "Peri","Peri.p", "SMC", 
                      "Neur", "Mela")


celltype.lv2.5 <- ref$celltype.lv2
celltype.lv4.5 <- ref$celltype.lv4

celltype.lv2.5[names(celltype.lv2.5) %in% p_Fib] <- "Fib.p"
celltype.lv2.5[names(celltype.lv2.5) %in% Fib.Kirrel3] <- "Fib.Krrel3"
celltype.lv2.5[names(celltype.lv2.5) %in% Fib.Gfpt2] <- "Fib.Gfpt2"
celltype.lv2.5[names(celltype.lv2.5) %in% Fib.Stng1] <- "Fib.Stng1"

celltype.lv2.5 <- celltype.lv2.5 %>% factor(levels= cell_annot_order) %>% droplevels()


celltype.lv4.5[names(celltype.lv4.5) %in% p_Fib] <- "Fib.p"
celltype.lv4.5[names(celltype.lv4.5) %in% Fib.Kirrel3] <- "Fib.Krrel3"
celltype.lv4.5[names(celltype.lv4.5) %in% Fib.Gfpt2] <- "Fib.Gfpt2"
celltype.lv4.5[names(celltype.lv4.5) %in% Fib.Stng1] <- "Fib.Stng1"

celltype.lv4.5 <- celltype.lv4.5 %>% factor(levels= cell_annot_order) %>% droplevels()

saveRDS(celltype.lv2.5, "./Data_Set_Creation/Manual_Annotation/Manual_Detection/celltype.lv2.5.rds")
saveRDS(celltype.lv4.5, "./Data_Set_Creation/Manual_Annotation/Manual_Detection/celltype.lv4.5.rds")
saveRDS(cell_annot_order, "./Data_Set_Creation/Manual_Annotation/cell_annot_order.rds")














common <- intersect(unique(p0$celltype.lv3), unique(p14$celltype.lv3))

all_in_common<- intersect(common,unique( p21$celltype.lv3))

sort(all_in_common)


alltypes.lv4 <- sort(union(union(p21_celltype.df$celltype.lv4, p14_celltype.df$celltype.lv4), p0_celltype.df$celltype.lv4))
alltypes.lv3 <- sort(union(union(p21_celltype.df$celltype.lv3, p14_celltype.df$celltype.lv3), p0_celltype.df$celltype.lv3))
alltypes.lv2 <- sort(union(union(p21_celltype.df$celltype.lv2, p14_celltype.df$celltype.lv2), p0_celltype.df$celltype.lv2))
alltypes.lv1 <- sort(union(union(p21_celltype.df$celltype.lv1, p14_celltype.df$celltype.lv1), p0_celltype.df$celltype.lv1))
alltypes <- union( union(alltypes.lv1, alltypes.lv2), union(alltypes.lv3, alltypes.lv4) )

alltypes <-sort(alltypes)

alltypes.pal <- hue_pal()(length(alltypes))
names(alltypes.pal) <- alltypes

cols <- alltypes.pal[as.character(unique(p14$"celltype.lv3"))]









celltype.marker <- list(      CM = c("Myocd","Myh7b","Myl7", "Mki67","Ankrd1", "Slit2", "Gm12381", "Tbx2",
                                     "Slit2", "Slit3", "Robo1", "Ddc", "Hmgcs2", "Zbtb16"),
                              Fib = c("Postn","Fbn1", "Col5a1"),
                              Endo = c("Pecam1","Flt1", "Tie1"),
                              "PC/SMC" =c("Myh11","Notch3","Cspg4"),
                              Epic = c("Wt1"),
                              Blood =  c("Mrc1","Fyb", "Cd86"),
                              Neur = c("Kcnh8"),
                              Mela = c("Oca2")
                              )

saveRDS(celltype.marker, "./rObjs/Man_annot/celltype.marker.rds")


ref$celltype.lv2 <- factor(ref$celltype.lv2,  levels= celltype.lv2.levels) %>% droplevels()

levels(ref$celltype.lv2)



levels(ref_CM$celltype.lv2)

pdf("./Manual_Annotation/ref_celltype_plots.pdf", width = 10, onefile = TRUE)
DimPlot(ref, group.by = "celltype.lv1", label = TRUE) 
DimPlot(ref, group.by = "celltype.lv2", label = TRUE) 
DimPlot(ref, group.by = "celltype.lv3", label = TRUE)
DimPlot(ref, group.by = "celltype.lv4", label = TRUE)
DimPlot(ref, group.by = "seurat_clusters", label = TRUE)

DimPlot(ref_CM, group.by = "celltype.lv1", label = TRUE) 
DimPlot(ref_CM, group.by = "celltype.lv2", label = TRUE) 
DimPlot(ref_CM, group.by = "celltype.lv3", label = TRUE)
DimPlot(ref_CM, group.by = "celltype.lv4", label = TRUE)
DimPlot(ref_CM, group.by = "seurat_clusters", label = TRUE)
dev.off()






pdf("./Manual_Annotation/Man_Annot_Markers_ref_summary.pdf", onefile = TRUE, width = 16)
Idents(ref)  <- "celltype.lv3"

DimPlot(ref, group.by = "celltype.lv1", label = TRUE ) + DimPlot(ref,  group.by = "celltype.lv2", label = TRUE) + DimPlot(ref, group.by = "celltype.lv3", label = TRUE)

for (gene in unlist(celltype.marker)){
  plot1 <- FeaturePlot(ref, gene, label=TRUE )
  plot2 <- VlnPlot(ref, gene, group.by = "celltype.lv3")

  
  plot <- plot_grid(plot1, plot2, nrow = 1 )
  print(plot)
  
}
dev.off()




DimPlot(p0, group.by = "celltype.lv1", label = TRUE ) + DimPlot(p0,  group.by = "celltype.lv2", label = TRUE) + DimPlot(p0, group.by = "celltype.lv3", label = TRUE)
plot_grid(DimPlot(p0_CM, group.by = "celltype.lv1", label = TRUE), 
          DimPlot(p0_CM, group.by = "celltype.lv2", label = TRUE),
          DimPlot(p0_CM, group.by = "celltype.lv3", label = TRUE),
          DimPlot(p0_CM, group.by = "celltype.lv4", label = TRUE, label.size = 2.5),
          nrow = 2)
DimPlot(p14, group.by = "celltype.lv1", label = TRUE ) + DimPlot(p14,  group.by = "celltype.lv2", label = TRUE) + DimPlot(p14, group.by = "celltype.lv3", label = TRUE)
plot_grid(DimPlot(p14_CM, group.by = "celltype.lv1", label = TRUE), 
          DimPlot(p14_CM, group.by = "celltype.lv2", label = TRUE),
          DimPlot(p14_CM, group.by = "celltype.lv3", label = TRUE),
          DimPlot(p14_CM, group.by = "celltype.lv4", label = TRUE, label.size = 2.5),
          nrow = 2)
DimPlot(p21, group.by = "celltype.lv1", label = TRUE ) + DimPlot(p21,  group.by = "celltype.lv2", label = TRUE) + DimPlot(p21, group.by = "celltype.lv3", label = TRUE)
DimPlot(p21_CM, group.by = "celltype.lv3", label = TRUE) + DimPlot(p21_CM, group.by = "celltype.lv4", label = TRUE)
plot_grid(DimPlot(p21_CM, group.by = "celltype.lv1", label = TRUE), 
          DimPlot(p21_CM, group.by = "celltype.lv2", label = TRUE),
          DimPlot(p21_CM, group.by = "celltype.lv3", label = TRUE),
          DimPlot(p21_CM, group.by = "celltype.lv4", label = TRUE, label.size = 2.5),
          nrow = 2)

pdf("./Manual_Annotation/celltype.lv4_comp.pdf", width = 15)
plot_grid(DimPlot(p0_CM, group.by = "celltype.lv4", label = TRUE,  label.size = 2.5), 
          DimPlot(p14_CM, group.by = "celltype.lv4", label = TRUE, label.size = 2.5),
          DimPlot(p21_CM, group.by = "celltype.lv4", label = TRUE, label.size = 2.5),
          nrow = 1)

plot_grid(DimPlot(p0, group.by =  "celltype.lv3", label = TRUE,  label.size = 2.5), 
          DimPlot(p14, group.by = "celltype.lv3", label = TRUE, label.size = 2.5),
          DimPlot(p21, group.by = "celltype.lv3", label = TRUE, label.size = 2.5),
          nrow = 1)


dev.off()

pdf("./Manual_Annotation/celltype.lv4_comp.pdf", width = 15, onefile = TRUE)
Idents(p0) <- "celltype.lv4"
plot_grid( FeaturePlot(p0, gene, label = TRUE),
           FeaturePlot(p14, gene, label = TRUE),
           FeaturePlot(p21, gene, label = TRUE),
          nrow = 1)
dev.off()

VlnPlot(p0_CM, )

pdf("./Manual_Annotation/AtriumvsVentricle.pdf", width = 5, onefile = TRUE, height = 10)
plot_grid(
  FeaturePlot(p0_CM, c("Nppa")),
  FeaturePlot(p0_CM, c("Myl4")),
  FeaturePlot(p0_CM, c("Myl7")),
ncol= 1)
 
plot_grid(
  FeaturePlot(p0_CM, c("Myh7")),
  FeaturePlot(p0_CM, c("Myl2")),
  FeaturePlot(p0_CM, c("Fhl2")),
  ncol= 1)

dev.off()


