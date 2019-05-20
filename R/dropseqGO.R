######################################################################
##### Gene ontology annotation for 10x datasets in Seurat object #####
#--------------------------------------------------------------------#
############ R script by Oliver Dietrich on 3 May 2019 ###############
######################################################################

# read input arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
DSname <- args[1]
arg2 <- args[2]
arg3 <- args[3]

# load dependencies
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library("clusterProfiler")
library("org.Mm.eg.db")
set.seed(154)

############################ Load dataset ############################
print(paste("Loading dataset", DSname))

# check to which project the dataset belongs, allows for usage of absolute path
HOME <- paste0("/home/", Sys.getenv("USERNAME"), "/Data/")
if(file.exists(paste0(HOME, "projects")) == TRUE) {
  projects <- read.table(file = paste0(HOME, "projects"), sep = "\t", header = TRUE)
  projectName <- as.character(projects$project[str_detect(projects$datasets, DSname)])
} else {
  print("projects file not available, aborting.")
  quit(status = 10)
}
# Load feature table
inputDir <- paste0(HOME, projectName, "/analysis/Input/")
features <- read.table(paste0(inputDir, "features.tsv"), sep = "\t", header = T)
# Load imageStats table
imageStats <- read.table("../Input/imageStats.csv", sep = ",", header = TRUE)
# Load marker table
markers <- read.table(paste0(DSname, "_markergenes.csv"), sep = ",", header = TRUE)

############################ Gene Ontology Annotation ############################
print(paste("Annotating Gene Ontology for", DSname))
dir.create("GO")



if(file.exists("GO/categories") & file.exists("GO/G6_GO_markers_pos.csv")) {
  
  categories <- readLines("GO/categories")
  markers <- read.csv("GO/G6_GO_markers_pos.csv", header = TRUE)
  
  
  markerTable <- markers[markers$Description %in% categories,]
  order <- markerTable$Description %>% as.character()
  
  a <- str_split(as.character(markerTable$GeneRatio), "/")
  for(i in seq(1, length(a), 1) ) {
    b <- a[i][[1]] %>% as.numeric()
    b <- b[1] / b[2]
    a[i] <- b
  }
  c <- unlist(a)
  markerTable$GeneRatio <- c
  
  plot <- ggplot(markerTable) + 
    geom_count(aes(
      x = Cluster, y = factor(Description, levels = categories), color = p.adjust, size = GeneRatio
    )) +
    scale_color_gradient(low = "blue2", high = "red") +
    labs(x = "", y = "", title = "")
  
  ggsave("GO_pos_assembled.eps", plot = plot, path = "GO/")
  
} else {
  
  markers <- markers[markers$avg_logFC > 0,]
  
  entrez <- bitr(markers$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  cluster <- markers$cluster[match(entrez$ENSEMBL, markers$gene)]
  
  entrez <- entrez %>% mutate("cluster" = cluster)
  
  clusters <- unique(entrez$cluster)
  
  EnsembleList <- list(clusters)
  
  for (i in clusters) {
    EnsembleList[[i]] <- entrez$ENTREZID[entrez$cluster == i]
  }
  
  ck1 <- compareCluster(geneCluster = EnsembleList, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
  write_csv(ck1@compareClusterResult, path = paste0("GO/", DSname, "_GO_markers_pos.csv"))
  
  p1 <- dotplot(ck1, showCategory = 3, font.size = 13, by = "geneRatio")
  ggsave(filename = paste0(DSname, "_GO_BP_pos.eps"), plot = p1, device = "eps", path = paste0("GO"), width = 10, height = 8)
  
  ### Only on short list of markergenes
  
  for(i in c(10,50)) {
  markers2plot <- markers %>% group_by(cluster) %>% top_n(i, avg_logFC) %>% as.data.frame()
  
  entrez <- bitr(markers2plot$gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  cluster <- markers2plot$cluster[match(entrez$ENSEMBL, markers2plot$gene)]
  
  entrez <- entrez %>% mutate("cluster" = cluster)
  
  clusters <- unique(entrez$cluster)
  
  EnsembleList <- list(clusters)
  
  for (i in clusters) {
    EnsembleList[[i]] <- entrez$ENTREZID[entrez$cluster == i]
  }
  
  ck1 <- compareCluster(geneCluster = EnsembleList, fun = "enrichGO", OrgDb = "org.Hs.eg.db", ont = "BP")
  write_csv(ck1@compareClusterResult, path = paste0("GO/", DSname, "_", i, "_GO_markers_pos.csv"))
  
  p1 <- dotplot(ck1, showCategory = 3, font.size = 13, by = "geneRatio")
  ggsave(filename = paste0(DSname, "_", i, "_GO_BP_pos.eps"), plot = p1, device = "eps", path = paste0("GO"), width = 10, height = 8)
  } # end of for loop
} # end of if statement

sessionInfo()
quit(status = 0)
