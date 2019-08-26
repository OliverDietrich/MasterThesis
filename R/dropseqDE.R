######################################################################
####### Gene set enrichment for 10x datasets in Seurat object ########
#--------------------------------------------------------------------#
############ R script by Oliver Dietrich on 2 May 2019 ###############
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
library(pheatmap)
library(destiny)
set.seed(154)

############################ Load dataset ############################
print(paste("Loading dataset", DSname))

# check to which project the dataset belongs, allows for usage of absolute path
HOME <- paste0("/home/", Sys.getenv("USER"), "/Data/")
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
imageStats <- read.table("Input/imageStats.csv", sep = ",", header = TRUE)
# Load dimensionalReduction table
labels <- read.table("Input/labels", sep = "\t", header = TRUE)
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

resolution <- colnames(labels)[1] ; resName <- str_replace_all(resolution, "\\.", "_")

############################ Gene Set Enrichment Analysis ############################
print("Analysing differential gene expression")

path <- paste0(resName, "/")

DS@meta.data$label <- as.character(DS@meta.data$label)
DS <- SetAllIdent(object = DS, id = "label")

markers <- FindAllMarkers(DS)
markers$feature <- features$feature[match(markers$gene, features$ensID)]
write_csv(markers, paste0(resName, "/", DSname, "_markergenes.csv"))

a <- t(data.frame(label = labels$label[order(labels$label)], color = labels$color[order(labels$label)]))
colnames(a) <- a["label",] ; a <- a["color",]

myColor <- list(
  label = a,
  cluster = a
)

for(i in c(5,10,15,30,50,100)) {
  ifelse(i < 30, labelRows <- TRUE, labelRows <- FALSE)
  markers2plot <- markers %>% group_by(cluster) %>% top_n(i, avg_logFC) %>% as.data.frame()
  data2plot <- DS@scale.data[markers2plot$gene,order(DS@meta.data$label)]
  rownames(data2plot) <- features$feature[match(rownames(data2plot), features$ensID)]
  data2plot <- unique(data2plot)
  
  rowClusters <- markers2plot[,c("feature", "cluster")]
  rowClusters <- rowClusters[match(unique(rowClusters$feature), rowClusters[,"feature"]),] %>% remove_rownames()
  rowClusters <- rowClusters %>% column_to_rownames(var = "feature")
  columnClusters <- DS@meta.data[,"label", drop = FALSE]
  
  plot <- pheatmap(data2plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, 
                   color = viridis(n = 100, option = "B"), labels_col = c(levels(DS@meta.data$label)), 
                   annotation_row = rowClusters, annotation_col = columnClusters, annotation_colors = myColor, 
                   show_rownames = labelRows)
  
  ggsave(paste0("DEheatmap_top", i, ".png"), plot = plot, path = resName, dpi = 300)
} # end of heatmap for loop

sessionInfo()
quit(status = 0)
