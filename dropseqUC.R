######################################################################
##### Unsupervised Clustering for 10x datasets in Seurat object ######
#--------------------------------------------------------------------#
########### R script by Oliver Dietrich on 29 April 2019 #############
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
imageStats <- read.table("Input/imageStats.csv", sep = ",", header = TRUE)
# Load dimensionalReduction table
clusterOptions <- read.table("Input/clusterOptions", sep = "\t", header = TRUE)
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

############################ Clustering ############################
print(paste("Loading dataset", DSname))

DS <- BuildSNN(DS, reduction.type = "pca", dims.use = DS@calc.params$RunUMAP$dims.use, prune.SNN = 0, 
               plot.SNN = FALSE, force.recalc = TRUE, save.SNN = TRUE)

stats <- function(parameter) {  as.character(clusterOptions$value[match(parameter, clusterOptions$parameter)])  }
numstats <- function(parameter) {  as.numeric(as.character(clusterOptions$value[match(parameter, clusterOptions$parameter)]))  }

if(stats("algorithm") %in% "Louvain") {algorithm <- 1} ; if(stats("algorithm") %in% "LMR") {algorithm <- 2}
if(stats("algorithm") %in% "SLM") {algorithm <- 3}

max <- numstats("max") ; min <- numstats("min") ; by <- numstats("by")

for(i in seq(min,max,by=by)) {
  print(paste("Resolution:",i))
DS <- FindClusters(DS, resolution = i, print.output = FALSE, reuse.SNN = TRUE, algorithm = algorithm)
}

############################ Visualization ############################
print(paste("Visualizing data in reduced dimension"))

# tSNE
if(DS@dr$tsne@key %in% "tSNE_") {
  path <- "DR/tSNE/clustering"
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$tSNE[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$tSNE[match(parameter, imageStats$parameters)]))  }
  if(imageStats$tSNE[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in colnames(DS@meta.data[,str_detect(colnames(DS@meta.data), "res")])) {
    
    if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) { 
      guide <- guides(alpha = FALSE) ; color <- scale_color_gradient(low = "grey95", high = "dodgerblue2")
    } else { 
      guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) 
      color <- geom_blank()
    }
    
    plot <- ggplot(as.data.frame(DS@dr$tsne@cell.embeddings), aes(x=DS@dr$tsne@cell.embeddings[,"tSNE_1"], 
                                                                  y=DS@dr$tsne@cell.embeddings[,"tSNE_2"])) + 
      geom_point(aes(
        alpha=numstats("alpha"), 
        color=DS@meta.data[,i]
      ), size = numstats("pointSize")) +
      geom_text(aes(label = DS@meta.data[,i]), size = numstats("pointSize")/2) +
      labs(
        x="tSNE_1", y="tSNE_2", title=i
      ) +
      theme(
        plot.title = element_text(size = numstats("titleSize")),
        axis.title = element_text(size = numstats("axisTitle")),
        axis.text.x = element_text(size = numstats("axisXsize"), angle = numstats("axisXangle"), 
                                   hjust = numstats("axisXhjust"), vjust = numstats("axisXvjust")),
        axis.text.y = element_text(size = numstats("axisYsize"), angle = numstats("axisYangle"), 
                                   hjust = numstats("axisYhjust"), vjust = numstats("axisYvjust")),
        axis.line.x = element_line(size = numstats("axisLineSize")),
        axis.line.y = element_line(size = numstats("axisLineSize")),
        axis.ticks = element_line(size = numstats("axisTicks")),
        legend.position = legend_position,
        legend.title = element_text(size = numstats("legendTitle")),
        legend.text = element_text(size = numstats("legendText")),
        legend.key.size = unit(numstats("legendKeySize"), "pt"),
        legend.background = element_rect(fill = stats("legendBgrdCol"))
      ) +
      guide +
      color
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_", i, ".", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot


# UMAP
if(DS@dr$umap@key %in% "UMAP") {
  path <- "DR/UMAP/clustering"
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$UMAP[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$UMAP[match(parameter, imageStats$parameters)]))  }
  if(imageStats$UMAP[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in colnames(DS@meta.data[,str_detect(colnames(DS@meta.data), "res")])) {
    
    if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) {
      guide <- guides(alpha = FALSE)
    } else { guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) }
    
    plot <- ggplot(as.data.frame(DS@dr$umap@cell.embeddings), aes(x=DS@dr$umap@cell.embeddings[,"UMAP1"], 
                                                                  y=DS@dr$umap@cell.embeddings[,"UMAP2"])) + 
      geom_point(aes(
        alpha=numstats("alpha"), 
        color=DS@meta.data[,i]
      ), size = numstats("pointSize")) +
      geom_text(aes(label = DS@meta.data[,i]), size = numstats("pointSize")/2) +
      labs(
        x="UMAP1", y="UMAP2", title=i
      ) +
      theme(
        plot.title = element_text(size = numstats("titleSize")),
        axis.title = element_text(size = numstats("axisTitle")),
        axis.text.x = element_text(size = numstats("axisXsize"), angle = numstats("axisXangle"), 
                                   hjust = numstats("axisXhjust"), vjust = numstats("axisXvjust")),
        axis.text.y = element_text(size = numstats("axisYsize"), angle = numstats("axisYangle"), 
                                   hjust = numstats("axisYhjust"), vjust = numstats("axisYvjust")),
        axis.line.x = element_line(size = numstats("axisLineSize")),
        axis.line.y = element_line(size = numstats("axisLineSize")),
        axis.ticks = element_line(size = numstats("axisTicks")),
        legend.position = legend_position,
        legend.title = element_text(size = numstats("legendTitle")),
        legend.text = element_text(size = numstats("legendText")),
        legend.key.size = unit(numstats("legendKeySize"), "pt"),
        legend.background = element_rect(fill = stats("legendBgrdCol"))
      ) +
      guide
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_", i, ".", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot

# Diffusion Map
if(DS@dr$dm@key %in% "DM") {
  path <- "DR/DM/clustering"
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$DM[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$DM[match(parameter, imageStats$parameters)]))  }
  if(imageStats$DM[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in colnames(DS@meta.data[,str_detect(colnames(DS@meta.data), "res")])) {
    
    if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) {
      guide <- guides(alpha = FALSE)
    } else { guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) }
    
    plot <- ggplot(as.data.frame(DS@dr$dm@cell.embeddings), aes(x=DS@dr$dm@cell.embeddings[,"DM1"], 
                                                                y=DS@dr$dm@cell.embeddings[,"DM2"])) + 
      geom_point(aes(
        alpha=numstats("alpha"), 
        color=DS@meta.data[,i]
      ), size = numstats("pointSize")) +
      geom_text(aes(label = DS@meta.data[,i]), size = numstats("pointSize")/2) +
      labs(
        x="DM1", y="DM2", title=i
      ) +
      theme(
        plot.title = element_text(size = numstats("titleSize")),
        axis.title = element_text(size = numstats("axisTitle")),
        axis.text.x = element_text(size = numstats("axisXsize"), angle = numstats("axisXangle"), 
                                   hjust = numstats("axisXhjust"), vjust = numstats("axisXvjust")),
        axis.text.y = element_text(size = numstats("axisYsize"), angle = numstats("axisYangle"), 
                                   hjust = numstats("axisYhjust"), vjust = numstats("axisYvjust")),
        axis.line.x = element_line(size = numstats("axisLineSize")),
        axis.line.y = element_line(size = numstats("axisLineSize")),
        axis.ticks = element_line(size = numstats("axisTicks")),
        legend.position = legend_position,
        legend.title = element_text(size = numstats("legendTitle")),
        legend.text = element_text(size = numstats("legendText")),
        legend.key.size = unit(numstats("legendKeySize"), "pt"),
        legend.background = element_rect(fill = stats("legendBgrdCol"))
      ) +
      guide
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_", i, ".", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot 

############################ Saving data ############################
print(paste("Saving dataset", DSname))

if(file.exists("labels") == FALSE & file.exists("Input/labels") == FALSE) {
  labels <- data.frame(
    "res.0.1" = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
    label = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U"),
    color = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U")
  )
  write.table(labels, "labels", sep = "\t", row.names = FALSE)
}

# save dataset
saveRDS(DS, paste0("RDS/", DSname, ".Rds"))

sessionInfo()
quit(status = 0)
