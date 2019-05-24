######################################################################
###### Print Expression Plots for 10x datasets in Seurat object ######
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
clusterOptions <- read.table("Input/clusterOptions", sep = "\t", header = TRUE)
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

############################ Expression Plots ############################
print(paste("Visualizing expression data in reduced dimension"))

markers <- as.character(clusterOptions$value[match("markers", clusterOptions$parameter)])
if( file.exists(paste0(HOME, projectName, "/markers/", markers)) == TRUE) {
  markergenes <- read.table(paste0(HOME, projectName, "/markers/", markers), sep = ",", header = TRUE)
  markers <- as.character(features$ensID[match(as.character(markergenes$gene), features$feature)])
  markers <- markers[markers %in% rownames(DS@data)]
} else { if(markers %in% "all") { markers <- rownames(DS@data)
  } else { markers <- DS@var.genes }}

# tSNE
if(DS@dr$tsne@key %in% "tSNE_") {
  print("tSNE")
  path <- "DR/tSNE/expression"
  dir.create(path, recursive = TRUE)
  
  stats <- function(parameter) {  as.character(imageStats$tSNE[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$tSNE[match(parameter, imageStats$parameters)]))  }
  if(imageStats$tSNE[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in markers) {
    j <- features$feature[match(i, features$ensID)]
    
    plot <- ggplot(as.data.frame(DS@dr$tsne@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$tsne@cell.embeddings[,"tSNE_1"], y=DS@dr$tsne@cell.embeddings[,"tSNE_2"], alpha=numstats("alpha"), 
        color=DS@data[i,]
      ), size = numstats("pointSize")) +
      labs(
        x="tSNE_1", y="tSNE_2", title=j
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
      guides(alpha = FALSE) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2")
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_",i,"_",j,".",stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot

# UMAP
if(DS@dr$umap@key %in% "UMAP") {
  print("UMAP")
  path <- "DR/UMAP/expression"
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$UMAP[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$UMAP[match(parameter, imageStats$parameters)]))  }
  if(imageStats$UMAP[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in markers) {
    j <- features$feature[match(i, features$ensID)]
    
    plot <- ggplot(as.data.frame(DS@dr$umap@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$umap@cell.embeddings[,"UMAP1"], y=DS@dr$umap@cell.embeddings[,"UMAP2"], alpha=numstats("alpha"), 
        color=DS@data[i,]
      ), size = numstats("pointSize")) +
      labs(
        x="UMAP1", y="UMAP2", title=j
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
      guides(alpha = FALSE) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2")
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_",i,"_",j,".",stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot

# Diffusion Map
if(DS@dr$dm@key %in% "DM") {
  print("Diffusion Map")
  path <- "DR/DM/expression"
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$DM[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$DM[match(parameter, imageStats$parameters)]))  }
  if(imageStats$DM[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
  for(i in markers) {
    j <- features$feature[match(i, features$ensID)]
    
    plot <- ggplot(as.data.frame(DS@dr$dm@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$dm@cell.embeddings[,"DM1"], y=DS@dr$dm@cell.embeddings[,"DM2"], alpha=numstats("alpha"), 
        color=DS@data[i,]
      ), size = numstats("pointSize")) +
      labs(
        x="DM1", y="DM2", title=j
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
      guides(alpha = FALSE) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2")
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_",i,"_",j,".",stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot 


sessionInfo()
quit(status = 0)
