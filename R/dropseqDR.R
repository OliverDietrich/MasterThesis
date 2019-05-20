######################################################################
###### Dimensional Reduction for 10x datasets in Seurat object #######
#--------------------------------------------------------------------#
########### R script by Oliver Dietrich on 26 April 2019 #############
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
dimensionalReduction <- read.table("Input/dimensionalReduction", sep = "\t", header = TRUE)
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

############################ Dimensional Reduction ############################
print(paste("Reducing dimensions of", DSname))

stats <- function(parameter) {as.character(dimensionalReduction$value[match(parameter, dimensionalReduction$parameter)])}
numstats <- function(parameter) {as.numeric(as.character(dimensionalReduction$value[match(parameter, dimensionalReduction$parameter)]))}

if(stats("tSNE") == TRUE) {DS <- RunTSNE(DS, reduction.use = "pca", dims.use = as.numeric(str_split(stats("PCs"), ",")[[1]])) }
if(stats("3DtSNE") == TRUE) {DS <- RunTSNE(DS, reduction.use = "pca", dims.use = as.numeric(str_split(stats("PCs"), ",")[[1]]),
                                           dim.embed = 3, reduction.key = "tSNE_", reduction.name = "tsne3")}
if(stats("UMAP") == TRUE) {DS <- RunUMAP(DS, reduction.use = "pca", dims.use = as.numeric(str_split(stats("PCs"), ",")[[1]])) }
if(stats("DM") == TRUE) {DS <- RunDiffusion(DS, reduction.use = "pca", dims.use = as.numeric(str_split(stats("PCs"), ",")[[1]])) }
if(stats("SPRING") == TRUE) { print("Do the SPRING plot") }

############################ Visualization ############################
print(paste("Visualizing data in reduced dimension"))

# tSNE
if(DS@dr$tsne@key %in% "tSNE_") {
path <- "DR/tSNE" ; print("tSNE")
dir.create(path, recursive = TRUE)

stats <- function(parameter) {  as.character(imageStats$tSNE[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$tSNE[match(parameter, imageStats$parameters)]))  }
if(imageStats$tSNE[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in colnames(DS@meta.data)) {
  
  if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) { 
    guide <- guides(alpha = FALSE) ; color <- scale_color_viridis()
    } else { 
    guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) 
    color <- geom_blank()
    }
  
  plot <- ggplot(as.data.frame(DS@dr$tsne@cell.embeddings)) + 
    geom_point(aes(
      x=DS@dr$tsne@cell.embeddings[,"tSNE_1"], y=DS@dr$tsne@cell.embeddings[,"tSNE_2"], alpha=numstats("alpha"), 
      color=DS@meta.data[,i]
    ), size = numstats("pointSize")) +
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
path <- "DR/UMAP" ; print("UMAP")
dir.create(path, recursive = FALSE)

stats <- function(parameter) {  as.character(imageStats$UMAP[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$UMAP[match(parameter, imageStats$parameters)]))  }
if(imageStats$UMAP[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in colnames(DS@meta.data)) {
  
  if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) {
    guide <- guides(alpha = FALSE) ; color <- scale_color_viridis()
  } else { 
    guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) 
    color <- geom_blank()
  }
  
  plot <- ggplot(as.data.frame(DS@dr$umap@cell.embeddings)) + 
    geom_point(aes(
      x=DS@dr$umap@cell.embeddings[,"UMAP1"], y=DS@dr$umap@cell.embeddings[,"UMAP2"], alpha=numstats("alpha"), 
      color=DS@meta.data[,i]
    ), size = numstats("pointSize")) +
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
    guide +
    color
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_tSNE_", i, ".", stats("fileType")),
         dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot

# Diffusion Map
if(DS@dr$dm@key %in% "DM") {
path <- "DR/DM" ; print("DM")
dir.create(path, recursive = FALSE)

stats <- function(parameter) {  as.character(imageStats$DM[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$DM[match(parameter, imageStats$parameters)]))  }
if(imageStats$DM[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in colnames(DS@meta.data)) {
  
  if(class(DS@meta.data[,i]) %in% c("numeric", "integer")) { 
    guide <- guides(alpha = FALSE) ; color <- scale_color_viridis()
  } else { 
    guide <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) 
    color <- geom_blank()
  }

  plot <- ggplot(as.data.frame(DS@dr$dm@cell.embeddings)) + 
    geom_point(aes(
      x=DS@dr$dm@cell.embeddings[,"DM1"], y=DS@dr$dm@cell.embeddings[,"DM2"], alpha=numstats("alpha"), 
      color=DS@meta.data[,i]
    ), size = numstats("pointSize")) +
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
    guide +
    color
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_tSNE_", i, ".", stats("fileType")),
         dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } } # end of plot 

# SPRING plot
# will hopefully be implemented soon ...

############################ Saving data ############################
print(paste("Saving dataset", DSname))

if(file.exists("clusterOptions") == FALSE & file.exists("Input/clusterOptions") == FALSE) {
clusterOptions <- data.frame(
  parameter = c("algorithm","min","max","by","markers"),
  value = c("SLM",0.1,3,0.1,"markergenes.csv")
)
write.table(clusterOptions, "clusterOptions", sep = "\t", row.names = FALSE)
}

# save dataset
saveRDS(DS, paste0("RDS/", DSname, ".Rds"))

sessionInfo()
quit(status = 0)
