######################################################################
######### Quality Control for 10x datasets in Seurat object ##########
#--------------------------------------------------------------------#
########### R script by Oliver Dietrich on 23 April 2019 #############
######################################################################

# read input arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
DSname <- args[1]

# load dependencies
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(ggpubr)
set.seed(154)

# make directories
dir.create("QC")
dir.create("QC/unfiltered")

# set trigger
status <- "unfiltered"

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
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

############################ Add metadata ############################
print(paste("Adding metadata to", DSname))

# automatically detect if mouse or human (other options not possible yet)
if (stringr::str_detect(features$ensID[1], "ENSMUSG") == TRUE) {organism <- "mouse"}
if (stringr::str_detect(features$ensID[1], "ENSG") == TRUE) {organism <- "human"}

# detect mitochondrial genes
if (organism == "human") {mitoGenes <- features$feature[grep(pattern = "^MT-", features$feature)]} 
if (organism == "mouse") {mitoGenes <- features$feature[grep(pattern = "^mt-", features$feature)]}
# convert gene names to ensemble IDs
mitoIDs <- features$ensID[match(mitoGenes, features$feature)] 
# calculate proportion of mitoGenes per cell
pMito <- Matrix::colSums(DS@raw.data[match(mitoIDs, rownames(DS@data)),]) / Matrix::colSums(DS@raw.data) * 100 
# add to DS@metadata
DS <- AddMetaData(DS, metadata = pMito, col.name = "pMito") 

# Add cell cycle scores if data is accessible
humanCCgenes <- paste0(inputDir, "regev_lab_cell_cycle_genes_GRCh38.txt")
mouseCCgenes <- paste0(inputDir, "regev_lab_cell_cycle_genes_mm10.txt")
if (file.exists(humanCCgenes) == TRUE | file.exists(mouseCCgenes) == TRUE ) {
  print(paste("Calculating cell cycle scores for", DSname))
  
  # load cell cycle genes (depending on organism) + write to file for subset analysis
  if (organism == "human") { 
    cc.genes <- readLines(humanCCgenes) } 
  if (organism == "mouse") { 
    cc.genes <- readLines(mouseCCgenes) }
  # split group of genes in two categories
  s.genes <- features$ensID[match(cc.genes[1:43], features$feature)]
  g2m.genes <- features$ensID[match(cc.genes[44:97], features$feature)]
  # add cell cycle scores and category to DS@meta.data
  DS <- CellCycleScoring(DS, s.genes = s.genes, g2m.genes = g2m.genes)
  DS@meta.data$CC.Difference <- DS@meta.data$S.Score - DS@meta.data$G2M.Score
} # end of cell cycle loop

############################ Filtering ############################
if(file.exists("Input/thresholds") == TRUE) {
  print(paste("Filtering genes and cells of dataset", DSname))
  dir.create("QC/filtered")
  
  # Load thresholds table
  thresholds <- read.table("Input/thresholds", sep = "\t", header = TRUE) %>% 
    column_to_rownames(var = "type")
  
  T0 <- dim(DS@raw.data)
  # Filter out zero-count genes
  DS@raw.data <- DS@raw.data[rowSums(DS@raw.data) > 0,]
  DS@data <- DS@data[rowSums(DS@data) > 0,]
  # Filter out cells based on provided thresholds
  a <- DS@meta.data[which(DS@meta.data$nGene > thresholds["low", "nGene"]),]
  b <- a[which(a$nGene < thresholds["high", "nGene"]),]
  c <- b[which(b$nUMI > thresholds["low", "nUMI"]),]
  d <- c[which(c$nUMI < thresholds["high", "nUMI"]),]
  e <- d[which(d$pMito > thresholds["low", "pMito"]),]
  f <- e[which(e$pMito < thresholds["high", "pMito"]),]
  
  DS <- SubsetData(DS, cells.use = rownames(f), subset.raw = TRUE)
  
  # Alternative
  #DS@raw.data <- DS@raw.data[,rownames(f)]
  #DS@data <- DS@data[,rownames(f)]
  #DS@meta.data <- DS@meta.data[rownames(f),]
  ##
  T1 <- dim(DS@raw.data)
  
  if(file.exists("RDS/dimensions") == TRUE) {
    dimensions <- read.table("RDS/dimensions")
    dimensions$after <- T1 ; dimensions$proportion <- round(dimensions$after/dimensions$before, digits = 3)
    write.table(dimensions, "RDS/dimensions")
  } else {
  dimensions <- data.frame(before = T0, after = T1) %>% mutate(proportion = round(after/before, digits = 3))
  rownames(dimensions) <- c("genes", "cells")
  write.table(dimensions, "RDS/dimensions")
  }
  

  
  status <- "filtered"
}

############################ Visualize data ############################
print(paste("Visualizing metadata for", DSname))

imageStats <- read.table("Input/imageStats.csv", sep = ",", header = TRUE)
path <- paste0("QC/", status)

stats <- function(plot,parameter) {as.character(imageStats[[plot]][match(parameter, imageStats$parameters)]) }
numstats <- function(plot,parameter) {as.numeric(as.character(imageStats[[plot]][match(parameter, imageStats$parameters)])) }

# Boxplots
class <- "boxplot"

if(stats(class,"legendPositionY") %in% "none") { legend_position <- stats(class, "legendPositionX") } else {
  legend_position <- c(numstats(class, "legendPositionX"), numstats(class, "legendPositionY")) }

for(i in c("nGene", "nUMI", "pMito")) {
  plot <- ggplot(DS@meta.data, aes(x=DS@meta.data[,"dataset"], y=DS@meta.data[,i])) + 
    #
    geom_violin(aes(color=DS@meta.data[,"dataset"])) +
    #
    geom_point(aes(alpha=numstats(class,"alpha")), position = stats(class,"position"), 
               size = numstats(class,"pointSize"), color=stats(class,"color")) +
    #
    geom_boxplot(aes(fill=DS@meta.data[,stats(class,"fill")])) +
    #
    labs(
      x="dataset", y=i, title=""
    ) +
    theme(
      plot.title = element_text(size = numstats(class,"titleSize")),
      axis.title = element_text(size = numstats(class,"axisTitle")),
      axis.text.x = element_text(size = numstats(class,"axisXsize"), angle = numstats(class,"axisXangle"), 
                                 hjust = numstats(class,"axisXhjust"), vjust = numstats(class,"axisXvjust")),
      axis.text.y = element_text(size = numstats(class,"axisYsize"), angle = numstats(class,"axisYangle"), 
                                 hjust = numstats(class,"axisYhjust"), vjust = numstats(class,"axisYvjust")),
      axis.line = element_line(size = numstats(class,"axisLineSize"), linetype = "solid"),
      legend.position = legend_position,
      legend.title = element_text(size = numstats(class,"legendTitle")),
      legend.text = element_text(size = numstats(class,"legendText")),
      legend.key.size = unit(numstats(class,"legendKeySize"), "pt"),
      legend.background = element_rect(fill = stats(class,"legendBgrdCol"))
    ) +
    guides(
      alpha = FALSE, color = FALSE
    )
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_", class, "_",i, ".", stats(class,"fileType")),
         dpi=numstats(class,"dpi"), width=numstats(class,"width"), height=numstats(class,"height")
  )   } # end of ggplot for loop

# Smooth density plots
class <- "smoothDensity"

if(stats(class,"legendPositionY") %in% "none") { legend_position <- stats(class, "legendPositionX") } else {
  legend_position <- c(numstats(class, "legendPositionX"), numstats(class, "legendPositionY")) } # end of legend if

for(i in c("nGene", "nUMI", "pMito")) {
  plot <- ggplot(DS@meta.data) + 
    #
    geom_density(aes( x=DS@meta.data[,i], fill=DS@meta.data[,stats(class, "fill")], 
                      alpha=numstats(class, "alpha"), color=DS@meta.data[,stats(class, "color")] )) +
    #
    labs( x=i, title="" ) +
    #
    theme(
      plot.title = element_text(size = numstats(class, "titleSize")),
      axis.title = element_text(size = numstats(class, "axisTitle")),
      axis.text.x = element_text(size = numstats(class, "axisXsize"), angle = numstats(class, "axisXangle"), 
                               hjust = numstats(class, "axisXhjust"), vjust = numstats(class, "axisXvjust")),
      axis.text.y = element_text(size = numstats(class, "axisYsize"), angle = numstats(class, "axisYangle"), 
                                 hjust = numstats(class, "axisYhjust"), vjust = numstats(class, "axisYvjust")),
      axis.line = element_line(size = numstats(class, "axisLineSize"), linetype = "solid"),
      legend.position = legend_position,
      legend.title = element_text(size = numstats(class, "legendTitle")),
      legend.text = element_text(size = numstats(class, "legendText")),
      legend.key.size = unit(numstats(class, "legendKeySize"), "pt"),
      legend.background = element_rect(fill = stats(class, "legendBgrdCol"))
    ) +
    guides(
      alpha = FALSE, color = FALSE
    )
  ggsave(plot = plot,
    path=path, filename=paste0(DSname,"_", class, "_",i, ".", stats(class, "fileType")),
    dpi=numstats(class, "dpi"), width=numstats(class, "width"), height=numstats(class, "height")
  )   } # end of ggplot for loop

# Scatterplots
class <- "scatter"

if(stats(class,"legendPositionY") %in% "none") { legend_position <- stats(class, "legendPositionX") } else {
  legend_position <- c(numstats(class, "legendPositionX"), numstats(class, "legendPositionY")) } # end of legend if

for(i in c("nUMI", "pMito")) {
  plot <- ggplot(DS@meta.data, aes(x=DS@meta.data[,"nGene"], y=DS@meta.data[,i])) + 
    #
    geom_point(aes( alpha=numstats(class, "alpha"), color=DS@meta.data[,stats(class, "color")] ), 
               size = numstats(class, "pointSize")) +
    #
    geom_smooth(method = lm) +
    stat_cor(method = "pearson", label.x = min(DS@meta.data$nGene), label.y = max(DS@meta.data[,i]), 
             size = numstats(class, "guideSize")) +
    #
    labs( x="nGene", y=i, title="" ) +
    #
    theme(
      plot.title = element_text(size = numstats(class, "titleSize")),
      axis.title = element_text(size = numstats(class, "axisTitle")),
      axis.text.x = element_text(size = numstats(class, "axisXsize"), angle = numstats(class, "axisXangle"), 
                                 hjust = numstats(class, "axisXhjust"), vjust = numstats(class, "axisXvjust")),
      axis.text.y = element_text(size = numstats(class, "axisYsize"), angle = numstats(class, "axisYangle"), 
                                 hjust = numstats(class, "axisYhjust"), vjust = numstats(class, "axisYvjust")),
      axis.line = element_line(size = numstats(class, "axisLineSize"), linetype = "solid"),
      legend.position = legend_position,
      legend.title = element_text(size = numstats(class, "legendTitle")),
      legend.text = element_text(size = numstats(class, "legendText")),
      legend.key.size = unit(numstats(class, "legendKeySize"), "pt"),
      legend.background = element_rect(fill = stats(class, "legendBgrdCol"))
    ) +
    #
    guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats(class, "guideSize"))))
  
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_", class, "_",i, ".", stats(class, "fileType")),
         dpi=numstats(class, "dpi"), width=numstats(class, "width"), height=numstats(class, "height")
  )   } # end of ggplot for loop

############################ Saving data ############################
print(paste("Saving dataset", DSname))

# print table to set filtering thresholds
if(file.exists("Input/thresholds") == FALSE) {
thresholds <- data.frame("type" = c("high", "low"), "nGene" = c(Inf, 0), "nUMI" = c(Inf, 0), "pMito" = c(Inf,-Inf))
write.table(thresholds, "thresholds", sep = "\t", row.names = FALSE) # print thresholds table 
} else {
  if(file.exists("Input/featureSelection") == FALSE) {
  featureSelection <- data.frame(
    "parameter" = c("normalization", "vars_2_regress", "featureSelection", "Xl_Xh_Yl_Yh", "PCs"),
    "value" = c("logNorm", "nUMI_pMito", "variance", "0.1_8_1_Inf", "50") )
  # write table with input options for clustering
  write.table(featureSelection, "featureSelection", sep = "\t", row.names = FALSE) }
}
# save dataset
saveRDS(DS, paste0("RDS/", DSname, ".Rds"))

sessionInfo()
quit(status = 0)
