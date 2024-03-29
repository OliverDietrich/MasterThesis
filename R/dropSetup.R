######################################################################
######### Setup Seurat object from Cell Ranger output files ##########
#--------------------------------------------------------------------#
########### R script by Oliver Dietrich on 17 April 2019 #############
######################################################################

# read input arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
DSname <- args[1]
arg2 <- args[2]
arg3 <- args[3]

# load dependencies
library(tidyverse)
library(Seurat)
library(hdf5r)
set.seed(154)

# check if the input is correct, otherwise abort script
if( arg2 %in% c("raw", "filtered")) { type <- arg2 ; print(paste("Type:", type)) } else {
  print("Error. Please specify data type ('raw' or 'filtered'). Usage: DS type format") ; quit(status = 10) }

if( arg3 %in% c("h5", "mtx") ) { format <- arg3 ; print(paste("Format:", format)) } else {
  print("Error. Please specify data format ('h5' or 'mtx'). Usage: DS type format") ; quit(status = 10) }

# check to which project the dataset belongs
HOME <- paste0("/home/", Sys.getenv("USER"), "/Data/")
if(file.exists(paste0(HOME, "projects")) == TRUE) {
  projects <- read.table(file = paste0("/home/", Sys.getenv("USER"), "/Data/", "projects"), sep = "\t", header = TRUE)
  projectName <- as.character(projects$project[str_detect(projects$datasets, DSname)])
} else { print("Error. Projects file not available or path corrupt.") ; quit(status = 12) }

dir.create(paste0(HOME, projectName, "/analysis/Input"), recursive = TRUE)

# load the data and convert to Seurat object, depends on the file format (mtx or hdf5)
if( format %in% "h5") {
  print("Error. The Hdf5 format cannot be used yet. Please refer to the matrix format instead.") ; quit(status = 13)
  # IN CONSTRUCTION ...
  # pathIN <- paste0(pathIN, projectName, "/datasets/", DSname, "/outs/", type, "_feature_bc_matrix.h5")
  # DS_data <- Read10X_h5(pathIN)
  # DS <- CreateSeuratObject(DS_data)
} else {
  # set path to data
  pathIN <- paste0(HOME, projectName, "/datasets/", DSname, "/outs/", type, "_feature_bc_matrix/")
  
  ### test the input, must lead from different inputs to same output
  print("Checking and sorting input")
  inputFolder <- paste0(HOME, projectName, "/analysis/Input/")
  if(file.exists(paste0(pathIN, "/genes.tsv")) == TRUE) {
    # load genes table, copy to features table 
    genes <- read.table(file = paste0(pathIN, "genes.tsv"), sep = "\t", header = FALSE)
    features <- genes
    # remove "/" and "_" from gene names and replaces with "-", necessary for printing out gene names
    levels(features$V2) <- gsub("/", "-", levels(features$V2)) ; levels(features$V2) <- gsub("_", "-", levels(features$V2))
    # change column names, does not depend on columns (only two)
    # change column names, depends on the number of columns
    if (length(colnames(features)) == 2) {colnames(features) <- c("ensID", "feature") } 
    if (length(colnames(features)) == 3) {colnames(features) <- c("ensID", "feature", "type") } 
    # save features table to file
    write.table(features, file = paste0(inputFolder, "features.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    # load features table
    features <- read.table(file = paste0(pathIN, "features.tsv"), sep = "\t", header = FALSE) 
    # remove "/" and "_" from gene names and replaces with "-", necessary for printing out gene names
    levels(features$V2) <- gsub("/", "-", levels(features$V2)) ; levels(features$V2) <- gsub("_", "-", levels(features$V2))
    # change column names, depends on the number of columns
    if (length(colnames(features)) == 2) {colnames(features) <- c("ensID", "feature") } 
    if (length(colnames(features)) == 3) {colnames(features) <- c("ensID", "feature", "type") } 
    # save features table to input files + copy as genes.tsv
    write.table(features, file = paste0(inputFolder, "features.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
    write.table(features, file = paste0(pathIN, "genes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  
  print("Converting matrix to Seurat object")
  # read matrix from file
  DS <- Read10X(data.dir = pathIN) 
  # change feature names to ensemble ID
  DS@Dimnames[[1]] <- as.character(features$ensID) 
  # convert to seurat object
  DS <- CreateSeuratObject(DS, project = DSname) 

  ### add metadata based on datasets file, necessary for distinguishing merged datasets and later visualizations
  print("Adding metadata")
  if(file.exists(paste0(HOME, projectName, "/datasets/datasets")) == TRUE) {
    datasets <- read.table(paste0(HOME, projectName, "/datasets/datasets"), sep = "\t", header = TRUE)
  } else { print("Datasets file not available, aborting.") ; quit(status = 14)}
  # add information for single dataset
  if (levels(DS@meta.data$orig.ident) %in% datasets$dataset) {
    for (i in colnames(datasets)) {
      DS@meta.data[i] <- datasets[match(DS@meta.data$orig.ident, datasets$dataset),i] } # end of for loop
  } # end of if statement
  # add information for multiple (merged) datasets
  if (str_detect(head(rownames(DS@meta.data), 1), "-") == TRUE) { 
    DS@meta.data$order <- as.factor(str_extract(string = rownames(DS@meta.data), "\\d"))
    # detect which datasets are being analyzed and subset the datasets table
    a <- which(str_detect(levels(DS@meta.data$orig.ident), as.character(datasets$dataset)) == TRUE)
    b <- datasets[a,]
    rownames(b) <- seq(1, length(rownames(b)))
    # add information to DS@meta.data (order is assumed by name, e.g. D6, E2, E12, F5, ...) 
    for (i in colnames(datasets)) {
      DS@meta.data[i] <- b[match(DS@meta.data$order, rownames(b)), i] } # end of for loop
  } # end of if statement
  
  # write R dataset to .Rds file
  print("Saving to file")
  pathOUT <- paste0(HOME, projectName, "/analysis/", DSname, "_", Sys.Date(), "/")
  dir.create(paste0(pathOUT, "RDS"), recursive = TRUE)
  info <- data.frame("parameter" = c("type", "format"), "value" = c(type, format))
  write.table(info, file = paste0(pathOUT, "RDS/setupInfo"), sep = "\t")
  saveRDS(DS, paste0(pathOUT, "RDS/", DSname, ".Rds"))
  
  # produce table to set values for image generation
  parameters <- c("pointSize","position","binSize","color","fill","shape","alpha","titleSize","axisTitle","axisXsize",
                  "axisXangle","axisXhjust","axisXvjust","axisYsize","axisYangle","axisYhjust","axisYvjust","axisLineSize",
                  "axisTicks","legendPositionX","legendPositionY","legendTitle","legendText","legendKeySize","legendBgrdCol",
                  "guideSize","fileType","dpi","width","height")
  
  boxplot <- c(1,"jitter",NA,"grey50","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  smoothDensity <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
                     0.7,0.7,0,35,35,"white",5,"png",350,10,10)
  scatter <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  hvgPlot <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "none","none",0,35,35,"white",5,"png",350,10,10)
  PCAplot <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  PCAloadings <- c(3,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
                   "right","none",0,35,35,"white",5,"png",350,10,10)
  tSNE <- c(2,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
            "right","none",0,35,35,"white",5,"png",350,10,10)
  UMAP <- c(2,"jitter",NA,"dataset","dataset",NA,0.5,35,0,0,0,0.5,0,0,0,0,0.5,1,1,
            "right","none",0,35,35,"white",5,"png",350,10,10)
  DM <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
          "right","none",0,35,35,"white",5,"png",350,10,10)
  barplot <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  sankey <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  DEheatmap <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  GO <- c(1,"jitter",NA,"dataset","dataset",NA,0.5,35,35,35,0,0.5,0,35,0,0,0.5,1,1,
               "right","none",0,35,35,"white",5,"png",350,10,10)
  
  
  imageStats <- data.frame(parameters = parameters, boxplot = boxplot, smoothDensity = smoothDensity, scatter = scatter,
                           hvgPlot=hvgPlot, PCAplot=PCAplot, PCAloadings=PCAloadings, tSNE=tSNE, UMAP=UMAP, DM=DM,
                           barplot=barplot, sankey=sankey, DEheatmap=DEheatmap, GO=GO)
 
  dir.create(paste0(pathOUT, "Input"))
  write.table(imageStats, file = paste0(pathOUT, "Input/imageStats.csv"), sep = ",", row.names = FALSE)
  print(paste("Finished. Seurat object for", DSname, "has been created and stored in", pathOUT))
  
} # end of else, mtx format

sessionInfo()
quit(status = 0)
