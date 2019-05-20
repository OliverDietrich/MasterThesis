######################################################################
######## Feature Selection for 10x datasets in Seurat object #########
#--------------------------------------------------------------------#
########### R script by Oliver Dietrich on 24 April 2019 #############
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
library(ggpubr)
library(GGally)
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
# Load featureSelection table
featureSelection <- read.table("Input/featureSelection", sep = "\t", header = TRUE)
# Load dimensions table
dimensions <- read.table("RDS/dimensions", sep = " ", header = TRUE)
# Load seurat object
DS <- readRDS(paste0("RDS/", DSname, ".Rds"))

############################ Normalization ############################
print("Normalization")

normalization <- as.character(featureSelection$value[match("normalization", featureSelection$parameter)])

if (normalization == "logNorm") { DS <- NormalizeData(DS) }
if (normalization == "bayNorm") { print("add bayesian normalization of DS@data") }

##### Find highly variable genes
print("Find highly variable genes")

cutoffs <- featureSelection$value[match("Xl_Xh_Yl_Yh", featureSelection$parameter)]
DS <- FindVariableGenes(DS, x.low.cutoff = as.numeric(str_split(cutoffs, "_")[[1]][1]), 
                        x.high.cutoff = as.numeric(str_split(cutoffs, "_")[[1]][2]), 
                        y.cutoff = as.numeric(str_split(cutoffs, "_")[[1]][3]), 
                        y.high.cutoff = as.numeric(str_split(cutoffs, "_")[[1]][4]), 
                        do.plot = FALSE, mean.function = "ExpMean", dispersion.function = "LogVMR", display.progress = FALSE)
if(length(colnames(dimensions)) < 4) {dimensions <- dimensions %>% add_column("hvg" = c(length(DS@var.genes), NA))} else {
  dimensions$hvg <- c(length(DS@var.genes), NA)
}
write.table(dimensions, "RDS/dimensions")

############################ Visualize data ############################
print(paste("Visualizing feature selection for", DSname))

path <- "FS/HVG/"
dir.create(path, recursive = TRUE)
varGenes <- rownames(DS@hvg.info) %in% DS@var.genes

# variable gene plot
stats <- function(parameter) {  as.character(imageStats$hvgPlot[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$hvgPlot[match(parameter, imageStats$parameters)]))  }
if(imageStats$hvgPlot[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in c("gene.dispersion", "gene.dispersion.scaled")) {
  plot <- ggplot(DS@hvg.info) + 
    geom_point(aes(
      x=gene.mean, y=DS@hvg.info[,i], alpha=numstats("alpha"), color=varGenes
    )) +
    scale_color_manual(values = c("black","red")) +
    labs(
      x="Average Expression", y="Dispersion", title=""
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
    )
    ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_", str_replace_all(i, "\\.", "_"), ".", stats("fileType")),
         width=numstats("width"), height=numstats("height"), dpi=numstats("dpi")) } # end of plot

############################ Scaling and Regression ############################
print("Data scaling and regression of confounding variation")

vars2regress <- featureSelection$value[match("vars_2_regress", featureSelection$parameter)]
if(vars2regress %in% c("", "none")) {  DS <- ScaleData(DS, display.progress = FALSE) } else {
  vars2regress <- str_split(vars2regress, "_")[[1]]
  DS <- ScaleData(DS, vars.to.regress = vars2regress, display.progress = FALSE) }

##### PCA - dimensional reduction
print("Calculating PCA")
DS <- RunPCA(DS, pc.genes = DS@var.genes, do.print = FALSE, 
             pcs.compute = as.numeric(as.character(featureSelection$value[match("PCs", featureSelection$parameter)])))
DS <- ProjectPCA(DS, do.print = FALSE)

############################ Visualize data ############################
print(paste("Visualizing PCA for", DSname))

# variable gene plot
path <- "FS/PCA/CellEmbedding"
dir.create(path, recursive = TRUE)

stats <- function(parameter) {  as.character(imageStats$PCAplot[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$PCAplot[match(parameter, imageStats$parameters)]))  }
if(imageStats$PCAplot[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in colnames(DS@dr$pca@cell.embeddings)) {
  plot <- ggplot(as.data.frame(DS@dr$pca@cell.embeddings)) + 
    geom_point(aes(
      x=DS@dr$pca@cell.embeddings[,"PC1"], y=DS@dr$pca@cell.embeddings[,i], alpha=numstats("alpha"), 
      color=DS@meta.data[,stats("color")]
    ), size = numstats("pointSize")) +
    labs(
      x="PC1", y=i, title=""
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
    guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize"))))
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_", i, ".", stats("fileType")),
         dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } # end of plot

# Gene Loadings 
path <- "FS/PCA/GeneLoading"
dir.create(path, recursive = FALSE)

stats <- function(parameter) {  as.character(imageStats$PCAloadings[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$PCAloadings[match(parameter, imageStats$parameters)]))  }
if(imageStats$PCAloadings[match("legendPositionY", imageStats$parameters)] %in% "none") {
  legend_position <- stats("legendPositionX") 
} else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }

for(i in colnames(DS@dr$pca@cell.embeddings)) {
  data <- sort(DS@dr$pca@gene.loadings[,i]) %>% head(30) %>% as.data.frame()
  colnames(data) <- i
  plot <- ggplot(data) + 
    geom_point(aes(
      x=data[,i], y=reorder(features$feature[match(rownames(data), features$ensID)],data[,i])
    ), color = "blue", size = numstats("pointSize")) +
    labs(
      x=i, y="", title=""
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
    guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize"))))
  ggsave(plot = plot,
         path=path, filename=paste0(DSname,"_", i, ".", stats("fileType")),
         dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } # end of plot


# PC Heatmap
path <- "FS/PCA/Heatmap"
dir.create(path, recursive = FALSE)

for(i in colnames(DS@dr$pca@gene.loadings)) {
data <- DS@scale.data[rownames(DS@dr$pca@gene.loadings),]
data <- data[order(DS@dr$pca@gene.loadings[,i]),order(DS@dr$pca@cell.embeddings[,i])]
rownames(data) <- features$feature[match(rownames(data), features$ensID)]
data <- rbind(head(data, 15),tail(data, 15))
#pheatmap(data, show_colnames = FALSE, main = i, cluster_rows = FALSE, cluster_cols = FALSE)
plot <- pheatmap(data, show_colnames = FALSE, main = i, cluster_rows = FALSE, cluster_cols = FALSE)
ggsave(plot = plot, path = path, filename = paste0(DSname, "_PCheatmap_", i, ".", stats("fileType")))
}

rownames(DS@scale.data) <- features$feature[match(rownames(DS@scale.data), features$ensID)]
rownames(DS@dr$pca@gene.loadings) <- features$feature[match(rownames(DS@dr$pca@gene.loadings), features$ensID)]
for(i in seq(1,length(colnames(DS@dr$pca@cell.embeddings)),by = 1)) {
  png(paste0("FS/PCA/Heatmap/", DSname, "_PCheatmap2_PC", i, ".png"), width = 1500, height = 1500, res = 300)
print(PCHeatmap(DS, pc.use = i, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE))
dev.off()
}
rownames(DS@dr$pca@gene.loadings) <- features$ensID[match(rownames(DS@dr$pca@gene.loadings), features$feature)]
rownames(DS@scale.data) <- features$ensID[match(rownames(DS@scale.data), features$feature)]

# Elbow Plot
png(paste0("FS/PCA/", DSname, "_elbow.png"), width = 1920, height = 1080, res = 300)
print(PCElbowPlot(DS, num.pc = 40))
dev.off()

# GGally
a <- c(1,2,3,4,5) ; b <- c(1,6,7,8,9) ; c <- c(1,10,11,12,13) ; d <- c(1,14,15,16,17)
e <- c(1,18,19,20,21) ; f <- c(1,22,23,24,25) 
dat <- data.frame(a=a, b=b, c=c, d=d, e=e, f=f)

for(i in c("a", "b", "c", "d", "e", "f")) {
  PCpairs <- ggpairs(as.data.frame(DS@dr$pca@cell.embeddings[,dat[,i]]), progress = FALSE)
  ggsave(filename = paste0("PCpairs_", i, ".png"), path = "FS/PCA/", plot = PCpairs)
}

############################ Saving data ############################
print(paste("Saving dataset", DSname))

dimensionalReduction <- data.frame(
  parameter = c("PCs", "UMAP", "nNeighbors", "tSNE", "3DtSNE", "perplexity", "DM", "qUse", "SPRING"),
  value = c("1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20", "TRUE", "30", "TRUE", "FALSE", "10", "TRUE", "0.01", "FALSE")
)
write.table(dimensionalReduction, "dimensionalReduction", sep = "\t", row.names = FALSE)

# save dataset
saveRDS(DS, paste0("RDS/", DSname, ".Rds"))

sessionInfo()
quit(status = 0)
