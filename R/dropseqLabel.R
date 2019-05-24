######################################################################
####### Annotating clusters for 10x datasets in Seurat object ########
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
library(ggalluvial)
library(gplots)

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

############################ Labelling ############################
print("Adding labels to metadata")

resolution <- colnames(labels)[1] ; resName <- str_replace_all(resolution, "\\.", "_")

DS@meta.data$label <- labels$label[match(DS@meta.data[,resolution], labels[,resolution])]
DS@meta.data$labelColor <- labels$color[match(DS@meta.data[,resolution], labels[,resolution])]

############################ Visualization ############################
print(paste("Visualizing data in reduced dimension"))

# tSNE
if(DS@dr$tsne@key %in% "tSNE_") {
  path <- paste0(resName)
  dir.create(path, recursive = TRUE)
  
  stats <- function(parameter) {  as.character(imageStats$tSNE[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$tSNE[match(parameter, imageStats$parameters)]))  }
  if(imageStats$tSNE[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
    plot <- ggplot(as.data.frame(DS@dr$tsne@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$tsne@cell.embeddings[,"tSNE_1"], y=DS@dr$tsne@cell.embeddings[,"tSNE_2"], alpha=numstats("alpha"), 
        color=DS@meta.data[,"label"]
      ), size = numstats("pointSize")) +
      labs(
        x="tSNE_1", y="tSNE_2", title=""
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) +
      scale_color_manual(values = col2hex(labels$color[order(labels$label)]))
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_tSNE_label.", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } # end of plot

# UMAP
if(DS@dr$umap@key %in% "UMAP") {
  path <- paste0(resName)
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$UMAP[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$UMAP[match(parameter, imageStats$parameters)]))  }
  if(imageStats$UMAP[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
    
    plot <- ggplot(as.data.frame(DS@dr$umap@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$umap@cell.embeddings[,"UMAP1"], y=DS@dr$umap@cell.embeddings[,"UMAP2"], alpha=numstats("alpha"), 
        color=DS@meta.data[,"label"]
      ), size = numstats("pointSize")) +
      labs(
        x="UMAP1", y="UMAP2", title=""
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) +
      scale_color_manual(values = col2hex(labels$color[order(labels$label)]))
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_UMAP_", "label", ".", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } # end of plot

# Diffusion Map
if(DS@dr$dm@key %in% "DM") {
  path <- paste0(resName)
  dir.create(path, recursive = FALSE)
  
  stats <- function(parameter) {  as.character(imageStats$DM[match(parameter, imageStats$parameters)])  }
  numstats <- function(parameter) {  as.numeric(as.character(imageStats$DM[match(parameter, imageStats$parameters)]))  }
  if(imageStats$DM[match("legendPositionY", imageStats$parameters)] %in% "none") {
    legend_position <- stats("legendPositionX") 
  } else { legend_position <- c(numstats("legendPositionX"),numstats("legendPositionY")) }
  
    plot <- ggplot(as.data.frame(DS@dr$dm@cell.embeddings)) + 
      geom_point(aes(
        x=DS@dr$dm@cell.embeddings[,"DM1"], y=DS@dr$dm@cell.embeddings[,"DM2"], alpha=numstats("alpha"), 
        color=DS@meta.data[,"label"]
      ), size = numstats("pointSize")) +
      labs(
        x="DM1", y="DM2", title=""
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = numstats("guideSize")))) +
      scale_color_manual(values = col2hex(labels$color[order(labels$label)]))
    ggsave(plot = plot,
           path=path, filename=paste0(DSname,"_DM_", "label", ".", stats("fileType")),
           dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  ) } # end of plot 

### Barplot
path <- paste0(resName)
stats <- function(parameter) {  as.character(imageStats$barplot[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$barplot[match(parameter, imageStats$parameters)]))  }

data <- table(DS@meta.data$label, factor(as.character(DS@meta.data$dataset))) %>% as.data.frame() ; colnames(data) <- c("label", "dataset", "Freq") 
data <- data %>% mutate(Percentage = round(Freq/sum(data$Freq)*100, 1))

plot <- ggplot(data, aes(x=dataset, y=Percentage, label=Percentage, fill=label)) +
  #
  geom_col(position = position_stack(vjust=0.5)) +
  #
  geom_text(position = position_stack(vjust=0.5), size = numstats("textSize")) +
  #
  scale_fill_manual(values = col2hex(labels$color[match(levels(data$label), labels$label)])) +
  #
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
# end of ggplot
ggsave(plot = plot,
       path=path, filename=paste0(DSname,"_barplot.", stats("fileType")),
       dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  )


### Sankey plot
path <- paste0(resName)
stats <- function(parameter) {  as.character(imageStats$sankey[match(parameter, imageStats$parameters)])  }
numstats <- function(parameter) {  as.numeric(as.character(imageStats$sankey[match(parameter, imageStats$parameters)]))  }

plot <- ggplot(data, aes(
  y=Freq, axis1=dataset, axis2=label
)) +
  geom_alluvium(aes(
    fill = label), width = 1/12
    ) +
  geom_stratum(
    fill = "black", color = "grey", width = 1/12
    ) +
  geom_label(
    stat="stratum", label.strata=TRUE
  ) +
  scale_fill_manual(values = col2hex(labels$color[order(labels$label)]))
ggsave(plot = plot,
       path=path, filename=paste0(DSname,"_sankey.", stats("fileType")),
       dpi=numstats("dpi"), width=numstats("width"), height=numstats("height")  )

############################ Saving data ############################
print(paste("Saving dataset", DSname))

subset1 <- c("cell type 1", "cell type 2")
subset2 <- "cell type 3"
subsets <- data.frame(subset1 = subset1, subset2 = subset2)

write.table(subsets, file = "subsets", sep = "\t")

# save dataset
saveRDS(DS, paste0("RDS/", DSname, ".Rds"))

sessionInfo()
quit(status = 0)
