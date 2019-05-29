library(Seurat)
library(tidyverse)

DS <- readRDS("/home/oliver/Data/Leo_PR00224/analysis/G3G4G5H4H5H6_2019-05-24/RDS/G3G4G5H4H5H6.Rds")
features <- read.table("/home/oliver/Data/Leo_PR00224/analysis/Input/features.tsv", header = TRUE, sep = "\t")

# New approach
datasets <- c("G3", "G4", "G5", "H4", "H5", "H6")
color <- "dataset"
gene <- "LAG3"
ens_id <- as.character(features$ensID[match(gene, features$feature)])
metadata <- DS@meta.data[which(DS@meta.data$dataset %in% datasets),]
cells <- rownames(metadata)

data <- as.data.frame(DS@dr$umap@cell.embeddings[cells,])
data$dataset <- as.character(DS@meta.data[rownames(data),"dataset"])
data$gene <- DS@data[ens_id, rownames(data)]

plot <- ggplot(data, aes(
  x = UMAP1, y = UMAP2, col = gene
), alpha = 0.5) +
  geom_point(size = .5) +
  facet_wrap(~ dataset) +
  scale_color_continuous(low = "grey90", high = "blue2")

ggsave(plot, path = "/home/oliver/Data/Leo_PR00224/analysis/G3G4G5H4H5H6_2019-05-24/", filename = paste0("expr_", gene, "_facet.png"), dpi = 300,
       width = 9, height = 6)


# calculate mean
data_plus <- DS@dr$umap@cell.embeddings %>% as.data.frame()
data_plus$dataset <- as.character(DS@meta.data[rownames(data_plus), "dataset"])
data_plus$label <- DS@meta.data[rownames(data_plus), "label"]
data_plus$LAG3 <- DS@data[features$ensID[match("LAG3", features$feature)], rownames(data_plus)]

plot2 <- data_plus %>% group_by(dataset, label) %>% summarise(mean = mean(LAG3)) %>% ggplot(aes(x = dataset, y = mean, color = label)) + geom_point(size = 4)
ggsave(plot2, path = "/home/oliver/Data/Leo_PR00224/analysis/G3G4G5H4H5H6_2019-05-24/", filename = "LAG3_statistics.png", dpi = 300,
       width = 8, height = 8)
