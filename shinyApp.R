############################ Shiny app programming ############################
###################### By Oliver Dietrich on 28 May 2019 ######################

library(shiny)

library(Seurat)
library(tidyverse)

DS <- readRDS("/home/od/Data/Cindrilla/analysis/G10_2019-03-02/G10_epithelial_2019-02-19/Rdatasets/G10_epithelial_labelled.Rds")
features <- read.table("/home/od/Data/Cindrilla/datasets/G12/outs/filtered_feature_bc_matrix/features.tsv", header = FALSE, sep = "\t")
colnames(features) <- c("ens_id", "gene", "type")

all_genes <- features$gene[match(rownames(DS@data), features$ens_id)]

# User Interface
ui <- fluidPage(
  # Input functions
  sliderInput(inputId = "alpha", label = "Choose a transparency level",
              value = 0.5, min = 0, max = 1),
  sliderInput(inputId = "size", label = "Choose dot size", value = 2, min = 0.1, max = 20),
  selectInput(inputId = "color", label = "Choose Gene", choices = all_genes),
  plotOutput("scatter")
  # Output functions
)

# Server Instructions
server <- function(input, output) {
  
  output$scatter <- renderPlot({
    
    ggplot(as.data.frame(DS@dr$umap@cell.embeddings), aes(
      x = UMAP1, y = UMAP2, col = DS@data[features$ens_id[match(input$color, features$gene)],]
    )) +
      geom_point(size = input$size, alpha = input$alpha)
  })
}

# execution
shinyApp(ui = ui, server = server)
