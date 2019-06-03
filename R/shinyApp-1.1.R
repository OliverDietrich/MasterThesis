############################ Shiny app programming ############################
###################### By Oliver Dietrich on 28 May 2019 ######################

library(shiny)

library(Seurat)
library(tidyverse)

DS <- readRDS("/home/od/Data/Cindrilla/analysis/G10_2019-03-02/G10_epithelial_2019-02-19/Rdatasets/G10_epithelial_labelled.Rds")
features <- read.table("/home/od/Data/Cindrilla/datasets/G12/outs/filtered_feature_bc_matrix/features.tsv", header = FALSE, sep = "\t")
colnames(features) <- c("ensID", "feature", "type")


# User Interface --------------------------------------------------------------
ui <- fluidPage(tabsetPanel(
  
  # Tab 1 - QC ------------------------
  tabPanel("Quality Control", 
           
           
           # Boxplot
           sidebarLayout(
             sidebarPanel(
               sliderInput(inputId = "boxplot_size", label = "Choose dot size", value = 2, min = 0.1, max = 20),
               sliderInput(inputId = "boxplot_alpha", label = "Change transparency",
                           value = 0.5, min = 0, max = 1),
               selectInput(inputId = "boxplot_selection", label = "Select category", choices = c("nGene", "nUMI", "pMito")),
               sliderInput(inputId = "boxplot_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
             ),
             mainPanel(
               plotOutput("boxplot")
             )
           ), # end of boxplot
           
           
           # Smooth Density
           sidebarLayout(
             sidebarPanel(
               sliderInput(inputId = "density_alpha", label = "Change transparency",
                           value = 0.5, min = 0, max = 1),
               selectInput(inputId = "density_selection", label = "Select category", choices = c("nGene", "nUMI", "pMito")),
               sliderInput(inputId = "density_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
             ),
             mainPanel(
               plotOutput("density")
             )
           ), # end of boxplot
           
           
           # Scatterplots
           sidebarLayout(
             sidebarPanel(
               sliderInput(inputId = "scatter_size", label = "Choose dot size", value = 2, min = 0.1, max = 20),
               sliderInput(inputId = "scatter_alpha", label = "Change transparency",
                           value = 0.5, min = 0, max = 1),
               selectInput(inputId = "scatter_selection", label = "Select category", choices = c("nUMI", "pMito")),
               sliderInput(inputId = "scatter_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
             ),
             mainPanel(
               plotOutput("scatter")
             )
           ) # end of scatterplots
           
           
  ), # end of tab 1
  
  # Tab2 - FS -------------------------
  tabPanel("Feature Selection",
  
           # HVG plot
           sidebarLayout(
             sidebarPanel(
               sliderInput(inputId = "hvg_plot_size", label = "Choose dot size", value = 2, min = 0.1, max = 20),
               selectInput(inputId = "hvg_plot_selection", label = "Select category", 
                           choices = c("gene.dispersion", "gene.dispersion.scaled")),
               sliderInput(inputId = "hvg_plot_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
             ),
             mainPanel(
               plotOutput("hvg_plot")
             )
           ) # end of hvg plot
           
           
          ), # end of tab 2


# Tab 3 - DR --------------------------
tabPanel("Dimensional Reduction", 
         
         "contents of tab 3"
         
         ) # end of tab 3


)) # end of input


# Server Instructions ---------------------------------------------------------
server <- function(input, output) {
  
  # Tab 1 - QC
  
    # boxplot
    output$boxplot <- renderPlot({
    as.data.frame(DS@meta.data) %>% ggplot(aes(
      x = dataset, y = DS@meta.data[,input$boxplot_selection]
    )) +
        geom_violin(aes(color = dataset)) +
        geom_point(size = input$boxplot_size, alpha = input$boxplot_alpha, position = "jitter") + 
        geom_boxplot(aes(fill = dataset), width = 0.8) +
        labs(y = input$boxplot_selection) + 
        theme(
          axis.title = element_text(size = input$boxplot_textsize),
          axis.text = element_text(size = input$boxplot_textsize),
          legend.title = element_text(size = input$boxplot_textsize),
          legend.text = element_text(size = input$boxplot_textsize)
          )
     }) # end of boxplot
  
    
    # smooth density
    output$density <- renderPlot({
      selection <- DS@meta.data[,input$density_selection]
      
      as.data.frame(DS@meta.data) %>% ggplot(aes(
        x = selection, fill = dataset, color = dataset)
      ) +
        geom_density(alpha = input$density_alpha) +
        labs(x = input$density_selection) + 
        theme(
          axis.title = element_text(size = input$density_textsize),
          axis.text = element_text(size = input$density_textsize),
          legend.text = element_text(size = input$density_textsize),
          legend.title = element_text(size = input$density_textsize)
        )
     }) # end of smooth density
    
    
    # scatterplots
    output$scatter <- renderPlot({
      selection <- DS@meta.data[,input$scatter_selection]
      
      as.data.frame(DS@meta.data) %>% ggplot(aes(
        x = selection, y = nGene, color = dataset)
      ) +
        geom_point(alpha = input$scatter_alpha, size = input$scatter_size) +
        geom_smooth(method = lm, color = "black") +
        stat_cor(method = "pearson", label.x = min(selection), label.y = max(DS@meta.data$nGene),
                 size = input$scatter_textsize/3, color = "black") +
        labs(y = "Number of Genes", x = input$scatter_selection) + 
        theme(
          axis.title = element_text(size = input$scatter_textsize),
          axis.text = element_text(size = input$scatter_textsize), 
          legend.text = element_text(size = input$scatter_textsize),
          legend.title = element_text(size = input$scatter_textsize)
        ) +
        guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$scatter_textsize/3)))
    }) # end of scatterplots
    
  
  # Tab 2 - FS
  
    
    # hvg plot
    output$hvg_plot <- renderPlot({
      varGenes <- rownames(DS@hvg.info) %in% DS@var.genes
      dispersion <- DS@hvg.info[,input$hvg_plot_selection]
      DS@hvg.info %>% ggplot(aes(
        x = gene.mean, y = dispersion, color = varGenes)
      ) +
        geom_point(size = input$hvg_plot_size) +
        scale_color_manual(values = c("black", "red")) +
        labs(x = input$hvg_plot_selection) + 
        theme(
          axis.title = element_text(size = input$hvg_plot_textsize),
          axis.text = element_text(size = input$hvg_plot_textsize),
          legend.text = element_text(size = input$hvg_plot_textsize),
          legend.title = element_text(size = input$hvg_plot_textsize)
        ) +
        guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$hvg_plot_textsize/3)))
    }) # end of hvg plot
    
    
  # Tab 3 - DR
  
} # end of server instructions

# execution -------------------------------------------------------------------
shinyApp(ui = ui, server = server)



