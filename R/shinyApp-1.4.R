###############################################################################
############################ Shiny app programming ############################
###################### By Oliver Dietrich on 28 May 2019 ######################
###############################################################################

####################
##### preamble #####
# library(shiny) ; library(shinydashboard) ; library(Seurat) ; library(tidyverse) ; library(ggpubr) ; library(pheatmap) ; library(gplots) ; library(ggalluvial) ; library(viridis)
# DS <- readRDS("/home/oliver/Data/Leo_PR00224/analysis/G3G4G5H4H5H6_2019-05-24/RDS/G3G4G5H4H5H6.Rds") ; features <- read.table("/home/oliver/Data/Leo_PR00224/analysis/Input/features.tsv", header = TRUE, sep = "\t")

##########################
##### User Interface #####
ui <- fluidPage(tabsetPanel(
  # Tab 1 - QC ----------------------------------------------------------------
  tabPanel("Quality Control", 
           fluidRow(
             # column 1
             column(4, # Boxplot
               plotOutput("boxplot"),
               # box 1
               box(title = "Image parameters", width = 6, background = "light-blue",
                   tabsetPanel(
                      # tab 1
                   tabPanel("1",
                            selectInput(inputId = "boxplot_selection", label = "Select category", choices = c("nGene", "nUMI", "pMito")),
                            sliderInput(inputId = "boxplot_size", label = "Choose dot size", value = 2, min = 0.1, max = 10)           
                            ), # end of tab 1
                     # tab 2
                   tabPanel("2",
                            sliderInput(inputId = "boxplot_alpha", label = "Change transparency", value = 0.5, min = 0, max = 1),
                            sliderInput(inputId = "boxplot_textsize", label = "Choose label size", value = 25, min = 1, max = 50)      
                   ) # end of tab 2
                   ) # end of tabset panel
                 ), # end of box 1
               # box 2
               box(title = "Download", width = 6,
                   tabsetPanel(
                     # tab 1
                     tabPanel("1",
                       textInput("boxplot_filename", label = "Enter filename", value = "Downloads/plot.png"),
                       actionButton("boxplot_download", "Save image")
                     ), # end of tab 1
                     # tab 2
                     tabPanel("2",
                       numericInput("boxplot_width", label = "Image width", value = 10, min = 5, max = 20),
                       numericInput("boxplot_height", label = "Image height", value = 10, min = 5, max = 20),
                       numericInput("boxplot_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                     ) # end of tab 2
                   ) # end of tabset panel
               ) # end of box 2
             ), # end of column 1
             # column 2
             column(4, # Smooth Density
                plotOutput("density"),
                # box 1
                box(title = "Image params", width = 6,
                    selectInput(inputId = "density_selection", label = "Select category", choices = c("nGene", "nUMI", "pMito")),
                    sliderInput(inputId = "density_alpha", label = "Change transparency",
                                value = 0.5, min = 0, max = 1),
                    sliderInput(inputId = "density_textsize", label = "Choose label size", value = 25, min = 1, max = 50)
                ), # end of box 1
                # box 2
                box(title = "Download", width = 6,
                    textInput("density_filename", label = "Enter filename", value = "Downloads/plot.png"),
                    actionButton("density_download", "Save image")      
                ) # end of box 2
                ), # end of column 2
             # column 3
             column(4, # Scatterplot
                plotOutput("scatter"),
                # box 1
                box(title = "Image params", width = 6,
                    selectInput(inputId = "scatter_selection", label = "Select category", choices = c("nUMI", "pMito")),
                    sliderInput(inputId = "scatter_size", label = "Choose dot size", value = 2, min = 0.1, max = 10),
                    sliderInput(inputId = "scatter_alpha", label = "Change transparency",
                                value = 0.5, min = 0, max = 1),
                    sliderInput(inputId = "scatter_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
                    ), # end of box 1
                # box 2
                box(title = "Download", width = 6,
                    textInput("scatter_filename", label = "Enter filename", value = "Downloads/plot.png"),
                    actionButton("scatter_download", "Save image")
                ) # end of box 2
                ) # end of column 3
           ) # end of fluid row
  ), # end of tab 1
  
  # Tab2 - FS -------------------------
  tabPanel("Feature Selection",
           tabsetPanel(
             # Part I
             tabPanel("Part I",
               # column 1
               column(4, # HVG
                      plotOutput("hvg_plot"),
                      # box 1
                      box(title = "Image params", width = 6,
                          selectInput(inputId = "hvg_plot_selection", label = "Select category", 
                                      choices = c("gene.dispersion", "gene.dispersion.scaled")),
                          sliderInput(inputId = "hvg_plot_size", label = "Choose dot size", value = 2, min = 0.1, max = 10),
                          sliderInput(inputId = "hvg_plot_textsize", label = "Choose label size", value = 25, min = 1, max = 50)
                          ), # end of box 1
                      # box 2
                      box(title = "Download", width = 6,
                          textInput("hvg_plot_filename", label = "Enter filename", value = "Downloads/plot.png"),
                          actionButton("hvg_plot_download", "Save image")
                          ) # end of box 2
               ), # end of column 1
               # column 2
               column(4, # elbow
                      plotOutput("PCelbow"),
                      # box 1
                      box(title = "Image params", width = 6,
                          sliderInput(inputId = "PCelbow_selection", label = "Select Range", value = 30, min = 10, max = 50),
                          sliderInput(inputId = "PCelbow_size", label = "Choose dot size", value = 4, min = 1, max = 10),
                          sliderInput(inputId = "PCelbow_textsize", label = "Choose label size", value = 25, min = 10, max = 40)
                      ), #end of box 1
                      # box 2
                      box(title = "Download", width = 6,
                          textInput("PCelbow_filename", label = "Enter filename", value = "Downloads/plot.png"),
                          actionButton("PCelbow_download", "Save image")
                      ) # end of box 2
               ), # end of column 2
               # column 3
               column(4, # embedding
                      plotOutput("PCembedding"),
                      # box 1
                      box(title = "Image params", width = 6,
                          selectInput(inputId = "PCembedding_selectionX", label = "Select category 1", 
                                      choices = colnames(DS@dr$pca@cell.embeddings)),
                          selectInput(inputId = "PCembedding_selectionY", label = "Select category 2", 
                                      choices = colnames(DS@dr$pca@cell.embeddings)),
                      sliderInput(inputId = "PCembedding_size", label = "Choose dot size", value = 2, min = 0.1, max = 10),
                      sliderInput(inputId = "PCembedding_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),

                      sliderInput(inputId = "PCembedding_textsize", label = "Choose label size", value = 25, min = 1, max = 50)
                      ), #end of box 1
                      # box 2
                      box(title = "Download", width = 6,
                      textInput("PCembedding_filename", label = "Enter filename", value = "Downloads/plot.png"),
                      actionButton("PCembedding_download", "Save image")
                      ) # end of box 2
               ) # end of column 3
             ), # end of Part I
             # Part II
             tabPanel("Part II",
               column(4, # heatmap
                      plotOutput("PCheatmap"),
                      box(title = "Image params", width = 6,
                      selectInput(inputId = "PCheatmap_selection", label = "Select category", 
                                  choices = colnames(DS@dr$pca@cell.embeddings)),
                      sliderInput(inputId = "PCheatmap_textsize", label = "Choose label size", value = 10, min = 0, max = 50)
                      ), #end of box 1
                      box(title = "Download", width = 6,
                      textInput("PCheatmap_filename", label = "Enter filename", value = "Downloads/plot.png"),
                      actionButton("PCheatmap_download", "Save image")
                      ) # end of box 2
               ), # end of column 1
               # column 2
               column(4, # loadings
                      plotOutput("PCloading"),
                      box(title = "Image params", width = 6,
                          selectInput(inputId = "PCloading_selection", label = "Select category", 
                                      choices = colnames(DS@dr$pca@cell.embeddings)),
                      sliderInput(inputId = "PCloading_size", label = "Choose dot size", value = 5, min = 0.1, max = 10),
                      sliderInput(inputId = "PCloading_range", label = "Choose range", value = 20, min = 5, max = 100),
                      sliderInput(inputId = "PCloading_textsize", label = "Choose label size", value = 20, min = 1, max = 50)
                      ), #end of box 1
                      box(title = "Download", width = 6,
                      textInput("PCloading_filename", label = "Enter filename", value = "Downloads/plot.png"),
                      actionButton("PCloading_download", "Save image")
                      ) # end of box 2
               ) # end of column 2
             ) # end of FS - tab 2
           ) # end of FS tabset panel
  ), # end of tab 2
  
  # Tab 3 - DR --------------------------
  tabPanel("Dimensional Reduction", 
           
           fluidRow( # 1
             # column 1
             column(4, # umap
                    plotOutput("umap"),
                    # box 1
                    box(title = "Image params", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                            selectInput(inputId = "dr_selection", label = "Select algorithm", choices = names(DS@dr)),
                            selectInput(inputId = "umap_selection", label = "Select category", 
                                        choices = colnames(DS@meta.data)),
                            checkboxGroupInput(inputId = "umap_datasets", label = "Choose datasets to display", inline = TRUE,
                                               choices = unique(as.character(DS@meta.data$dataset)), 
                                               selected = unique(as.character(DS@meta.data$dataset)))
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                            sliderInput(inputId = "umap_size", label = "Choose dot size", value = 4, min = 0.1, max = 10),
                            sliderInput(inputId = "umap_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),
                            sliderInput(inputId = "umap_labelsize", label = "Adjust point labels", value = 0, min = 0, max = 5.5),
                            sliderInput(inputId = "umap_textsize", label = "Choose label size", value = 0, min = 0, max = 40),
                            checkboxInput(inputId = "umap_axes", label = "Display axes", FALSE),
                            sliderInput(inputId = "umap_legendsize", label = "Choose legend size", value = 25, min = 0, max = 40)
                          ), # end of tab 2
                          # tab 3
                          tabPanel("3",
                            numericInput(inputId = "umap_dim1", label = "Choose dimension 1", value = 1, min = 1, max = 50, step = 1),
                            numericInput(inputId = "umap_dim2", label = "Choose dimension 2", value = 2, min = 1, max = 50, step = 1)
                          ) # end of tab 3
                        ) # end of tabset panel
                    ), #end of box 1
                    # box 2
                    box(title = "Download", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                                   textInput("umap_filename", label = "Enter filename", value = "Downloads/plot.png"),
                                   actionButton("umap_download", "Save image")
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                                   numericInput("umap_width", label = "Image width", value = 8, min = 2, max = 20),
                                   numericInput("umap_height", label = "Image height", value = 6, min = 2, max = 20),
                                   numericInput("umap_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                          ) # end of tab 2
                        ) # end of tabset panel
                    ) # end of box 2
             ) # end of column 1
           ) # end of fluid row           
  ), # end of tab 3

  # Tab 4 - Expression -------------------------
  tabPanel("Expression",
           # row 1
           fluidRow(
             # column 1
             column(4, # umap_expr
                    plotOutput("umap_expr"),
                    # box 1
                    box(title = "Image params", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                                   selectInput(inputId = "dr_expr_selection", label = "Select algorithm", choices = names(DS@dr)),
                                   selectInput(inputId = "umap_expr_selection", label = "Select gene", 
                                               choices = levels(features$feature[match(rownames(DS@data), features$ensID)])),
                                   checkboxGroupInput(inputId = "umap_expr_datasets", label = "Choose datasets to display", inline = TRUE,
                                                      choices = unique(as.character(DS@meta.data$dataset)), 
                                                      selected = unique(as.character(DS@meta.data$dataset)))
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                                   sliderInput(inputId = "umap_expr_size", label = "Choose dot size", value = 4, min = 1, max = 10),
                                   sliderInput(inputId = "umap_expr_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),
                                   sliderInput(inputId = "umap_expr_labelsize", label = "Adjust point labels", value = 0, min = 0, max = 5.5),
                                   sliderInput(inputId = "umap_expr_textsize", label = "Choose label size", value = 0, min = 0, max = 40),
                                   checkboxInput(inputId = "umap_expr_axes", label = "Display axes", FALSE),
                                   sliderInput(inputId = "umap_expr_legendsize", label = "Choose legend size", value = 25, min = 0, max = 40)
                                   ), # end of tab 2
                          # tab 3
                          tabPanel("3",
                                   numericInput(inputId = "umap_expr_dim1", label = "Choose dimension 1", value = 1, min = 1, max = 50, step = 1),
                                   numericInput(inputId = "umap_expr_dim2", label = "Choose dimension 2", value = 2, min = 1, max = 50, step = 1)
                          ) # end of tab 3
                        ) # end of tabset panel
                    ), #end of box 1
                    # box 2
                    box(title = "Download", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                                   textInput("umap_expr_filename", label = "Enter filename", value = "Downloads/plot.png"),
                                   actionButton("umap_expr_download", "Save image")
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                                   numericInput("umap_expr_width", label = "Image width", value = 8, min = 2, max = 20),
                                   numericInput("umap_expr_height", label = "Image height", value = 6, min = 2, max = 20),
                                   numericInput("umap_expr_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                          ) # end of tab 2
                        ) # end of tabset panel
                    ) # end of box 2
             ), # end of column 1
             # column 2
             column(4, # dimplot facet
               plotOutput("expr_scatter"),
               # box 1
               box(title = "Image params", width = 6,
                 tabsetPanel(
                   # tab 1
                   tabPanel("1",
                            selectInput(inputId = "expr_scatter_selectionX", label = "Select gene 1", 
                                        choices = levels(features$feature[match(rownames(DS@data), features$ensID)])),
                            selectInput(inputId = "expr_scatter_selectionY", label = "Select gene 2", 
                                        choices = levels(features$feature[match(rownames(DS@data), features$ensID)]),
                                        selected = levels(features$feature[match(rownames(DS@data), features$ensID)])[2]),
                            selectInput(inputId = "expr_scatter_selection", label = "Select category", 
                                        choices = colnames(DS@meta.data))
                   ), # end of tab 1
                   # tab 2
                   tabPanel("2",
                            sliderInput(inputId = "expr_scatter_size", label = "Choose dot size", value = 4, min = 1, max = 10),
                            sliderInput(inputId = "expr_scatter_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),
                            sliderInput(inputId = "expr_scatter_textsize", label = "Choose label size", value = 25, min = 0, max = 40),
                            checkboxInput(inputId = "expr_scatter_axes", label = "Display axes", TRUE),
                            sliderInput(inputId = "expr_scatter_legendsize", label = "Choose legend size", value = 25, min = 0, max = 40)
                    ) # end of tab 2
                 ) # end of tabset panel
               ), # end of box 1
               # box 2
               box(title = "Download", width = 6,
                 tabsetPanel(
                   # tab 1
                   tabPanel("1",
                            textInput("expr_scatter_filename", label = "Enter filename", value = "Downloads/plot.png"),
                            actionButton("expr_scatter_download", "Save image")
                   ), # end of tab 1
                   # tab 2
                   tabPanel("2",
                            numericInput("expr_scatter_width", label = "Image width", value = 8, min = 2, max = 20),
                            numericInput("expr_scatter_height", label = "Image height", value = 6, min = 2, max = 20),
                            numericInput("expr_scatter_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                   ) # end of tab 2
                 ) # end of tabset panel
               ) # end of box 2
             ), # end of column 2
             # column 3
             column(4, # expr boxplot
                    plotOutput("expr_boxplot"),
                    # box 1
                    box(title = "Image params", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                                   selectInput(inputId = "expr_boxplot_gene", label = "Select gene", 
                                               choices = levels(features$feature[match(rownames(DS@data), features$ensID)])),
                                   selectInput(inputId = "expr_boxplot_selection", label = "Select category", 
                                               choices = colnames(DS@meta.data))
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                                   sliderInput(inputId = "expr_boxplot_size", label = "Choose dot size", value = 4, min = 1, max = 10),
                                   sliderInput(inputId = "expr_boxplot_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),
                                   sliderInput(inputId = "expr_boxplot_textsize", label = "Choose label size", value = 25, min = 0, max = 40),
                                   checkboxInput(inputId = "expr_boxplot_axes", label = "Display axes", TRUE),
                                   sliderInput(inputId = "expr_boxplot_legendsize", label = "Choose legend size", value = 25, min = 0, max = 40)
                          ) # end of tab 2
                        ) # end of tabset panel
                    ), # end of box 1
                    # box 2
                    box(title = "Download", width = 6,
                        tabsetPanel(
                          # tab 1
                          tabPanel("1",
                                   textInput("expr_boxplot_filename", label = "Enter filename", value = "Downloads/plot.png"),
                                   actionButton("expr_boxplot_download", "Save image")
                          ), # end of tab 1
                          # tab 2
                          tabPanel("2",
                                   numericInput("expr_boxplot_width", label = "Image width", value = 8, min = 2, max = 20),
                                   numericInput("expr_boxplot_height", label = "Image height", value = 6, min = 2, max = 20),
                                   numericInput("expr_boxplot_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                          ) # end of tab 2
                        ) # end of tabset panel
                    ) # end of box 2
             ) # end of column 3
           ), # end of fluid row  1
           # row 2
           fluidRow(
             # column 1
             column(4, # dimplot facet
                    plotOutput("expr_dr_facet"),
                    # box 1
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                selectInput(inputId = "expr_dr_facet_algorithm", label = "Select algorithm", choices = names(DS@dr)),
                                selectInput(inputId = "expr_dr_facet_gene", label = "Select gene", 
                                      choices = levels(features$feature[match(rownames(DS@data), features$ensID)]))
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 sliderInput(inputId = "expr_dr_facet_size", label = "Choose dot size", value = 4, min = 1, max = 10),
                                 sliderInput(inputId = "expr_dr_facet_transparency", label = "Change Transparency", value = 0.5, min = 0, max = 1),
                                 sliderInput(inputId = "expr_dr_facet_textsize", label = "Choose label size", value = 15, min = 0, max = 40),
                                 checkboxInput(inputId = "expr_dr_facet_axes", label = "Display axes", TRUE),
                                 sliderInput(inputId = "expr_dr_facet_legendsize", label = "Choose legend size", value = 25, min = 0, max = 40)
                        ), # end of tab 2
                        # tab 3
                        tabPanel("3",
                                 numericInput(inputId = "expr_dr_facet_dim1", label = "Choose dimension 1", value = 1, min = 1, max = 50, step = 1),
                                 numericInput(inputId = "expr_dr_facet_dim2", label = "Choose dimension 2", value = 2, min = 1, max = 50, step = 1)
                        ) # end of tab 3
                      ) # end of tabset panel
                    ), # end of box 1
                    # box 2
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                 textInput("expr_dr_facet_filename", label = "Enter filename", value = "Downloads/plot.png"),
                                 actionButton("expr_dr_facet_download", "Save image")
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 numericInput("expr_dr_facet_width", label = "Image width", value = 8, min = 2, max = 20),
                                 numericInput("expr_dr_facet_height", label = "Image height", value = 6, min = 2, max = 20),
                                 numericInput("expr_dr_facet_dpi", label = "Image resolution", value = 300, min = 100, max = 1000)
                        ) # end of tab 2
                      ) # end of tabset panel
                    ) # end of box 2
             ), # end of column 1
             # column 2
             column(4, # some plot
                    #plotOutput("some_plot"),
                    # box 1
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                 "tab 1"
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 "tab 1"
                        ) # end of tab 2
                      ) # end of tabset panel
                    ), # end of box 1
                    # box 2
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                 "tab 1"
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 "tab 1"
                        ) # end of tab 2
                      ) # end of tabset panel
                    ) # end of box 2
             ), # end of column 2
             # column 3
             column(4, # some plot
                    #plotOutput("some_plot"),
                    # box 1
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                 "tab 1"
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 "tab 1"
                        ) # end of tab 2
                      ) # end of tabset panel
                    ), # end of box 1
                    # box 2
                    box(
                      tabsetPanel(
                        # tab 1
                        tabPanel("1",
                                 "tab 1"
                        ), # end of tab 1
                        # tab 2
                        tabPanel("2",
                                 "tab 1"
                        ) # end of tab 2
                      ) # end of tabset panel
                    ) # end of box 2
             ) # end of column 3
           ) # end of row 2
  ), # end of tab 4
  
  # Tab 5 - Label --------------------------
  tabPanel("Labelled Dataset", 
           
           fluidRow( # 1
             # column 1
             column(4, # Sankey plot
                    plotOutput("sankey"),
                    # box 1
                    box(title = "Image params", width = 6,
                        sliderInput(inputId = "sankey_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
                    ), #end of box 1
                    # box 2
                    box(title = "Download", width = 6,
                        textInput("sankey_filename", label = "Enter filename", value = "Downloads/plot.png"),
                        actionButton("sankey_download", "Save image")
                    ) # end of box 2
             ), # end of column 1
             # column 2
             column(4, # Composition barplot
                    plotOutput("comp_barplot"),
                    # box 1
                    box(title = "Image params", width = 6,
                        selectInput(inputId = "comp_barplot", label = "Select category", 
                                    choices = colnames(DS@dr$pca@cell.embeddings)),
                        sliderInput(inputId = "comp_barplot_size", label = "Choose text size", value = 5, min = 1, max = 10),
                        sliderInput(inputId = "comp_barplot_textsize", label = "Choose label size", value = 20, min = 1, max = 100)
                    ), #end of box 1
                    # box 2
                    box(title = "Download", width = 6,
                      textInput("comp_barplot_filename", label = "Enter filename", value = "Downloads/plot.png"),
                      actionButton("comp_barplot_download", "Save image")
                    ) # end of box 2
             ) # end of column 2
           ) # end of fluid row
  ) # end of tab 5
  
)) # end of input

###############################
##### Server Instructions #####
server <- function(input, output) {
  # Tab 1 - QC ----------------------------------------------------------------
  
  # boxplot
  output$boxplot <- renderPlot({
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
      x = dataset, y = DS@meta.data[,input$boxplot_selection]
    )) +
      geom_point(size = input$boxplot_size, alpha = input$boxplot_alpha, position = "jitter") + 
      geom_violin(aes(fill = dataset, alpha = input$boxplot_alpha)) +
      geom_boxplot(aes(fill = dataset), width = 0.8) +
      labs(y = input$boxplot_selection) + 
      theme(
        axis.title = element_text(size = input$boxplot_textsize),
        axis.text = element_text(size = input$boxplot_textsize),
        legend.title = element_text(size = input$boxplot_textsize),
        legend.text = element_text(size = input$boxplot_textsize)
      ) +
      guides(alpha = FALSE)
    print(plot)
  })
  observeEvent(input$boxplot_download, {
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
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
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$boxplot_filename,
           width = input$boxplot_width, height = input$boxplot_height, dpi = input$boxplot_dpi)
  }) # end of boxplot
  
  
  # smooth density
  output$density <- renderPlot({
    selection <- DS@meta.data[,input$density_selection]
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
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
    print(plot)
  })
  observeEvent(input$density_download, {
    selection <- DS@meta.data[,input$density_selection]
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
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
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$density_filename)
  }) # end of smooth density
  
  
  # scatterplots
  output$scatter <- renderPlot({
    selection <- DS@meta.data[,input$scatter_selection]
    
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$scatter_textsize/3, alpha = 1)))
    print(plot)
  })
  observeEvent(input$scatter_download, {
    selection <- DS@meta.data[,input$scatter_selection]
    
    plot <- as.data.frame(DS@meta.data) %>% ggplot(aes(
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$scatter_textsize/3, alpha = 1)))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$scatter_filename)
  }) # end of scatterplots
  
  # Tab 2 - FS --------------------
  
  # hvg plot
  output$hvg_plot <- renderPlot({
    varGenes <- rownames(DS@hvg.info) %in% DS@var.genes
    dispersion <- DS@hvg.info[,input$hvg_plot_selection]
    plot <- DS@hvg.info %>% ggplot(aes(
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$hvg_plot_textsize/3, alpha = 1)))
    print(plot)
  })
  observeEvent(input$hvg_plot_download, {
    varGenes <- rownames(DS@hvg.info) %in% DS@var.genes
    dispersion <- DS@hvg.info[,input$hvg_plot_selection]
    plot <- DS@hvg.info %>% ggplot(aes(
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
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$hvg_plot_textsize/3, alpha = 1)))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$hvg_plot_filename)
  }) # end of hvg plot
  
  # PCelbow
  output$PCelbow <- renderPlot({
    data <- DS@dr$pca@sdev[1:input$PCelbow_selection] %>% as.data.frame() %>% rownames_to_column() ; colnames(data) <- c("PC", "sdev")
    data$PC <- as.numeric(data$PC)
      plot <- ggplot(data, aes(x = PC, y = sdev)) + 
      geom_point(size = input$PCelbow_size) +
        labs(x = "PC", y = "Standard Deviation") +
      theme(
        axis.title = element_text(size = input$PCelbow_textsize),
        axis.text = element_text(size = input$PCelbow_textsize),
        legend.text = element_text(size = input$PCelbow_textsize),
        legend.title = element_text(size = input$PCelbow_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCelbow_textsize/3, alpha = 1)))
      print(plot)
  })
  observeEvent(input$PCelbow_download, {
    data <- DS@dr$pca@sdev[1:input$PCelbow_selection] %>% as.data.frame() %>% rownames_to_column() ; colnames(data) <- c("PC", "sdev")
    data$PC <- as.numeric(data$PC)
    plot <- ggplot(data, aes(x = PC, y = sdev)) + 
      geom_point(size = input$PCelbow_size) +
      labs(x = "PC", y = "Standard Deviation") +
      theme(
        axis.title = element_text(size = input$PCelbow_textsize),
        axis.text = element_text(size = input$PCelbow_textsize),
        legend.text = element_text(size = input$PCelbow_textsize),
        legend.title = element_text(size = input$PCelbow_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCelbow_textsize/3, alpha = 1)))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$PCelbow_filename)
  }) # end of PCelbow
  
  # PCembedding
  output$PCembedding <- renderPlot({
    data <- DS@dr$pca@cell.embeddings %>% as.data.frame()
    dataset <- DS@meta.data$dataset
    plot <- ggplot(data, aes(
      y = data[,input$PCembedding_selectionY], x = data[,input$PCembedding_selectionX], color = dataset)) +
      geom_point(size = input$PCembedding_size, alpha = input$PCembedding_transparency) +
      labs(x = input$PCembedding_selectionX, y = input$PCembedding_selectionY) +
      theme(
        axis.title = element_text(size = input$PCembedding_textsize),
        axis.text = element_text(size = input$PCembedding_textsize),
        legend.text = element_text(size = input$PCembedding_textsize),
        legend.title = element_text(size = input$PCembedding_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCembedding_textsize/3), alpha = 1))
    print(plot)
  }) 
  observeEvent(input$PCembedding_download, {
    data <- DS@dr$pca@cell.embeddings %>% as.data.frame()
    dataset <- DS@meta.data$dataset
    plot <- ggplot(data, aes(
      y = data[,input$PCembedding_selectionY], x = data[,input$PCembedding_selectionX], color = dataset)) +
      geom_point(size = input$PCembedding_size, alpha = input$PCembedding_transparency) +
      labs(x = input$PCembedding_selectionX, y = input$PCembedding_selectionY) +
      theme(
        axis.title = element_text(size = input$PCembedding_textsize),
        axis.text = element_text(size = input$PCembedding_textsize),
        legend.text = element_text(size = input$PCembedding_textsize),
        legend.title = element_text(size = input$PCembedding_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCembedding_textsize/3, alpha = 1)))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$PCembedding_filename)
  })
  # end of PCembedding
  
  # PCloading
  output$PCloading <- renderPlot({
    data <- sort(DS@dr$pca@gene.loadings[,input$PCloading_selection]) %>% head(input$PCloading_range) %>% as.data.frame()
    rownames(data) <- features$feature[match(rownames(data), features$ensID)] ; colnames(data) <- "score"
    plot <- ggplot(data, aes(
      x = score, y = reorder(rownames(data), rev(score)))) +
      geom_point(size = input$PCloading_size, color = "blue2") +
      labs(x = input$PCloading_selection, y = "") +
      theme(
        axis.title = element_text(size = input$PCloading_textsize),
        axis.text = element_text(size = input$PCloading_textsize),
        legend.text = element_text(size = input$PCloading_textsize),
        legend.title = element_text(size = input$PCloading_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCloading_textsize/3, alpha = 1)))
    print(plot)
  })
  observeEvent(input$PCloading_download, {
    data <- sort(DS@dr$pca@gene.loadings[,input$PCloading_selection]) %>% head(input$PCloading_range) %>% as.data.frame()
    rownames(data) <- features$feature[match(rownames(data), features$ensID)] ; colnames(data) <- "score"
    plot <- ggplot(data, aes(
      x = score, y = reorder(rownames(data), rev(score)))) +
      geom_point(size = input$PCloading_size, color = "blue2") +
      labs(x = input$PCloading_selection, y = "") +
      theme(
        axis.title = element_text(size = input$PCloading_textsize),
        axis.text = element_text(size = input$PCloading_textsize),
        legend.text = element_text(size = input$PCloading_textsize),
        legend.title = element_text(size = input$PCloading_textsize)
      ) +
      guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$PCloading_textsize/3, alpha = 1)))
    print(plot)
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$PCloading_filename)
  }) # end of PCembedding
  
  # PCheatmap
  output$PCheatmap <- renderPlot({
    data <- DS@scale.data[rownames(DS@dr$pca@gene.loadings),]
    data <- data[order(DS@dr$pca@gene.loadings[,input$PCheatmap_selection]),order(DS@dr$pca@cell.embeddings[,input$PCheatmap_selection])]
    rownames(data) <- features$feature[match(rownames(data), features$ensID)]
    data <- rbind(head(data, 15),tail(data, 15))
    plot <-pheatmap(data, show_colnames = FALSE, main = input$PCheatmap_selection, cluster_rows = FALSE, cluster_cols = FALSE, silent = TRUE, 
                    fontsize = input$PCheatmap_textsize)
    print(plot)
  })
  observeEvent(input$PCheatmap_download, {
    data <- DS@scale.data[rownames(DS@dr$pca@gene.loadings),]
    data <- data[order(DS@dr$pca@gene.loadings[,input$PCheatmap_selection]),order(DS@dr$pca@cell.embeddings[,input$PCheatmap_selection])]
    rownames(data) <- features$feature[match(rownames(data), features$ensID)]
    data <- rbind(head(data, 15),tail(data, 15))
    plot <-pheatmap(data, show_colnames = FALSE, main = input$PCheatmap_selection, cluster_rows = FALSE, cluster_cols = FALSE, silent = TRUE, 
                    fontsize = input$PCheatmap_textsize)
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$PCheatmap_filename, dpi = 300)
  }) # end of PCheatmap (aborts after printing one heatmap ...)
  
  # Tab 3 - DR --------------------
  
  # umap
  output$umap <- renderPlot({
    cells <- which(DS@meta.data$dataset %in% input$umap_datasets)
    data <- DS@dr[[input$dr_selection]]@cell.embeddings[,] %>% as.data.frame() %>% mutate(selection = DS@meta.data[, input$umap_selection])
    if(class(DS@meta.data[,input$umap_selection]) %in% c("numeric", "integer")) {
           guides <- guides(alpha = FALSE, color = guide_colorbar(override.aes = list(size = input$umap_legendsize/2, alpha = 1), 
                                                                  barwidth = input$umap_legendsize/15, 
                                                                  barheight = input$umap_legendsize*0.75))
           color <- scale_color_viridis() } else {
           guides <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$umap_legendsize/3, alpha = 1)))
           color <- geom_blank()
           }
    if (input$umap_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    plot <- data[cells,] %>% ggplot(aes(
      x = data[cells,colnames(data)[input$umap_dim1]], y = data[cells,colnames(data)[input$umap_dim2]] 
    )) +
      geom_point(aes(color = selection), size = input$umap_size, alpha = input$umap_transparency) +
      geom_text(aes(label = selection), size = input$umap_labelsize) +
      labs(color = input$umap_selection, x = colnames(data)[input$umap_dim1], y = colnames(data)[input$umap_dim2]) +
      theme(
        axis.title = element_text(size = input$umap_textsize),
        axis.text = element_text(size = input$umap_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$umap_legendsize),
        legend.title = element_text(size = input$umap_legendsize)
      ) +
      guides +
      color +
      geom_blank(data = data, aes(
        x = data[,colnames(data)[input$umap_dim1]], y = data[,colnames(data)[input$umap_dim2]], inherit.aes = FALSE))
    print(plot)
  }) 
  observeEvent(input$umap_download, {
    cells <- which(DS@meta.data$dataset %in% input$umap_datasets)
    data <- DS@dr[[input$dr_selection]]@cell.embeddings[,] %>% as.data.frame() %>% mutate(selection = DS@meta.data[, input$umap_selection])
    if(class(DS@meta.data[,input$umap_selection]) %in% c("numeric", "integer")) {
      guides <- guides(alpha = FALSE, color = guide_colorbar(override.aes = list(size = input$umap_legendsize/2, alpha = 1), 
                                                             barwidth = input$umap_legendsize/15, 
                                                             barheight = input$umap_legendsize*0.75))
      color <- scale_color_viridis() } else {
        guides <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$umap_legendsize/3, alpha = 1)))
        color <- geom_blank()
      }
    if (input$umap_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    plot <- data[cells,] %>% ggplot(aes(
      x = data[cells,colnames(data)[input$umap_dim1]], y = data[cells,colnames(data)[input$umap_dim2]] 
    )) +
      geom_point(aes(color = selection), size = input$umap_size, alpha = input$umap_transparency) +
      geom_text(aes(label = selection), size = input$umap_labelsize) +
      labs(color = input$umap_selection, x = colnames(data)[input$umap_dim1], y = colnames(data)[input$umap_dim2]) +
      theme(
        axis.title = element_text(size = input$umap_textsize),
        axis.text = element_text(size = input$umap_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$umap_legendsize),
        legend.title = element_text(size = input$umap_legendsize)
      ) +
      guides +
      color +
      geom_blank(data = data, aes(
        x = data[,colnames(data)[input$umap_dim1]], y = data[,colnames(data)[input$umap_dim2]], inherit.aes = FALSE))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$umap_filename,
           width = input$umap_width, height = input$umap_height, dpi = input$umap_dpi)
  }) # end of umap
  
  # Tab 4 - Expression ------------------
  
  # expression dimplot
  output$umap_expr <- renderPlot({
    gene <- features$ensID[match(input$umap_expr_selection, features$feature)]
    cells <- which(DS@meta.data$dataset %in% input$umap_expr_datasets)
    if (input$umap_expr_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- DS@dr[[input$dr_expr_selection]]@cell.embeddings[,] %>% as.data.frame() %>% mutate(selection = DS@data[gene,])
    plot <- data[cells,] %>% ggplot(aes(
      x = data[cells,colnames(data)[input$umap_expr_dim1]], y = data[cells,colnames(data)[input$umap_expr_dim2]], color = selection)) +
      geom_point(size = input$umap_expr_size, alpha = input$umap_expr_transparency) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2") +
      labs(color = input$umap_expr_selection) +
      theme(
        plot.title = element_text(size = input$umap_expr_legendsize),
        axis.title = element_text(size = input$umap_expr_textsize),
        axis.text = element_text(size = input$umap_expr_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$umap_expr_legendsize),
        legend.title = element_text(size = input$umap_expr_legendsize)
      ) +
      guides(color = guide_colorbar(barwidth = input$umap_expr_legendsize/15, 
                                    barheight = input$umap_legendsize*0.75)) +
      geom_blank(data = data, aes(
        x = data[,colnames(data)[input$umap_expr_dim1]], y = data[,colnames(data)[input$umap_expr_dim2]], inherit.aes = FALSE))
    print(plot)
  }) 
  observeEvent(input$umap_expr_download, {
    gene <- features$ensID[match(input$umap_expr_selection, features$feature)]
    cells <- which(DS@meta.data$dataset %in% input$umap_expr_datasets)
    if (input$umap_expr_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- DS@dr[[input$dr_expr_selection]]@cell.embeddings[,] %>% as.data.frame() %>% mutate(selection = DS@data[gene,])
    plot <- data[cells,] %>% ggplot(aes(
      x = data[cells,colnames(data)[input$umap_expr_dim1]], y = data[cells,colnames(data)[input$umap_expr_dim2]], color = selection)) +
      geom_point(size = input$umap_expr_size, alpha = input$umap_expr_transparency) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2") +
      labs(color = input$umap_expr_selection) +
      theme(
        plot.title = element_text(size = input$umap_expr_legendsize),
        axis.title = element_text(size = input$umap_expr_textsize),
        axis.text = element_text(size = input$umap_expr_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$umap_expr_legendsize),
        legend.title = element_text(size = input$umap_expr_legendsize)
      ) +
      guides(color = guide_colorbar(barwidth = input$umap_expr_legendsize/15, 
                                    barheight = input$umap_legendsize*0.75)) +
      geom_blank(data = data, aes(
        x = data[,colnames(data)[input$umap_expr_dim1]], y = data[,colnames(data)[input$umap_expr_dim2]], inherit.aes = FALSE))
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$umap_expr_filename,
           width = input$umap_expr_width, height = input$umap_expr_height, dpi = input$umap_expr_dpi)
  }) # end of dimplot of expression
  
  # expression correlation scatter
  output$expr_scatter <- renderPlot({
    if(class(DS@meta.data[,input$expr_scatter_selection]) %in% c("numeric", "integer")) {
      guides <- guides(alpha = FALSE, color = guide_colorbar(override.aes = list(size = input$expr_scatter_legendsize/2, alpha = 1), 
                                                             barwidth = input$expr_scatter_legendsize/15, 
                                                             barheight = input$expr_scatter_legendsize*0.75))
      color <- scale_color_viridis() } else {
        guides <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$expr_scatter_legendsize/3, alpha = 1)))
        color <- geom_blank()
      }
    gene_x <- as.character(features$ensID[match(input$expr_scatter_selectionX, features$feature)])
    gene_y <- as.character(features$ensID[match(input$expr_scatter_selectionY, features$feature)])
    if (input$expr_scatter_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- t(DS@data[c(gene_x, gene_y),]) %>% as.matrix() %>% as.data.frame() %>% mutate(selection = DS@meta.data[,input$expr_scatter_selection],
                                                                                          dataset = DS@meta.data$dataset,
                                                                                          label = DS@meta.data$label)
    data %>% ggplot(aes(
      x = data[,gene_x], y = data[,gene_y], color = selection
    )) +
      geom_point() +
      geom_smooth(method = lm, color = "black") +
      stat_cor(method = "pearson", label.x = min(data[, gene_x]), label.y = max(data[, gene_y]),
               size = input$expr_scatter_textsize/3, color = "black") +
      labs(x = input$expr_scatter_selectionX, y = input$expr_scatter_selectionY, color = input$expr_scatter_selection) +
      theme(
        plot.title = element_text(size = input$expr_scatter_legendsize),
        axis.title = element_text(size = input$expr_scatter_textsize),
        axis.text = element_text(size = input$expr_scatter_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$expr_scatter_legendsize),
        legend.title = element_text(size = input$expr_scatter_legendsize)
      ) +
      guides +
      color + facet_wrap(~ dataset)
  })
  observeEvent(input$expr_scatter_download, {
    if(class(DS@meta.data[,input$expr_scatter_selection]) %in% c("numeric", "integer")) {
      guides <- guides(alpha = FALSE, color = guide_colorbar(override.aes = list(size = input$expr_scatter_legendsize/2, alpha = 1), 
                                                             barwidth = input$expr_scatter_legendsize/15, 
                                                             barheight = input$expr_scatter_legendsize*0.75))
      color <- scale_color_viridis() } else {
        guides <- guides(alpha = FALSE, color = guide_legend(override.aes = list(size = input$expr_scatter_legendsize/3, alpha = 1)))
        color <- geom_blank()
      }
    gene_x <- as.character(features$ensID[match(input$expr_scatter_selectionX, features$feature)])
    gene_y <- as.character(features$ensID[match(input$expr_scatter_selectionY, features$feature)])
    if (input$expr_scatter_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- t(DS@data[c(gene_x, gene_y),]) %>% as.matrix() %>% as.data.frame() %>% mutate(selection = DS@meta.data[,input$expr_scatter_selection])
    plot <- data %>% ggplot(aes(
      x = data[,gene_x], y = data[,gene_y], color = selection
    )) +
      geom_point() +
      geom_smooth(method = lm, color = "black") +
      stat_cor(method = "pearson", label.x = min(data[, gene_x]), label.y = max(data[, gene_y]),
               size = input$expr_scatter_textsize/3, color = "black") +
      labs(x = input$expr_scatter_selectionX, y = input$expr_scatter_selectionY, color = input$expr_scatter_selection) +
      theme(
        plot.title = element_text(size = input$expr_scatter_legendsize),
        axis.title = element_text(size = input$expr_scatter_textsize),
        axis.text = element_text(size = input$expr_scatter_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$expr_scatter_legendsize),
        legend.title = element_text(size = input$expr_scatter_legendsize)
      ) +
      guides +
      color
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$expr_scatter_filename,
           width = input$expr_scatter_width, height = input$expr_scatter_height, dpi = input$expr_scatter_dpi)
  })# end of expression summary
  
  # expression boxplot
  output$expr_boxplot <- renderPlot({
    gene <- features$ensID[match(input$expr_boxplot_gene, features$feature)]
    data <- DS@data[gene,] %>% as.data.frame() %>% mutate(selection = DS@meta.data[,input$expr_boxplot_selection])
    plot <- data %>% ggplot(aes(
      x = selection, y = data[,1]
    )) +
      geom_boxplot(aes(fill = selection)) +
      labs(x = "", y = "", fill = input$expr_boxplot_selection, title = input$expr_boxplot_gene) +
      theme(
        plot.title = element_text(size = input$expr_boxplot_legendsize),
        axis.title = element_text(size = input$expr_boxplot_textsize),
        axis.text = element_text(size = input$expr_boxplot_textsize),
        legend.text = element_text(size = input$expr_boxplot_legendsize),
        legend.title = element_text(size = input$expr_boxplot_legendsize)
      )
    print(plot)
  }) # end of expression summary
  
  # expr_dr_facet
  output$expr_dr_facet <- renderPlot({
    gene <- features$ensID[match(input$expr_dr_facet_gene, features$feature)]
    if (input$expr_dr_facet_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- DS@dr[[input$expr_dr_facet_algorithm]]@cell.embeddings[,] %>% as.data.frame() %>% 
      mutate(gene = DS@data[gene,], selection = DS@meta.data$dataset)
    plot <- data[,] %>% ggplot(aes(
      x = data[,colnames(data)[input$expr_dr_facet_dim1]], y = data[,colnames(data)[input$expr_dr_facet_dim2]], color = gene)) +
      geom_point(size = input$expr_dr_facet_size, alpha = input$expr_dr_facet_transparency) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2") +
      labs(color = input$expr_dr_facet_gene, x = colnames(data)[input$expr_dr_facet_dim1], 
           y = colnames(data)[input$expr_dr_facet_dim2]) +
      theme(
        plot.title = element_text(size = input$expr_dr_facet_legendsize),
        axis.title = element_text(size = input$expr_dr_facet_textsize),
        axis.text = element_text(size = input$expr_dr_facet_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$expr_dr_facet_legendsize),
        legend.title = element_text(size = input$expr_dr_facet_legendsize)
      ) +
      guides(color = guide_colorbar(barwidth = input$expr_dr_facet_legendsize/15, 
                                    barheight = input$expr_dr_facet_legendsize*0.75)) +
      facet_wrap(~ selection)
    print(plot) 
  })
  # download
  observeEvent(input$expr_dr_facet_download, {
    gene <- features$ensID[match(input$expr_dr_facet_gene, features$feature)]
    if (input$expr_dr_facet_axes == "TRUE") {umap_axes <- element_line()} else {umap_axes <- element_blank()}
    data <- DS@dr[[input$expr_dr_facet_algorithm]]@cell.embeddings[,] %>% as.data.frame() %>% 
      mutate(gene = DS@data[gene,], selection = DS@meta.data$dataset)
    plot <- data[,] %>% ggplot(aes(
      x = data[,colnames(data)[input$expr_dr_facet_dim1]], y = data[,colnames(data)[input$expr_dr_facet_dim2]], color = gene)) +
      geom_point(size = input$expr_dr_facet_size, alpha = input$expr_dr_facet_transparency) +
      scale_color_gradient(low = "grey90", high = "dodgerblue2") +
      labs(color = input$expr_dr_facet_gene, x = colnames(data)[input$expr_dr_facet_dim1], 
           y = colnames(data)[input$expr_dr_facet_dim2]) +
      theme(
        plot.title = element_text(size = input$expr_dr_facet_legendsize),
        axis.title = element_text(size = input$expr_dr_facet_textsize),
        axis.text = element_text(size = input$expr_dr_facet_textsize),
        axis.line = umap_axes,
        axis.ticks = umap_axes,
        legend.text = element_text(size = input$expr_dr_facet_legendsize),
        legend.title = element_text(size = input$expr_dr_facet_legendsize)
      ) +
      guides(color = guide_colorbar(barwidth = input$expr_dr_facet_legendsize/15, 
                                    barheight = input$expr_dr_facet_legendsize*0.75)) +
      facet_wrap(~ selection)
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$expr_dr_facet_filename,
           width = input$expr_dr_facet_width, height = input$expr_dr_facet_height, dpi = input$expr_dr_facet_dpi)
  })# end of expr_dr_facet
  
  # Tab 5 - Label --------------------
  
  # Composition barplot
  output$comp_barplot <- renderPlot({
    data <- table(as.character(DS@meta.data$label), as.character(DS@meta.data$dataset)) %>% as.matrix() %>% as.data.frame() %>% as_tibble()
    colnames(data) <- c("Label", "Dataset", "Frequency") 
    data <- data %>% group_by(Dataset) %>% mutate(Proportion = round(Frequency / sum(Frequency) * 100, digits = 2)) %>% ungroup(Dataset)
    
    colors <- DS@meta.data %>% select("label", "labelColor") %>% unique() %>% remove_rownames() %>% arrange(label)
    
    plot <- data %>% ggplot(aes(x = Dataset, y = Proportion, fill = Label, label = Proportion)) +
      geom_col(position = position_stack(vjust = 0.5)) + 
      geom_text(position = position_stack(vjust = 0.5), check_overlap = TRUE, size = input$comp_barplot_size) +
      scale_fill_manual(values = col2hex(colors$labelColor)) +
      theme(
        axis.title = element_text(size = input$comp_barplot_textsize),
        axis.text = element_text(size = input$comp_barplot_textsize),
        legend.text = element_text(size = input$comp_barplot_textsize),
        legend.title = element_text(size = input$comp_barplot_textsize)
      )
    print(plot)
  }) 
  observeEvent(input$comp_barplot_download, {
    data <- table(as.character(DS@meta.data$label), as.character(DS@meta.data$dataset)) %>% as.matrix() %>% as.data.frame() %>% as_tibble()
    colnames(data) <- c("Label", "Dataset", "Frequency") 
    data <- data %>% group_by(Dataset) %>% mutate(Proportion = round(Frequency / sum(Frequency) * 100, digits = 2)) %>% ungroup(Dataset)
    
    colors <- DS@meta.data %>% select("label", "labelColor") %>% unique() %>% remove_rownames() %>% arrange(label)
    
    plot <- data %>% ggplot(aes(x = Dataset, y = Proportion, fill = Label, label = Proportion)) +
      geom_col(position = position_stack(vjust = 0.5)) + 
      geom_text(position = position_stack(vjust = 0.5), check_overlap = TRUE, size = input$comp_barplot_size) +
      scale_fill_manual(values = col2hex(colors$labelColor)) +
      theme(
        axis.title = element_text(size = input$comp_barplot_textsize),
        axis.text = element_text(size = input$comp_barplot_textsize),
        legend.text = element_text(size = input$comp_barplot_textsize),
        legend.title = element_text(size = input$comp_barplot_textsize)
      )
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$comp_barplot_filename)
  })# end of comp_barplot
  
  # Sankey plot
  output$sankey <- renderPlot({
    data <- table(as.character(DS@meta.data$label), as.character(DS@meta.data$dataset)) %>% as.matrix() %>% as.data.frame() %>% as_tibble()
    colnames(data) <- c("Label", "Dataset", "Frequency") 
    data <- data %>% group_by(Dataset) %>% mutate(Proportion = round(Frequency / sum(Frequency) * 100, digits = 2)) %>% ungroup(Dataset)
    
    colors <- DS@meta.data %>% select("label", "labelColor") %>% unique() %>% remove_rownames() %>% arrange(label)
    
    strata_colors <- rev(rainbow(length(c(levels(data$Dataset), length(data$Label)))))
    alluvium_colors <- rainbow(length(DS@meta.data$label))
    
    colors <- col2hex(colors$labelColor)
    alluvium_colors <- rep(colors, 12)
    stratum_colors <- c(rev(rainbow(length(levels(data$Dataset)), alpha = 0.5)), rev(colors))
    
    plot <- ggplot(data, aes(y = Proportion, axis1 = Dataset, axis2 = Label)) +
      geom_alluvium(fill = alluvium_colors, width = 1/12) +
      geom_stratum(fill = stratum_colors, width = 1/12) +
      geom_label(stat="stratum", label.strata=TRUE) +
      theme(
        axis.title = element_text(size = input$sankey_textsize),
        axis.text = element_text(size = input$sankey_textsize),
        legend.text = element_text(size = input$sankey_textsize),
        legend.title = element_text(size = input$sankey_textsize)
      )
    print(plot)
    }) 
  observeEvent(input$sankey_download, {
    data <- table(as.character(DS@meta.data$label), as.character(DS@meta.data$dataset)) %>% as.matrix() %>% as.data.frame() %>% as_tibble()
    colnames(data) <- c("Label", "Dataset", "Frequency") 
    data <- data %>% group_by(Dataset) %>% mutate(Proportion = round(Frequency / sum(Frequency) * 100, digits = 2)) %>% ungroup(Dataset)
    
    colors <- DS@meta.data %>% select("label", "labelColor") %>% unique() %>% remove_rownames() %>% arrange(label)
    
    strata_colors <- rev(rainbow(length(c(levels(data$Dataset), length(data$Label)))))
    alluvium_colors <- rainbow(length(DS@meta.data$label))
    
    colors <- col2hex(colors$labelColor)
    alluvium_colors <- rep(colors, 12)
    stratum_colors <- c(rev(rainbow(length(levels(data$Dataset)), alpha = 0.5)), rev(colors))
    
    plot <- ggplot(data, aes(y = Proportion, axis1 = Dataset, axis2 = Label)) +
      geom_alluvium(fill = alluvium_colors, width = 1/12) +
      geom_stratum(fill = stratum_colors, width = 1/12) +
      geom_label(stat="stratum", label.strata=TRUE) +
      theme(
        axis.title = element_text(size = input$sankey_textsize),
        axis.text = element_text(size = input$sankey_textsize),
        legend.text = element_text(size = input$sankey_textsize),
        legend.title = element_text(size = input$sankey_textsize)
      )
    ggsave(plot = plot, path = paste0("/home/", Sys.getenv("USER"), "/"), filename = input$sankey_filename)
  })# end of sankey
  
} # end of server instructions

#####################
##### execution #####
shinyApp(ui = ui, server = server)



