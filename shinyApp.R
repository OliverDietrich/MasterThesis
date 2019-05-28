############################ Shiny app programming ############################
###################### By Oliver Dietrich on 28 May 2019 ######################

library(shiny)

library(seurat)
library(tidyverse)

# User Interface
ui <- fluidPage(
  # Input functions
  sliderInput(inputId = "num", label = "Choose the binwidth",
              value = 5, min = 1, max = 100),
  selectInput(inputId = "categoryX", label = "Choose a category for the x-axis", choices = colnames(murders)),
  selectInput(inputId = "categoryY", label = "Choose a category for the y-axis", choices = colnames(murders)),
  selectInput(inputId = "color", label = "Choose a category for the color", choices = colnames(murders)),
  plotOutput("scatter")
  # Output functions
)

# Server Instructions
server <- function(input, output) {
  
  readRDS("/home/od/Data/Cindrilla/analysis/G10_2019-03-02/G10_epithelial_2019-02-19/Rdatasets/G10_epithelial_labelled.Rds")
  
  output$scatter <- renderPlot({
    ggplot(murders, aes(
      x = murders[,input$categoryX], 
      y = murders[,input$categoryY],
      col = murders[,input$color]
      )) + 
      geom_point(size = input$num)
  })
}

# execution
shinyApp(ui = ui, server = server)
