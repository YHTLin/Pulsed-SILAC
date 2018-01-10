# ShinyApp accompanying pulsed_SILAC_ozlem.R

library(shiny)

ui = fluidPage(
  # Input for missing value filter: DMSO
  sliderInput("DMSO_valid", label = h3("Slider"), min = 0, max = 5, value = 5),
  # Input for missing value filter: M1071
  sliderInput("M1071_valid", label = h3("Slider"), min = 0, max = 4, value = 4),
  # Input for missing value filter: Rapa
  sliderInput("Rapa_valid", label = h3("Slider"), min = 0, max = 4, value = 4),
  
  # Outputs histogram after applying filter and normalization
  plotOutput("hist"),
  
  # Download button
  actionButton(inputId = "download", label = "Download")
  
  
)

server = function(input, output) {
  output$value <- renderPlot({ input$slider1 })
  
  # Update button
  
  # Download data table as CSV
  observeEvent(input$download, {
    #write.csv()
  })
}


shinyApp(ui = ui, server = server)