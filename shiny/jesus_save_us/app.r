library(shiny)
library(gridlayout)
library(ggplot2)

ui <- navbarPage(
  title = "Chick Weights",
  selected = "Line Plots",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "My Shiny App",
    textInput(
      inputId = "myTextInput",
      label = "Text Input",
      value = ""
    )
  )
)

server <- function(input, output) {
  output$linePlots <- renderPlot({
    obs_to_include <- as.integer(ChickWeight$Chick) <= input$numChicks
    chicks <- ChickWeight[obs_to_include, ]

    ggplot(
      chicks,
      aes(
        x = Time,
        y = weight,
        group = Chick
      )
    ) +
      geom_line(alpha = 0.5) +
      ggtitle("Chick weights over time")
  })

  output$dists <- renderPlot({
    ggplot(
      ChickWeight,
      aes(x = weight)
    ) +
      facet_wrap(input$distFacet) +
      geom_density(fill = "#fa551b", color = "#ee6331") +
      ggtitle("Distribution of weights by diet")
  })
}

shinyApp(ui, server)
