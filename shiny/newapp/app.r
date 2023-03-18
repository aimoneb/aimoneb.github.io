library(shiny)
library(gridlayout)
library(DT)
library(plotly)
library(tidyverse)
library(BiocManager)
library(ShortRead)
library(Biostrings)

ui <- grid_page(
  layout = c(
    "area2 area2 area2",
    "area3 area3 table",
    "area3 area3 table",
    "area3 area3 table"
  ),
  row_sizes = c(
    "95px",
    "1.66fr",
    "0.34fr",
    "1fr"
  ),
  col_sizes = c(
    "250px",
    "0.59fr",
    "1.41fr"
  ),
  gap_size = "1rem",
  grid_card(
    area = "table",
    title = "Step 2: determine amplified sequence",
    scrollable = TRUE,
    item_gap = "12px",
    textInput(
      inputId = "ForwardPrimerInput",
      label = "Forward Primer",
      value = ""
    ),
    textInput(
      inputId = "ReversePrimerInput",
      label = "Reverse Primer",
      value = ""
    ),
    textInput(
      inputId = "sequence",
      label = "Sequence Input",
      value = ""
    ),
    actionButton(
      inputId = "ampButton",
      label = "GO",
      width = "25%"
    ),
    textOutput(outputId = "amplifiedSequence"),
    actionButton(
      inputId = "copyAmpSeq",
      label = "copy",
      width = "25%"
    )
  ),
  grid_card(
    area = "area3",
    title = "Step 1: determine compliment primer",
    textInput("text", "text (not isolated):", "input text"),
      inputId = "PI")
    )
    actionButton(
      inputId = "complimentButton",
      label = "GO",
      width = "25%"
  
  grid_card_text(
    content = "Determine compliment primer and get amplified sequence",
    alignment = "center",
    area = "area2"
  )



server <- function(input, output){
  function(x){
  x<-input$PI
  x<-str_squish(x)
  cdna<-as.character(x)
  dna<-DNAString(cdna)
  L<-as.character(dna)
  output$PO<-renderText(input$PI)
  
  output$PO<-renderText(input$PI,paste0(input$PI,"fghjk"))
}}

# This will be used for step 2 
  #set<-input$sequence%>% DNAString()
  #y<-input$ReversePrimerInput %>% DNAString()
  #complement(y)
  #matchPattern(set,y)
 

shinyApp(ui, server)
