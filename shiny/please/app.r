library(shiny)
library(gridlayout)
library(DT)
library(plotly)
library(tidyverse)
library(BiocManager)
library(ShortRead)
library(Biobase)
library(Biostrings)
library(stringr)

reverseprimer<-function(x){
  x<-str_squish(x)
  # make sure the entered primer as a character set
  cdna<-as.character(x)
  # then convert that set to a DNA sting to be used 
  dna<-DNAString(cdna)
  # then get the reverse compliment 
  reverseComplement(dna)
}
chunkasDNA<-function(x){
  # store the chunk as something
  ci<-x
  # concatenate the strings together 
  ci<-c("ci") 
  # turn it inot a DNAString to use later 
  cunkstr<-ci %>% DNAString()
}
CompRPrime<- function(x){
  # store the reverse primer as something 
  rprim<-x
  # get rid of any weird spaces that might be in there 
  rprim1<-str_squish(rprim)
  # make sure it is a character 
  rprim2<-as.character(rprim1)
  # turn it inot a DNA string
  rprim3<-DNAString(rprim2)
  # find the complement to use to look for a stop position in the chunk DNA string 
  rprim4<-complement(rprim3)
}
forward<-function(x){
  # store the forward primer as something 
  fprim<-x
  # get rid of any weird spaces that might be in there 
  fprim1<-str_squish(fprim)
  # make sure that it is classified as a character 
  fprim2<-as.character(fprim1)
  # makes the forward primer into a DNA string to look for a start potision in the chunk DNA string 
  fprim3<-DNAString(fprim2)
}


ampseq<-function(ci,ri,fi){
  # finds the start and stop position for 
  startpos<-matchPattern(fi,ci)
  stoppos<-matchPattern(ri,ci)
  startpos1<-startpos@ranges@start + startpos@ranges@width
  stoppos1<-stoppos@ranges@start
  st <- startpos@ranges@start + length(fi)
  en <- st + startpos@ranges@width
  string_norev <- Biostrings::mask(ci,en,length(stringd))
  string_norev <- Biostrings::mask(stringd,en,length(stringd))
  string_noprimers <- string_norev[[1]] %>% 
    Biostrings::mask(1,length(dp))
  return(string_noprimers@subject)
}


# fix this bitch 
mamamia<-function(ci,ri,fi){
  # this is the chunk function 
  s<-ci
  ci<-c("s") 
  stringd<-ci %>% DNAString()
  # this is the reverse primer function 
  rprim<-ri
  rprim1<-str_squish(rprim)
  rprim2<-as.character(rprim1)
  rprim3<-DNAString(rprim2)
  rprim4<-complement(rprim3)
  # this is the firward primer function 
  fprim<-fi
  fprim1<-str_squish(fprim)
  fprim2<-as.character(fprim1)
  fprim3<-DNAString(fprim2)
  # this is where the magic happens 
  startpos<-matchPattern(fi,ci)
  stoppos<-matchPattern(ri,ci)
  startpos1<-startpos@ranges@start + startpos@ranges@width
  stoppos1<-stoppos@ranges@start
  st <- startpos@ranges@start + length(fi)
  en <- st + startpos@ranges@width
  string_norev <- Biostrings::mask(ci,en,length(stringd))
  string_norev <- Biostrings::mask(stringd,en,length(stringd))
  string_noprimers <- string_norev[[1]] %>% 
    Biostrings::mask(1,length(fprim3))
  return(string_noprimers)
}


ui <- grid_page(
  layout = c(
    "header        header      header     ",
    "PrimerInArea1 chunkarea   revprimarea",
    "revcomparea   chunkoutput dnacarea   ",
    ". . ."
  ),
  row_sizes = c(
    "65px",
    "1.25fr",
    "1.58fr",
    "0.17fr"
  ),
  col_sizes = c(
    "250px",
    "0.59fr",
    "1.41fr"
  ),
  gap_size = "1rem",
  grid_card_text(
    area = "header",
    content = "This made me want to die! ",
    alignment = "start",
    is_title = FALSE
  ),
  grid_card(
    area = "PrimerInArea1",
    textInput(
      inputId = "pin",
      label = "get reverse compliment ",
      value = ""
    )
  ),
  grid_card(
    area = "chunkarea",
    textInput(
      inputId = "fin",
      label = "forward primer",
      value = ""
    )
  ),
  grid_card(
    area = "chunkoutput",
    textOutput(outputId = "chunkout")
  ),
  grid_card(
    area = "revprimarea",
    textInput(
      inputId = "rin",
      label = "reverse primer",
      value = ""
    )
  ),
  grid_card(
    area = "dnacarea",
    textInput(
      inputId = "chunkin",
      label = "DNA chunk",
      value = ""
    )
  ),
  grid_card(
    area = "revcomparea",
    textOutput(outputId = "pout")
  )
)

server <- function(input, output) {
  
  
  output$pout<-renderText({
    rc <- input$pin %>% reverseprimer()
    rc %>% as.character()
  })
 
  
 output$chunkout<-renderText({
   g<-input$chunkin
   y<-input$rin
   b<-input$fin
   
   fprime<-c(b)
   rprime<-c(y)
   chunk<-c(g)
   
   Dfprime<-fprime %>% DNAString()
   Drprime<-rprime %>% DNAString()
   Dchunk<-chunk %>% DNAString()
   
   CDrprime<-Drprime %>% complement()
   
   revpos<-matchPattern(CDrprime, Dchunk)
   forpos<-matchPattern(Dfprime, Dchunk)
   
   start_p<-forpos@ranges@start+length(Dfprime)
   # we need the start of the reverse primer 
   stop<-revpos@ranges@start
   # then we need the start position of the stop primer
   string_no_reverse_primer<-mask(Dchunk,stop,length(Dchunk))
   # then we can mask the parts we don't want to get the amplified sequence 
   amplified_sequence<-string_no_reverse_primer[[1]] %>% mask(start=1,end=start_p-1)
   character_out<-as.character(amplified_sequence)
   return(character_out)
 

   
 })
}

shinyApp(ui, server)
