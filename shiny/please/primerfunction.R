library(tidyverse)
library(BiocManager)
library(ShortRead)
library(Biostrings)
library(stringr)



blah<-c("ACTGCTGACACATCGCATA")
as.list(blah)
as.DNAbin(blah)



# my reverse primer function 
reverseprimer<-function(x){
  x<-str_squish(x)
  # make sure the entered primer as a character set
  cdna<-as.character(x)
  # then convert that set to a DNA sting to be used 
  dna<-DNAString(cdna)
  # then get the reverse compliment 
  reverseComplement(dna)
}

reverseprimer(blah)
dna<-DNAString(blah)

reverseComplement(dna)


reverseprimer(blah)




# Now lets try to figure out how to find the patter Bitch 

# First we have to change everything that is input to DNA strings (the chunk and the primer) 
# then we just feed it to match pattern which will give us the start and stop date 
    # for the forward
primer<-c("ACTGTCT")
dp<-primer %>% DNAString()
ds<-c("ACTGTTCGTTGACGCTACTGTCTTTTTGTCCCATATCTAGTCTAGTCATGTCTGTCAT")
stringd<-ds %>% DNAString()
matchedpatt <- matchPattern(dp,stringd)

st <- matchedpatt@ranges@start + length(dp)
en <- st + matchedpatt@ranges@width

# we have to get rid of the reverse first 
string_norev <- Biostrings::mask(stringd,en,length(stringd))
string_norev
# then we can get rid of the forward then convert it to character and it will get spit back out 
string_noprimers <- string_norev[[1]] %>% 
Biostrings::mask(1,length(dp))
string_noprimers@subject %>% as.character

# for the reverse we have to get the compliment of the reverse primer because that is how it lives in the DNA strand
rprime<-reverseprimer(primer)
searchb<-complement(rprime)
matchPattern(searchb,stringd)
rprime

# now we have the start and stop positions for the primers in the DNA chunk, so we have to do something with them 
# we need the stop position of the forward primer and the start of the reverse 
b<-matchPattern(dp,stringd) 
b

searchb

chunkasDNA<-function(x){
  ci<-x
  ci<-c("ci") 
  cunkstr<-ci %>% DNAString()
}

CompRPrime<- function(x){
  rprim<-x
  rprim1<-str_squish(rprim)
  rprim2<-as.character(rprim1)
  rprim3<-DNAString(rprim2)
  rprim4<-complement(rprim3)
}

forward<-function(x){
  fprim<-x
  fprim1<-str_squish(fprim)
  fprim2<-as.character(fprim1)
  fprim3<-DNAString(fprim2)
}
dp
ampseq<-function(ci,ri,fi){
  # this is the chunk function 
  s<-ci
  ci<-c("s") 
  cunkstr<-ci %>% DNAString()
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
    Biostrings::mask(1,length(dp))
  return(string_noprimers)
}
ampseq("ACTGTTCGTTGACGCTACTGTCTTTTTGTCCCATATCTAGTCTAGTCATGTCTGTCAT",searchb,dp) 

#
#
#
#
#

# Store the inputs 
fprime<-c("AAA")
rprime<-c("GGG")
chunk<-c("TAAAGTGTGCTCTATGTAGTCATGCCCGAG")

# convert the inputs to DNA 
Dfprime<-fprime %>% DNAString()
Drprime<-rprime %>% DNAString()
Dchunk<-chunk %>% DNAString()

#get the complement of the reverse primer 
CDrprime<-Drprime %>% complement()

# get the start and stop positions for the forward and reverse primers 
revpos<-matchPattern(CDrprime, Dchunk)
forpos<-matchPattern(Dfprime, Dchunk)
# now we have to extract the numbers for the positions we want 
  # we need to stop position of the forward primer (forpos)
    # to do this we need to get at the start position and then add the length of the primer 
start_p<-forpos@ranges@start+length(Dfprime)
  # we need the start ofthe reverse primer 
stop<-revpos@ranges@start

# now we need to mask the regions of DNA that we do not want using the start and stop positions we identified 
string_no_reverse_primer<-mask(Dchunk,stop,length(Dchunk))
amplified_sequence<-string_no_reverse_primer[[1]] %>% mask(start=1,end=start_p-1)
as.character(amplified_sequence)

#
#
#
#
#

daddy<-function(b,y,f){
  # enter inputs 
  fprime<-c("b")
  rprime<-c("y")
  chunk<-c("f")
  
  # changed them to DNA 
  Dfprime<-fprime %>% DNAString()
  Drprime<-rprime %>% DNAString()
  Dchunk<-chunk %>% DNAString()
  
  # get the complement of the reverse primer 
  CDrprime<-Drprime %>% complement()
  
  # get the start and stop positions for the forward and reverse primers 
  revpos<-matchPattern(CDrprime, Dchunk)
  forpos<-matchPattern(Dfprime, Dchunk)
  
  # now we have to extract the numbers for the positions we want 
  # we need to stop position of the forward primer (forpos)
  # to do this we need to get at the start position and then add the length of the primer 
  start_p<-forpos@ranges@start+length(Dfprime)
  # we need the start ofthe reverse primer 
  stop<-revpos@ranges@start
  
  # now we need to mask the regions of DNA that we do not want using the start and stop positions we identified 
  string_no_reverse_primer<-mask(Dchunk,stop,length(Dchunk))
  amplified_sequence<-string_no_reverse_primer[[1]] %>% mask(start=1,end=start_p-1)
  character_out<-as.character(amplified_sequence)
  return(character_out)
}

daddy("AAA","GGG","TAAAGTGTGCTCTATGTAGTCATGCCCGAG")
