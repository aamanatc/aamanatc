
###Project C: Geography and Evolutionary Diversification

#I would encourage you to organize your file and divide each section by adding 3 or more hashtags at the beginning and the end of you title.

####Introduction#### 

###The Danio genus, a group of freshwater fish that originate from South and Southeast Asia, like India has captivated scientists, fish enthusiasts and geneticists alike due, to its diversity and ability to thrive in different aquatic environments(1). These small charming fish, including the known zebrafish (Danio rerio) have not only become popular additions to our aquariums but have also played a significant role in scientific research. One of the tools used to unravel the evolutionary mysteries surrounding Danio species is the Cytochrome c Oxidase subunit I (COI) gene(2). COI, which is found in mitochondria is highly regarded for its function as both a powerhouse and a molecular barcode(3). This makes it an invaluable resource for studying phylogenetics and evolutionary biology(3). With this gene, as our compass we embark on a journey to comprehend the geographical factors that have shaped the evolutionary fate of Danio species. Through this exploration we hope to shed light on their tapestry of life.

###Our main objective, for this project is to investigate how genetic data and geographic factors interact in the diversification of Danio species. Specifically, we plan to use the COI gene as a tool to create trees that reveal the evolutionary relationships among different Danio species. These phylogenetic trees will provide insights into the history of Danio. By combining data with information, we aim to answer important questions about how geography influences diversification and whether closely related species tend to inhabit the same geographic areas or if their evolution has been driven by large scale separation. Through this study we hope to uncover the connections, within the Danio genus and contribute to our understanding of speciation and evolutionary diversification on a scale.

#How does it sound to you if I suggest you to set up your directory just after your introduction?
##### Set up directory ####
#set up specific directory in porder to save my data into.
setwd("/Users/alireza/Desktop/Bioinformatic/Assignments /#2/Assignment 2")
#investigating current work directory path 
getwd()

##### Loading necessary packages ####

#Hey there! I just wanted to share a suggestion that might make your code a bit more accurate and organized. In my opinion, it would be helpful to allocate one section of your code for loading and adding necessary packages before diving into your main code. What do you think? Therefore, in the following lines, I will replace the package code from the main part to this beginning part:

#Install the rentrez package
install.packages("rentrez")
#Install the seqinr package
install.packages("seqinr")
#Install the Biostrings package
install.packages("Biostrings")
# Install the tidyverse package (which includes multiple packages)
install.packages("tidyverse")
#loading all acquired packages via library()
library("rentrez")
library("seqinr")
library("Biostrings")
library("tidyverse")
#loading dplyr package 
library(dplyr)
#installing required packages 
install.packages("ape")
#loading ape() package 
library(ape)
#installing required packages for DECIPHER
install.packages("RSQLite") 
#loading the package 
library(RSQLite)
#install Biocmanager to do tasks related to sequences 
install.packages("BiocManager")
#library the package 
library(BiocManager)
#Then, install any needed packages. for alignment and clustering 
BiocManager::install(c("Biostrings", "muscle", "msa", "DECIPHER"))
#Load libraries
library(Biostrings)
library(muscle)
library(DECIPHER)
#install required package 
install.packages("phangorn")
#load the library
library(phangorn)
#install required package
install.packages("ggplot2")
#load the library 
library("ggplot2")
#install required package
BiocManager::install("ggtree")
#circular phylogram 
library(ggtree)
#install packages 
install.packages("phytools")
#load the package
library(phytools)
#install the package
install.packages("maps")
#load the package 
library(maps)
#loading by library 
library(mapdata)

#Also, you can merge all the aforementioned packages installation lines into one line of code:
install.packages(c("rentrez", "seqinr", "Biostrings", "tidyverse", "ape", "RSQLite", "BiocManager", "phangorn", "ggplot2", "ggtree", "phytools", "maps"))

##### Main Codes ####
#performing the search in pubmed by searching the phylum 
search.res.danio <- entrez_search(db = "pubmed", term = "Danio")
#now lets see how the classification of teh result (the result)
class(search.res.danio)
#see the hits resulted 
search.res.danio
#see the sample id relating to our search 
search.res.danio$ids
#classifying the variables which is PMID from pubmed
class(search.res.danio$ids)
#let see how many results we got from our search
length(search.res.danio$ids)
#as there are just 20 which is because of default retmax rate so we need to change the retmax since it limit the number of records returned by by the search 
#using to see what are the availbale contents to be searched under "nuccore" to identify the abbreviation using for search.

#I just wanted to suggest you, maybe it would be good idea if you explain why you used pubmed as the data base, because I myself, unfortunately, cannot understand why you start with pubmed, and then you continued with "nuccore"

entrez_db_searchable("nuccore")
#performing the search in pubmed based on nucleotide data base and using phoronis in terms of genus searching 
danio_search <- entrez_search(db="nuccore", term = "Danio[ORGN]")
#see the quantity of hits 
danio_search

#as COI gene sequences are mostly between 400 to 700, the search needs to be narrow down by the sequence length as well. Also, I put retmax number 100 as I need more than 20 data 
danio_search.COI <- entrez_search(db="nuccore", term = "Danio[ORGN] AND COI[gene] AND 400:700[SLEN]", retmax=100)
#see the quantity of the results 
length(danio_search.COI$ids)
#review the summary of the result which is assigned to COI_summary 
COI_summary <- entrez_summary(db = "nuccore", id = danio_search.COI$ids)
COI_summary
#see the classification 
class(COI_summary)

#Helper function called extract_from_esummary takes specific elements from each of your list elements.(take a quick look at organisms file)
extract_from_esummary(COI_summary,"organism")

#Now we have to use entrez_fetch() function from "rentrez" package to solicit the modified data from NCBI
COI_fetch <- entrez_fetch(db = "nuccore", id = danio_search.COI$ids, rettype = "fasta")
#see what is the classification of acquired file from NCBI
class(COI_fetch)
#quick data inspection to get a sense of what the data looks like without displaying the entire dataset.
head(COI_fetch)

#write my file in specific working directory as set up earlier, and keep copy of that!
write(COI_fetch, "COI_fetch.fasta", sep = "\n")
#After checking the length of the data and making sure that we haven't downloaded WGS by mistake. we need to read data!
COI.stringset <- readDNAStringSet("COI_fetch.fasta")
#now lets see what is our file classification and quick overview on 
class(COI.stringset)
#head() to display few row from columns that is provided by command names()
head(names(COI.stringset))
#now we need to make a data frame for further analysis!
dfCOI <- data.frame(COI_title = names(COI.stringset), COI_sequence = paste(COI.stringset))
#lets see how it looks like!
view(dfCOI)

#the first column from our data frame indicate that the names are not very clear and nice! So, I decided to add another column to my dfCOI, named Species_name 
#to achieve this, we can use dplyr package to run pipe line

#Use the pipe operator (%) to apply a sequence of operations
dfCOI <- dfCOI %>%
# First, use the mutate function to create a new column "Species_Name"
# This new column extracts words 2 to 3 from the "COI_Title" column
  mutate(species_name = word(COI_title, 2L, 3L)) %>%
# Finally, use the select function to rearrange the columns as specified
  select(COI_title, species_name, COI_sequence) %>%
# view the ultimate data frame 
  view()
# I would like to suggest you to use a function for the five previous lines, as this process will be iterate in part c, for the data from NCBI:
processDataFrame <- function(input_df) {
  result_df <- input_df %>%
    # First, use the mutate function to create a new column "Species_Name"
    # This new column extracts words 2 to 3 from the "COI_Title" column
    mutate(species_name = word(COI_title, 2L, 3L)) %>%
    # Finally, use the select function to rearrange the columns as specified
    select(COI_title, species_name, COI_sequence)
  
  return(result_df)
}
dfCOI <- processDataFrame(dfCOI)


unique(dfCOI$species_name)
#Afterward, the data needs to be cleared up regarding alignment; so some packages needs to be installed, followed by loading them.

#lets take a looka at to the summmary of the sequence lengths
summary(nchar(dfCOI$COI_sequence))
#using histogram to show the frequency of their length 
x <- nchar(dfCOI$COI_sequence)
hist(x)
#the frequency shows good result as most of them are more than 650bp
#we need to also put accession ID as a column as it is likely for species_name to be repeated thorugh the data.
dfCOI <- dfCOI %>%
## This command extracts first words from the "COI_Title" column then name the column as asccession.id
  mutate(accession.id = word(COI_title, 1L)) %>%
  #selecting the useful columns 
  select(accession.id, COI_title, species_name, COI_sequence) %>%
  #view the data frame
  view()
#as there are more than one sequences for each species so we one random sample from them. Indeed, this step will be helpful when I want to merge 2 data frames as there is no lon and lat data from NCBI. 
dfCOI <- dfCOI%>%
#grouping by the species_name column 
  group_by(species_name)%>%
#a random sample for each 
  sample_n(1)%>%
  view()

#see the classification 
class(dfCOI)
#see the classification of sequences 
class(dfCOI$COI_sequence)
dfCOI <- as.data.frame(dfCOI)
#we need to change it to the format to do DNA sequences task 
dfCOI$COI_sequence <- DNAStringSet(dfCOI$COI_sequence)
#see the classification again!
class(dfCOI$COI_sequence)
#lets put species_name as a individual column otherwise we get our result based on accession.id 
names(dfCOI$COI_sequence) <- dfCOI$species_name
#see the columns 
names(dfCOI$COI_sequence)
#take a look at sequences in R
dfCOI$COI_sequence
#using BrowsSeq to look at them in HTML format which is easier to read through 
BrowseSeqs(dfCOI$COI_sequence)
#double-check whether ther are any N/A data
sum(is.na(dfCOI$COI_sequence))

#### ALIGNMENT 
#alignment needs to be done by muscle and the gap open set up at -300 at first try 
dfCOI.alignment <- DNAStringSet(muscle::muscle(dfCOI$COI_sequence, gapopen= -300), use.names=TRUE)
BrowseSeqs(dfCOI.alignment)
#now run the alignment again with default -600 as there are not many gaps 
dfCOI.alignment <- DNAStringSet(muscle::muscle(dfCOI$COI_sequence, gapopen= -600), use.names=TRUE)
#lets see the results in HTML version 
BrowseSeqs(dfCOI.alignment)
#lets take a look at first sequence length 
length(dfCOI.alignment[[1]])

#let see the quantity of gaps exist in our alignment
mean(unlist(lapply(as.character(dfCOI.alignment), str_count, "-")))
#also getting a summary 
summary(unlist(lapply(as.character(dfCOI.alignment), str_count, "-")))

#now it looks like ultimately we reached reasonable rates for the gaps by making comparison of various gapopen rates between -300to -600. 
#now we look for the frequency of the gaps based on the sequence length as well to give us a insight how is our result.
gap.Freq <- unlist(lapply(as.character(dfCOI.alignment), str_count, "-"))
hist(gap.Freq)
#then we can save the alignment file as we might use it later in other softwares just in case.
writeXStringSet(dfCOI.alignment, file = "dfCOI.alignment.fasta")

###clustering and phylogenetic analysis

#our alignment needs to be set up as DNAbin to do clustering and phylogeny 
dna.bin.COI.alignment <- as.DNAbin(dfCOI.alignment)
#see the class 
class(dfCOI.alignment)
#for clustering 3 percent is assigned 
threshhold <- 0.03
#TN93 model applied for analysing the distance between the sequences 
taken.model <- "TN93"
#The given R code calculates a genetic distance matrix from a DNA sequence alignment. It uses a specified modelto determine genetic distances and handles missing data with pairwise deletion. The resulting distance matrix, stored in distanceMatrix, is useful for tasks such as phylogenetic tree construction and genetic diversity assessment in DNA sequence analysis.
distanceMatrix <- dist.dna(dna.bin.COI.alignment, model = taken.model, as.matrix = TRUE, pairwise.deletion = TRUE)
#view first rows and columns data about the distance rate 
head(distanceMatrix)
#This R code uses the DECIPHER package to perform hierarchical clustering on a genetic distance matrix (distanceMatrix). It employs the "single" linkage method to group sequences into clusters based on their genetic distances. The threshold variable determines when to stop forming clusters. If showPlot is set to TRUE, it displays a dendrogram plot illustrating the clustering results. The outcome is stored in the dfCOI.cluster variable, providing insights into the relationships among the DNA sequences.
dfCOI.cluster <- DECIPHER::TreeLine(myDistMatrix = distanceMatrix,
                                    method = "single",
                                    cutoff = threshhold,
                                    showPlot = TRUE,
                                    type = "dendrogram",
                                  verbose = TRUE)

#figure 2: circular phylogram
#make a matrix and assign it matrix.COI.alignment
matrix.COI.alignment <- as.matrix(distanceMatrix)
#making a tree with neighbor joining method 
tree.ape <- nj(matrix.COI.alignment)
#see the classification 
class(dna.bin.COI.alignment)
#review the tip labels and internals nodes of the tree
tree.ape

#The provided R code using the ggtree package creates a circular phylogenetic tree plot from the tree.ape object. It suppresses the legend with theme_tree2(legend.position = "none"), and labels the tree tips with geom_tiplab(align = TRUE, size = 2). The resulting tree_plot variable stores the tree plot with these customization for visualization and analysis of phylogenetic relationships in a circular layout.
tree_plot <- ggtree(tree.ape, layout = "circular") + 
  theme_tree2(legend.position = "none") +
  geom_tiplab(align = TRUE, size = 2)
# view the tree 
tree_plot

#Next step would be merging data frames from BOLD and NCBI as NCBI doesn't provide data related to GPS(longitude & Latitude )
#Organizing the data mined by Bold as it is including GPS data 
#getting data from Bold via API 
dfpho <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Danio&format=tsv")
#saving the file 
write_tsv(dfpho,"Phoro.Bold")
#see the classification 
class(dfpho)
#viewing the data 
view(dfpho)
#see the including columns 
names(dfpho)
#Using pipe line to organizing data from bold 
dfpho.COI.geo <- dfpho %>%
  # Filter rows with marker code "COI-5P"
  filter(markercode == "COI-5P") %>%  
  # Filter out rows with missing lat or lon values
  filter(!is.na(lat) & !is.na(lon)) %>% 
  # Filter out rows with missing species_name
  filter(!is.na(species_name)) %>%
  #selecting columns contains of data related to species_name, latitude, and longitude 
  select(species_name, lat, lon) %>%
  #review the resulted data frame 
  view()
#In this part of you analysis, I think you should add add one more filter to exclude rows with zero latitude and longitude. In this way, this filter ensures that only raws with validate geographic coordinates are are included in further analysis:
dfpho.COI.geo <- dfpho %>%
  filter(markercode=="COI-5P" & !is.na(lat) & !is.na(lon) & lat!=0 & lon!=0) %>%
  select(species_name,lat,lon) %>%
  view()
#lets see the unique species in their column
unique(dfpho.COI.geo$species_name)
#see the sequences classification 
class(dfCOI$COI_sequence)
#as it is biostring it is needed to be character for further manipulation 
dfCOI$COI_sequence <- as.character(dfCOI$COI_sequence)

#preparing data from NCBI regarding merging to data from BOLD
#c the data from, NCBI
  dfCOI.V2 <- dfCOI %>%
  # First, use the mutate function to create a new column "Species_Name"
  # This command extracts second words   from the "COI_Title" column
  mutate(species_name = word(COI_title, 2L, 3L)) %>%
  # Finally, use the select function to rearrange the columns as specified
  select(species_name, COI_sequence) %>%
  filter(!is.na(species_name))%>%
  # view the ultimate data frame 
  view()
#I think here you can use the function:
dfCOI.V2 <- processDataFrame(dfCOI)
view(dfCOI.V2)
#lets merge these 2 data frame together.
dfmerge <- merge(dfCOI.V2, dfpho.COI.geo, by="species_name", all=F)
#view the data frame 
view(dfmerge)
#lets distinct combination of these data frames in species_name, lat, lon
dfmerge2 <- dfmerge %>%
  distinct(species_name, lat, lon) %>%
  view()
#lets reconstruct the data by removing a column and convert it into a matrix to make sure there is no duplicate 
row_names <- dfmerge2$species_name
dfmerge2$species_name <- NULL
dfmerge2 <- as.matrix(dfmerge2)
row.names(dfmerge2) <- row_names
# lets convert our tree to hclus which is proper choice for using clusters in hierarchical shape
dfCOI.cluster <- as.hclust(dfCOI.cluster)
#to use our tree in phytools we need them in phytool format 
dfCOI.cluster.phylo <- as.phylo(dfCOI.cluster)
#view the file and tips 
dfCOI.cluster.phylo

#Figure 3 Phylogeny tree on mapping plot 
## now we have to merge the data from NCBI and BOLD as there are no data related to GPS or geographical analysis

#by using keep tip we just want to use tip labels in our figure to match to the spots 
tree <- keep.tip(dfCOI.cluster.phylo, unique(rownames(dfmerge2)))
#This R code employs the phylo.to.map function to create a map visualization with a phylogenetic tree (phylogram) overlaid on it. The tree variable represents the phylogenetic tree, and dfmerge2 contains geographic data. It sets the type to "phylogram," avoids tree rotation (rotate = FALSE), and doesn't display the map immediately (plot = FALSE). The resulting map visualization is stored in the objective variable for further examination or plotting if needed.
objective<- phylo.to.map(tree, dfmerge2, type = "phylogram", rotate = F, plot = F)
#This R code uses the plot function to create a plot of the map visualization stored in the objective variable. It customizes various aspects of the plot, including the panel split, font size, font type, aspect ratio, line style, background color. 
plot(objective, split = c(0.5, 0.5), fsize = 0.75, ftype = "i", asp = 1, from.tip = F, lty = "dotted", map.bg = "purple", map.fill = "lightblue", lwd = 1, pts = F, cex.points = 1, delimit_map = T)









