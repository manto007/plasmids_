# load packages
install.packages("phytools")
library(phytools)
install.packages("ape")
library(ape)
library(ggplot2)

library(dplyr)

#setwd("~/Downloads/Ecoli_Kallonen_Subtrees/")

phylogroup = "D"

# setting up tree
#Drug = "Amoxicillin"
Plasmid = "FIA"
#Drug = "Cefuroxime"
#Drug = "PipTaz"
#Drug = "Trimethoprim"
#Drug = "Gentamicin"
#Drug = "PipTaz"

# FOR NWK FILES UNCOMMENT LINES BELOW
#phylogroup = "A"
#Resfile = paste0("Data/E_coli_phylogroup_",phylogroup,".tsv")
#treefile = paste0("Data/E_coli_phylogroup_",phylogroup,".nwk")
#myTree <- ape::read.tree(treefile)
#ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE)
#aacColumns = which(substr(names(ResData),1,3)=="aac")
#ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings=c("","NA"))[,c(1:42,aacColumns)]
#names(ResData)[1]<-"Label" # rename first column to "Label"

# FOR CONTREE FILE UNCOMMENT LINES BELOW
treefile = paste0("Data/cleaned_core_alignments.contree") # read tree from contree file
myTree <- ape::read.tree(treefile) # load tree

#make.resdata <- function(phylogroup) { # function to make a ResData df for any phylogroup
  #Resfile = paste0("Data/E_coli_phylogroup_",phylogroup,".tsv") 
  #Resfile = paste0("results_",phylogroup,".csv") 
  #ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE)
  #ResData <- read.csv("results_", phylogroup, ".csv")
  #aacColumns = which(substr(names(ResData),1,3)=="aac")
  #ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings=c("","NA"))[,c(1:72,aacColumns)]
  #names(ResData)[1]<-"Label" # rename first column to "Label"
  #return(ResData)
#}

make.resdata <- function(phylogroup) {
  Resfile <- paste0("Data/results_", phylogroup, ".csv")
  columns_to_read <- c("Accession", "FIA", "FIB", "FIC", "PipTaz", "AmoxiClav",	"Amoxicillin",	"Ceftazidime",	"Ciprofloxacin",	"Cefotaxime",	"Cefuroxime",	"Gentamicin")
  ResData <- read.csv(Resfile, header = TRUE, stringsAsFactors = TRUE)[, columns_to_read]
  #names(ResData)[1] <- "Label"
  return(ResData)
}

A_ResData <- make.resdata("A")#read.csv("results_A.csv") # call function to create phylogroup A ResData
#A_ResData[A_ResData=='0'] <- "No FIA"
B1_ResData <- make.resdata("B1") # repeat for B1
B2_ResData <- make.resdata("B2") # repeat for B2
D_ResData <- make.resdata("D") # repeat for D
#F_ResData <- make.resdata("F") # repeat for F

add.outgroup <- function(ResData, phylogroup) {
  D_ResData <- make.resdata("D") # make df for phylogroup D to grab an outgroup row from
  A_ResData <- make.resdata("A") # make df for phylogroup A to grab an outgroup row from
  D_row <- D_ResData[55,] # take a row from phylogroup D df to make as the outgroup for the other phylogroups
  A_row <- A_ResData[90,] # take a row from phylogroup A df to make as the outgroup for phylogroup D
  if (phylogroup == "D"){ResData <- rbind(ResData, A_row)} # if phylogroup is D add outgroup from A df
  else {ResData <- rbind(ResData, D_row)} # else for other phylogroups use the outgroup from D df
  return(ResData)
}

A_ResData <- add.outgroup(A_ResData, "A")
B1_ResData <- add.outgroup(B1_ResData, "B1")
B2_ResData <- add.outgroup(B2_ResData, "B2")
D_ResData <- add.outgroup(D_ResData, "D")
#F_ResData <- add.outgroup(F_ResData, "F")

make.subtree <- function(df, phylogroup, myTree) { # function that filters a subtree for any phylogroup
  err_tips <- df$Label #isolate label vector to get ERR accessions
  err_tips <- paste(err_tips) #convert elements in the vector to strings
  sub_tree <- keep.tip(myTree, err_tips) #only keep the tips
  #if (phylogroup == "D"){sub_tree <- root(sub_tree,  outgroup='ERR435342_abricate_results.txt', resolve.root = TRUE)}
  #else {sub_tree <- root(sub_tree,  outgroup='ERR435048', resolve.root = TRUE)}
  sub_tree <- multi2di(sub_tree, tol=1e-10) # needed to avoid Error in .check.tree(tree, ...) : 'tree must be bifurcating (no polytomies or unbranched nodes)'
  return(sub_tree)
}

make.subtree.no_outgroup <- function(ResData, phylogroup, CompleteTree) { # function that filters a subtree for any phylogroup
  err_tips <- ResData$Accession #isolate label vector to get ERR accessions
  err_tips <- paste(err_tips) #convert elements in the vector to strings
  sub_tree <- keep.tip(CompleteTree, err_tips) #only keep the tips
  #if (phylogroup == "D"){sub_tree <- root(sub_tree,  outgroup='ERR434558', resolve.root = TRUE)}
  #else {sub_tree <- root(sub_tree,  outgroup='ERR435048', resolve.root = TRUE)}
  sub_tree <- multi2di(sub_tree, tol=1e-10) # needed to avoid Error in .check.tree(tree, ...) : 'tree must be bifurcating (no polytomies or unbranched nodes)'
  return(sub_tree)
}

A_ResData <- subset(A_ResData, Label != "ERR434624")
D_ResData <- subset(D_ResData, Accession != "ERR435348")


A_tree <- make.subtree(A_ResData,"A", myTree)
B1_tree <- make.subtree(B1_ResData,"B1", myTree)
B2_tree <- make.subtree(B2_ResData,"B2", myTree)
D_tree <- make.subtree(D_ResData,"D", myTree)
F_tree <- make.subtree(F_ResData,"F", myTree)

make.y.plasmid <- function(ResData, myTree, Plasmid) { # function that creates y variable
  ResistanceArray <-c()
  for (l in myTree$tip.label){
    #find right place in ResDataDetails
    #r = which(ResData$Accession == l) - #look up in ResData and look at right line whether it is the right way 
    if (length(which(ResData$Accession == l)) == 1){
      if (is.na(ResData[r, Plasmid])) ResistanceArray <- c(ResistanceArray, 1) ###Missing data coded as susceptible 
      if (!is.na(ResData[r, Plasmid])){
        if (ResData[r, Plasmid] == "FIA") ResistanceArray <- c(ResistanceArray, 2)
        #if (ResData[r, Plasmids] == "S") ResistanceArray <- c(ResistanceArray, 1)
        if (ResData[r, Plasmid] == "No_FIA") ResistanceArray <- c(ResistanceArray, 1) ##NOTE: for now missing data coded as susceptible. Should be removed. 
        #if (ResData[r, Plasmids] == "I") ResistanceArray <- c(ResistanceArray, 3) ##NOTE: intermediate is coded as intermediate here!
      }
    }
  }
  y<-setNames(ResistanceArray,myTree$tip.label)
  return(y)
}


make.y.plasmid <- function(ResData, myTree, Plasmid) { 
  ResistanceArray <- c()
  
  for (l in myTree$tip.label) {
    # Find the right place in ResDataDetails
    r <- which(ResData$Accession == l)  # Look up in ResData and get the correct line
    
    if (length(r) == 1) {  # Ensure there is exactly one match
      if (is.na(ResData[r, Plasmid])) {
        ResistanceArray <- c(ResistanceArray, 1)  # Missing data coded as susceptible
      } else {
        if (ResData[r, Plasmid] == "FIC") {
          ResistanceArray <- c(ResistanceArray, 2)
        } else if (ResData[r, Plasmid] == "No_FIC") {
          ResistanceArray <- c(ResistanceArray, 1)  # Missing data coded as susceptible
        }
      }
    } else {
      print(paste("No matching index found for tip label:", l))
    }
    
    print(paste("Current ResistanceArray:", paste(ResistanceArray, collapse = ", ")))
  }
  
  if (length(ResistanceArray) == 0) {
    stop("ResistanceArray is empty. Check if the conditions are being met.")
  }
  
  y <- setNames(ResistanceArray, myTree$tip.label)
  return(y)
}

# Example call to the function
y_list <- make.y.plasmid(ResData, myTree, Plasmid)




A_y <- make.y.plasmid(A_ResData, A_tree, Plasmid)
B1_y <- make.y.plasmid(B1_ResData, B1_tree, Plasmid)
B2_y <- make.y.plasmid(B2_ResData, B2_tree, Plasmid)
D_y <- make.y(D_ResData, D_tree, Drug)
F_y <- make.y(F_ResData, F_tree, Drug)

#we are going to adjust the code 
plot.tree <- function(myTree,y) { # function for plotting trees 
  if (sum(y==3)>=1){ # checks if there is an intermediate in the y variable, if TRUE then plot with pink
    plotTree(myTree,type="fan",fsize=0.8,ftype="i") # uncomment if you want fans style tree
    #plotTree(myTree,fsize=0.5,ftype="off")
    tips_palette <- c("white", "red", "pink") # palette for the tip labels
    mycols<-setNames(tips_palette[1:length(unique(y))],sort(unique(y))) # match palette colors to phenotypes
    print(mycols)
    tiplabels(pie=to.matrix(y,sort(unique(y))),piecol=mycols,cex=0.1)
    legend_palette <- c("white", "pink", "red") # palette for the legend column 
    my_pheno <- c("S", "I", "R") # phenotypes for the legend column
    names(legend_palette) <-my_pheno # assign palette colors to the respective phenotypes
    legend_col <- legend_palette # set this to a new variable name
    add.simmap.legend(colors=legend_col,prompt=FALSE, vertical=FALSE, shape="circle") # add the legend
  } # TIP: Set prompt=TRUE above if you want to choose where the legend goes (need to click area on plot for it to appear)
  else { # if no intermediate plot without pink 
    plotTree(myTree,type="fan",fsize=0.8,ftype="i") # uncomment if you want fans style tree
    #plotTree(myTree,fsize=0.5,ftype="off") # plot initial tree
    tips_palette <- c("white", "red") # palette for the tip labels 
    mycols<-setNames(tips_palette[1:length(unique(y))],sort(unique(y))) # match palette colors to phenotypes
    print(mycols)
    tiplabels(pie=to.matrix(y,sort(unique(y))),piecol=mycols,cex=0.1)
    legend_palette <- c("white","red") # palette for the legend column
    my_pheno <- c("S","R") # phenotypes for the legend column
    names(legend_palette) <-my_pheno # assign palette colors to the respective phenotypes
    legend_col <- legend_palette # set this to a new variable name
    add.simmap.legend(colors=legend_col,prompt=FALSE, vertical=FALSE, shape="circle") # add the legend 
  } # TIP: Set prompt=TRUE above if you want to choose where the legend goes (need to click area on plot for it to appear)
  
  #add.simmap.legend(colors=mycols,prompt=FALSE,x=0.9*par()$usr[1],
  #y=-max(nodeHeights(myTree)),fsize=0.8)
}
# make plots for each phylogroup
A_plot <- plot.tree(A_tree, A_y)
B1_plot <- plot.tree(B1_tree, B1_y)
B2_plot <- plot.tree(B2_tree, B2_y)
D_plot <- plot.tree(D_tree, D_y)
F_plot <- plot.tree(F_tree, F_y)

# check tip colors
check_df <- F_ResData %>% 
  select(Label, Gentamicin)

# if you want to extract a certain clade in the tree:
nodelabels(font = 0.1, width = 0.1, height = 0.1) # use to check nodes, width and height parameters don't really work tho :(
ex.clade<-extract.clade(F_tree,node=158) # extract clade at node 158
plotTree(ex.clade,ftype="i") # plot it

#installing packages needed for table visualisation
install.packages("data.table")       # Install & load data.table
library("data.table")



#creating a dataframe for the mutation data > change working directory
mutation_data <- read.csv("Supplemental_Data_S1.csv")

#trying this for our data
mutation_data <- subset (mutation_data, select = -ERR.accession)
mutation_data
View(mutation_data)

#Create dataframe for mutation
FIA <- data.frame("Label" = A_ResData[, c("Label")], "FIA" = A_ResData[,c("FIA")] )
#rownames(gyrA) <- B2_ResData$ERS
#colnames(gyrA)[0] ="ERS"
#colnames(gyrA)[1] ="mutation"
View(FIA)

#create dataframe for cipro resistance
cipR <- data.frame("Ciprofloxacin" = B2_ResData[,c("Ciprofloxacin")])
rownames(cipR) <- B2_ResData$ERS

library(heatmap)

h1 <-  gheatmap(Label, FIA,                                 # we add a heatmap layer of the gender dataframe to our tree plot
                offset = 10,                               # offset shifts the heatmap to the right,
                width = 0.10,                              # width defines the width of the heatmap column,
                color = NULL,                              # color defines the boarder of the heatmap columns
                colnames = FALSE) +                               # hides column names for the heatmap
  scale_fill_manual(name = "Mutation",                       # define the coloring scheme and legend for gender
                    values = c("green", "purple", "red", "yellow"),
                    breaks = c("FIA", "Not FIA"),
                    labels = c("FIA", "Not FIA")) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.box = "vertical", legend.margin = margin())

data(A_ResData)
cols<-setNames(palette()[1:6],mapped.states(A_ResData))
plot(A_ResData,cols,type="fan",fsize=0.8,lwd=3,ftype="i")
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE,fsize=0.9)

#code to run ggtree
# Bioconductor version
library(BiocManager)
BiocManager::install("ggtree")

# github version
devtools::install_github("YuLab-SMU/ggtree")

getwd()
tryout_tree <- ape::read.tree("Shigella_tree.txt")

pacman::p_install_gh("appliedepi/epirhandbook")

plotBranchbyTrait(B2_tree, x, mode=c("edges","tips","nodes"), palette="rainbow", 
                  legend=TRUE, xlims=NULL, ...)

install.packages(ggtree)
library(ggtree)

#plotting with ggtree
ggtree(A_tree)
ggtree(A_tree) + theme_tree2()
ggtree(A_tree, branch.length="none") #this looks weird

ggtree(A_tree) + geom_text(aes(label=node), hjust=-.3)

#clean and inspect
head(A_tree$tip.label) #looks good we see the sample identifiers
colnames(A_ResData) 

#visualize the tree
p <- ggtree(B1_tree) 
p

#visualise trees
D <- ggtree(D_tree) 
D
F1 <- ggtree(F_tree) 
F1

#try this out
p <- p %<+% B1_ResData + geom_tippoint(aes(color=Ciprofloxacin))

#trait matrix
#dataframe$varname
trait_df = data.frame(
  species =  B1_tree$tip.label,
  trait1 = rnorm(B1_ResData$Ciprofloxacin)) # Example of continuous trait

#plotting
plot_tree = ggtree(A_tree, layout = "fan", right = TRUE, size = 0.1)


# Add trait information into tree
plot_tree %<+% trait_df + geom_tippoint(aes(color = trait1))#shows trait in blue color

#another way of plotting > not working
p <- p %<+% trait_df + geom_tippoint(aes(color=B1_ResData$Ciprofloxacin))

#trying out the code from epidemiology handbook 
colnames(B2_ResData)

ggtree(A_tree, layout="circular", branch.length = 'none')  
ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_ResData + 
  aes(color = I(FIA), fill = I(FIA)) +
  scale_color_manual(name = "FIA", 
                     breaks = c("FIA", "Not FIA"), 
                     labels = c("FIA", "Not FIA"), 
                     values = c("blue", "white")) +
  scale_fill_manual(name = "FIA", 
                    breaks = c("FIA", "Not FIA"), 
                    labels = c("FIA", "Not FIA"), 
                    values = c("blue", "white"))


#trees
ggtree(D_tree, layout="circular", branch.length = "none")
ggtree(F_tree, layout="circular", branch.length = "none")

bp <- ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(FIA))+                       # color the branches according to a variable in your dataframe
  #branches <- list(A=1,B=2,C=3)
#tree <- groupOTU(tree,branches)
#  scale_color_manual(
#    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
#    breaks = c("I", "R", "S"),                     # the different options in your variable
#    labels = c("Intermediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
#    values = c("pink", "red", "green"),                  # the color you want to assign to the variable 
#    na.value = "black") +                        # color NA values in black as well
#  new_scale_color()+                             # allows to add an additional color scheme for another variable
geom_tippoint(
  mapping = aes(color = gyrA),          # tip color by continent. You may change shape adding "shape = "
  size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "gyrA",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "black") +                    # for the NA values we choose the color grey
  #geom_tiplab(                             # adds name of sample to tip of its branch 
  #  color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
  #  offset = 1,
  #  size = 1,
  #  geom = "text",
  #  align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 



#for B2 > parC

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(FIA))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("FIA", "Not FIA"),                     # the different options in your variable
    labels = c("FIA", "Not FIA"),        # how you want the different options named in your legend, allows for formatting
    values = c("blue", "white"),                  # the color you want to assign to the variable 
    na.value = "grey") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "black") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup A in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

#for B2 > parC

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "gyrA",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "black") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "gyrA",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 


#A ESBL
ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "esbl",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup A in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

#B1 ESBL
ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% B1_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "esbl",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B1 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

#F ESBL
ggtree(B1_tree, layout = "circular", branch.length = 'none') %<+% B1_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "esbl",                    # name of your color scheme (will show up in the legend like this)
    #    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B1 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

#B2 esbl > PipTaz
#F ESBL
ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Cefuroxime))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Intermediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("yellow", "red", "green"),                  # the color you want to assign to the variable 
    na.value = "black") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = esbl),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    #  name = "esbl",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set3",                      # we choose a set of colors coming with the brewer package
    na.value = "grey") +                    # for the NA values we choose the color grey
  # geom_tiplab(                             # adds name of sample to tip of its branch 
  #    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
  #   offset = 1,
  #   size = 1,
  #  geom = "text",
  #    align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

library(ggtree)
library(treeio)
library(ggplot2)
library(ggtreeExtra)
library(tidytree)
library(ggnewscale)
library(ggimage)


#tryout
ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  # scale_color_manual(
  # name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
  #breaks = c("I", "R", "S"),                     # the different options in your variable
  #labels = c("Intermediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
  # values = c("yellow", "green", "blue"),                  # the color you want to assign to the variable 
  # na.value = "grey") +                        # color NA values in black as well
  
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = parC),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "SNPs in parC",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "black") +                    # for the NA values we choose the color grey
  # geom_tiplab(                             # adds name of sample to tip of its branch 
  # color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
  #  offset = 1,
  # size = 1,
  #geom = "text",
  #align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Intermediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("yellow", "red", "green"),                  # the color you want to assign to the variable 
    na.value = "black")                        # color NA values in black as well
#  new_scale_color()+                             # allows to add an additional color scheme for another variable
geom_tippoint(
  mapping = aes(color = gyrA),          # tip color by continent. You may change shape adding "shape = "
  size = 1.5)+                               # define the size of the point at the tip
  scale_color_brewer(
    name = "SNPs in parC",                    # name of your color scheme (will show up in the legend like this)
    palette = "Set1",                      # we choose a set of colors coming with the brewer package
    na.value = "black") +                    # for the NA values we choose the color grey
  # geom_tiplab(                             # adds name of sample to tip of its branch 
  # color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
  #  offset = 1,
  # size = 1,
  #geom = "text",
  #align = TRUE)+    
  ggtitle("Phylogenetic tree of phylogroup B2 in E.coli")+       # title of your graph
  theme(
    axis.title.x = element_blank(), # removes x-axis title
    axis.title.y = element_blank(), # removes y-axis title
    legend.title = element_text(    # defines font size and format of the legend title
      face = "bold",
      size = 12),   
    legend.text=element_text(       # defines font size and format of the legend text
      face = "bold",
      size = 10),  
    plot.title = element_text(      # defines font size and format of the plot title
      size = 12,
      face = "bold"),  
    legend.position = "bottom",     # defines placement of the legend
    legend.box = "vertical",        # defines placement of the legend
    legend.margin = margin()) 

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_ResData + # %<+% adds dataframe with sample data to tree
  aes(color = I(Cefuroxime))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Intermediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("yellow", "red", "green"),                  # the color you want to assign to the variable 
    na.value = "black")

