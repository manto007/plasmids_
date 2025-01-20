# load packages
library(phytools)
library(ape)
library(ggplot2)
library(dplyr)

setwd("~/Downloads/Ecoli_Kallonen_Subtrees/")

# setting up tree
PhyloGroupList = c("A", "B1", "B2", "D") #REMOVE F BECAUSE IT IS PART OF D
PlasmidsList = c("FIA", "FIB","FIC")

#pdf("No_Tips_May2023_FIB", width = 20, height = 10)#for everything else 
pdf("No_Tips_May2023_B2", width = 20, height = 40) #for B2 since there are so much 

for(Plasmid in PlasmidsList) #run individually PLASMID THEN PHYLOGROUP
  #print("Plasmid:")
  Plasmid = "FIA"
print(Plasmid)
for (phylogroup in PhyloGroupList)
  phylogroup = "B1"
  phylogroup = "B2"
#print("phylogroup:")
print(phylogroup)
#Drug = "Amoxicillin"
#Drug = "Ciprofloxacin"
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
#pdf_filename = paste0("Plots_E_coli_Plasmid_",Plasmid,".pdf")
# FOR CONTREE FILE UNCOMMENT LINES BELOW
treefile = paste0("Data/cleaned_core_alignments.contree") # read tree from contree file
myTree <- ape::read.tree(treefile) # load tree

make.resdata <- function(phylogroup) { # function to make a ResData df for any phylogroup
  Resfile = paste0("Data/E_coli_phylogroup_",phylogroup,".tsv") 
  ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE)
  aacColumns = which(substr(names(ResData),1,3)=="aac")
  ResData <- read.csv(Resfile, header = TRUE, sep="\t", stringsAsFactors = TRUE, na.strings=c("","NA"))[,c(1:72,aacColumns)]
  names(ResData)[1]<-"Label" # rename first column to "Label"
  return(ResData)
}

A_ResData <- read.csv("results_A.csv") #make.resdata("A") # call function to create phylogroup A ResData
A_ResData[A_ResData==''] <- "No FIA"
B1_ResData <- make.resdata("B1") # repeat for B1
B2_ResData <- make.resdata("B2") # repeat for B2
D_ResData <- make.resdata("D") # repeat for D
#F_ResData <- make.resdata("F") # repeat for F


make.subtree <- function(df, phylogroup, myTree) { # function that filters a subtree for any phylogroup
  err_tips <- df$Accession #isolate label vector to get ERR accessions
  err_tips <- paste(err_tips) #convert elements in the vector to strings
  sub_tree <- keep.tip(myTree, err_tips) #only keep the tips
  #if (phylogroup == "D"){sub_tree <- root(sub_tree,  outgroup='ERR434558', resolve.root = TRUE)}
  #else {sub_tree <- root(sub_tree,  outgroup='ERR435048', resolve.root = TRUE)}
  sub_tree <- multi2di(sub_tree, tol=1e-10) # needed to avoid Error in .check.tree(tree, ...) : 'tree must be bifurcating (no polytomies or unbranched nodes)'
  return(sub_tree)
}

A_tree <- make.subtree(A_ResData,"A", myTree)
B1_tree <- make.subtree(B1_ResData,"B1", myTree)
B2_tree <- make.subtree(B2_ResData,"B2", myTree)
D_tree <- make.subtree(D_ResData,"D", myTree)
#F_tree <- make.subtree(F_ResData,"F", myTree)

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
        if (ResData[r, Plasmid] == "Not FIA") ResistanceArray <- c(ResistanceArray, 1) ##NOTE: for now missing data coded as susceptible. Should be removed. 
        #if (ResData[r, Plasmids] == "I") ResistanceArray <- c(ResistanceArray, 3) ##NOTE: intermediate is coded as intermediate here!
      }
    }
  }
  y<-setNames(ResistanceArray,myTree$tip.label)
  return(y)
}

A_y <- make.y.plasmid(A_ResData, A_tree, Plasmid)
B1_y <- make.y.plasmid(B1_ResData, B1_tree, Plasmid)
B2_y <- make.y.plasmid(B2_ResData, B2_tree, Plasmid)
D_y <- make.y.plasmid(D_ResData, D_tree, Plasmid)
#F_y <- make.y.plasmid(F_ResData, F_tree, Plasmid)

#we are going to adjust the code 
plot.tree <- function(myTree,y) { # function for plotting trees 
  #if (sum(y==3)>=1){ # checks if there is an intermediate in the y variable, if TRUE then plot with pink
  #plotTree(myTree,type="fan",fsize=0.8,ftype="i", show.tip.label = FALSE) # uncomment if you want fans style tree
  plotTree(myTree,fsize=0.5,ftype="off", show.tip.label = FALSE)
  tips_palette <- c("white", "blue") # palette for the tip labels
  mycols<-setNames(tips_palette[1:length(unique(y))],sort(unique(y))) # match palette colors to phenotypes
  #print(mycols)
  tiplabels(pie=to.matrix(y,sort(unique(y))),piecol=mycols,cex=0.1)
  legend_palette <- c("white", "blue") # palette for the legend column 
  #my_pheno <- c("A", "P") # phenotypes for the legend column
  #names(legend_palette) <-my_pheno # assign palette colors to the respective phenotypes
  #legend_col <- legend_palette # set this to a new variable name
  #add.simmap.legend(colors=legend_col,prompt=FALSE, vertical=FALSE, shape="circle") # add the legend
  title(main=paste0(phylogroup," ", Plasmid),font.main=3,
        line=-2)
} # TIP: Set prompt=TRUE above if you want to choose where the legend goes (need to click area on plot for it to appear)
# TIP: Set prompt=TRUE above if you want to choose where the legend goes (need to click area on plot for it to appear)

#add.simmap.legend(colors=mycols,prompt=FALSE,x=0.9*par()$usr[1],
#y=-max(nodeHeights(myTree)),fsize=0.8)

#pdf("No_Tips_May2023_", width = 20, height = 10)

# make plots for each phylogroup
A_plot <- plot.tree(A_tree, A_y)
B1_plot <- plot.tree(B1_tree, B1_y)
B2_plot <- plot.tree(B2_tree, B2_y)
D_plot <- plot.tree(D_tree, D_y)
#F_plot <- plot.tree(F_tree, F_y)
}
}
#pdf("No_Tips_May2023_", width = 20, height = 10)

#png(paste("Output/No_Tips_May2023_B2_F1A", phylogroup, "_",Plasmid,".png"), width = 20, height = 10, units = "in", res = 200)
dev.off()
