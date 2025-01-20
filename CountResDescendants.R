
#Here is my plan:
#determine which branch the switch from sus to res happens.
#determine the next node from that branch.
#use getDescendants(tree, node, curr=NULL) to get all the descendents.
#determine how many of the descencents are resistant.
library(phytools)
#setwd("~/Dropbox/Team-Phylo/")
#ClusterSizeDF<-data.frame(Drug = c(), phylogroup = c(), clustersize = c())
ClusterSizeDF<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Plasmid", "phylogroup", "clustersize", "lastresistantparent"))

# FOR CONTREE FILE UNCOMMENT LINES BELOW
#treefile = paste0("Data/cleaned_core_alignments_copy.contree") # read tree from contree file
CompleteTree <- ape::read.tree("Data/cleaned_core_alignments.contree") # load tree
PlasmidList = c("FIA")
PhyloGroupList = c("A", "B1", "B2", "D") #REMOVE F BECAUSE IT IS PART OF D

source("TreeFunctions_Plasmids.R")

for (Plasmid in PlasmidList){
  #Drug = "Ceftazidime"
  print(Plasmid)
  Plasmid="FIA"
  for (phylogroup in PhyloGroupList){
    #phylogroup = "B2"
    print("phylogroup:")
    print(phylogroup)
    
    ResData <- make.resdata(phylogroup) # call function to create phylogroup A ResData
    #ResData <- add.outgroup(ResData,phylogroup) # can it work without outgroup?
    myTree <- make.subtree.no_outgroup(ResData,phylogroup, CompleteTree)
    y_list <- make.y.plasmid(ResData, myTree, Plasmid)
    plot.tree(myTree, y_list, phylogroup, Plasmid)
    
    pi<-setNames(c(1,0),c("1","2"))
    if (length(which(y_list == 3))>0)pi<-setNames(c(1,0,0),c("1","2","3"))
    ####OK, in the B2 tree there are a few (2) branches with length 0, and I wonder if that is what is causing the problem. 
    #Let's see what happens if I make them none 0'
    myTree$edge.length[which(myTree$edge.length==0)]= min(myTree$edge.length[myTree$edge.length>0])
    
    ###LOOK AT ONE POSSIBLE HISTORY
    #An alternative approach to the one outline above is to use an MCMC approach to sample character histories from their posterior probability distribution. 
    #This is called stochastic character mapping (Huelsenbeck et al. 2003). 
    #The model is the same but in this case we get a sample of unambiguous histories for our discrete character's evolution on the tree - 
    #rather than a probability distribution for the character at nodes.
    
    #For instance, given the data simulated above - 
    #we can generate the stochastic character map as follows:
    
    ## simulate single stochastic character map using empirical Bayes method
    mymtree<-make.simmap(myTree,y_list,model="ER", pi=pi) #Added pi = pi to fix root node
    #mymtree<-make.simmap(myTree,y_list,model="ER") #Added pi = pi to fix root node
    
    #I want the tree to show the tip numbers, so I change the tip labels into the tip numbers
    mymtree$tip.label<-1:length(mymtree$tip.label)
    
    #a summary of what happened on the tree. 
    #countSimmap(mymtree, states=NULL, message=TRUE)
    
    
    #for each resistant leave, determine if its parent is resistant. 
    #If not: then singleton, clustersize = 1
    #If parent is resistant, find last resistant parent going up in tree and then determine clade size
    #Also, all descendants should be on ignore list to avoid double counting .
    
    #get the resistant leaves 
    resleaves = names(which(getStates(mymtree,"tips")=="2"))
    clustersizes = c(); listToIgnore = c(); resnodesclusters = c(); lastresistantparents = c()
    for (resleave in  resleaves){
      #print(paste("leave: ", resleave))
      parent = getParent(mymtree, resleave) 
      #find the last resistant parent, if singleton, then last resistant parent is itself
      #if parent not resistant, then singleton clade 
      if (getStates(mymtree,"nodes")[as.character(parent)] %in% c("1", "3")) { ### NOV 2022 added "3" as option () = Intermediate)
        #print("singleton clade") #getStates is a named vector
        lastresistantparents <-c(lastresistantparents, resleave)
      }
      #if parent resistant, count clade size 
      if (getStates(mymtree,"nodes")[as.character(parent)] == "2") {
        #print("not singleton clade") #if not a singleton, get the last resistant parent (going back in time)
        while (getStates(mymtree,"nodes")[as.character(parent)] == "2"){
          previousparent = parent
          parent = getParent(mymtree, parent) #I got to remember what the previous parent was though! 
          #print(paste(previousparent, parent))
        }
        lastresistantparent =   previousparent
        lastresistantparents <-c(lastresistantparents, lastresistantparent)
      }
    }
    
    clustersizes = as.data.frame(table(lastresistantparents))$Freq
    resnodesclusters = as.character(as.data.frame(table(lastresistantparents))$lastresistantparents)
    
    print("clustersizes:")
    print(clustersizes)
    for (c in 1: length(clustersizes)){
      ClusterSizeDF[nrow(ClusterSizeDF)+1,]<-list(Plasmid, phylogroup, clustersizes[c], resnodesclusters[c])
    }
  }
}

write.csv(x = ClusterSizeDF, file = "Output/ClusterSizeDF_June2024.csv", row.names = FALSE)


make.y.plasmid <- function(ResData, myTree, Plasmid) {
  # Extract unique tip labels from the phylogenetic tree
  tip_labels <- unique(myTree$tip.label)
  
  # Create the ResistanceArray based on the Plasmid column
  ResistanceArray <- ifelse(tip_labels %in% ResData$Accession,
                            ifelse(is.na(ResData$Plasmid[ResData$Accession %in% tip_labels]) |
                                     ResData$Plasmid[ResData$Accession %in% tip_labels] == "", 1,
                                   ifelse(ResData$Plasmid[ResData$Accession %in% tip_labels] == "FIA", 2, 1)),
                            NA)
  
  # Create the named vector y
  y <- setNames(ResistanceArray, tip_labels)
  
  return(y)
}

