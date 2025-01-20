setwd("~/Downloads/Plasmid_Thesis")

library(ggtree)

A_ResData <- read.csv("results_A.csv")
A_ResData$FIA[A_ResData$FIA == ''] <- "No FIA"
FIA_PipTaz <- data.frame("Accession" = A_ResData[, c("Accession")], "FIA" = A_ResData[,c("FIA")], "PipTaz" = A_ResData[,c("PipTaz")] )
FIA_PipTaz$FIA[FIA_PipTaz$FIA == ''] <- "No FIA"
FIA_PipTaz$PipTaz[FIA_PipTaz$PipTaz == ''] <- "S"

# Assuming 'FIA_PipTaz' already has the 'FIA' and 'PipTaz' columns, you can use the following code:

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% FIA_PipTaz + # %<+% adds dataframe with sample data to tree
  aes(color = I(PipTaz))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "white"),                  # the color you want to assign to the variable 
    na.value = "") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  scale_color_brewer(
    name = "",                   
    palette = "Set1",     
    drop=TRUE,
    na.value = "black")+
  ggtitle("Phylogenetic tree of phylogroup A_PipTaz in E.coli")+       # title of your graph
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

A_ResData <- read.csv("results_A.csv")
FIA_AmoxiClav <- data.frame("Accession" = A_ResData[, c("Accession")], "FIA" = A_ResData[,c("FIA")], "AmoxiClav" = A_ResData[,c("AmoxiClav")] )
FIA_AmoxiClav$FIA[FIA_AmoxiClav$FIA == ''] <- "No FIA"
FIA_AmoxiClav$AmoxiClav[FIA_AmoxiClav$AmoxiClav == ''] <- "S"

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% FIA_AmoxiClav + # %<+% adds dataframe with sample data to tree
  aes(color = I(AmoxiClav))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "white"),                  # the color you want to assign to the variable 
    na.value = "") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  scale_color_brewer(
    name = "",                   
    palette = "Set1",     
    drop=TRUE,
    na.value = "black")+
  ggtitle("Phylogenetic tree of phylogroup A_AmoxiClav in E.coli")+       # title of your graph
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

A_ResData <- read.csv("results_A.csv")
FIA_Ciprofloxacin <- data.frame("Accession" = A_ResData[, c("Accession")], "FIA" = A_ResData[,c("FIA")], "Ciprofloxacin" = A_ResData[,c("Ciprofloxacin")] )
FIA_Ciprofloxacin$FIA[FIA_Ciprofloxacin$FIA == ''] <- "No FIA"
FIA_AmoxiClav$AmoxiClav[FIA_AmoxiClav$AmoxiClav == ''] <- "S"

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% FIA_Ciprofloxacin + # %<+% adds dataframe with sample data to tree
  aes(color = I(Ciprofloxacin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "white"),                  # the color you want to assign to the variable 
    na.value = "") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  scale_color_brewer(
    name = "",                   
    palette = "Set1",     
    drop=TRUE,
    na.value = "black")+
  ggtitle("Phylogenetic tree of phylogroup A_Ciprofloxacin in E.coli")+       # title of your graph
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

A_ResData <- read.csv("results_A.csv")
FIA_Gentamicin <- data.frame("Accession" = A_ResData[, c("Accession")], "FIA" = A_ResData[,c("FIA")], "Gentamicin" = A_ResData[,c("Gentamicin")] )
FIA_Gentamicin$FIA[FIA_Gentamicin$FIA == ''] <- "No FIA"
FIA_AmoxiClav$AmoxiClav[FIA_AmoxiClav$AmoxiClav == ''] <- "S"

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% FIA_Gentamicin + # %<+% adds dataframe with sample data to tree
  aes(color = I(Gentamicin))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "white"),                  # the color you want to assign to the variable 
    na.value = "") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  scale_color_brewer(
    name = "",                   
    palette = "Set1",     
    drop=TRUE,
    na.value = "black")+
  ggtitle("Phylogenetic tree of phylogroup A_Gentamicin in E.coli")+       # title of your graph
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

B1_ResData <- read.csv("results_B1.csv")
B1_PipTaz <- data.frame("Accession" = B1_ResData[, c("Accession")], "FIA" = B1_ResData[,c("FIA")], "PipTaz" = B1_ResData[,c("PipTaz")] )
B1_PipTaz$FIA[B1_PipTaz$FIA == ''] <- "No FIA"
FIA_PipTaz$PipTaz[FIA_PipTaz$PipTaz == ''] <- "S"

ggtree(B1_tree, layout = "circular", branch.length = 'none') %<+% B1_PipTaz + # %<+% adds dataframe with sample data to tree
  aes(color = I(PipTaz))+                       # color the branches according to a variable in your dataframe
  scale_color_manual(
    name = "Resistant",                      # name of your color scheme (will show up in the legend like this)
    breaks = c("I", "R", "S"),                     # the different options in your variable
    labels = c("Interemediate", "Resistant", "Susceptible"),        # how you want the different options named in your legend, allows for formatting
    values = c("pink", "red", "blue"),                  # the color you want to assign to the variable 
    na.value = "") +                        # color NA values in black as well
  #  new_scale_color()+                             # allows to add an additional color scheme for another variable
  geom_tippoint(
    mapping = aes(color = FIA),          # tip color by continent. You may change shape adding "shape = "
    size = 1.5)+                               # define the size of the point at the tip
  geom_tiplab(                             # adds name of sample to tip of its branch 
    color = 'black',                       # (add as many text lines as you wish with + , but you may need to adjust offset value to place them next to each other)
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE)+    
  scale_color_brewer(
    name = "",                   
    palette = "Set2",     
    drop=TRUE,
    na.value = "black")+
  ggtitle("Phylogenetic tree of phylogroup B1_PipTaz in E.coli")+       # title of your graph
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
