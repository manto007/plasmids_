library(ggtree)
library(ggplot2)
library(readr)

# Define the plasmid and drug lists
plasmidlist <- c("FIA", "FIB", "FIC")
druglist <- c("PipTaz", "AmoxiClav", "Amoxicillin", "Ceftazidime", "Ciprofloxacin", "Cefotaxime", "Cefuroxime", "Gentamicin")

# Load your data
A_ResData <- read.csv("Data/results_A.csv")
B1_ResData <- read.csv("Data/results_B1.csv")
B2_ResData <- read.csv("Data/results_B2.csv")
D_ResData <- read.csv("Data/results_D.csv")


# Loop through each plasmid and drug
for (plasmid in plasmidlist) {
  for (drug in druglist) {
    # Construct the dataframe dynamically
    df_name <- paste0(plasmid, "_", drug, "_D")
    df <- data.frame("Accession" = D_ResData[, "Accession"], 
                     "plasmid" = D_ResData[, plasmid], 
                     "drug" = D_ResData[, drug])
    
    assign(df_name, df)
    
    # Create the plot
    plot <- ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% get(df_name) + 
      aes(color = drug) +                       
      scale_color_manual(
        name = "Resistant",                      
        breaks = c("I", "R", "S"),                     
        labels = c("Intermediate", "Resistant", "Susceptible"),        
        values = c("pink", "red", "white"),                  
        na.value = "") +                        
      geom_tippoint(
        mapping = aes(color = plasmid),          
        size = 1.5) +                               
      geom_tiplab(                             
        color = 'black',                       
        offset = 1,
        size = 1,
        geom = "text",
        align = TRUE) +    
      scale_color_brewer(
        name = "",                   
        palette = "Set1",     
        drop = TRUE,
        na.value = "black") +
      ggtitle(paste("Phylogenetic tree of phylogroup", plasmid, drug, "in E.coli")) +       
      theme(
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.title = element_text(    
          face = "bold",
          size = 12),   
        legend.text = element_text(       
          face = "bold",
          size = 10),  
        plot.title = element_text(      
          size = 12,
          face = "bold"),  
        legend.position = "bottom",     
        legend.box = "vertical",        
        legend.margin = margin()
      )
    
    # Save the plot as a PNG file
    file_name <- paste0("plot_", plasmid, "_", drug, "_Phylogroup_D.png")
    ggsave(filename = file_name, plot = plot, width = 10, height = 8, dpi = 300)
  }
}
