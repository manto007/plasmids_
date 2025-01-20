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
    df <- data.frame("Accession" = A_ResData[, "Accession"], 
                     "plasmid" = A_ResData[, plasmid], 
                     "drug" = A_ResData[, drug])
    
    assign(df_name, df)
    
    # Create the plot
    plot <- ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% get(df_name) + 
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
      ggtitle(paste("Phylogenetic tree of phylogroup A", plasmid, drug, "in E.coli")) +       
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
    file_name <- paste0("plot_", plasmid, "_", drug, "_Phylogroup_A.png")
    ggsave(filename = file_name, plot = plot, width = 10, height = 8, dpi = 300)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)

# Load your data
A_ResData <- read.csv("Data/results_A.csv")
B1_ResData <- read.csv("Data/results_B1.csv")
B2_ResData <- read.csv("Data/results_B2.csv")
D_ResData <- read.csv("Data/results_D.csv")


# Reshape the data to long format for plasmids
long_plasmid_data <- B2_ResData %>%
  pivot_longer(cols = c("FIA", "FIB", "FIC"), names_to = "plasmid", values_to = "plasmid_value")

# Reshape the data to long format for drugs
long_drug_data <- long_plasmid_data %>%
  pivot_longer(cols = c("PipTaz", "AmoxiClav", "Amoxicillin", "Ceftazidime", "Ciprofloxacin", "Cefotaxime", "Cefuroxime", "Gentamicin"), names_to = "drug", values_to = "resistance")

# Filter out the rows where plasmid_value is "No_FIA", "No_FIB", or "No_FIC"
filtered_data <- long_drug_data %>%
  filter(plasmid_value != "No_FIA" & plasmid_value != "No_FIB" & plasmid_value != "No_FIC")

# Summarize the data
summary_table <- filtered_data %>%
  group_by(plasmid, drug, resistance) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = resistance, values_from = count, values_fill = 0)

print(summary_table)

# Perform Chi-Square test (if appropriate)
perform_chisq_test <- function(data) {
  table <- xtabs(~ resistance + drug, data = data)
  test <- chisq.test(table)
  return(test)
}

chisq_results <- lapply(unique(filtered_data$plasmid), function(plasmid) {
  perform_chisq_test(filter(filtered_data, plasmid == plasmid))
})

names(chisq_results) <- unique(filtered_data$plasmid)
print(chisq_results)

# Reshape the summarized data for plotting
summary_long <- summary_table %>%
  pivot_longer(cols = c(S, I, R), names_to = "Resistance_Type", values_to = "Count")

# Plot the summarized data
ggplot(summary_long, aes(x = drug, y = Count, fill = Resistance_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ plasmid) +
  labs(title = "Resistance Types for Each Plasmid-Drug Combination in D",
       x = "Drug",
       y = "Count",
       fill = "Resistance Type") +
  theme_minimal()
