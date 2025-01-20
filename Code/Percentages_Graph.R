# Load required library
library(ggplot2)

# Create a data frame with the provided data
data <- data.frame(
  Phylogroup = rep(c("A", "B1", "B2", "D"), each = 18),
  Drug = rep(rep(c("Ciprofloxacin", "Gentamicin", "AmoxiClav"), each = 6), 4),
  Plasmid = rep(c("FIA", "No FIA", "FIB", "No FIB", "FIC", "No FIC"), 12),
  Percentage = c(
    # Phylogroup A
    60, 10, 24, 12, 96, 12, 50, 1, 17, 0, 25, 5, 57, 26, 43, 12, 48, 26,
    # Phylogroup B1
    27, 14, 22, 7, 23, 14, 9, 5, 10, 0, 15, 4, 27, 21, 25, 17, 23, 21,
    # Phylogroup B2
    68, 3, 13, 21, 11, 17, 22, 2, 6, 6, 2, 7, 66, 19, 26, 33, 18, 32,
    # Phylogroup D
    70, 7, 16, 38, 6, 25, 20, 2, 9, 4, 3, 8, 46, 27, 31, 36, 21, 34
  )
)

# Filter data for the individual phylogroups
phylogroups <- unique(data$Phylogroup)

# Create the plots for each phylogroup
for (phylogroup in phylogroups) {
  plot_data <- subset(data, Phylogroup == phylogroup)
  
  # Create a scatter plot with drug symbols and plasmid key
  p <- ggplot(plot_data, aes(x = Plasmid, y = Percentage, color = Drug, shape = Plasmid)) +
    geom_point(size = 4) +
    labs(
      title = paste("Percentage of Plasmid Types for Phylogroup", phylogroup),
      x = "Plasmid Type",
      y = "Percentage",
      color = "Drug",
      shape = "Plasmid Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot for each phylogroup
  print(p)
}
