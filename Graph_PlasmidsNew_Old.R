# Install and load required packages
install.packages("tidyverse") # if not installed
library(tidyverse)


# Your dataset, assuming you read it from a CSV file (replace 'your_data_file.csv' with the actual file path)
# Assuming the structure matches the table you provided.
data <- read.csv("Plasmids_Full.csv")

# Ensure empty strings are treated as NA
data$FIA <- ifelse(data$FIA == "" | data$FIA == " ", NA, data$FIA)
data$FIB <- ifelse(data$FIB == "" | data$FIB == " ", NA, data$FIB)
data$FIC <- ifelse(data$FIC == "" | data$FIC == " ", NA, data$FIC)

# Count the non-empty entries for each category
fia_count_old <- sum(!is.na(data$FIA))
fib_count_old <- sum(!is.na(data$FIB))
fic_count_old <- sum(!is.na(data$FIC))

# Create a summary table
category_count <- data.frame(
  Category = c("FIA", "FIB", "FIC"),
  Count = c(fia_count_old, fib_count_old, fic_count_old)
)

# Calculate percentages
category_count <- category_count %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

# Plot a more aesthetic pie chart
ggplot(category_count, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.7) +  # Add a black outline
  coord_polar("y", start = 0) +  # Align start of the chart
  theme_void() +  # Remove axes
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#7CAE00")) +  # Aesthetic color palette
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 5, 
            fontface = "bold") +  # Adjust label size and color
  labs(title = "Proportion of FIA, FIB, and FIC_old") +  # Chart title
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style title
    legend.position = "right",  # Adjust legend position
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 12)  # Increase legend text size
  )

#WITH THE NEW DATASET 

data_new <- read.csv("Plasmids_New.csv")

# Ensure empty strings are treated as NA
data_new$FIA <- ifelse(data_new$FIA == "" | data_new$FIA == " ", NA, data_new$FIA)
data_new$FIB <- ifelse(data_new$FIB == "" | data_new$FIB == " ", NA, data_new$FIB)
data_new$FIC <- ifelse(data_new$FIC == "" | data_new$FIC == " ", NA, data_new$FIC)

# Count the non-empty entries for each category
fia_count_new <- sum(!is.na(data_new$FIA))
fib_count_new <- sum(!is.na(data_new$FIB))
fic_count_new <- sum(!is.na(data_new$FIC))

# Create a summary table
category_count <- data.frame(
  Category = c("FIA", "FIB", "FIC"),
  Count = c(fia_count_new, fib_count_new, fic_count_new)
)

# Calculate percentages
category_count <- category_count %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

# Plot a more aesthetic pie chart
ggplot(category_count, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.7) +  # Add a black outline
  coord_polar("y", start = 0) +  # Align start of the chart
  theme_void() +  # Remove axes
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#7CAE00")) +  # Aesthetic color palette
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white", 
            size = 5, 
            fontface = "bold") +  # Adjust label size and color
  labs(title = "Proportion of FIA, FIB, and FIC_new") +  # Chart title
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and style title
    legend.position = "right",  # Adjust legend position
    legend.title = element_blank(),  # Remove legend title for simplicity
    legend.text = element_text(size = 12)  # Increase legend text size
  )

