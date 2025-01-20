library(dplyr)
library(tidyr)

# Load your data
A_ResData <- read.csv("Data/results_A.csv")
B1_ResData <- read.csv("Data/results_B1.csv")
B2_ResData <- read.csv("Data/results_B2.csv") 
D_ResData <- read.csv("Data/results_D.csv")


# Reshape the data to long format for drugs
long_data_A <- A_ResData %>%
  pivot_longer(cols = starts_with("PipTaz"):starts_with("Gentamicin"), names_to = "drug", values_to = "resistance")
long_data_B1 <- B1_ResData %>%
  pivot_longer(cols = starts_with("AmoxiClav"):starts_with("PipTaz"), names_to = "drug", values_to = "resistance")
long_data_B2 <- B2_ResData %>%
  pivot_longer(cols = starts_with("AmoxiClav"):starts_with("PipTaz"), names_to = "drug", values_to = "resistance")
long_data_D <- D_ResData %>%
  pivot_longer(cols = starts_with("AmoxiClav"):starts_with("PipTaz"), names_to = "drug", values_to = "resistance")

# Filter for strains with FIC
with_FIC_A <- long_data_A %>%
  filter(FIC == "FIC")
with_FIC_B1 <- long_data_B1 %>%
  filter(FIC == "FIC")
with_FIC_B2 <- long_data_B2 %>%
  filter(FIC == "FIC")
with_FIC_D <- long_data_D %>%
  filter(FIC == "FIC")

# Filter for strains without FIC
without_FIC_A <- long_data_A %>%
  filter(FIC == "No_FIC")
without_FIC_B1 <- long_data_B1 %>%
  filter(FIC == "No_FIC")
without_FIC_B2 <- long_data_B2 %>%
  filter(FIC == "No_FIC")
without_FIC_D <- long_data_D %>%
  filter(FIC == "No_FIC")

# Calculate the fraction of resistant strains for each drug with FIC
summary_with_FIC_A <- with_FIC_A %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_with_FIC_B1 <- with_FIC_B1 %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_with_FIC_B2 <- with_FIC_B2 %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_with_FIC_D <- with_FIC_D %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)

# Calculate the fraction of resistant strains for each drug without FIC
summary_without_FIC_A <- without_FIC_A %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_without_FIC_B1 <- without_FIC_B1 %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_without_FIC_B2 <- without_FIC_B2 %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)
summary_without_FIC_D <- without_FIC_D %>%
  group_by(drug) %>%
  summarise(total = n(),
            resistant = sum(resistance == "R"),
            fraction_resistant = resistant / total)

# Combine the results
summary_combined_A <- summary_with_FIC_A %>%
  rename(total_with_FIC_A = total,
         resistant_with_FIC_A = resistant,
         fraction_resistant_with_FIC_A = fraction_resistant) %>%
  full_join(summary_without_FIC_A %>%
              rename(total_without_FIC_A = total,
                     resistant_without_FIC_A = resistant,
                     fraction_resistant_without_FIC_A = fraction_resistant),
            by = "drug")

print(summary_combined_A)

summary_combined_B1 <- summary_with_FIC_B1 %>%
  rename(total_with_FIC_B1 = total,
         resistant_with_FIC_B1 = resistant,
         fraction_resistant_with_FIC_B1 = fraction_resistant) %>%
  full_join(summary_without_FIC_B1 %>%
              rename(total_without_FIC_B1 = total,
                     resistant_without_FIC_B1 = resistant,
                     fraction_resistant_without_FIC_B1 = fraction_resistant),
            by = "drug")

print(summary_combined_B1)

summary_combined_B2 <- summary_with_FIC_B2 %>%
  rename(total_with_FIC_B2 = total,
         resistant_with_FIC_B2 = resistant,
         fraction_resistant_with_FIC_B2 = fraction_resistant) %>%
  full_join(summary_without_FIC_B2 %>%
              rename(total_without_FIC_B2 = total,
                     resistant_without_FIC_B2 = resistant,
                     fraction_resistant_without_FIC_B2 = fraction_resistant),
            by = "drug")

print(summary_combined_B2)

summary_combined_D <- summary_with_FIC_D %>%
  rename(total_with_FIC_D = total,
         resistant_with_FIC_D = resistant,
         fraction_resistant_with_FIC_D = fraction_resistant) %>%
  full_join(summary_without_FIC_D %>%
              rename(total_without_FIC_D = total,
                     resistant_without_FIC_D = resistant,
                     fraction_resistant_without_FIC_D = fraction_resistant),
            by = "drug")

print(summary_combined_D)

# Reshape data for plotting
plot_data_A <- summary_combined_A %>%
  pivot_longer(cols = c(fraction_resistant_with_FIC_A, fraction_resistant_without_FIC_A), 
               names_to = "Group", 
               values_to = "Fraction_Resistant")
plot_data_B1 <- summary_combined_B1 %>%
  pivot_longer(cols = c(fraction_resistant_with_FIC_B1, fraction_resistant_without_FIC_B1), 
               names_to = "Group", 
               values_to = "Fraction_Resistant")
plot_data_B2 <- summary_combined_B2 %>%
  pivot_longer(cols = c(fraction_resistant_with_FIC_B2, fraction_resistant_without_FIC_B2), 
               names_to = "Group", 
               values_to = "Fraction_Resistant")
plot_data_D <- summary_combined_D %>%
  pivot_longer(cols = c(fraction_resistant_with_FIC_D, fraction_resistant_without_FIC_D), 
               names_to = "Group", 
               values_to = "Fraction_Resistant")
# Rename groups for better readability
plot_data_A$Group <- recode(plot_data_A$Group,
                          "fraction_resistant_with_FIC_A" = "With FIC",
                          "fraction_resistant_without_FIC_A" = "Without FIC")
plot_data_B1$Group <- recode(plot_data_B1$Group,
                            "fraction_resistant_with_FIC_B1" = "With FIC",
                            "fraction_resistant_without_FIC_B1" = "Without FIC")
plot_data_B2$Group <- recode(plot_data_B2$Group,
                             "fraction_resistant_with_FIC_B2" = "With FIC",
                             "fraction_resistant_without_FIC_B2" = "Without FIC")
plot_data_D$Group <- recode(plot_data_D$Group,
                             "fraction_resistant_with_FIC_D" = "With FIC",
                             "fraction_resistant_without_FIC_D" = "Without FIC")

pdf("AllPlasmid_Plots.pdf")

# Plot the results
ggplot(plot_data_A, aes(x = drug, y = Fraction_Resistant, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fraction of Resistant Strains by Drug and FIC Presence in Phylogroup A",
       x = "Drug",
       y = "Fraction Resistant",
       fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_data_B1, aes(x = drug, y = Fraction_Resistant, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fraction of Resistant Strains by Drug and FIC Presence in Phylogroup B1",
       x = "Drug",
       y = "Fraction Resistant",
       fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_data_B2, aes(x = drug, y = Fraction_Resistant, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fraction of Resistant Strains by Drug and FIC Presence in Phylogroup B2",
       x = "Drug",
       y = "Fraction Resistant",
       fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_data_D, aes(x = drug, y = Fraction_Resistant, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fraction of Resistant Strains by Drug and FIC Presence in Phylogroup D",
       x = "Drug",
       y = "Fraction Resistant",
       fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


library(dplyr)
library(tidyr)
library(ggplot2)

# Define the plasmid list
plasmid_list <- c("FIC", "FIC", "FIC")

# Load your data
A_ResData <- read.csv("Data/results_A.csv")

# Reshape the data to long format for drugs
long_data <- A_ResData %>%
  pivot_longer(cols = c("PipTaz", "AmoxiClav", "Amoxicillin", "Ceftazidime", "Ciprofloxacin", "Cefotaxime", "Cefuroxime", "Gentamicin"), names_to = "drug", values_to = "resistance")

# Initialize an empty data frame to store the combined results
combined_results <- data.frame()

# Loop through each plasmid
for (plasmid in plasmid_list) {
  # Filter for strains with the current plasmid
  with_plasmid <- long_data %>%
    filter(get(plasmid) == plasmid)
  
  # Filter for strains without the current plasmid
  without_plasmid <- long_data %>%
    filter(get(plasmid) == paste("No_", plasmid, sep = ""))
  
  # Calculate the fraction of resistant strains for each drug with the plasmid
  summary_with_plasmid <- with_plasmid %>%
    group_by(drug) %>%
    summarise(total = n(),
              resistant = sum(resistance == "R"),
              fraction_resistant = resistant / total) %>%
    mutate(plasmid_group = paste(plasmid, "With"))
  
  # Calculate the fraction of resistant strains for each drug without the plasmid
  summary_without_plasmid <- without_plasmid %>%
    group_by(drug) %>%
    summarise(total = n(),
              resistant = sum(resistance == "R"),
              fraction_resistant = resistant / total) %>%
    mutate(plasmid_group = paste(plasmid, "Without"))
  
  # Combine the results for the current plasmid
  combined_plasmid_results <- bind_rows(summary_with_plasmid, summary_without_plasmid)
  
  # Add the results to the combined data frame
  combined_results <- bind_rows(combined_results, combined_plasmid_results)
}

# Reshape the combined results for plotting
plot_data <- combined_results %>%
  select(drug, fraction_resistant, plasmid_group)

# Plot the results
ggplot(plot_data, aes(x = drug, y = fraction_resistant, fill = plasmid_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Fraction of Resistant Strains by Drug and Plasmid Presence",
       x = "Drug",
       y = "Fraction Resistant",
       fill = "Plasmid Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


