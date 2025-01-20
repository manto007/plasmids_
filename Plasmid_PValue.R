# Load necessary libraries
library(dplyr)

# Read your dataset
df <- read.csv("Data/results_A.csv")

# Load necessary libraries
library(dplyr)

# Read your dataset
#df <- read.csv("Data/results_A.csv")

# Define a function to calculate the percentages for one FIA/FIB/FIC vs. one resistance column
calculate_percentages <- function(fia_col, resistance_col) {
  # Select relevant columns
  selected_df <- df %>%
    select(Accession, {{fia_col}}, {{resistance_col}})
  
  # Count occurrences for each combination and calculate percentages
  count_df <- selected_df %>%
    group_by({{fia_col}}, {{resistance_col}}) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(percentage = count / sum(count) * 100) %>%
    rename(FIA_FIB_FIC = {{fia_col}}, Resistance_Drug = {{resistance_col}})
  
  return(count_df)
}

# Define the columns to loop through
fia_fib_fic_cols <- c("FIA", "FIB", "FIC")
resistance_cols <- c("Ciprofloxacin", "Gentamicin", "AmoxiClav")

# Initialize an empty data frame to store results
results_df <- data.frame()

# Loop through each FIA/FIB/FIC and each resistance column
for (fia_col in fia_fib_fic_cols) {
  for (res_col in resistance_cols) {
    # Print the current columns being processed
    cat("Percentage for", fia_col, "vs.", res_col, ":\n")
    
    # Calculate the percentages
    percentage_df <- calculate_percentages(!!sym(fia_col), !!sym(res_col))
    
    # Add the Drug column with the actual drug name
    percentage_df <- percentage_df %>%
      mutate(Drug = res_col)  # Assign drug name here
    
    # Print the result in tibble format
    print(percentage_df)
    
    # Append the results to the results_df using bind_rows
    results_df <- bind_rows(results_df, percentage_df)
  }
}

# Optionally, save the results to a CSV file
write.csv(results_df, "results_percentages_combined_D.csv", row.names = FALSE) #change for every phylogroup when saving 

# Display the combined results data frame
print(results_df)


# Now, create the contingency table
contingency_table <- xtabs(count ~ FIA_FIB_FIC + Resistance_Drug, data = ciprofloxacin_df)
contingency_table <- xtabs(count ~ FIA_FIB_FIC + Resistance_Drug, data = gentamicin_df)
contingency_table <- xtabs(count ~ FIA_FIB_FIC + Resistance_Drug, data = amoxiclav_df)

# Print the contingency table
print(contingency_table)

# Calculate p-value using chi-squared test
chisq_test_result <- chisq.test(contingency_table)

# Print the p-value
print(chisq_test_result$p.value)

# Print the contingency table
print(contingency_table)

# Print the p-value
print(chisq_test_result$p.value)





