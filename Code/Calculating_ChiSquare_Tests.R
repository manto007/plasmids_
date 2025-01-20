# Resistance and susceptibility data setup for multiple plasmids
data <- list(
  A = list(
    Ciprofloxacin = list(
      Res = c(F1A = 16, `No F1A` = 10, F1B = 21, `No F1B` = 5, F1C = 28, `No F1C` = 10),
      Sus = c(F1A = 11, `No F1A` = 92, F1B = 67, `No F1B` = 36, F1C = 1, `No F1C` = 75)
    ),
    Gentamicin = list(
      Res = c(F1A = 14, `No F1A` = 1, F1B = 15, `No F1B` = 0, F1C = 11, `No F1C` = 4),
      Sus = c(F1A = 14, `No F1A` = 99, F1B = 73, `No F1B` = 40, F1C = 33, `No F1C` = 80)
    ),
    Amoxicillin = list(
      Res = c(F1A = 16, `No F1A` = 27, F1B = 38, `No F1B` = 5, F1C = 21, `No F1C` = 22),
      Sus = c(F1A = 12, `No F1A` = 75, F1B = 51, `No F1B` = 36, F1C = 23, `No F1C` = 64)
    )
  ),
  B1 = list(
    Ciprofloxacin = list(
      Res = c(F1A = 3, `No F1A` = 8, F1B = 9, `No F1B` = 2, F1C = 3, `No F1C` = 8),
      Sus = c(F1A = 8, `No F1A` = 50, F1B = 31, `No F1B` = 27, F1C = 10, `No F1C` = 48)
    ),
    Gentamicin = list(
      Res = c(F1A = 1, `No F1A` = 3, F1B = 4, `No F1B` = 0, F1C = 1, `No F1C` = 2),
      Sus = c(F1A = 10, `No F1A` = 55, F1B = 36, `No F1B` = 29, F1C = 11, `No F1C` = 54)
    ),
    Amoxicillin = list(
      Res = c(F1A = 3, `No F1A` = 12, F1B = 10, `No F1B` = 5, F1C = 3, `No F1C` = 12),
      Sus = c(F1A = 8, `No F1A` = 46, F1B = 30, `No F1B` = 24, F1C = 10, `No F1C` = 44)
    )
  ),
  B2 = list(
    Ciprofloxacin = list(
      Res = c(F1A = 142, `No F1A` = 22, F1B = 77, `No F1B` = 87, F1C = 22, `No F1C` = 142),
      Sus = c(F1A = 65, `No F1A` = 783, F1B = 522, `No F1B` = 326, F1C = 174, `No F1C` = 674)
    ),
    Gentamicin = list(
      Res = c(F1A = 45, `No F1A` = 15, F1B = 37, `No F1B` = 23, F1C = 3, `No F1C` = 57),
      Sus = c(F1A = 161, `No F1A` = 783, F1B = 559, `No F1B` = 385, F1C = 190, `No F1C` = 754)
    ),
    Amoxicillin = list(
      Res = c(F1A = 139, `No F1A` = 155, F1B = 157, `No F1B` = 137, F1C = 35, `No F1C` = 259),
      Sus = c(F1A = 72, `No F1A` = 652, F1B = 448, `No F1B` = 276, F1C = 162, `No F1C` = 562)
    )
  ),
  D = list(
    Ciprofloxacin = list(
      Res = c(F1A = 52, `No F1A` = 14, F1B = 32, `No F1B` = 34, F1C = 2, `No F1C` = 64),
      Sus = c(F1A = 22, `No F1A` = 198, F1B = 164, `No F1B` = 56, F1C = 32, `No F1C` = 188)
    ),
    Gentamicin = list(
      Res = c(F1A = 17, `No F1A` = 5, F1B = 18, `No F1B` = 4, F1C = 1, `No F1C` = 21),
      Sus = c(F1A = 68, `No F1A` = 286, F1B = 180, `No F1B` = 86, F1C = 33, `No F1C` = 233)
    ),
    Amoxicillin = list(
      Res = c(F1A = 36, `No F1A` = 58, F1B = 61, `No F1B` = 33, F1C = 7, `No F1C` = 87),
      Sus = c(F1A = 42, `No F1A` = 155, F1B = 139, `No F1B` = 58, F1C = 27, `No F1C` = 170)
    )
  )
)

# Function to perform chi-square tests
perform_chi_square_tests <- function(phylogroup, drug, res, sus) {
  cat("\nPhylogroup:", phylogroup, "| Drug:", drug, "\n")
  
# Loop through each plasmid
for (plasmid in c("F1A", "F1B", "F1C")) {
    # Create 2x2 table for the selected plasmid and drug
    plasmid_table <- matrix(c(res[plasmid], res[paste0("No ", plasmid)],  # Resistant
                              sus[plasmid], sus[paste0("No ", plasmid)]), # Susceptible
                            nrow = 2, byrow = TRUE)
    
    # Naming the rows and columns
    rownames(plasmid_table) <- c("Resistant", "Susceptible")
    colnames(plasmid_table) <- c(paste(plasmid, "Present"), paste(plasmid, "Absent"))
    
    # Perform chi-square test
    test_result <- chisq.test(plasmid_table)
    
    # Print results for each plasmid
    print(plasmid_table)
    print(test_result)
  }
}

# Loop through each phylogroup and drug
for (phylogroup in names(data)) {
  for (drug in names(data[[phylogroup]])) {
    res <- data[[phylogroup]][[drug]]$Res
    sus <- data[[phylogroup]][[drug]]$Sus
    perform_chi_square_tests(phylogroup, drug, res, sus)
  }
}

# Initialize a list to store results
results_list <- list()

# Function to perform chi-square tests and store results
perform_chi_square_tests <- function(phylogroup, drug, res, sus) {
  # Loop through each plasmid
  for (plasmid in c("F1A", "F1B", "F1C")) {
    # Create 2x2 table for the selected plasmid and drug
    plasmid_table <- matrix(c(res[plasmid], res[paste0("No ", plasmid)],  # Resistant
                              sus[plasmid], sus[paste0("No ", plasmid)]), # Susceptible
                            nrow = 2, byrow = TRUE)
    
    # Perform chi-square test
    test_result <- chisq.test(plasmid_table)
    
    # Store results in a list
    results_list[[length(results_list) + 1]] <- data.frame(
      Phylogroup = phylogroup,
      Drug = drug,
      Plasmid = plasmid,
      Chi_Square = test_result$statistic,
      P_Value = test_result$p.value,
      Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  }
}

# Loop through each phylogroup and drug
for (phylogroup in names(data)) {
  for (drug in names(data[[phylogroup]])) {
    res <- data[[phylogroup]][[drug]]$Res
    sus <- data[[phylogroup]][[drug]]$Sus
    perform_chi_square_tests(phylogroup, drug, res, sus)
  }
}

# Combine results into a single data frame
results_table <- do.call(rbind, results_list)

# Print the results table
print(results_table)

# Optionally, write to a CSV file
write.csv(results_table, "chi_square_results.csv", row.names = FALSE)
