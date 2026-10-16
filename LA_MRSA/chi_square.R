

# Load data from CSV as a matrix (assuming first column is row names)
# data_matrix <- as.matrix(read.csv("~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/nov5_431sample_reordered_freq_by_amr.csv", row.names = 1))
# data_matrix <- as.matrix(read.csv("~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/nov5_431sample_reordered_alldata_LA_freq_by_amr.csv", row.names = 1))
data_matrix <- as.matrix(read.csv("~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/nov5_431sample_freq_by_amr_LA264_by_host.csv", row.names = 1))
# Function to perform Chi-square test for independence on each row
# This assumes each row has counts in the order: Human_Y, Human_N, Animal_Y, Animal_N
test_row <- function(row) {
  # Reshape the row into a 2x2 contingency table
  contingency_table <- matrix(row, nrow = 2, byrow = TRUE,
                              dimnames = list(c("Human", "Animal"), c("Y", "N")))
  print(contingency_table)  
  # Perform Chi-square test
  test_result <- chisq.test(contingency_table)
  
  # Return p-value (or you could return other values like test_result$statistic)
  return(test_result$p.value)
}

# Apply the test to each row and store the p-values
p_values <- apply(data_matrix, 1, test_row)

# Create a data frame with the results, including row names for easy reference
results <- data.frame(Row = rownames(data_matrix), P_Value = p_values)

# Print results
print(results)

write.csv(results, "~/Output_Nov5ccir_2rpt/pt04_pathgonwatch/amr_LA264_by_host_chi_square.csv", row.names=FALSE, quote=FALSE) 
# Confirmation message
print("Results saved to 'chi_square_all_gene_results.csv'")

