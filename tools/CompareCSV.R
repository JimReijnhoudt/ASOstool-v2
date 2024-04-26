library(dplyr)
library(readr)

# Read the CSV files
csv1 <- read_csv("Results_output_2024-04-05 13:02:11.csv")
csv2 <- read_csv("Results_output_2024-04-05 14:02:24.csv")

# Compare the CSV files and find differences
# Find rows in csv2 that are not in csv1
differences_second_to_first <- csv2 %>%
  anti_join(csv1)

# Find rows in csv1 that are not in csv2
differences_first_to_second <- csv1 %>%
  anti_join(csv2)

# Function to print detailed differences
print_detailed_differences <- function(differences, csv_source, csv_comparison) {
  if (nrow(differences) > 0) {
    print(paste("Differences found in", csv_source, "not in", csv_comparison, ":"))
    for (i in 1:nrow(differences)) {
      different_row <- differences[i, ]
      print(paste("Row", i, "in", csv_source, "does not exist in", csv_comparison, ". Detailed differences:"))
      for (col in names(differences)) {
        if (!is.na(different_row[[col]])) {
          print(paste(" - Column", col, "value:", different_row[[col]]))
        }
      }
    }
  } else {
    print(paste("No differences found in", csv_source, "not in", csv_comparison, "."))
  }
}

# Print the differences along with the line numbers
print_detailed_differences(differences_second_to_first, "second CSV", "first CSV")
print_detailed_differences(differences_first_to_second, "first CSV", "second CSV")
