library(lintr)

# 1. Run the linter and save the results to an object
my_lints <- lint_package()

# 2. Convert the results to a data frame
lint_df <- as.data.frame(my_lints)

# 3. Export the data frame to a CSV file in your working directory
write.csv(lint_df, file = "lint_results.csv", row.names = FALSE)
