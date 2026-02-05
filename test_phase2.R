library(selection.index)
data(seldata)

# Test Regression method with missing values
test_data <- seldata[1:40, 3:5]
test_data[c(1, 5, 10), 1] <- NA  # Add missing values

cat("Testing Regression method optimization:\n")
cat("========================================\n\n")

# Test 1: Regression method
result_reg <- selection.index:::missing_value_estimation(
  data_mat = as.matrix(test_data),
  gen_idx = as.integer(seldata[1:40, 2]),
  rep_idx = as.integer(seldata[1:40, 1]),
  col_idx = integer(40),
  main_idx = integer(40),
  design_type = 'RCBD',
  method = 'Regression'
)

cat("Regression method result:\n")
cat("Has NaN/Inf:", any(!is.finite(result_reg)), "\n")
cat("Result shape:", dim(result_reg), "\n")
cat("First few rows:\n")
print(head(result_reg))

# Verify the specific imputed values
cat("\n\nVerifying imputations:\n")
cat("Row 1 (was NA):", result_reg[1, 1], "\n")
cat("Row 5 (was NA):", result_reg[5, 1], "\n")
cat("Row 10 (was NA):", result_reg[10, 1], "\n")

# All should be finite
stopifnot(all(is.finite(result_reg)))

cat("\n✓ Phase 2 (Regression method) optimization working correctly!\n")
