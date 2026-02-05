library(selection.index)
data(seldata)

# Exact test case from test-mean_performance.R:169-178
test_data <- seldata[1:40, 3:5]
test_data[c(1, 5, 10), 1] <- NA
test_data[c(2, 8), 2] <- NA

cat("Test data NAs:\n")
cat("Column 1 NAs at rows:", which(is.na(test_data[,1])), "\n")
cat("Column 2 NAs at rows:", which(is.na(test_data[,2])), "\n")
cat("Column 3 NAs at rows:", which(is.na(test_data[,3])), "\n")

tryCatch({
  performance <- mean_performance(data = test_data, 
                                genotypes = seldata[1:40, 2], 
                                replications = seldata[1:40, 1])
  
  cat("\nPerformance result:\n")
  print(performance)
  
  cat("\nResult subset [,-1]:\n")
  result_subset <- as.matrix(performance[, -1])
  cat("Has NaN/Inf:", any(!is.finite(result_subset)), "\n")
  
  if (any(!is.finite(result_subset))) {
    cat("Locations of non-finite values:\n")
    print(which(!is.finite(result_subset), arr.ind=TRUE))
  }
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  print(traceback())
})
