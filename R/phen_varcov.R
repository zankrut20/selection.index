#' Phenotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments (sub-plot treatments in SPD)
#' @param replication vector containing replication/blocks (RCBD) or rows (LSD)
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param main_plots vector containing main plot treatments (required for Split Plot Design only)
#' @param design_type experimental design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Phenotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' # RCBD example
#' phen_varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
#' 
#' # Latin Square Design example (requires columns parameter)
#' # phen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#' #            replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
#' 
#' # Split Plot Design example (requires main_plots parameter)
#' # phen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot, 
#' #            replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
phen_varcov<- function (data, genotypes, replication, columns = NULL, main_plots = NULL, 
                        design_type = c("RCBD", "LSD", "SPD"), 
                        method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"))
{
  design_type <- match.arg(design_type)
  
  # Convert to numeric matrix once - avoids repeated data.frame/list conversions
  # storage.mode assignment is faster than as.numeric() on each column
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  headings <- colnames(data)
  
  # Convert factors once - eliminates colnumber² redundant conversions
  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)
  repli <- nlevels(replication)
  genotype <- nlevels(genotypes)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }
  
  # Validate Split Plot Design requirements
  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }
  if (design_type == "SPD") {
    main_plots <- as.factor(main_plots)
  }
  
  # MISSING VALUE HANDLING: Use modular engine for imputation
  # Only process method parameter if missing values are detected
  if (any(!is.finite(data_mat))) {
    # Check if user explicitly provided a method
    method_provided <- !missing(method)
    method <- match.arg(method)
    
    # Warn user if they have missing values but didn't explicitly specify a method
    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
              "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
              call. = FALSE)
    }
    
    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }
  
  # Convert to integer indices - design engine uses integer grouping
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
  
  # Pre-allocate result matrix - avoids growing vector (memory reallocation)
  phenotypic.cov <- matrix(0, nrow = colnumber, ncol = colnumber,
                           dimnames = list(headings, headings))
  
  # DESIGN ENGINE: Use modular design calculations (RCBD, LSD, or SPD)
  # Eliminates ~60 lines of repeated design calculation code
  # Centralizes correction factor, sums of products, mean products computations
  for (i in seq_len(colnumber)) {
    trait1 <- data_mat[, i]
    
    for (j in seq_len(colnumber)) {
      trait2 <- data_mat[, j]
      
      # Single call to design engine replaces ~15 lines of manual calculations
      if (design_type == "RCBD") {
        design_stats <- design_stats(trait1, trait2, gen_idx, rep_idx, 
                                     design_type = "RCBD", calc_type = "mean_products")
        # Genotypic covariance = (GMP - EMP) / r
        GCov <- (design_stats$GMP - design_stats$EMP) / design_stats$n_replications
      } else if (design_type == "LSD") {
        design_stats <- design_stats(trait1, trait2, gen_idx, rep_idx, col_idx,
                                     design_type = "LSD", calc_type = "mean_products")
        # For LSD: Genotypic covariance = (GMP - EMP) / t
        GCov <- (design_stats$GMP - design_stats$EMP) / design_stats$n_genotypes
      } else {
        # SPD
        design_stats <- design_stats(trait1, trait2, gen_idx, rep_idx, main_plots = main_idx,
                                     design_type = "SPD", calc_type = "mean_products")
        # For SPD: Genotypic covariance = (GMP - EMP) / (r * a)
        GCov <- (design_stats$GMP - design_stats$EMP) / (design_stats$n_replications * design_stats$n_main_plots)
      }
      
      # Phenotypic covariance = genotypic + environmental
      phenotypic.cov[i, j] <- GCov + design_stats$EMP
    }
  }
  
  return(phenotypic.cov)
}
