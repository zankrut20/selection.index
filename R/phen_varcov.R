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
  
  # C++ OPTIMIZATION: Vectorized variance-covariance computation
  # Replaces nested R loops with single C++ call using Eigen linear algebra
  # Processes all trait pairs simultaneously with optimized grouped sums
  # Expected speedup: 5-20x for 7-30 traits
  design_code <- switch(design_type, "RCBD" = 1L, "LSD" = 2L, "SPD" = 3L)
  
  phenotypic.cov <- cpp_varcov_iterator(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code,
    cov_type = 2L  # 2 = phenotypic
  )
  
  # Restore dimension names
  dimnames(phenotypic.cov) <- list(headings, headings)
  
  return(phenotypic.cov)
}
