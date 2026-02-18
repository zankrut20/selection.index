#' Missing Value Imputation for Experimental Designs
#'
#' @description
#' Exported wrapper for missing value estimation in RCBD, Latin Square, and 
#' Split Plot designs. Calls the centralized design_stats engine for ANOVA
#' computations.
#'
#' @name missing-value
NULL


#' Impute Missing Values in Experimental Data
#'
#' @description
#' Estimates and imputes missing values in randomized complete block design (RCBD),
#' Latin square design (LSD), or split plot design (SPD) experimental data.
#' 
#' Uses one of six methods: REML, Yates, Healy, Regression, Mean, or Bartlett.
#'
#' @param data Matrix or data.frame with observations (rows) by traits (columns).
#'   May contain missing values (NA, NaN, Inf).
#' @param genotypes Vector indicating genotype/treatment for each observation 
#'   (sub-plot treatments in SPD).
#' @param replications Vector indicating replication/block (RCBD) or row (LSD) 
#'   for each observation.
#' @param columns Vector indicating column index for each observation 
#'   (required for Latin Square Design only).
#' @param main_plots Vector indicating main plot treatment for each observation 
#'   (required for Split Plot Design only).
#' @param design Character string specifying experimental design:
#'   \itemize{
#'     \item "RCBD" - Randomized Complete Block Design (default)
#'     \item "LSD" - Latin Square Design
#'     \item "SPD" - Split Plot Design
#'   }
#' @param method Character string specifying the estimation method:
#'   \itemize{
#'     \item \strong{REML} - Restricted Maximum Likelihood with BLUP. Most robust 
#'       for complex missing patterns. (RCBD, LSD only)
#'     \item \strong{Yates} - Traditional iterative formula. Simple and fast. 
#'       Good for simple missing patterns. (RCBD, LSD only)
#'     \item \strong{Healy} - Healy & Westmacott weighted adjustment method. 
#'       More stable than Yates for multiple missing values. (RCBD, LSD only)
#'     \item \strong{Regression} - Linear regression with QR decomposition. 
#'       Non-iterative, fast and stable. (RCBD, LSD only)
#'     \item \strong{Mean} - Mean substitution using treatment and block effects. 
#'       Non-iterative, fastest. (RCBD, LSD, SPD)
#'     \item \strong{Bartlett} - ANCOVA using other traits as covariates. 
#'       Best when traits are correlated. (RCBD, LSD only)
#'   }
#' @param tolerance Numeric convergence criterion for iterative methods. 
#'   Iteration stops when maximum change in estimated values falls below this 
#'   threshold. Default: 1e-6.
#'
#' @return Matrix of the same dimensions as \code{data} with all missing values
#'   replaced by estimates.
#'
#' @details
#' The function handles missing values by iteratively estimating them based on
#' the experimental design structure:
#' 
#' **RCBD:** 2-way blocking (genotypes × blocks)  
#' **LSD:** 3-way blocking (genotypes × rows × columns)  
#' **SPD:** Nested structure (blocks > main plots > sub-plots)
#' 
#' \strong{Method Availability:}
#' \itemize{
#'   \item RCBD: All methods (REML, Yates, Healy, Regression, Mean, Bartlett)
#'   \item LSD: All methods (REML, Yates, Healy, Regression, Mean, Bartlett)
#'   \item SPD: Mean only (other methods fall back to Mean)
#' }
#' 
#' \strong{Method Selection Guide:}
#' \itemize{
#'   \item Use \strong{REML} for complex missing patterns or when precision is critical
#'   \item Use \strong{Yates} for balanced designs with few missing values
#'   \item Use \strong{Healy} when multiple values missing from same treatment/block
#'   \item Use \strong{Regression} for fast, deterministic estimation
#'   \item Use \strong{Mean} for quick estimation when precision is less critical
#'   \item Use \strong{Bartlett} when traits are highly correlated
#' }
#' 
#' The function uses the centralized design_stats engine for all ANOVA computations,
#' ensuring consistency with gen_varcov(), phen_varcov(), and mean_performance().
#'
#' @references
#' Yates, F. (1933). The analysis of replicated experiments when the field results
#' are incomplete. \emph{Empire Journal of Experimental Agriculture}, 1, 129-142.
#'
#' Healy, M. J. R., & Westmacott, M. (1956). Missing values in experiments analysed
#' on automatic computers. \emph{Applied Statistics}, 5(3), 203-206.
#'
#' Bartlett, M. S. (1937). Some examples of statistical methods of research in
#' agriculture and applied biology. \emph{Supplement to the Journal of the Royal
#' Statistical Society}, 4(2), 137-183.
#'
#' @keywords internal
#' @noRd
#' @examples
#' # RCBD example with missing values
#' data(seldata)
#' test_data <- seldata[, 3:5]
#' test_data[c(1, 10, 25), 1] <- NA
#' test_data[c(5, 15), 2] <- NA
#' 
#' # Impute using Yates method
#' imputed <- impute_missing(test_data, seldata$treat, seldata$rep, method = "Yates")
#' 
#' # Check that no NA remain
#' anyNA(imputed)  # Should be FALSE
#' 
#' \dontrun{
#' # Latin Square Design example
#' # lsd_data should have genotypes, rows, and columns
#' imputed_lsd <- impute_missing(
#'   data = lsd_data[, 3:7],
#'   genotypes = lsd_data$treat,
#'   replications = lsd_data$row,
#'   columns = lsd_data$col,
#'   design = "LSD",
#'   method = "REML"
#' )
#' 
#' # Split Plot Design example
#' # spd_data should have sub-plots, blocks, and main plots
#' imputed_spd <- impute_missing(
#'   data = spd_data[, 3:7],
#'   genotypes = spd_data$subplot,
#'   replications = spd_data$block,
#'   main_plots = spd_data$mainplot,
#'   design = "SPD",
#'   method = "Mean"
#' )
#' }
impute_missing <- function(data, genotypes, replications, 
                          columns = NULL, main_plots = NULL,
                          design = c("RCBD", "LSD", "SPD"),
                          method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
                          tolerance = 1e-6) {
  
  # Match arguments
  design <- match.arg(design)
  method <- match.arg(method)
  
  # Convert data to matrix
  data_mat <- as.matrix(data)
  
  # Convert indices to integer vectors
  # Simply convert to factor then integer - works for numeric, character, factor
  gen_idx <- as.integer(factor(genotypes))
  rep_idx <- as.integer(factor(replications))
  col_idx <- if (!is.null(columns)) as.integer(factor(columns)) else NULL
  main_idx <- if (!is.null(main_plots)) as.integer(factor(main_plots)) else NULL
  
  # Validate design-specific requirements
  design_code <- switch(design,
                        "RCBD" = DESIGN_RCBD,
                        "LSD" = DESIGN_LSD,
                        "SPD" = DESIGN_SPD)
  
  validate_design_args(design_code, col_idx, main_idx)
  
  # Validate indices match data dimensions
  validate_indices(nrow(data_mat), gen_idx, rep_idx, col_idx, main_idx, "data")
  
  # Call internal estimation function
  result <- missing_value_estimation(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design,
    method = method,
    tolerance = tolerance
  )
  
  # Preserve column names
  if (!is.null(colnames(data))) {
    colnames(result) <- colnames(data)
  }
  
  return(result)
}
