#' Calculate Variance-Covariance Components (Internal Helper)
#'
#' @description
#' Internal function to compute genotypic or phenotypic variance-covariance
#' matrices using math primitives. Used by gen_varcov() and phen_varcov().
#'
#' @param data_mat Numeric matrix of trait data (n_obs x n_traits)
#' @param gen_idx Integer vector of genotype indices
#' @param rep_idx Integer vector of replication indices
#' @param col_idx Integer vector of column indices (for LSD, optional)
#' @param main_idx Integer vector of main plot indices (for SPD, optional)
#' @param design_type Integer design code: 1=RCBD, 2=LSD, 3=SPD
#' @param cov_type Integer: 1=genotypic, 2=phenotypic
#'
#' @return Symmetric variance-covariance matrix (n_traits x n_traits)
#'
#' @details
#' MSG and MSE are mean square matrices computed by design_stats_api(),
#' the centralized engine for experimental design statistics.
#' For phenotypic variance, returns MSG directly.
#' For genotypic variance, returns (MSG - MSE) / r where r is the replication factor.
#'
#' @keywords internal
#' @noRd
.calculate_varcov <- function(data_mat, gen_idx, rep_idx,
                              col_idx = NULL, main_idx = NULL,
                              design_type = 1L, cov_type = 1L) {
  if (!design_type %in% c(1L, 2L, 3L)) {
    stop("design_type must be 1 (RCBD), 2 (LSD), or 3 (SPD).")
  }
  if (!cov_type %in% c(1L, 2L)) {
    stop("cov_type must be 1 (genotypic) or 2 (phenotypic).")
  }
  if (design_type == 2L && is.null(col_idx)) {
    stop("col_idx is required for Latin Square Design (design_type = 2).")
  }
  if (design_type == 3L && is.null(main_idx)) {
    stop("main_idx is required for Split Plot Design (design_type = 3).")
  }


  anova_result <- design_stats_api(
    data_mat, gen_idx, rep_idx,
    col_idx, main_idx, design_type
  )

  MSG <- anova_result$MSG # Mean square for genotypes (n_traits x n_traits matrix)
  MSE <- anova_result$MSE # Mean square for error (n_traits x n_traits matrix)

  if (cov_type == 2L) {
    return(MSG)
  }

  n_rep <- length(unique(rep_idx))

  if (design_type == 1L || design_type == 2L) {
    Vg <- (MSG - MSE) / n_rep
  } else if (design_type == 3L) {
    n_main <- length(unique(main_idx))
    Vg <- (MSG - MSE) / (n_rep * n_main)
  }

  Vg
}


#' Genotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments (sub-plot treatments in SPD)
#' @param replication vector containing replication/blocks (RCBD) or rows (LSD)
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param main_plots vector containing main plot treatments (required for Split Plot Design only)
#' @param design_type experimental design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Genotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' # RCBD example
#' gen_varcov(data = seldata[, 3:9], genotypes = seldata$treat, replication = seldata$rep)
#'
#' # Latin Square Design example (requires columns parameter)
#' # gen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat,
#' #           replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
#'
#' # Split Plot Design example (requires main_plots parameter)
#' # gen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot,
#' #           replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
gen_varcov <- function(data, genotypes, replication, columns = NULL, main_plots = NULL,
                       design_type = c("RCBD", "LSD", "SPD"),
                       method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")) {
  design_type <- match.arg(design_type)

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"


  headings <- colnames(data)

  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)


  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }

  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }
  if (design_type == "SPD") {
    main_plots <- as.factor(main_plots)
  }

  if (any(!is.finite(data_mat))) {
    method_provided <- !missing(method)
    method <- match.arg(method)

    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
        "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
        call. = FALSE
      )
    }

    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }

  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL

  design_code <- switch(design_type,
    "RCBD" = 1L,
    "LSD" = 2L,
    "SPD" = 3L
  )

  genetic.cov <- .calculate_varcov(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code,
    cov_type = 1L
  )

  dimnames(genetic.cov) <- list(headings, headings)

  genetic.cov
}


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
#' phen_varcov(data = seldata[, 3:9], genotypes = seldata$treat, replication = seldata$rep)
#'
#' # Latin Square Design example (requires columns parameter)
#' # phen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat,
#' #            replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
#'
#' # Split Plot Design example (requires main_plots parameter)
#' # phen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot,
#' #            replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
phen_varcov <- function(data, genotypes, replication, columns = NULL, main_plots = NULL,
                        design_type = c("RCBD", "LSD", "SPD"),
                        method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")) {
  design_type <- match.arg(design_type)

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"


  headings <- colnames(data)

  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)


  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }

  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }
  if (design_type == "SPD") {
    main_plots <- as.factor(main_plots)
  }

  if (any(!is.finite(data_mat))) {
    method_provided <- !missing(method)
    method <- match.arg(method)

    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
        "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
        call. = FALSE
      )
    }

    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }

  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL

  design_code <- switch(design_type,
    "RCBD" = 1L,
    "LSD" = 2L,
    "SPD" = 3L
  )

  phenotypic.cov <- .calculate_varcov(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code,
    cov_type = 2L
  )

  dimnames(phenotypic.cov) <- list(headings, headings)

  phenotypic.cov
}
