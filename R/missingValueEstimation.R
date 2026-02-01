#' Missing Value Estimation for RCBD Design
#'
#' Estimates missing values in a Randomized Complete Block Design (RCBD) using
#' one of three methods: REML, Yates, or Bartlett (Covariate).
#'
#' @param data_mat Numeric matrix with observations (rows) by traits (columns).
#'   May contain missing values (NA, NaN, Inf).
#' @param gen_idx Integer vector indicating genotype/treatment index for each observation.
#' @param rep_idx Integer vector indicating replicate/block index for each observation.
#' @param method Character string specifying the estimation method. One of:
#'   \describe{
#'     \item{REML}{Restricted Maximum Likelihood with variance components and BLUP
#'       (Best Linear Unbiased Prediction). Most robust for complex missing patterns.
#'       Default method. Iterations: 100.}
#'     \item{Yates}{Traditional iterative formula: (r*T + t*B - G) / ((r-1)(t-1)).
#'       Simple and fast. Good for simple missing patterns. Iterations: 50.}
#'     \item{Healy}{Healy & Westmacott method using weighted adjustment of treatment
#'       and block means. More stable than Yates for multiple missing values. Iterations: 50.}
#'     \item{Regression}{Linear regression on treatment and block effects with QR
#'       decomposition. Non-iterative, single-pass estimation. Fast and stable.}
#'     \item{Mean}{Simple mean substitution using treatment and block effects.
#'       Non-iterative, fastest method. Good for quick estimation.}
#'     \item{Bartlett}{ANCOVA using other traits as covariates with QR decomposition.
#'       Best when traits are correlated. Iterations: 30.}
#'   }
#' @param tolerance Numeric convergence criterion. Iteration stops when maximum
#'   change in estimated values is below this threshold. Default: 1e-6.
#'
#' @return A numeric matrix of the same dimensions as \code{data_mat} with all
#'   missing values replaced by estimates.
#'
#' @details
#' The function handles missing values in RCBD experiments by iteratively
#' estimating them until convergence. Each method has different strengths:
#'
#' \strong{REML Method:}
#' Uses variance component estimation with restricted maximum likelihood.
#' Applies shrinkage to predictions based on estimated variance ratios.
#' Most computationally intensive but handles complex patterns well.
#'
#' \strong{Yates Method:}
#' Uses the classical formula accounting for treatment, block, and grand totals.
#' Fast and simple, suitable for balanced designs with few missing values.
#'
#' \strong{Healy & Westmacott Method:}
#' Uses weighted adjustment based on residuals from treatment and block effects.
#' More stable than Yates when multiple values are missing from same treatment or block.
#' Particularly effective for unbalanced missing patterns.
#'
#' \strong{Regression Method:}
#' Fits linear model using complete observations with treatment and block as factors.
#' Non-iterative single-pass estimation using QR decomposition for numerical stability.
#' Fast and deterministic - no convergence iterations needed.
#'
#' \strong{Mean Substitution Method:}
#' Replaces missing values with treatment mean + block mean - grand mean.
#' Non-iterative, fastest method. Simple additive model without regression.
#' Best for quick estimation when precision is less critical.
#'
#' \strong{Bartlett Method:}
#' Performs ANCOVA using other traits as covariates to predict missing values.
#' Leverages trait correlations for better predictions when traits are related.
#' Uses QR decomposition for numerical stability.
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
#' @examples
#' # Create sample data with missing values
#' data <- matrix(rnorm(100), ncol = 5)
#' data[c(3, 7, 15)] <- NA
#' genotypes <- rep(1:10, each = 2)
#' replicates <- rep(1:2, times = 10)
#'
#' # Estimate using REML (default)
#' data_reml <- missingValueEstimation(data, genotypes, replicates)
#'
#' # Estimate using Yates method
#' data_yates <- missingValueEstimation(data, genotypes, replicates, method = "Yates")
#'
#' # Estimate using Healy & Westmacott method
#' data_healy <- missingValueEstimation(data, genotypes, replicates, method = "Healy")
#'
#' # Estimate using Regression method
#' data_regression <- missingValueEstimation(data, genotypes, replicates, method = "Regression")
#'
#' # Estimate using Mean substitution method
#' data_mean <- missingValueEstimation(data, genotypes, replicates, method = "Mean")
#'
#' # Estimate using Bartlett method
#' data_bartlett <- missingValueEstimation(data, genotypes, replicates, method = "Bartlett")
#'
#' @export
missingValueEstimation <- function(data_mat, gen_idx, rep_idx,
                                   method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
                                   tolerance = 1e-6) {
  # Validate and match method argument
  method <- match.arg(method)
  
  # Get dimensions
  nobs <- nrow(data_mat)
  ncols <- ncol(data_mat)
  genotype <- length(unique(gen_idx))
  repli <- length(unique(rep_idx))
  
  # Check if there are any missing values
  if (!any(!is.finite(data_mat))) {
    return(data_mat)  # No missing values, return as is
  }
  
  # Create working copy for imputation
  data_imputed <- data_mat
  
  # Set maximum iterations based on method
  # Note: Regression and Mean methods are non-iterative (single pass)
  max_iter <- if (method == "REML") {
    100
  } else if (method == "Yates" || method == "Healy") {
    50
  } else if (method == "Regression" || method == "Mean") {
    1  # Non-iterative
  } else {
    30
  }
  
  # Apply selected method
  if (method == "Yates") {
    # YATES METHOD: Traditional iterative formula
    # Missing value = (r*T + t*B - G) / ((r-1)(t-1))
    # where T=treatment total, B=block total, G=grand total
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        # Initial estimate: use column mean
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean
        
        # Yates iterative algorithm
        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]
          
          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]
            
            # Calculate totals excluding current missing value
            current_val <- data_imputed[idx, col]
            data_imputed[idx, col] <- 0
            
            T_total <- sum(data_imputed[gen_idx == g, col], na.rm = TRUE)
            B_total <- sum(data_imputed[rep_idx == r, col], na.rm = TRUE)
            G_total <- sum(data_imputed[, col], na.rm = TRUE)
            
            # Yates formula
            numerator <- repli * T_total + genotype * B_total - G_total
            denominator <- (repli - 1) * (genotype - 1)
            data_imputed[idx, col] <- numerator / denominator
          }
          
          # Check convergence
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
    
  } else if (method == "Healy") {
    # HEALY & WESTMACOTT METHOD: Weighted adjustment method
    # More stable than Yates for multiple missing values
    # Uses residual-based weights to account for missing pattern
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        # Initial estimate: use column mean
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean
        
        # Healy & Westmacott iterative algorithm
        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]
          
          # Calculate treatment and block means from available data
          treat_means <- numeric(genotype)
          block_means <- numeric(repli)
          treat_counts <- numeric(genotype)
          block_counts <- numeric(repli)
          
          for (i in seq_len(nobs)) {
            if (is.finite(data_mat[, col][i])) {
              g <- gen_idx[i]
              r <- rep_idx[i]
              treat_means[g] <- treat_means[g] + data_imputed[i, col]
              block_means[r] <- block_means[r] + data_imputed[i, col]
              treat_counts[g] <- treat_counts[g] + 1
              block_counts[r] <- block_counts[r] + 1
            }
          }
          
          # Complete means (including current estimates)
          for (i in seq_len(nobs)) {
            if (!is.finite(data_mat[, col][i])) {
              g <- gen_idx[i]
              r <- rep_idx[i]
              treat_means[g] <- treat_means[g] + data_imputed[i, col]
              block_means[r] <- block_means[r] + data_imputed[i, col]
              treat_counts[g] <- treat_counts[g] + 1
              block_counts[r] <- block_counts[r] + 1
            }
          }
          
          # Calculate means
          treat_means <- treat_means / pmax(treat_counts, 1)
          block_means <- block_means / pmax(block_counts, 1)
          grand_mean <- mean(data_imputed[, col])
          
          # Healy & Westmacott estimation with weights
          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]
            
            # Number of missing in same treatment and block
            n_miss_treat <- sum(!is.finite(data_mat[gen_idx == g, col]))
            n_miss_block <- sum(!is.finite(data_mat[rep_idx == r, col]))
            
            # Weights based on number of observations
            w_treat <- (repli - n_miss_treat) / repli
            w_block <- (genotype - n_miss_block) / genotype
            
            # Normalize weights
            w_sum <- w_treat + w_block
            if (w_sum > 0) {
              w_treat <- w_treat / w_sum
              w_block <- w_block / w_sum
            } else {
              w_treat <- 0.5
              w_block <- 0.5
            }
            
            # Weighted estimate
            treat_effect <- treat_means[g] - grand_mean
            block_effect <- block_means[r] - grand_mean
            
            data_imputed[idx, col] <- grand_mean + 
                                     w_treat * treat_effect + 
                                     w_block * block_effect
          }
          
          # Check convergence
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
    
  } else if (method == "Regression") {
    # REGRESSION METHOD: Linear model with treatment and block factors
    # Non-iterative single-pass estimation using complete observations
    # Fast and deterministic - uses QR decomposition for stability
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        # Identify complete observations
        complete_idx <- which(is.finite(data_mat[, col]))
        
        if (length(complete_idx) > (genotype + repli)) {
          # Build design matrix from complete observations
          n_complete <- length(complete_idx)
          
          # Extract complete data
          y_complete <- data_mat[complete_idx, col]
          treat_complete <- gen_idx[complete_idx]
          block_complete <- rep_idx[complete_idx]
          
          # Create treatment effects matrix (contrast coding)
          X_treat <- matrix(0, n_complete, genotype - 1)
          for (i in seq_len(n_complete)) {
            if (treat_complete[i] < genotype) {
              X_treat[i, treat_complete[i]] <- 1
            } else {
              X_treat[i, ] <- -1
            }
          }
          
          # Create block effects matrix (contrast coding)
          X_block <- matrix(0, n_complete, repli - 1)
          for (i in seq_len(n_complete)) {
            if (block_complete[i] < repli) {
              X_block[i, block_complete[i]] <- 1
            } else {
              X_block[i, ] <- -1
            }
          }
          
          # Combine: intercept + treatment + block effects
          X <- cbind(1, X_treat, X_block)
          
          # Fit linear model using QR decomposition
          qr_fit <- qr(X)
          
          if (qr_fit$rank == ncol(X)) {
            # Estimate coefficients
            beta <- qr.coef(qr_fit, y_complete)
            
            # Predict missing values
            for (idx in missing_idx) {
              # Build predictor vector for missing observation
              x_new <- numeric(ncol(X))
              x_new[1] <- 1  # Intercept
              
              # Treatment effect
              treat_idx <- gen_idx[idx]
              if (treat_idx < genotype) {
                x_new[1 + treat_idx] <- 1
              } else {
                x_new[2:genotype] <- -1
              }
              
              # Block effect
              block_idx <- rep_idx[idx]
              block_start <- genotype
              if (block_idx < repli) {
                x_new[block_start + block_idx] <- 1
              } else {
                x_new[(block_start + 1):(block_start + repli - 1)] <- -1
              }
              
              # Predicted value from regression
              data_imputed[idx, col] <- sum(x_new * beta, na.rm = TRUE)
            }
          } else {
            # Fallback if design matrix is rank deficient
            # Use simple mean imputation
            treat_means <- tapply(data_mat[complete_idx, col], 
                                 gen_idx[complete_idx], mean)
            block_means <- tapply(data_mat[complete_idx, col], 
                                 rep_idx[complete_idx], mean)
            grand_mean <- mean(data_mat[complete_idx, col])
            
            for (idx in missing_idx) {
              g <- gen_idx[idx]
              r <- rep_idx[idx]
              treat_eff <- if (g %in% names(treat_means)) treat_means[g] - grand_mean else 0
              block_eff <- if (r %in% names(block_means)) block_means[r] - grand_mean else 0
              data_imputed[idx, col] <- grand_mean + treat_eff + block_eff
            }
          }
        } else {
          # Not enough complete observations - use simple mean
          data_imputed[missing_idx, col] <- mean(data_mat[, col], na.rm = TRUE)
        }
      }
    }
    
  } else if (method == "Mean") {
    # MEAN SUBSTITUTION METHOD: Simple additive model
    # Non-iterative single-pass estimation using treatment and block means
    # Fastest method: missing = grand_mean + treatment_effect + block_effect
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        # Calculate means from available (non-missing) observations only
        complete_idx <- which(is.finite(data_mat[, col]))
        
        if (length(complete_idx) > 0) {
          y_complete <- data_mat[complete_idx, col]
          treat_complete <- gen_idx[complete_idx]
          block_complete <- rep_idx[complete_idx]
          
          # Calculate treatment means from complete observations
          treat_means <- tapply(y_complete, treat_complete, mean)
          
          # Calculate block means from complete observations
          block_means <- tapply(y_complete, block_complete, mean)
          
          # Grand mean from complete observations
          grand_mean <- mean(y_complete)
          
          # Estimate missing values using additive model
          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]
            
            # Get treatment effect (or 0 if no data for this treatment)
            treat_effect <- if (g %in% as.integer(names(treat_means))) {
              treat_means[as.character(g)] - grand_mean
            } else {
              0
            }
            
            # Get block effect (or 0 if no data for this block)
            block_effect <- if (r %in% as.integer(names(block_means))) {
              block_means[as.character(r)] - grand_mean
            } else {
              0
            }
            
            # Additive model: grand mean + treatment effect + block effect
            data_imputed[idx, col] <- grand_mean + treat_effect + block_effect
          }
        } else {
          # No complete observations - should not happen but handle gracefully
          data_imputed[missing_idx, col] <- 0
        }
      }
    }
    
  } else if (method == "REML") {
    # REML METHOD: Variance components with BLUP
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean
        
        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]
          
          # Estimate variance components
          complete_idx <- which(is.finite(data_mat[, col]))
          if (length(complete_idx) > 3) {
            y_complete <- data_imputed[complete_idx, col]
            treat_complete <- gen_idx[complete_idx]
            block_complete <- rep_idx[complete_idx]
            
            grand_mean <- mean(y_complete)
            treat_means <- tapply(y_complete, treat_complete, mean)
            block_means <- tapply(y_complete, block_complete, mean)
            
            ss_treat <- sum(tapply(y_complete, treat_complete, 
                                   function(x) length(x) * (mean(x) - grand_mean)^2))
            ss_block <- sum(tapply(y_complete, block_complete, 
                                   function(x) length(x) * (mean(x) - grand_mean)^2))
            
            residuals <- y_complete - treat_means[treat_complete] - 
                        block_means[block_complete] + grand_mean
            ss_error <- sum(residuals^2)
            
            df_treat <- genotype - 1
            df_block <- repli - 1
            df_error <- length(complete_idx) - genotype - repli + 1
            
            if (df_error > 0) {
              ms_treat <- ss_treat / df_treat
              ms_block <- ss_block / df_block
              ms_error <- ss_error / df_error
              
              var_error <- ms_error
              var_treat <- max(0, (ms_treat - ms_error) / repli)
              var_block <- max(0, (ms_block - ms_error) / genotype)
              
              total_var <- var_treat + var_block + var_error
              if (total_var > 0) {
                shrinkage <- var_treat / total_var
              } else {
                shrinkage <- 0.5
              }
            } else {
              shrinkage <- 0.5
            }
            
            # Predict missing values with shrinkage
            for (idx in missing_idx) {
              g <- gen_idx[idx]
              r <- rep_idx[idx]
              
              treat_mean_val <- if (g %in% names(treat_means)) treat_means[g] else grand_mean
              block_mean_val <- if (r %in% names(block_means)) block_means[r] else grand_mean
              
              predicted <- treat_mean_val + block_mean_val - grand_mean
              data_imputed[idx, col] <- shrinkage * predicted + (1 - shrinkage) * grand_mean
            }
          }
          
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
    
  } else {
    # BARTLETT (COVARIATE) METHOD: Uses ANCOVA with other traits as covariates
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean
        
        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]
          
          complete_idx <- which(is.finite(data_mat[, col]))
          
          if (length(complete_idx) > 3) {
            covariate_cols <- setdiff(seq_len(ncols), col)
            
            valid_covariates <- integer(0)
            for (cov_col in covariate_cols) {
              if (all(is.finite(data_imputed[complete_idx, cov_col]))) {
                valid_covariates <- c(valid_covariates, cov_col)
              }
            }
            
            if (length(valid_covariates) > 0) {
              treat_complete <- gen_idx[complete_idx]
              block_complete <- rep_idx[complete_idx]
              y_complete <- data_imputed[complete_idx, col]
              
              n_complete <- length(complete_idx)
              
              # Treatment effects (contrast coding)
              X_treat <- matrix(0, n_complete, genotype - 1)
              for (i in seq_len(n_complete)) {
                if (treat_complete[i] < genotype) {
                  X_treat[i, treat_complete[i]] <- 1
                } else {
                  X_treat[i, ] <- -1
                }
              }
              
              # Block effects (contrast coding)
              X_block <- matrix(0, n_complete, repli - 1)
              for (i in seq_len(n_complete)) {
                if (block_complete[i] < repli) {
                  X_block[i, block_complete[i]] <- 1
                } else {
                  X_block[i, ] <- -1
                }
              }
              
              # Covariate effects (centered)
              X_cov <- matrix(0, n_complete, length(valid_covariates))
              cov_means <- numeric(length(valid_covariates))
              for (k in seq_along(valid_covariates)) {
                cov_col <- valid_covariates[k]
                cov_values <- data_imputed[complete_idx, cov_col]
                cov_means[k] <- mean(cov_values)
                X_cov[, k] <- cov_values - cov_means[k]
              }
              
              X <- cbind(1, X_treat, X_block, X_cov)
              
              # Fit using QR decomposition
              qr_fit <- qr(X)
              if (qr_fit$rank == ncol(X)) {
                beta <- qr.coef(qr_fit, y_complete)
                
                # Predict missing values
                for (idx in missing_idx) {
                  x_new <- numeric(ncol(X))
                  x_new[1] <- 1
                  
                  treat_idx <- gen_idx[idx]
                  if (treat_idx < genotype) {
                    x_new[1 + treat_idx] <- 1
                  } else {
                    x_new[2:genotype] <- -1
                  }
                  
                  block_idx <- rep_idx[idx]
                  block_start <- genotype
                  if (block_idx < repli) {
                    x_new[block_start + block_idx] <- 1
                  } else {
                    x_new[(block_start + 1):(block_start + repli - 1)] <- -1
                  }
                  
                  cov_start <- genotype + repli - 1
                  for (k in seq_along(valid_covariates)) {
                    cov_col <- valid_covariates[k]
                    if (is.finite(data_imputed[idx, cov_col])) {
                      x_new[cov_start + k] <- data_imputed[idx, cov_col] - cov_means[k]
                    }
                  }
                  
                  data_imputed[idx, col] <- sum(x_new * beta, na.rm = TRUE)
                }
              }
            } else {
              # No valid covariates: fall back to simple mean adjustment
              treat_means <- tapply(data_imputed[complete_idx, col], 
                                   gen_idx[complete_idx], mean)
              block_means <- tapply(data_imputed[complete_idx, col], 
                                   rep_idx[complete_idx], mean)
              grand_mean <- mean(data_imputed[complete_idx, col])
              
              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                treat_eff <- treat_means[g] - grand_mean
                block_eff <- block_means[r] - grand_mean
                data_imputed[idx, col] <- grand_mean + treat_eff + block_eff
              }
            }
          }
          
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
  }
  
  return(data_imputed)
}
