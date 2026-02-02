#' Missing Value Estimation for RCBD and Latin Square Designs
#'
#' Estimates missing values in Randomized Complete Block Design (RCBD) or 
#' Latin Square Design (LSD) using one of six methods: REML, Yates, Healy,
#' Regression, Mean, or Bartlett.
#'
#' @param data_mat Numeric matrix with observations (rows) by traits (columns).
#'   May contain missing values (NA, NaN, Inf).
#' @param gen_idx Integer vector indicating genotype/treatment index for each observation.
#' @param rep_idx Integer vector indicating replicate/block (RCBD) or row (LSD) index for each observation.
#' @param col_idx Integer vector indicating column index for each observation (required for LSD only).
#' @param design_type Character string specifying experimental design: "RCBD" (default) or "LSD" (Latin Square).
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
#' The function handles missing values in RCBD or LSD experiments by iteratively
#' estimating them until convergence. For RCBD, uses 2-way blocking (genotypes × blocks);
#' for LSD, uses 3-way blocking (genotypes × rows × columns). Each method has different strengths:
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
#' @keywords internal
missingValueEstimation <- function(data_mat, gen_idx, rep_idx, col_idx = NULL,
                                   design_type = c("RCBD", "LSD"),
                                   method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
                                   tolerance = 1e-6) {
  # Validate and match arguments
  design_type <- match.arg(design_type)
  method <- match.arg(method)
  
  # Validate LSD requirements
  if (design_type == "LSD" && is.null(col_idx)) {
    stop("Latin Square Design requires 'col_idx' parameter")
  }
  
  # Get dimensions
  nobs <- nrow(data_mat)
  ncols <- ncol(data_mat)
  genotype <- length(unique(gen_idx))
  repli <- length(unique(rep_idx))
  ncol_blocks <- if (!is.null(col_idx)) length(unique(col_idx)) else 0
  
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
    # RCBD: Missing value = (r*T + t*B - G) / ((r-1)(t-1))
    # LSD:  Missing value = (t*T + t*R + t*C - 2*G) / ((t-1)(t-2))
    # where T=treatment total, B/R=block/row total, C=column total, G=grand total
    
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
            R_total <- sum(data_imputed[rep_idx == r, col], na.rm = TRUE)
            G_total <- sum(data_imputed[, col], na.rm = TRUE)
            
            if (design_type == "RCBD") {
              # RCBD Yates formula
              numerator <- repli * T_total + genotype * R_total - G_total
              denominator <- (repli - 1) * (genotype - 1)
            } else {
              # LSD Yates formula
              c <- col_idx[idx]
              C_total <- sum(data_imputed[col_idx == c, col], na.rm = TRUE)
              numerator <- genotype * T_total + genotype * R_total + genotype * C_total - 2 * G_total
              denominator <- (genotype - 1) * (genotype - 2)
            }
            
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
          
          # Calculate treatment, row/block, and column means from available data
          treat_means <- numeric(genotype)
          block_means <- numeric(repli)
          treat_counts <- numeric(genotype)
          block_counts <- numeric(repli)
          
          # For LSD, also track column means
          col_means <- if (design_type == "LSD") numeric(ncol_blocks) else NULL
          col_counts <- if (design_type == "LSD") numeric(ncol_blocks) else NULL
          
          for (i in seq_len(nobs)) {
            if (is.finite(data_mat[, col][i])) {
              g <- gen_idx[i]
              r <- rep_idx[i]
              treat_means[g] <- treat_means[g] + data_imputed[i, col]
              block_means[r] <- block_means[r] + data_imputed[i, col]
              treat_counts[g] <- treat_counts[g] + 1
              block_counts[r] <- block_counts[r] + 1
              
              if (design_type == "LSD") {
                c <- col_idx[i]
                col_means[c] <- col_means[c] + data_imputed[i, col]
                col_counts[c] <- col_counts[c] + 1
              }
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
              
              if (design_type == "LSD") {
                c <- col_idx[i]
                col_means[c] <- col_means[c] + data_imputed[i, col]
                col_counts[c] <- col_counts[c] + 1
              }
            }
          }
          
          # Calculate means
          treat_means <- treat_means / pmax(treat_counts, 1)
          block_means <- block_means / pmax(block_counts, 1)
          if (design_type == "LSD") {
            col_means <- col_means / pmax(col_counts, 1)
          }
          grand_mean <- mean(data_imputed[, col])
          
          # Healy & Westmacott estimation with weights
          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]
            
            # Number of missing in same treatment and block/row
            n_miss_treat <- sum(!is.finite(data_mat[gen_idx == g, col]))
            n_miss_block <- sum(!is.finite(data_mat[rep_idx == r, col]))
            
            if (design_type == "RCBD") {
              # RCBD: 2-way weighting
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
            } else {
              # LSD: 3-way weighting
              c <- col_idx[idx]
              n_miss_col <- sum(!is.finite(data_mat[col_idx == c, col]))
              
              w_treat <- (genotype - n_miss_treat) / genotype
              w_row <- (genotype - n_miss_block) / genotype
              w_col <- (genotype - n_miss_col) / genotype
              
              # Normalize weights
              w_sum <- w_treat + w_row + w_col
              if (w_sum > 0) {
                w_treat <- w_treat / w_sum
                w_row <- w_row / w_sum
                w_col <- w_col / w_sum
              } else {
                w_treat <- 1/3
                w_row <- 1/3
                w_col <- 1/3
              }
              
              # Weighted estimate
              treat_effect <- treat_means[g] - grand_mean
              row_effect <- block_means[r] - grand_mean
              col_effect <- col_means[c] - grand_mean
              
              data_imputed[idx, col] <- grand_mean + 
                                       w_treat * treat_effect + 
                                       w_row * row_effect +
                                       w_col * col_effect
            }
          }
          
          # Check convergence
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
    
  } else if (method == "Regression") {
    # REGRESSION METHOD: Linear model with treatment and block factors
    # RCBD: Y ~ Treatment + Block (2-way)
    # LSD:  Y ~ Treatment + Row + Column (3-way)
    # Non-iterative single-pass estimation using complete observations
    # Fast and deterministic - uses QR decomposition for stability
    
    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      
      if (length(missing_idx) > 0) {
        # Identify complete observations
        complete_idx <- which(is.finite(data_mat[, col]))
        
        # Determine minimum required observations based on design
        min_obs <- if (design_type == "RCBD") {
          genotype + repli
        } else {
          genotype + repli + ncol_blocks
        }
        
        if (length(complete_idx) > min_obs) {
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
          
          # Create row/block effects matrix (contrast coding)
          X_block <- matrix(0, n_complete, repli - 1)
          for (i in seq_len(n_complete)) {
            if (block_complete[i] < repli) {
              X_block[i, block_complete[i]] <- 1
            } else {
              X_block[i, ] <- -1
            }
          }
          
          # For LSD, add column effects
          if (design_type == "LSD") {
            col_complete <- col_idx[complete_idx]
            X_col <- matrix(0, n_complete, ncol_blocks - 1)
            for (i in seq_len(n_complete)) {
              if (col_complete[i] < ncol_blocks) {
                X_col[i, col_complete[i]] <- 1
              } else {
                X_col[i, ] <- -1
              }
            }
            # Combine: intercept + treatment + row + column effects
            X <- cbind(1, X_treat, X_block, X_col)
          } else {
            # RCBD: Combine intercept + treatment + block effects
            X <- cbind(1, X_treat, X_block)
          }
          
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
              
              # Row/Block effect
              block_idx <- rep_idx[idx]
              block_start <- genotype
              if (block_idx < repli) {
                x_new[block_start + block_idx] <- 1
              } else {
                x_new[(block_start + 1):(block_start + repli - 1)] <- -1
              }
              
              # Column effect (LSD only)
              if (design_type == "LSD") {
                col_idx_val <- col_idx[idx]
                col_start <- genotype + repli - 1
                if (col_idx_val < ncol_blocks) {
                  x_new[col_start + col_idx_val] <- 1
                } else {
                  x_new[(col_start + 1):(col_start + ncol_blocks - 1)] <- -1
                }
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
            
            if (design_type == "LSD") {
              col_means <- tapply(data_mat[complete_idx, col], 
                                 col_idx[complete_idx], mean)
              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                c <- col_idx[idx]
                treat_eff <- if (g %in% names(treat_means)) treat_means[g] - grand_mean else 0
                row_eff <- if (r %in% names(block_means)) block_means[r] - grand_mean else 0
                col_eff <- if (c %in% names(col_means)) col_means[c] - grand_mean else 0
                data_imputed[idx, col] <- grand_mean + treat_eff + row_eff + col_eff
              }
            } else {
              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                treat_eff <- if (g %in% names(treat_means)) treat_means[g] - grand_mean else 0
                block_eff <- if (r %in% names(block_means)) block_means[r] - grand_mean else 0
                data_imputed[idx, col] <- grand_mean + treat_eff + block_eff
              }
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
    # RCBD: missing = grand_mean + treatment_effect + block_effect
    # LSD:  missing = grand_mean + treatment_effect + row_effect + column_effect
    # Fastest method
    
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
          
          # Calculate row/block means from complete observations
          block_means <- tapply(y_complete, block_complete, mean)
          
          # Grand mean from complete observations
          grand_mean <- mean(y_complete)
          
          if (design_type == "LSD") {
            # For LSD, also calculate column means
            col_complete <- col_idx[complete_idx]
            col_means <- tapply(y_complete, col_complete, mean)
            
            # Estimate missing values using 3-way additive model
            for (idx in missing_idx) {
              g <- gen_idx[idx]
              r <- rep_idx[idx]
              c <- col_idx[idx]
              
              # Get treatment effect (or 0 if no data for this treatment)
              treat_effect <- if (g %in% as.integer(names(treat_means))) {
                treat_means[as.character(g)] - grand_mean
              } else {
                0
              }
              
              # Get row effect (or 0 if no data for this row)
              row_effect <- if (r %in% as.integer(names(block_means))) {
                block_means[as.character(r)] - grand_mean
              } else {
                0
              }
              
              # Get column effect (or 0 if no data for this column)
              col_effect <- if (c %in% as.integer(names(col_means))) {
                col_means[as.character(c)] - grand_mean
              } else {
                0
              }
              
              # Additive model: grand mean + treatment effect + row effect + column effect
              data_imputed[idx, col] <- grand_mean + treat_effect + row_effect + col_effect
            }
          } else {
            # RCBD: 2-way additive model
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
          }
        } else {
          # No complete observations - should not happen but handle gracefully
          data_imputed[missing_idx, col] <- 0
        }
      }
    }
    
  } else if (method == "REML") {
    # REML METHOD: Variance components with BLUP
    # RCBD: Estimates treatment and block variance components
    # LSD:  Estimates treatment, row, and column variance components
    
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
            
            if (design_type == "RCBD") {
              # RCBD variance components
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
            } else {
              # LSD variance components
              col_complete <- col_idx[complete_idx]
              col_means <- tapply(y_complete, col_complete, mean)
              
              ss_col <- sum(tapply(y_complete, col_complete, 
                                   function(x) length(x) * (mean(x) - grand_mean)^2))
              
              residuals <- y_complete - treat_means[treat_complete] - 
                          block_means[block_complete] - col_means[col_complete] + 2 * grand_mean
              ss_error <- sum(residuals^2)
              
              df_treat <- genotype - 1
              df_row <- genotype - 1
              df_col <- genotype - 1
              df_error <- (genotype - 1) * (genotype - 2)
              
              if (df_error > 0 && length(complete_idx) > (3 * genotype - 2)) {
                ms_treat <- ss_treat / df_treat
                ms_row <- ss_block / df_row
                ms_col <- ss_col / df_col
                ms_error <- ss_error / df_error
                
                var_error <- ms_error
                var_treat <- max(0, (ms_treat - ms_error) / genotype)
                var_row <- max(0, (ms_row - ms_error) / genotype)
                var_col <- max(0, (ms_col - ms_error) / genotype)
                
                total_var <- var_treat + var_row + var_col + var_error
                if (total_var > 0) {
                  shrinkage <- var_treat / total_var
                } else {
                  shrinkage <- 0.33
                }
              } else {
                shrinkage <- 0.33
              }
              
              # Predict missing values with shrinkage (3-way)
              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                c <- col_idx[idx]
                
                treat_mean_val <- if (g %in% names(treat_means)) treat_means[g] else grand_mean
                row_mean_val <- if (r %in% names(block_means)) block_means[r] else grand_mean
                col_mean_val <- if (c %in% names(col_means)) col_means[c] else grand_mean
                
                predicted <- treat_mean_val + row_mean_val + col_mean_val - 2 * grand_mean
                data_imputed[idx, col] <- shrinkage * predicted + (1 - shrinkage) * grand_mean
              }
            }
          }
          
          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
    
  } else {
    # BARTLETT (COVARIATE) METHOD: Uses ANCOVA with other traits as covariates
    # RCBD: Y ~ Treatment + Block + Covariates
    # LSD:  Y ~ Treatment + Row + Column + Covariates
    
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
              
              # Row/Block effects (contrast coding)
              X_block <- matrix(0, n_complete, repli - 1)
              for (i in seq_len(n_complete)) {
                if (block_complete[i] < repli) {
                  X_block[i, block_complete[i]] <- 1
                } else {
                  X_block[i, ] <- -1
                }
              }
              
              # For LSD, add column effects
              if (design_type == "LSD") {
                col_complete <- col_idx[complete_idx]
                X_col <- matrix(0, n_complete, ncol_blocks - 1)
                for (i in seq_len(n_complete)) {
                  if (col_complete[i] < ncol_blocks) {
                    X_col[i, col_complete[i]] <- 1
                  } else {
                    X_col[i, ] <- -1
                  }
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
              
              # Combine design matrix based on design type
              if (design_type == "LSD") {
                X <- cbind(1, X_treat, X_block, X_col, X_cov)
              } else {
                X <- cbind(1, X_treat, X_block, X_cov)
              }
              
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
                  
                  # Column effects for LSD
                  if (design_type == "LSD") {
                    col_idx_val <- col_idx[idx]
                    col_start <- genotype + repli - 1
                    if (col_idx_val < ncol_blocks) {
                      x_new[col_start + col_idx_val] <- 1
                    } else {
                      x_new[(col_start + 1):(col_start + ncol_blocks - 1)] <- -1
                    }
                    cov_start <- genotype + repli + ncol_blocks - 2
                  } else {
                    cov_start <- genotype + repli - 1
                  }
                  
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
              
              if (design_type == "LSD") {
                col_means <- tapply(data_imputed[complete_idx, col], 
                                   col_idx[complete_idx], mean)
                for (idx in missing_idx) {
                  g <- gen_idx[idx]
                  r <- rep_idx[idx]
                  c <- col_idx[idx]
                  treat_eff <- treat_means[g] - grand_mean
                  row_eff <- block_means[r] - grand_mean
                  col_eff <- col_means[c] - grand_mean
                  data_imputed[idx, col] <- grand_mean + treat_eff + row_eff + col_eff
                }
              } else {
                for (idx in missing_idx) {
                  g <- gen_idx[idx]
                  r <- rep_idx[idx]
                  treat_eff <- treat_means[g] - grand_mean
                  block_eff <- block_means[r] - grand_mean
                  data_imputed[idx, col] <- grand_mean + treat_eff + block_eff
                }
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
