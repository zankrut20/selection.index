#' Stochastic Simulation of Selection Indices (Chapter 10)
#' @name stochastic_simulation
#'
#' @description
#' Implements stochastic simulation methods for evaluating the long-term
#' performance of selection indices over multiple breeding cycles.
#' These methods model linkage between loci using Haldane's mapping function
#' and track genetic gain across generations.
#'
#' Methods included:
#' - Haldane's mapping function for recombination
#' - Multi-cycle simulation framework
#' - Comparisons of LPSI, ESIM, RLPSI, and RESIM over time
#'
#' @references
#' Haldane, J. B. S. (1919). The combination of linkage values and the calculation
#' of distances between the loci of linked factors. Journal of Genetics, 8(29), 299-309.
#'
#' Smith, H. F. (1936). A discriminant function for plant selection.
#' Annals of Eugenics, 7(3), 240-250.
#'
#' Hazel, L. N. (1943). The genetic basis for constructing selection indexes.
#' Genetics, 28(6), 476.
#'
#' Kempthorne, O., & Nordskog, A. W. (1959). Restricted selection indices.
#' Biometrics, 15(1), 10-19.
#'
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapter 7 & 10.
#'
#' @keywords internal
#' @importFrom stats rnorm runif rbinom var
#' @importFrom graphics lines legend
NULL

# ==============================================================================
# HALDANE'S MAPPING FUNCTION
# ==============================================================================

#' Haldane's Mapping Function
#'
#' @description
#' Converts genetic distance (in Morgans) to recombination fraction using
#' Haldane's mapping function. This function models the relationship between
#' genetic distance and the probability of recombination between loci.
#'
#' @param distance Genetic distance in Morgans (scalar or vector).
#'   One Morgan corresponds to a 50\% recombination frequency.
#'
#' @return Recombination fraction (r) ranging from 0 to 0.5.
#'   - r = 0 indicates complete linkage (no recombination)
#'   - r = 0.5 indicates independent assortment (unlinked loci)
#'
#' @details
#' \strong{Mathematical Formula (Chapter 10, Section 10.1):}
#'
#' The relationship between recombination fraction (r) and genetic distance (d):
#' \deqn{r = \frac{1}{2}(1 - e^{-2d})}
#'
#' Where:
#' - d = Genetic distance in Morgans
#' - r = Recombination fraction (probability of recombination per meiosis)
#'
#' This function assumes no crossover interference beyond that implied
#' by the mapping function itself.
#'
#' @export
#' @examples
#' # Zero distance means complete linkage (no recombination)
#' haldane_mapping(0) # Returns 0
#'
#' # 1 Morgan distance
#' haldane_mapping(1) # Returns ~0.43
#'
#' # Large distance approaches 0.5 (independent assortment)
#' haldane_mapping(10) # Returns ~0.5
#'
#' # Vector of distances
#' distances <- c(0, 0.1, 0.5, 1.0, 2.0)
#' haldane_mapping(distances)
haldane_mapping <- function(distance) {
  if (!is.numeric(distance) || any(distance < 0)) {
    stop("distance must be non-negative numeric value(s)")
  }

  # Haldane's formula: r = (1/2)(1 - e^(-2d))
  recombination_fraction <- 0.5 * (1 - exp(-2 * distance))

  return(recombination_fraction)
}


#' Inverse Haldane Mapping Function
#'
#' @description
#' Converts recombination fraction back to genetic distance (in Morgans).
#' This is the inverse of Haldane's mapping function.
#'
#' @param recombination_fraction Recombination fraction (r) between 0 and 0.5.
#'
#' @return Genetic distance in Morgans.
#'
#' @details
#' \strong{Mathematical Formula:}
#'
#' Solving Haldane's equation for d:
#' \deqn{d = -\frac{1}{2} \ln(1 - 2r)}
#'
#' @export
#' @examples
#' # Convert recombination fraction to distance
#' inverse_haldane_mapping(0.25) # Returns ~0.347 Morgans
#' inverse_haldane_mapping(0.5) # Returns Inf (unlinked)
inverse_haldane_mapping <- function(recombination_fraction) {
  if (!is.numeric(recombination_fraction) ||
    any(recombination_fraction < 0) ||
    any(recombination_fraction > 0.5)) {
    stop("recombination_fraction must be between 0 and 0.5")
  }

  # Handle special case of r = 0.5 (infinite distance)
  if (any(recombination_fraction == 0.5)) {
    warning("recombination_fraction = 0.5 corresponds to infinite genetic distance")
  }

  # Inverse Haldane: d = -(1/2) * ln(1 - 2r)
  distance <- -0.5 * log(1 - 2 * recombination_fraction)

  return(distance)
}


# ==============================================================================
# SIMULATION HELPER FUNCTIONS
# ==============================================================================

#' Create Initial Population for Simulation
#'
#' @description
#' Creates a simulated diploid population with quantitative trait loci (QTL)
#' and marker information. Each individual has two haplotypes with additive
#' genetic effects.
#'
#' @param n_individuals Number of individuals in the population
#' @param n_loci Number of QTL per trait
#' @param n_traits Number of quantitative traits
#' @param qtl_effects Optional matrix of QTL effects (n_loci x n_traits).
#'   If NULL, effects are sampled from standard normal distribution.
#' @param genetic_distances Optional vector of genetic distances between adjacent loci (in Morgans).
#'   If NULL, assumes 0.1 Morgan spacing (roughly 10 cM).
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{haplotype1} - First haplotype matrix (n_individuals x n_loci)
#'     \item \code{haplotype2} - Second haplotype matrix (n_individuals x n_loci)
#'     \item \code{qtl_effects} - QTL effect matrix (n_loci x n_traits)
#'     \item \code{genetic_distances} - Distances between loci (in Morgans)
#'     \item \code{recombination_fractions} - Recombination fractions between adjacent loci
#'   }
#'
#' @keywords internal
#' @noRd
.create_initial_population <- function(n_individuals, n_loci, n_traits,
                                       qtl_effects = NULL,
                                       genetic_distances = NULL) {
  # Set default QTL effects (standardized)
  if (is.null(qtl_effects)) {
    qtl_effects <- matrix(rnorm(n_loci * n_traits, mean = 0, sd = 1),
      nrow = n_loci, ncol = n_traits
    )
  }

  # Set default genetic distances (0.1 Morgan = 10 cM between adjacent loci)
  if (is.null(genetic_distances)) {
    genetic_distances <- rep(0.1, n_loci - 1)
  }

  # Calculate recombination fractions using Haldane's function
  recombination_fractions <- haldane_mapping(genetic_distances)

  # Initialize haplotypes (0 or 1 alleles, representing allele dosage)
  haplotype1 <- matrix(rbinom(n_individuals * n_loci, 1, 0.5),
    nrow = n_individuals, ncol = n_loci
  )
  haplotype2 <- matrix(rbinom(n_individuals * n_loci, 1, 0.5),
    nrow = n_individuals, ncol = n_loci
  )

  list(
    haplotype1 = haplotype1,
    haplotype2 = haplotype2,
    qtl_effects = qtl_effects,
    genetic_distances = genetic_distances,
    recombination_fractions = recombination_fractions
  )
}


#' Recombine Two Haplotypes
#'
#' @description
#' Simulates meiotic recombination between two parental haplotypes using
#' the recombination fractions derived from Haldane's mapping function.
#'
#' @param haplotype1 First parental haplotype (vector of length n_loci)
#' @param haplotype2 Second parental haplotype (vector of length n_loci)
#' @param recombination_fractions Vector of recombination fractions between
#'   adjacent loci (length n_loci - 1)
#'
#' @return Recombinant haplotype (vector of length n_loci)
#'
#' @keywords internal
#' @noRd
.recombine_haplotypes <- function(haplotype1, haplotype2, recombination_fractions) {
  n_loci <- length(haplotype1)

  # Start with first haplotype
  gamete <- numeric(n_loci)
  current_parent <- 1 # Start with haplotype 1

  gamete[1] <- if (current_parent == 1) haplotype1[1] else haplotype2[1]

  # Process each locus
  for (i in 2:n_loci) {
    # Check if recombination occurs
    if (runif(1) < recombination_fractions[i - 1]) {
      # Recombination: switch parent
      current_parent <- 3 - current_parent # Toggle between 1 and 2
    }

    # Assign allele from current parent
    gamete[i] <- if (current_parent == 1) haplotype1[i] else haplotype2[i]
  }

  return(gamete)
}


#' Compute Genetic Values for Population
#'
#' @description
#' Calculates breeding values (genetic values) for all individuals in the
#' population based on their diploid genotypes and QTL effects.
#'
#' @param haplotype1 First haplotype matrix (n_individuals x n_loci)
#' @param haplotype2 Second haplotype matrix (n_individuals x n_loci)
#' @param qtl_effects QTL effect matrix (n_loci x n_traits)
#'
#' @return Matrix of genetic values (n_individuals x n_traits)
#'
#' @keywords internal
#' @noRd
.compute_genetic_values <- function(haplotype1, haplotype2, qtl_effects) {
  # Diploid genotype (allele dosage: 0, 1, or 2)
  genotype <- haplotype1 + haplotype2

  # Genetic value: G = X * beta (genotype matrix x effect matrix)
  genetic_values <- genotype %*% qtl_effects

  return(genetic_values)
}


#' Compute Phenotypic Values with Environmental Noise
#'
#' @description
#' Adds environmental variation to genetic values to create phenotypic observations.
#'
#' **IMPORTANT**: Environmental variance should be calculated ONCE at generation 0
#' and held constant across all cycles. Recalculating var_e each generation
#' artificially maintains heritability and overestimates genetic gain (fails to
#' model the Bulmer effect correctly).
#'
#' @param genetic_values Matrix of genetic values (n_individuals x n_traits)
#' @param environmental_variance Vector of environmental variances (constant across cycles)
#'
#' @return Matrix of phenotypic values (n_individuals x n_traits)
#'
#' @keywords internal
#' @noRd
.compute_phenotypes <- function(genetic_values, environmental_variance) {
  n_individuals <- nrow(genetic_values)
  n_traits <- ncol(genetic_values)

  # Ensure environmental_variance is a vector matching n_traits
  if (length(environmental_variance) == 1) {
    environmental_variance <- rep(environmental_variance, n_traits)
  }

  # Add constant environmental noise
  # Environmental variance remains fixed; genetic variance decreases due to selection
  # This naturally causes heritability to decline (Bulmer effect)
  phenotypes <- genetic_values +
    matrix(rnorm(n_individuals * n_traits), nrow = n_individuals) %*% diag(sqrt(environmental_variance))

  return(phenotypes)
}


#' Generate Offspring Population Through Recombination
#'
#' @description
#' Creates the next generation by randomly mating selected parents and
#' simulating meiotic recombination.
#'
#' @param selected_parents List containing haplotypes and metadata of selected individuals
#' @param n_offspring Number of offspring to generate
#'
#' @return List with offspring haplotypes
#'
#' @keywords internal
#' @noRd
.generate_offspring <- function(selected_parents, n_offspring) {
  n_parents <- nrow(selected_parents$haplotype1)
  n_loci <- ncol(selected_parents$haplotype1)
  recomb_frac <- selected_parents$recombination_fractions

  # Pre-allocate offspring matrices
  offspring_hap1 <- matrix(0, nrow = n_offspring, ncol = n_loci)
  offspring_hap2 <- matrix(0, nrow = n_offspring, ncol = n_loci)

  for (i in 1:n_offspring) {
    # Randomly select two parents
    parent1_idx <- sample(n_parents, 1)
    parent2_idx <- sample(n_parents, 1)

    # Parent 1 contributes one gamete (recombined)
    gamete1 <- .recombine_haplotypes(
      selected_parents$haplotype1[parent1_idx, ],
      selected_parents$haplotype2[parent1_idx, ],
      recomb_frac
    )

    # Parent 2 contributes one gamete (recombined)
    gamete2 <- .recombine_haplotypes(
      selected_parents$haplotype1[parent2_idx, ],
      selected_parents$haplotype2[parent2_idx, ],
      recomb_frac
    )

    offspring_hap1[i, ] <- gamete1
    offspring_hap2[i, ] <- gamete2
  }

  list(
    haplotype1 = offspring_hap1,
    haplotype2 = offspring_hap2,
    qtl_effects = selected_parents$qtl_effects,
    genetic_distances = selected_parents$genetic_distances,
    recombination_fractions = recomb_frac
  )
}


# ==============================================================================
# MAIN SIMULATION FUNCTION
# ==============================================================================

#' Simulate Multi-Cycle Selection Using Different Indices
#'
#' @description
#' Performs a stochastic simulation comparing the long-term performance of
#' four selection indices (LPSI, ESIM, RLPSI, RESIM) over multiple breeding cycles.
#' Models linkage between loci using Haldane's mapping function.
#'
#' @param n_cycles Number of selection cycles to simulate (default: 50)
#' @param n_individuals Initial population size (default: 1000)
#' @param n_loci Number of QTL per trait (default: 100)
#' @param n_traits Number of quantitative traits (default: 3)
#' @param heritability Heritability for all traits (scalar or vector, default: 0.5)
#' @param selection_proportion Proportion of individuals selected (default: 0.1)
#' @param economic_weights Economic weights for LPSI (default: equal weights)
#' @param restricted_traits Trait indices to restrict to zero gain for RLPSI/RESIM (default: NULL)
#' @param genetic_distances Vector of genetic distances between loci in Morgans (default: 0.1)
#' @param qtl_effects Optional matrix of QTL effects (n_loci x n_traits)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{lpsi_gain} - Matrix of genetic gains per cycle for LPSI (n_cycles x n_traits)
#'     \item \code{esim_gain} - Matrix of genetic gains per cycle for ESIM
#'     \item \code{rlpsi_gain} - Matrix of genetic gains per cycle for RLPSI
#'     \item \code{resim_gain} - Matrix of genetic gains per cycle for RESIM
#'     \item \code{lpsi_mean} - Mean genetic value per cycle for LPSI (n_cycles x n_traits)
#'     \item \code{esim_mean} - Mean genetic value per cycle for ESIM
#'     \item \code{rlpsi_mean} - Mean genetic value per cycle for RLPSI
#'     \item \code{resim_mean} - Mean genetic value per cycle for RESIM
#'     \item \code{parameters} - List of simulation parameters used
#'   }
#'
#' @details
#' \strong{Simulation Procedure (Chapter 10):}
#'
#' 1. Initialize population with random QTL alleles at linked loci
#' 2. For each cycle:
#'    - Compute genetic values from diploid genotypes
#'    - Add environmental noise to create phenotypes
#'    - Calculate variance-covariance matrices
#'    - Apply each selection index (LPSI, ESIM, RLPSI, RESIM)
#'    - Select top individuals based on index scores
#'    - Generate offspring through recombination (using Haldane's function)
#' 3. Track genetic gain and population mean across cycles
#'
#' \strong{The Four Indices (Chapter 10, Section 10.2):}
#'
#' \itemize{
#'   \item \strong{LPSI}: Unrestricted index maximizing genetic gain: b = P^(-1)Gw
#'   \item \strong{ESIM}: Unrestricted eigen index maximizing accuracy: (P^(-1)C - lambda^2 I)b = 0
#'   \item \strong{RLPSI}: Restricted LPSI with constraints: U'Gb = 0
#'   \item \strong{RESIM}: Restricted ESIM with constraints: U'Cb = 0
#' }
#'
#' \strong{Important Simulation Assumptions:}
#'
#' \strong{Note:} Environmental variance is calculated once at generation 0
#' and held constant across all cycles. Heritability will naturally decline
#' over cycles as genetic variance is depleted (Bulmer effect). This models
#' the biological reality that selection exhausts additive genetic variance
#' while environmental variation remains stable.
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic simulation with 3 traits over 20 cycles
#' results <- simulate_selection_cycles(
#'   n_cycles = 20,
#'   n_individuals = 500,
#'   n_loci = 50,
#'   n_traits = 3,
#'   heritability = 0.5,
#'   selection_proportion = 0.1,
#'   economic_weights = c(10, 5, 3),
#'   seed = 123
#' )
#'
#' # Plot genetic gains
#' plot(1:20, results$lpsi_mean[, 1],
#'   type = "l", col = "blue",
#'   ylab = "Mean Genetic Value", xlab = "Cycle",
#'   main = "Genetic Gain - Trait 1"
#' )
#' lines(1:20, results$esim_mean[, 1], col = "red")
#' lines(1:20, results$rlpsi_mean[, 1], col = "green")
#' lines(1:20, results$resim_mean[, 1], col = "orange")
#' legend("topleft", c("LPSI", "ESIM", "RLPSI", "RESIM"),
#'   col = c("blue", "red", "green", "orange"), lty = 1
#' )
#'
#' # Restrict trait 2 to zero gain
#' results_restricted <- simulate_selection_cycles(
#'   n_cycles = 20,
#'   n_traits = 3,
#'   restricted_traits = 2,
#'   seed = 456
#' )
#' }
simulate_selection_cycles <- function(n_cycles = 50,
                                      n_individuals = 1000,
                                      n_loci = 100,
                                      n_traits = 3,
                                      heritability = 0.5,
                                      selection_proportion = 0.1,
                                      economic_weights = NULL,
                                      restricted_traits = NULL,
                                      genetic_distances = NULL,
                                      qtl_effects = NULL,
                                      seed = NULL) {
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (n_cycles < 1) stop("n_cycles must be at least 1")
  if (n_individuals < 10) stop("n_individuals must be at least 10")
  if (n_loci < 1) stop("n_loci must be at least 1")
  if (n_traits < 1) stop("n_traits must be at least 1")
  if (any(heritability <= 0) || any(heritability >= 1)) {
    stop("heritability must be between 0 and 1")
  }
  if (selection_proportion <= 0 || selection_proportion >= 1) {
    stop("selection_proportion must be between 0 and 1")
  }

  # Set default economic weights
  if (is.null(economic_weights)) {
    economic_weights <- rep(1, n_traits)
  }

  # Calculate number of parents to select
  n_selected <- max(2, round(n_individuals * selection_proportion))

  # Initialize storage for results
  lpsi_mean <- matrix(0, nrow = n_cycles, ncol = n_traits)
  esim_mean <- matrix(0, nrow = n_cycles, ncol = n_traits)
  rlpsi_mean <- matrix(0, nrow = n_cycles, ncol = n_traits)
  resim_mean <- matrix(0, nrow = n_cycles, ncol = n_traits)

  lpsi_gain <- matrix(0, nrow = n_cycles, ncol = n_traits)
  esim_gain <- matrix(0, nrow = n_cycles, ncol = n_traits)
  rlpsi_gain <- matrix(0, nrow = n_cycles, ncol = n_traits)
  resim_gain <- matrix(0, nrow = n_cycles, ncol = n_traits)

  # Create initial population (shared across all index methods)
  pop_base <- .create_initial_population(
    n_individuals = n_individuals,
    n_loci = n_loci,
    n_traits = n_traits,
    qtl_effects = qtl_effects,
    genetic_distances = genetic_distances
  )

  # Initialize separate populations for each index method
  pop_lpsi <- pop_base
  pop_esim <- pop_base
  pop_rlpsi <- pop_base
  pop_resim <- pop_base

  # ============================================================================
  # CRITICAL: Calculate environmental variance ONCE at Generation 0
  # ============================================================================
  # Environmental variance stays constant across cycles. Genetic variance
  # decreases due to selection (Bulmer effect), causing heritability to naturally
  # decline. Do NOT recalculate var_e each cycle.

  # Calculate initial genetic values
  g0 <- .compute_genetic_values(
    pop_base$haplotype1, pop_base$haplotype2,
    pop_base$qtl_effects
  )

  # Ensure heritability is a vector
  if (length(heritability) == 1) {
    heritability <- rep(heritability, n_traits)
  }

  # Calculate initial genetic variance
  var_g0 <- apply(g0, 2, var)

  # Calculate CONSTANT environmental variance from h² = σ²_G / (σ²_G + σ²_E)
  # Rearranging: σ²_E = σ²_G * (1 - h²) / h²
  environmental_variance <- var_g0 * (1 - heritability) / heritability

  # Run simulation cycles
  for (cycle in 1:n_cycles) {
    # ========================================================================
    # LPSI Selection
    # ========================================================================
    g_lpsi <- .compute_genetic_values(
      pop_lpsi$haplotype1, pop_lpsi$haplotype2,
      pop_lpsi$qtl_effects
    )
    p_lpsi <- .compute_phenotypes(g_lpsi, environmental_variance)

    # Calculate variance-covariance matrices
    G_lpsi <- cov(g_lpsi)
    P_lpsi <- cov(p_lpsi)

    # Calculate LPSI coefficients
    tryCatch(
      {
        b_lpsi <- cpp_symmetric_solve(P_lpsi, G_lpsi %*% economic_weights)

        # Calculate index scores
        scores_lpsi <- p_lpsi %*% b_lpsi

        # Select top individuals
        selected_idx <- order(scores_lpsi, decreasing = TRUE)[1:n_selected]

        # Store results
        lpsi_mean[cycle, ] <- colMeans(g_lpsi)
        if (cycle > 1) {
          lpsi_gain[cycle, ] <- lpsi_mean[cycle, ] - lpsi_mean[cycle - 1, ]
        }

        # Create next generation
        pop_lpsi_selected <- list(
          haplotype1 = pop_lpsi$haplotype1[selected_idx, , drop = FALSE],
          haplotype2 = pop_lpsi$haplotype2[selected_idx, , drop = FALSE],
          qtl_effects = pop_lpsi$qtl_effects,
          genetic_distances = pop_lpsi$genetic_distances,
          recombination_fractions = pop_lpsi$recombination_fractions
        )
        pop_lpsi <- .generate_offspring(pop_lpsi_selected, n_individuals)
      },
      error = function(e) {
        message("LPSI failed at cycle ", cycle, ": ", e$message)
        lpsi_mean[cycle, ] <- if (cycle > 1) lpsi_mean[cycle - 1, ] else colMeans(g_lpsi)
      }
    )

    # ========================================================================
    # ESIM Selection
    # ========================================================================
    g_esim <- .compute_genetic_values(
      pop_esim$haplotype1, pop_esim$haplotype2,
      pop_esim$qtl_effects
    )
    p_esim <- .compute_phenotypes(g_esim, environmental_variance)

    G_esim <- cov(g_esim)
    P_esim <- cov(p_esim)

    tryCatch(
      {
        # Solve eigenproblem: (P^-1 G - lambda^2 I)b = 0
        # Use ginv() for robustness to near-singular P as variances decline
        PinvG <- suppressWarnings(MASS::ginv(P_esim)) %*% G_esim
        eigen_result <- eigen(PinvG, symmetric = FALSE)

        # Check for significant imaginary parts (indicates numerical issues)
        max_imag <- max(abs(Im(eigen_result$values)))
        if (max_imag > 1e-5) {
          warning(sprintf(
            "ESIM cycle %d: Eigenvalues have significant imaginary parts (max = %.2e). Check covariance matrices.",
            cycle, max_imag
          ))
        }

        # Use leading eigenvector
        b_esim <- Re(eigen_result$vectors[, 1])

        # Calculate index scores
        scores_esim <- p_esim %*% b_esim

        # Select top individuals
        selected_idx <- order(scores_esim, decreasing = TRUE)[1:n_selected]

        # Store results
        esim_mean[cycle, ] <- colMeans(g_esim)
        if (cycle > 1) {
          esim_gain[cycle, ] <- esim_mean[cycle, ] - esim_mean[cycle - 1, ]
        }

        # Create next generation
        pop_esim_selected <- list(
          haplotype1 = pop_esim$haplotype1[selected_idx, , drop = FALSE],
          haplotype2 = pop_esim$haplotype2[selected_idx, , drop = FALSE],
          qtl_effects = pop_esim$qtl_effects,
          genetic_distances = pop_esim$genetic_distances,
          recombination_fractions = pop_esim$recombination_fractions
        )
        pop_esim <- .generate_offspring(pop_esim_selected, n_individuals)
      },
      error = function(e) {
        message("ESIM failed at cycle ", cycle, ": ", e$message)
        esim_mean[cycle, ] <- if (cycle > 1) esim_mean[cycle - 1, ] else colMeans(g_esim)
      }
    )

    # ========================================================================
    # RLPSI Selection (if restrictions specified)
    # ========================================================================
    g_rlpsi <- .compute_genetic_values(
      pop_rlpsi$haplotype1, pop_rlpsi$haplotype2,
      pop_rlpsi$qtl_effects
    )
    p_rlpsi <- .compute_phenotypes(g_rlpsi, environmental_variance)

    G_rlpsi <- cov(g_rlpsi)
    P_rlpsi <- cov(p_rlpsi)

    tryCatch(
      {
        if (!is.null(restricted_traits)) {
          # Create constraint matrix
          U <- diag(n_traits)[, restricted_traits, drop = FALSE]

          # Calculate RLPSI coefficients
          # b_R = [I - P^-1 G U (U' G P^-1 G U)^-1 U' G] P^-1 G w
          # Use ginv() for robustness to near-singular matrices
          PinvG <- suppressWarnings(MASS::ginv(P_rlpsi)) %*% G_rlpsi
          GPinvG <- G_rlpsi %*% PinvG

          UtGPinvGU <- t(U) %*% GPinvG %*% U
          middle_term <- U %*% suppressWarnings(MASS::ginv(UtGPinvGU)) %*% t(U) %*% G_rlpsi

          K <- diag(n_traits) - PinvG %*% middle_term
          b_rlpsi <- K %*% PinvG %*% economic_weights
        } else {
          # No restrictions: RLPSI = LPSI
          b_rlpsi <- cpp_symmetric_solve(P_rlpsi, G_rlpsi %*% economic_weights)
        }

        # Calculate index scores
        scores_rlpsi <- p_rlpsi %*% b_rlpsi

        # Select top individuals
        selected_idx <- order(scores_rlpsi, decreasing = TRUE)[1:n_selected]

        # Store results
        rlpsi_mean[cycle, ] <- colMeans(g_rlpsi)
        if (cycle > 1) {
          rlpsi_gain[cycle, ] <- rlpsi_mean[cycle, ] - rlpsi_mean[cycle - 1, ]
        }

        # Create next generation
        pop_rlpsi_selected <- list(
          haplotype1 = pop_rlpsi$haplotype1[selected_idx, , drop = FALSE],
          haplotype2 = pop_rlpsi$haplotype2[selected_idx, , drop = FALSE],
          qtl_effects = pop_rlpsi$qtl_effects,
          genetic_distances = pop_rlpsi$genetic_distances,
          recombination_fractions = pop_rlpsi$recombination_fractions
        )
        pop_rlpsi <- .generate_offspring(pop_rlpsi_selected, n_individuals)
      },
      error = function(e) {
        message("RLPSI failed at cycle ", cycle, ": ", e$message)
        rlpsi_mean[cycle, ] <- if (cycle > 1) rlpsi_mean[cycle - 1, ] else colMeans(g_rlpsi)
      }
    )

    # ========================================================================
    # RESIM Selection (if restrictions specified)
    # ========================================================================
    g_resim <- .compute_genetic_values(
      pop_resim$haplotype1, pop_resim$haplotype2,
      pop_resim$qtl_effects
    )
    p_resim <- .compute_phenotypes(g_resim, environmental_variance)

    G_resim <- cov(g_resim)
    P_resim <- cov(p_resim)

    tryCatch(
      {
        if (!is.null(restricted_traits)) {
          # Create constraint matrix
          U <- as.matrix(diag(n_traits)[, restricted_traits, drop = FALSE])

          # Calculate projection matrix K (Chapter 7.2, Section 7.2.3)
          # K = I - P^-1 G U (U' G P^-1 G U)^-1 U' G
          # This ensures the restriction U' C b = 0 is enforced
          # Use ginv() for robustness to near-singular matrices
          PinvG <- suppressWarnings(MASS::ginv(P_resim)) %*% G_resim
          GPinvG <- G_resim %*% PinvG

          # Middle term: U (U' G P^-1 G U)^-1 U' G
          UtGPinvGU <- t(U) %*% GPinvG %*% U
          middle_term <- U %*% suppressWarnings(MASS::ginv(UtGPinvGU)) %*% t(U) %*% G_resim

          # Projection matrix (same as RLPSI but used differently)
          K <- diag(n_traits) - PinvG %*% middle_term

          # Solve restricted eigenproblem: (K P^-1 C - lambda^2 I)b = 0
          # where C = genotypic covariance matrix
          KPinvC <- K %*% PinvG
          eigen_result <- eigen(KPinvC, symmetric = FALSE)

          # Check for significant imaginary parts (indicates numerical issues)
          max_imag <- max(abs(Im(eigen_result$values)))
          if (max_imag > 1e-5) {
            warning(sprintf(
              "RESIM (restricted) cycle %d: Eigenvalues have significant imaginary parts (max = %.2e). Check covariance matrices.",
              cycle, max_imag
            ))
          }

          # Use leading eigenvector
          b_resim <- Re(eigen_result$vectors[, 1])
        } else {
          # No restrictions: RESIM = ESIM
          # Use ginv() for robustness to near-singular P
          PinvG <- suppressWarnings(MASS::ginv(P_resim)) %*% G_resim
          eigen_result <- eigen(PinvG, symmetric = FALSE)

          # Check for significant imaginary parts (indicates numerical issues)
          max_imag <- max(abs(Im(eigen_result$values)))
          if (max_imag > 1e-5) {
            warning(sprintf(
              "RESIM (unrestricted) cycle %d: Eigenvalues have significant imaginary parts (max = %.2e). Check covariance matrices.",
              cycle, max_imag
            ))
          }

          b_resim <- Re(eigen_result$vectors[, 1])
        }

        # Calculate index scores
        scores_resim <- p_resim %*% b_resim

        # Select top individuals
        selected_idx <- order(scores_resim, decreasing = TRUE)[1:n_selected]

        # Store results
        resim_mean[cycle, ] <- colMeans(g_resim)
        if (cycle > 1) {
          resim_gain[cycle, ] <- resim_mean[cycle, ] - resim_mean[cycle - 1, ]
        }

        # Create next generation
        pop_resim_selected <- list(
          haplotype1 = pop_resim$haplotype1[selected_idx, , drop = FALSE],
          haplotype2 = pop_resim$haplotype2[selected_idx, , drop = FALSE],
          qtl_effects = pop_resim$qtl_effects,
          genetic_distances = pop_resim$genetic_distances,
          recombination_fractions = pop_resim$recombination_fractions
        )
        pop_resim <- .generate_offspring(pop_resim_selected, n_individuals)
      },
      error = function(e) {
        message("RESIM failed at cycle ", cycle, ": ", e$message)
        resim_mean[cycle, ] <- if (cycle > 1) resim_mean[cycle - 1, ] else colMeans(g_resim)
      }
    )
  } # End of cycles loop

  # Return results
  structure(
    list(
      lpsi_gain = lpsi_gain,
      esim_gain = esim_gain,
      rlpsi_gain = rlpsi_gain,
      resim_gain = resim_gain,
      lpsi_mean = lpsi_mean,
      esim_mean = esim_mean,
      rlpsi_mean = rlpsi_mean,
      resim_mean = resim_mean,
      parameters = list(
        n_cycles = n_cycles,
        n_individuals = n_individuals,
        n_loci = n_loci,
        n_traits = n_traits,
        heritability = heritability,
        selection_proportion = selection_proportion,
        economic_weights = economic_weights,
        restricted_traits = restricted_traits,
        genetic_distances = if (is.null(genetic_distances)) rep(0.1, n_loci - 1) else genetic_distances
      )
    ),
    class = "selection_simulation"
  )
}


# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print Method for Selection Simulation Results
#'
#' @param x Object of class \code{selection_simulation}
#' @param ... Additional arguments (not used)
#'
#' @export
print.selection_simulation <- function(x, ...) {
  cat("Stochastic Selection Index Simulation Results\n")
  cat("==============================================\n\n")

  cat("Simulation Parameters:\n")
  cat("  Cycles:", x$parameters$n_cycles, "\n")
  cat("  Population size:", x$parameters$n_individuals, "\n")
  cat("  QTL per trait:", x$parameters$n_loci, "\n")
  cat("  Number of traits:", x$parameters$n_traits, "\n")
  cat("  Heritability:", paste(round(x$parameters$heritability, 2), collapse = ", "), "\n")
  cat("  Selection proportion:", x$parameters$selection_proportion, "\n")

  if (!is.null(x$parameters$restricted_traits)) {
    cat("  Restricted traits:", paste(x$parameters$restricted_traits, collapse = ", "), "\n")
  }

  cat("\nFinal Mean Genetic Values (last cycle):\n")
  final_cycle <- x$parameters$n_cycles

  cat("  LPSI: ", paste(round(x$lpsi_mean[final_cycle, ], 3), collapse = ", "), "\n")
  cat("  ESIM: ", paste(round(x$esim_mean[final_cycle, ], 3), collapse = ", "), "\n")
  cat("  RLPSI:", paste(round(x$rlpsi_mean[final_cycle, ], 3), collapse = ", "), "\n")
  cat("  RESIM:", paste(round(x$resim_mean[final_cycle, ], 3), collapse = ", "), "\n")

  invisible(x)
}


#' Summary Method for Selection Simulation Results
#'
#' @param object Object of class \code{selection_simulation}
#' @param ... Additional arguments (not used)
#'
#' @export
summary.selection_simulation <- function(object, ...) {
  cat("Stochastic Selection Index Simulation Summary\n")
  cat("==============================================\n\n")

  print(object)

  cat("\nTotal Genetic Gain (sum across all cycles):\n")
  cat("  LPSI: ", paste(round(colSums(object$lpsi_gain), 3), collapse = ", "), "\n")
  cat("  ESIM: ", paste(round(colSums(object$esim_gain), 3), collapse = ", "), "\n")
  cat("  RLPSI:", paste(round(colSums(object$rlpsi_gain), 3), collapse = ", "), "\n")
  cat("  RESIM:", paste(round(colSums(object$resim_gain), 3), collapse = ", "), "\n")

  cat("\nMean Genetic Gain per Cycle:\n")
  cat("  LPSI: ", paste(round(colMeans(object$lpsi_gain[-1, , drop = FALSE]), 3), collapse = ", "), "\n")
  cat("  ESIM: ", paste(round(colMeans(object$esim_gain[-1, , drop = FALSE]), 3), collapse = ", "), "\n")
  cat("  RLPSI:", paste(round(colMeans(object$rlpsi_gain[-1, , drop = FALSE]), 3), collapse = ", "), "\n")
  cat("  RESIM:", paste(round(colMeans(object$resim_gain[-1, , drop = FALSE]), 3), collapse = ", "), "\n")

  invisible(object)
}


#' Plot Method for Selection Simulation Results
#'
#' @param x Object of class \code{selection_simulation}
#' @param trait_index Which trait to plot (default: 1)
#' @param type Type of plot: "mean" for mean genetic value, "gain" for genetic gain per cycle
#' @param ... Additional arguments passed to plot
#'
#' @export
plot.selection_simulation <- function(x, trait_index = 1, type = c("mean", "gain"), ...) {
  type <- match.arg(type)
  n_cycles <- x$parameters$n_cycles
  cycles <- 1:n_cycles

  if (type == "mean") {
    # Plot mean genetic values
    ylab <- paste("Mean Genetic Value - Trait", trait_index)
    ylim <- range(c(
      x$lpsi_mean[, trait_index], x$esim_mean[, trait_index],
      x$rlpsi_mean[, trait_index], x$resim_mean[, trait_index]
    ))

    plot(cycles, x$lpsi_mean[, trait_index],
      type = "l", col = "blue", lwd = 2,
      xlab = "Selection Cycle", ylab = ylab, ylim = ylim,
      main = "Genetic Progress Across Cycles", ...
    )
    lines(cycles, x$esim_mean[, trait_index], col = "red", lwd = 2)
    lines(cycles, x$rlpsi_mean[, trait_index], col = "green", lwd = 2)
    lines(cycles, x$resim_mean[, trait_index], col = "orange", lwd = 2)
  } else {
    # Plot genetic gains
    ylab <- paste("Genetic Gain - Trait", trait_index)
    ylim <- range(c(
      x$lpsi_gain[, trait_index], x$esim_gain[, trait_index],
      x$rlpsi_gain[, trait_index], x$resim_gain[, trait_index]
    ))

    plot(cycles, x$lpsi_gain[, trait_index],
      type = "l", col = "blue", lwd = 2,
      xlab = "Selection Cycle", ylab = ylab, ylim = ylim,
      main = "Genetic Gain per Cycle", ...
    )
    lines(cycles, x$esim_gain[, trait_index], col = "red", lwd = 2)
    lines(cycles, x$rlpsi_gain[, trait_index], col = "green", lwd = 2)
    lines(cycles, x$resim_gain[, trait_index], col = "orange", lwd = 2)
  }

  legend("topleft", c("LPSI", "ESIM", "RLPSI", "RESIM"),
    col = c("blue", "red", "green", "orange"), lwd = 2, bty = "n"
  )

  invisible(x)
}
