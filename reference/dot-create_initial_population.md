# Create Initial Population for Simulation

Creates a simulated diploid population with quantitative trait loci
(QTL) and marker information. Each individual has two haplotypes with
additive genetic effects.

## Usage

``` r
.create_initial_population(
  n_individuals,
  n_loci,
  n_traits,
  qtl_effects = NULL,
  genetic_distances = NULL
)
```

## Arguments

- n_individuals:

  Number of individuals in the population

- n_loci:

  Number of QTL per trait

- n_traits:

  Number of quantitative traits

- qtl_effects:

  Optional matrix of QTL effects (n_loci x n_traits). If NULL, effects
  are sampled from standard normal distribution.

- genetic_distances:

  Optional vector of genetic distances between adjacent loci (in
  Morgans). If NULL, assumes 0.1 Morgan spacing (roughly 10 cM).

## Value

List with components:

- `haplotype1` - First haplotype matrix (n_individuals x n_loci)

- `haplotype2` - Second haplotype matrix (n_individuals x n_loci)

- `qtl_effects` - QTL effect matrix (n_loci x n_traits)

- `genetic_distances` - Distances between loci (in Morgans)

- `recombination_fractions` - Recombination fractions between adjacent
  loci
