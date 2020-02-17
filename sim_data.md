## Reproduce simulated data

### Please read:
> This script reproduces simulated data used in *Deciphering sex-specific
genetic architectures using local Bayesian regressions*. The simulation used real
genotypes from the UK Biobank, which we cannot provide in this repository. Therefore it
**is assumed the user of this script has access to their own UKBB array genotypes**.
To work with this script, the array genotypes are assumed merged (not split across
chromosomes) and in Plink (bed/fam/bim) format.

First install and load the required R package, BGData. BGData is a suite of
packages used for memory-mapping large datasets.

```r
install.packages("BGData")
library(BGData)
```

Set the parameters used for the simulation including the number of QTL, the prop.
of variance explained by the QTL, and the proportion of QTL with GxS interacrtions.

```r
nQTL <- 150 # The total number of QTL
h2 <- 0.05 # Prop. var. explained by the QTL
p_sexdiff <- 0.4 # Proportion of QTL with sex-specific effects
nRep <- 30 # The number of MC replicates.
```

Using BGData, memory-map your UKBB genotypes. This line will depend on where
your plink file is stored on the name of the prefix.

```r
X <- BEDMatrix("path/to/plink_file_prefix")
```

Load IDs of individuals and SNPs used in the simulation. Both data is provided
here in this repository. Note that not all SNPs genome-wide were used in the simulation.
Instead, roughly 1/10 of the genome was used (60,000 SNPs; 6,000 SNPs taken from
10 different chromosomes).

```r
load("simulation_full/LABELS.RData") # Contains ids and sex of each individual
load("simulation_full/SNPS.RData") # Contains SNPs used in simulation
```

Match up the rows and columns of your memory-mapped genotype matrix with the IDs
and SNPs used in the simulation.

```r
colsX <- colnames(X) %in% SNPS
rowsX <- rownames(X) %in% LABELS$ID
IDs <- LABELS$ID[match(rownames(X)[rowsX], LABELS$ID)]
stopifnot(all(rownames(X)[rowsX] == IDs))
n <- length(rownames(X)[rowsX])
nSNPs <- length(colnames(X)[colsX])
```

Prepare the outputs. Outputs are as follows:
- `Y`: a matrix of simulated phenotypes for individuals (in rows) for each replicate
of the simulation (in columns).
- `E`: a matrix of residuals, matching the order of individuals and replicates as in `Y`.
- `SIGNAL`: a matrix of simulated additive genetic values, matching the order of individuals
and replicates as in `Y`.
- `B_MEN`: a matrix of male QTL effects with QTL in rows and replicates in columns.
- `B_WOMEN`: a matrix of female QTL effects with QTL in rows and replicates in columns.
- `QTL`: a matrix of QTL IDs (in rows) for each replicate (in columns).


```r
Y <- matrix(nrow = n, ncol = nRep)
rownames(Y) <- rownames(X)[rowsX]
E <- Y
SIGNAL <- Y
B_MEN <- matrix(nrow = nQTL, ncol = nRep)
B_WOMEN <- matrix(nrow = nQTL, ncol = nRep)
QTL <- matrix(nrow = nQTL, ncol = nRep)
```

For each simulation replicate, follow this basic strategy:
1. Sample QTL
2. Sample which QTL will have identical effects between men and women. Sample these QTL
effects from a gamma distribution.
3. Sample which QTL will have different nonzero effects between men and women.
Sample these effects from a gamma distribution. Which sex has a larger effect is
also randomized.
4. The remaining QTL have a non-null effect in one sex and a null effect in the other.
Sample these effects from a gamma distribution. Which sex has a non-null effect
is also randomized.
5. Calculate additive genetic values, and obtain additive genetic variance for each sex.
6. Set the residual variance for each sex, according to the additive genetic variance
for each sex and the heritability, which is assumed the same for both sexes.
7. Sample residuals for males and females, given male and female residual variances.
8. Simulate the phenotype, given additive genetic values and residuals for each individual.

```r
set.seed(195021)
for (i in 1:nRep) {
  # Sample QTN
  QTL_full <- sample(1:nSNPs, size = nQTL, replace = FALSE)

  # Select 60% of QTL to have common effects between men and women
  QTL_set1 <- sample(QTL_full, size = nQTL * (1 - p_sexdiff), replace = FALSE)
  b_set1 <- rgamma(n = length(QTL_set1),
                   shape = 10,
                   scale = 1) * ifelse(runif(n = length(QTL_set1)) > .5, 1, -1)

  # Select 20% QTL to have different nonzero effects between men and women
  QTL_set2 <- sample(QTL_full[!QTL_full %in% QTL_set1],
                     size = nQTL * (p_sexdiff / 2),
                     replace = FALSE)
  sex_idx <- runif(n = length(QTL_set2)) > 0.5
  b_set2_men <- rgamma(n = length(QTL_set2),
                       shape = ifelse(sex_idx, 5, 20),
                       scale = 1)
  b_set2_women <- rgamma(n = length(QTL_set2),
                         shape = ifelse(1 - sex_idx, 5, 20),
                         scale = 1)

  # Select 20% QTL to have nonzero effects in one sex and zero effects in the other
  QTL_set3 <- sample(QTL_full[!QTL_full %in% c(QTL_set1, QTL_set2)],
                     size = nQTL * (p_sexdiff / 2),
                     replace = FALSE)
  sex_idx <- runif(n = length(QTL_set3)) > 0.5
  b_set3_men <- rgamma(n = length(QTL_set3),
                       shape = 10,
                       scale = 1) * ifelse(sex_idx, 1, 0)
  b_set3_women <- rgamma(n = length(QTL_set3),
                         shape = 10,
                         scale = 1) * ifelse(1 - sex_idx, 1, 0)

  # Extract QTL, center and scale QTL matrix
  QTL_idx <- colsX[c(QTL_set1, QTL_set2, QTL_set3)]
  Z <- X[rowsX, QTL_idx]
  Z <- scale(Z, center = TRUE, scale = TRUE)
  Z[is.na(Z)] <- 0
  Z <- scale(Z, center = TRUE, scale = TRUE)

  # Calculate polygenic scores for men and women
  b_men <- c(b_set1, b_set2_men, b_set3_men)
  b_women <- c(b_set1, b_set2_women, b_set3_women)
  men <- LABELS$ID[LABELS$sex == 1]
  women  <- LABELS$ID[LABELS$sex == 0]
  signal_men <- Z[rownames(Z) %in% men, ] %*% b_men
  signal_women <- Z[rownames(Z) %in% women, ] %*% b_women
  varU_men <- var(signal_men)
  varU_women <- var(signal_women)

  # Scale error variance based on calculated additive variance and heritability
  # and simulate trait.
  varE_men <- varU_men / h2 - varU_men
  varE_women <- varU_women / h2 - varU_women
  e_men <- rnorm(length(signal_men), mean = 0, sd = sqrt(varE_men))
  y_men <- signal_men + e_men
  e_women <- rnorm(length(signal_women), mean = 0, sd = sqrt(varE_women))
  y_women <- signal_women + e_women
  y <- rbind(y_men, y_women)
  error <- c(e_men, e_women)
  signal <- rbind(signal_men, signal_women)

  # Comibine results for men and women, in the order that they are in Y, E, and SIGNAL
  idx <- match(rownames(Y), rownames(y))
  y <- y[idx, ]
  error <- error[idx]
  signal <- signal[idx, ]
  QTL[, i] <- colnames(Z)
  Y[, i] <- y
  E[, i] <- error
  SIGNAL[, i] <- signal
  B_MEN[, i] <- b_men
  B_WOMEN[, i] <- b_women
}
```

Now the outputs are available to save.
