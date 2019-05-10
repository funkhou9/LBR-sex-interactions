# Function to process posterior samples obtained by a LBR GxS model.

# Requires:
# 1. Posterior samples from three parameters: a main genetic effect, female-specific
# interaction effect, and male-specific interaction effect.
# 2. getWindows() function to define LD-based windows.
# 3. readBinMat from the BGLR package.

# Arguments:
# X--genotype matrix that was used to fit LBR model
# main--path to main genetic effect posterior samples
# male_int--path to male interaction effect posterior samples
# female_int--path to female interaction effect posterior samples
# var_threshold--a proportion p, such that when using window variances to declare
# a GxS interaction, the maginitude of the difference must be p * mean of all male and female
# samples.
# rSq--allele dosage R2 threshold used to define LD.
# maxGaps--maximum number of consecutive SNPs that do not meet rSq threshold.

# Output:
# A tibble with columns: SNP,
process_samples <- function(X,
                            main,
                            male_int,
                            female_int,
                            var_threshold = 0.1,
                            rSq = 0.1,
                            maxGaps = 2) {
  BMain <- readBinMat(main)
  BM <- readBinMat(male_int)
  BF <- readBinMat(female_int)

  # Estimating and inferring sex-specific SNP effects and differences
  BMale <- BMain+BM
  bMale <- colMeans(BMale)
  pNeg <- colMeans(BMale < 0)
  pPos <- colMeans(BMale > 0)
  pMale <- apply(X = cbind(pNeg, pPos),
                 MARGIN = 1,
                 FUN = max)

  BFemale <- BMain + BF
  bFemale <- colMeans(BFemale)
  pNeg <- colMeans(BFemale < 0)
  pPos <- colMeans(BFemale > 0)
  pFemale <- apply(X = cbind(pNeg, pPos),
                   MARGIN = 1,
                   FUN = max)

  BDif <- BFemale - BMale
  pNeg <- colMeans(BDif < 0)
  pPos <- colMeans(BDif > 0)
  pSexDif <- apply(X = cbind(pNeg, pPos),
                   MARGIN = 1,
                   FUN = max)

  # Estimating and inferring sex-specific window variances and differences
  windows <- getWindows(X,
                        center_snp = 1:ncol(X),
                        rSq = rSq,
                        maxGaps = maxGaps)
  bMale_window <- c()
  bFemale_window <- c()
  pMale_window <- c()
  pFemale_window <- c()
  pSexDif_window <- c()

  for (i in 1:length(windows)) {
    cols <- windows[[i]]
    XtX <- crossprod(X[, cols], X[, cols])
    varA_men <- c()
    varA_women <- c()
    for (j in 1:nrow(BMain)) {
      b_men <- BMale[j, cols]
      b_women <- BFemale[j, cols]
      varA_men <- c(varA_men, (t(b_men) %*% XtX %*% b_men) / (nrow(X) - 1))
      varA_women <- c(varA_women, (t(b_women) %*% XtX %*% b_women) / (nrow(X) - 1))
    }
    pM <- mean(varA_men > 0)
    pF <- mean(varA_women > 0)
    threshold <- mean(c(varA_men, varA_women)) * 0.1
    pdiff <- max(mean((varA_men - varA_women) > threshold),
                 mean((varA_women - varA_men) > threshold))
    bMale_window <- c(bMale_window, mean(varA_men))
    bFemale_window <- c(bFemale_window, mean(varA_women))
    pMale_window <- c(pMale_window, pM)
    pFemale_window <- c(pFemale_window, pF)
    pSexDif_window <- c(pSexDif_window, pdiff)
  }

  if (is.null(colnames(X))) {
    snp_names <- 1:ncol(X)
  } else {
    snp_names <- colnames(X)
  }

  results <- tibble("SNP" = snp_names,
                    "b_male" = bMale,
                    "b_female" = bFemale,
                    "PP_male" = pMale,
                    "PP_female" = pFemale,
                    "PP_diff" = pSexDif,
                    "b_male_window" = bMale_window,
                    "b_female_window" = bFemale_window,
                    "PP_male_window" = pMale_window,
                    "PP_female_window" = pFemale_window,
                    "PP_diff_window" = pSexDif_window)
  return(results)
}
