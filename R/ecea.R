getECI <- function(smd1, smd2, p1, p2) {
  #' Generates equivalent change index, given two vectors of
  #' effect sizes.
  #'
  #' NOTE: This function will subset the data to include only the
  #' entries that are in common to both vectors.
  #'
  #' @param smd1 Named vector of effect sizes for first treatment.
  #' @param smd2 Named vector of effect sizes for second treatment.
  #' @param p1 Named vector of unadjusted p-values for first treatment.
  #' @param p2 Named vector of unadjusted p-values for second treatment.

  # Find the common entries in both vectors.
  common_vars <- intersect(names(smd1), names(smd2))

  # If the number of common entries are less than the max in either
  # vector, generate warning.
  if(length(common_vars) < max(length(smd1), length(smd2))) {
    message('Number of common variables less than max length.',
            appendLF = T)
  }

  if(length(common_vars) == 0) {
    warning('No variables in common between vectors.')
    return(list(eci=NA, base= NA))
  }

  # Subset and reorder to the common variables.
  smd1 <- smd1[common_vars]
  smd2 <- smd2[common_vars]

  # Calculate the ECI
  eci <- list()
  base <- list()
  for(i in 1:length(common_vars)) {
    eci[[i]] <- sign(smd1[i] * -smd2[i]) * (max(abs(smd1[i]), abs(smd2[i])) -
                                              (abs(smd1[i]) + abs(smd2[i])))
    eci[[i]] <- eci[[i]] / max(abs(smd1[i]), abs(smd2[i]))
  }

  eci <- unlist(eci)
  eci <- setNames(eci, common_vars)

  # Find the max p-value for each variable.
  p1 <- p1[common_vars]
  p2 <- p2[common_vars]
  comb_p = pmax(p1, p2)

  # Adjust the ECI by p-values.
  eci <- eci * (1 - comb_p)

  return(eci)
}

ECEA <- function(gene_sets, eci, min_size=10, max_size=500, num_perm=1000, fdr_cutoff=.2, seed=1) {
  #' Performs Equivalent Change Enrichment Analysis
  #'
  #' @param gene_sets list of gene sets, as in fgsea
  #' @param eci Named vector of equivalent change indices from getECI()
  #' @param min_size Minimum size of pathways
  #' @param max_size Maximum size of pathways
  #' @param num_perm Number of permutations for permutation test
  #' @param fdr_cutoff FDR threshold for returned pathways
  #' @param seed Random seed for permutation test
  #'
  require(fgsea)

  set.seed(seed)
  ecea_res <- fgsea(gene_sets, eci, minSize=min_size, maxSize=max_size, nperm=num_perm)
  return(ecea_res[ecea_res$padj < fdr_cutoff,])
}
