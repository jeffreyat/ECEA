getECI <- function(smd1, smd2, p1=NULL, p2=NULL) {
  #' Generates equivalent change index, given two vectors of
  #' effect sizes.
  #'
  #' If you do not want to adjust the ECI by p-value, then
  #' do not pass p-values.
  #'
  #' NOTE: This function will subset the data to include only the
  #' entries that are in common to both vectors.
  #'
  #' @param smd1 Named vector of effect sizes for first treatment.
  #' @param smd2 Named vector of effect sizes for second treatment.
  #' @param p1 Named vector of unadjusted p-values for first treatment.
  #' @param p2 Named vector of unadjusted p-values for second treatment.
  #'
  #' @export

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

  if(!is.null(p1) & !is.null(p2)) {
    # Find the max p-value for each variable.
    p1 <- p1[common_vars]
    p2 <- p2[common_vars]
    comb_p = pmax(p1, p2)

    # Adjust the ECI by p-values.
    eci <- eci * (1 - comb_p)
  }

  return(eci)
}

doECEA <- function(gene_sets, eci, min_size=10, max_size=500, num_perm=1000, fdr_cutoff=.2, seed=1) {
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
  #' @export
  #'
  require(fgsea)

  set.seed(seed)
  ecea_res <- fgsea(gene_sets, eci, minSize=min_size, maxSize=max_size, nperm=num_perm)
  return(ecea_res[ecea_res$padj < fdr_cutoff,])
}

getReactome <- function(species = 'human', progress=TRUE) {
  #' Gets gene sets for use with ECEA.
  #'
  #' @param species character, currently accepts 'human' or 'mouse'
  #' @param progress Boolean value indicating if progress bar should be shown.
  #'
  #' @export
  #'
  require(reactome.db)
  require(annotate)

  db = ''

  if(species == 'human') {
    require(org.Hs.eg.db)
    db <- 'org.Hs.eg'
  } else if(species == 'mouse') {
    require(org.Mm.eg.db)
    db <- 'org.Mm.eg'
  } else {
    stop(paste0('Species ', species, ' not supported.'))
  }

  reactome_sets_full <- as.list(reactomePATHID2EXTID)

  pb <- txtProgressBar(min = 0, max = length(reactome_sets_full), style = 3)

  reactome_sets <- list()

  for(i in 1:length(reactome_sets_full)) {
    reactome_sets[[i]] <- as.vector(na.omit(getSYMBOL(reactome_sets_full[[i]], data=db)))
    if(progress) {
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)
  names(reactome_sets) <- names(reactome_sets_full)

  xx = as.list(reactomePATHID2NAME)

  names(reactome_sets) <- xx[names(reactome_sets)]

  reactome_sets <- reactome_sets[lapply(reactome_sets, length)>0]

  return(reactome_sets)
}
