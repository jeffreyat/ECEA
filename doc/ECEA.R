## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(ECEA)
library(GEOquery)
library(mgu74av2.db)
library(limma)
library(fgsea)
library(knitr)

## -----------------------------------------------------------------------------
mouse_probes_to_symbols <- function(dat) {
  x <- mgu74av2SYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  
  genes <- unlist(xx[rownames(dat)])
  
  common <- intersect(names(genes), rownames(dat))
  
  genes <- genes[common]
  dat <- dat[common,]
  
  message("Aggregating...", appendLF = FALSE)
  dat <- aggregate(dat, by=list(genes), FUN=median)
  rownames(dat) = dat[,1]
  dat <- dat[,-1]
  message("done!", appendLF=TRUE)
  
  return(dat)
}

## ---- warning = FALSE, message = FALSE----------------------------------------
gse <- getGEO('GSE35378', GSEMatrix = T)
dat <- exprs(gse[[1]])
dat <- log2(dat)
dat <- mouse_probes_to_symbols(dat)

## ---- warning = FALSE, message = FALSE----------------------------------------
des <- model.matrix(~c(1,1,1,0,0,0))

mod1 <- lmFit(dat[,1:6], des)
ebmod1 <- eBayes(mod1)
mod1tt <- topTable(ebmod1, coef=2, number=100000)

mod2 <- lmFit(dat[,7:12], des)
ebmod2 <- eBayes(mod2)
mod2tt <- topTable(ebmod2, coef=2, number=100000)

## ---- warning = FALSE, message = FALSE----------------------------------------
mod1_fc <- mod1tt$logFC
names(mod1_fc) <- rownames(mod1tt)

mod2_fc <- mod2tt$logFC
names(mod2_fc) <- rownames(mod2tt)

mod1_p <- mod1tt$P.Value
names(mod1_p) <- rownames(mod1tt)

mod2_p <- mod2tt$P.Value
names(mod2_p) <- rownames(mod2tt)

## ---- warning = FALSE, message = FALSE----------------------------------------
eci = getECI(mod1_fc, mod2_fc, mod1_p, mod2_p)

## ---- warning = FALSE, message = FALSE----------------------------------------
reactome_sets = getReactome(species = 'mouse', progress=FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
ecea_res = doECEA(reactome_sets, eci)

## ---- warning = FALSE, message = FALSE----------------------------------------
kable(ecea_res[ecea_res$padj < .25,1:7])

