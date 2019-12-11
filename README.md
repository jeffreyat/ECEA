---
title: "README"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Equivalent Change Enrichment Analysis is a new approach to *in silico* functional genomic analysis that can help you to determine if biological pathways are changed in similar (equivalent) or opposing (inverse) ways across two genomic experiments. For example, it can help you to determine if there are similar functional impacts from two drugs, or potentially if a new drug reverses some of the changes associated with a disease. Although it is possible to do similar analyses using *ad hoc* approaches, ECEA uses a statistically robust methodology, based on GSEA but using a new local statistic that indicates equivalent or inverse change among genes in a pathway.

To install the package in R, you can simply do the following, but not that building the vignette will take some time and is optional:

```{r}
library(devtools)

install_github('jeffreyat/ECEA', build_vignettes=TRUE)
```

As noted, a vignette is available to help you get started. Once the package is installed, assuming you built the vignette, you can do the following:

```{r}
vignette('ECEA')
```
