
DepLink_MultiGenes <- function(geneList, direction = c(1,0), fileDropList, drugs, permDistr) {
  
  # check if gene is included in Dropbox and download perturbations signatures of individual genes
  
  pertSigs <- as.data.frame(matrix(data = NA, nrow = length(geneList), ncol = dim(drugs)[2]))
  pertSigs_p <- as.data.frame(matrix(data = NA, nrow = length(geneList), ncol = dim(drugs)[2]))
  for (i in 1:length(geneList)) {
      ix <- which(fileDropList[,1]==geneList[i])
      tmp <- read.table(url(fileDropList[ix,2]), na.strings = "NaN")
      pertSigs[i,] <- tmp[1,]
      pertSigs_p[i,] <- tmp[2,]
      row.names(pertSigs)[i] <- geneList[i]
      row.names(pertSigs_p)[i] <- paste0(geneList[i]," pVal")
  }
  pertSigs <- pertSigs + 10^-5*matrix(rnorm(dim(pertSigs)[1]*dim(pertSigs)[2]), nrow = dim(pertSigs)[1])
  meanPerDrug <- as.data.frame(t(colMeans(pertSigs)))
  row.names(meanPerDrug) <- "Mean Perturbation Score"
  
  # comparison against zero
  if (direction==1) { # positive regulation of gene expression
    pVal_t <- apply(X = pertSigs, MARGIN = 2, function(x) {t.test(x, alternative = "greater")$p.value})
  } else { # negative regulation
    pVal_t <- apply(X = pertSigs, MARGIN = 2, function(x) {t.test(x, alternative = "less")$p.value})
  }
  pVal_t <- t(as.data.frame(pVal_t))
  row.names(pVal_t) <- "T-test pVal"
  
  # comparison against random gene lists
  pVal_perm <- as.data.frame(matrix(data = "NA", nrow = 1, ncol = dim(pertSigs)[2]))
  row.names(pVal_perm) <- "Permutation pVal"
  if (direction==1) { # positive regulation of gene expression
    for (i in 1:dim(pertSigs)[2]) {
      pVal_perm[i] <- pnorm(q = as.numeric(meanPerDrug[i]), mean = permDistr[1, i], sd = permDistr[2, i], lower.tail = FALSE)}
  } else { # negative regulation
    for (i in 1:dim(pertSigs)[2]) {
      pVal_perm[i] <- pnorm(q = as.numeric(meanPerDrug[i]), mean = permDistr[1, i], sd = permDistr[2, i], lower.tail = TRUE)}
  }
  
  return(list("pertSigs" = pertSigs, "pertSigs_p" = pertSigs_p, "meanPerDrug" = meanPerDrug, "permDistr" = permDistr, 
              "drugs" = drugs, "pVal_t" = pVal_t, "pVal_perm" = pVal_perm))
}

vline <- function(x = 0,  color = "black") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color)
  )
}

annotation <- function(x= 0,y,color = "black"){
  list(
  x = x,
  y = y,
  text = "Observed in query genes",
  xref = "x",
  yref = "y",
  showarrow = TRUE,
  arrowhead = 7,
  ax = -20,
  ay = -40,
  size = 10)
}

