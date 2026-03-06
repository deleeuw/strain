


library(RSpectra)
library(nnls)

strainSSMissing <- function(theData, # MDS data structure with non-missing elements
                            ndim = 2,
                            xinit = NULL,
                            itmax = 1000,
                            eps = 1e-10,
                            verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delta <- theData$delta
  if (is.null(xinit)) {
    # make xold
  } else {
    xold <- xinit
  }
  dmat <- dnul <- matrix(0, nobj, nobj) 
  indi <- 1 - diag(nobj) 
  mdel <- mean(theData$delta)
  mdel <- 0
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- delta[k]
    dnul[i, j] <- dnul[j, i] <- delta[k]
    indi[i, j] <- indi[j, i] <- 0
  }
  nmis <- sum(indi) / 2
  if (nmis > 0) {
    kind <- matrix(0, nmis, 2)
    k <- 1
    for (i in 2:nobj) {
      for (j in 1:(i - 1)) {
        if (indi[i, j] == 1) {
          kind[k, 1] <- i
          kind[k, 2] <- j
          k <- k + 1
        }
      }
    }
    amat <- matrix(0, nobj^2, nmis)
    for (k in 1:nmis) {
      ii <- ei(kind[k, 1], nobj)
      jj <- ei(kind[k, 2], nobj)
      eij <- outer(ei(kind[k, 1], nobj), ei(kind[k, 2], nobj))
      amat[, k] <- as.vector(doubleCenterSS(eij + t(eij)))
    }
  }
  itel <- 1
  sold <- Inf
  repeat {
    bmat <- -doubleCenterSS(dmat^2) / 2
    smat <- sum(bmat^2)
    ev <- eigs_sym(bmat, ndim, which = "LA")
    fv <- pmax(0, ev$values)
    smid <- smat - sum(fv^2)
    evev <- diag(sqrt(fv))
    x <- ev$vectors %*% evev
    d <- as.matrix(dist(x))
    h <- as.vector(doubleCenterSS(d^2 - dnul^2))
    if (nmis == 0) {
      snew <- sum(h^2) / 4
      return(list(
        conf = x,
        dmat = dmat,
        itel = 0,
        loss = snew
      ))
    }
    tn <- nnls(amat, h)
    snew <- tn$deviance / 4
    if (verbose) {
      cat(
        " itel ",
        formatC(itel, format = "d", width = 6),
        " sold ",
        formatC(
          sold,
          format = "f",
          digits = 10,
          width = 15
        ),
        " smid ",
        formatC(
          smid,
          format = "f",
          digits = 10,
          width = 15
        ),
        " snew ",
        formatC(
          snew,
          format = "f",
          digits = 10,
          width = 15
        ),
        "\n"
      )
    }
    th <- tn$x
    for (k in 1:nmis) {
      i <- kind[k, 1]
      j <- kind[k, 2]
      dmat[i, j] <- dmat[j, i] <- sqrt(th[k])
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
  }
  return(list(
    conf = x,
    dmat = dmat,
    itel = itel,
    loss = snew
  ))
}

strainFillMissing <- function(theData) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delta <- theData$delta
  dmat <- dnul <- matrix(0, nobj, nobj) 
  indi <- 1 - diag(nobj) 
  mdel <- mean(theData$delta)
  mdel <- 0
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- delta[k]
    dnul[i, j] <- dnul[j, i] <- delta[k]
    indi[i, j] <- indi[j, i] <- 0
  }
  nmis <- sum(indi) / 2
  if (nmis > 0) {
    kind <- matrix(0, nmis, 2)
    k <- 1
    for (i in 2:nobj) {
      for (j in 1:(i - 1)) {
        if (indi[i, j] == 1) {
          kind[k, 1] <- i
          kind[k, 2] <- j
          k <- k + 1
        }
      }
    }
    amat <- matrix(0, nobj^2, nmis)
    for (k in 1:nmis) {
      ii <- ei(kind[k, 1], nobj)
      jj <- ei(kind[k, 2], nobj)
      eij <- outer(ei(kind[k, 1], nobj), ei(kind[k, 2], nobj))
      amat[, k] <- as.vector(doubleCenterSS(eij + t(eij)))
    }
  }
  return(dnul = dnul, dmat = dmat, amat = amat, nmis = nmis)
}

doubleCenterSS <- function(x) {
  r <- apply(x, 1, mean)
  m <- mean(x)
  return(x - outer(r, r, "+") + m)
}

ei <- function(i, n) {
  ifelse(i == 1:n, 1, 0)
}