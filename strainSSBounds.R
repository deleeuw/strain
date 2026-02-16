
library(RSpectra)

strainSSBounds <- function(theData,
                           ndim = 2,
                           ties = 1,
                           itmax = 1000,
                           eps = 1e-6,
                           verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  blks <- theData$blocks
  wght <- theData$weights
  iind <- theData$iind
  jind <- theData$jind
  lower <- theData$lower
  upper <- theData$upper
  dhat <- matrix(0, nobj, nobj)
  lower <
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dhat[i, j] <- dhat[j, i] <- theData$delta[k]
  }
  bmat <- -doubleCenterSS(dhat^2) / 2
  ev <- eigs_sym(bmat, ndim, which = "LA")
  fv <- pmax(0, ev$values)
  evev <- diag(sqrt(fv))
  xold <- ev$vectors %*% evev
  dold <- as.matrix(dist(xold))
  hold <- doubleCenterSS(dhat - dold^2)
  sold <- sum(hold^2) / 4
  daux <- rep(0, ndat)
  itel <- 1
  repeat {
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      daux[k] <- dhat[i, j] - hold[i, j]
      if (daux[k] < lower[k]) {
        daux[k] <- lower[k]
      }
      if (daux[k] > upper[k]) {
        daux[k] <- upper[k]
      }
      dhat[i, j] <- dhat[j, i] <- daux[k]
    }
    hmid <- doubleCenterSS(dhat - dold^2)
    smid <- sum(hmid^2) / 4
    bmat <- -doubleCenterSS(dhat) / 2
    ev <- eigs_sym(bmat, ndim, which = "LA")
    fv <- pmax(0, ev$values)
    evev <- diag(sqrt(fv))
    xnew <- ev$vectors %*% evev
    dnew <- as.matrix(dist(xnew))
    hnew <- doubleCenterSS(dhat - dnew^2)
    snew <- sum(hnew^2) / 4
    if (verbose) {
      cat(
        " itel ",
        formatC(itel, format = "d", width = 6),
        " sold ",
        formatC(
          sold,
          format = "f",
          digits = 6,
          width = 10
        ),
        " smid ",
        formatC(
          smid,
          format = "f",
          digits = 6,
          width = 10
        ),
        " snew ",
        formatC(
          snew,
          format = "f",
          digits = 6,
          width = 10
        ),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
    dold <- dnew
    hold <- hnew
  }
  ddis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    ddis[k] <- dnew[i, j]
  }
  return(list(
    conf = xnew,
    delt = theData$delta,
    daux = daux,
    ddis = ddis,
    itel = itel,
    loss = snew
  ))
}

doubleCenterSS <- function(x) {
  r <- apply(x, 1, mean)
  m <- mean(x)
  return(x - outer(r, r, "+") + m)
}

ei <- function(i, n) {
  ifelse(i == 1:n, 1, 0)
}
