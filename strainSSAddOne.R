library(RSpectra)
library(polynom)
dyn.load("jeffrey.so")

strainSSAddOne <- function(theData,
                           ndim = 2,
                           addc = 0,
                           pos = TRUE,
                           itmax = 100,
                           eps = 1e-6,
                           verbose = TRUE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delt <- theData$delta
  bond <- -min(delt)
  dnul <- matrix(0, nobj, nobj) # order n matrix of dissimilarities
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dnul[i, j] <- dnul[j, i] <- delt[k]
  }
  dmat <- dnul + addc * (1 - diag(nobj))
  bmat <- -doubleCenterSS(dmat^2) / 2
  ev <- eigs_sym(bmat, ndim, which = "LA")
  fv <- pmax(0, ev$values)
  evev <- diag(sqrt(fv))
  xold <- ev$vectors %*% evev
  dold <- as.matrix(dist(xold))
  sold <- sum(bmat^2) - sum(fv^2)
  itel <- 1
  repeat {
    h <- doubleCenterSS(dnul^2 - dold^2)
    g <- doubleCenterSS(dnul)
    a0 <- sum(h^2)
    a1 <- 4 * sum(h * g)
    a2 <- 4 * sum(g^2) - 2 * sum(diag(h))
    a3 <- -4 * sum(diag(g))
    a4 <- (nobj - 1)
    aa <- c(a0, a1, a2, a3, a4)
    u <- .C("jeffrey",
            a = as.double(aa),
            mw = as.double(0),
            mv = as.double(0))
    if (pos) {
      addc <- max(bond, Re(u$mw))
    } else {
      addc <- Re(u$mw)
    }
    dmat <- dnul + addc * (1 - diag(nobj))
    smid <- u$mv / 4.0
    bmat <- -doubleCenterSS(dmat^2) / 2
    ev <- eigs_sym(bmat, ndim, which = "LA")
    fv <- pmax(0, ev$values)
    evev <- diag(sqrt(fv))
    xnew <- ev$vectors %*% evev
    dnew <- as.matrix(dist(xnew))
    h <- doubleCenterSS(dmat^2 - dnew^2)
    snew <- sum(bmat ^ 2) - sum(fv ^ 2)
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
      " addc ",
      formatC(
        addc,
        format = "f",
        digits = 10,
        width = 15
      ),
      "\n"
    )
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
    xold <- xnew
    dold <- dnew
  }
  return(list(
    conf = xnew,
    addc = addc,
    dmat = dmat
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
