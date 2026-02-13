procrustus <- function(x,
                       itmax = 100,
                       eps = 1e-6,
                       verbose = TRUE) {
  z <- x[[1]]
  m <- length(x)
  nr <- dim(z)[1]
  nc <- dim(z)[2]
  itel <- 1
  sold <- 0.0
  for (k in 1:m) {
    sold <- sold + sum((z - x[[k]])^2)
  }
  repeat {
    cm <- matrix(0, nr, nc)
    smid <- 0.0
    for (k in 1:m) {
      s <- svd(crossprod(x[[k]], z))
      r <- tcrossprod(s$u, s$v)
      x[[k]] <- x[[k]] %*% r
      cm <- cm + x[[k]]
      smid <- smid + sum((z - x[[k]])^2)
    }
    z <- cm / m
    snew <- 0.0
    for (k in 1:m) {
      snew <- snew + sum((z - x[[k]])^2)
    }
    if (verbose) {
      print(c(sold, smid, snew))
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
  }
  return(list(
    x = x,
    itel = itel,
    z = z,
    loss = snew
  ))
}
