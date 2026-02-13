tryMe <- function(i, j, k, l, n) {
  e <- rep(1, n)
  ei <- ifelse(i == 1:n, 1, 0)
  ej <- ifelse(j == 1:n, 1, 0)
  ek <- ifelse(k == 1:n, 1, 0)
  el <- ifelse(l == 1:n, 1, 0)
  dik <- dki <- ifelse(i == k, 1, 0)
  dil <- dli <- ifelse(i == l, 1, 0)
  djk <- dkj <- ifelse(j == k, 1, 0)
  djl <- dlj <- ifelse(j == l, 1, 0)
  jj <- diag(n) - 1/n
  xij <- outer(ei, ej) + outer(ej, ei)
  xkl <- outer(ek, el) + outer(el, ek)
  mij <- outer(ei + ej, e) + outer(e, ei + ej)
  mkl <- outer(ek + el, e) + outer(e, ek + el)
  nn <- matrix(1, n, n)
  alpa <- (dik * djl) + (dil * djk)
  beta <- dik + djl + djk + dil
  print(c(alpa, beta))
  print(cxx <- c(sum(diag(xij %*% xkl)), 2 * alpa))
  print(cmm <- c(sum(diag(mij %*% mkl)), 8 + 2 * n * beta))
  print(cxm <- c(sum(diag(xij %*% mkl)), 2 * beta))
  print(cxn <- c(sum(diag(xij %*% nn)), 2))
  print(cmn <- c(sum(diag(mij %*% nn)), 4 * n))
  print(cnn <- c(sum(diag(nn %*% nn)), n ^ 2))
  hij <- jj %*% (outer(ei, ej) + outer(ej, ei)) %*% jj
  hkl <- jj %*% (outer(ek, el) + outer(el, ek)) %*% jj
  hhh <- 2 * alpa - (4 * beta) / n + (4 + 2 * n * beta) / n^2 
  return(c(sum(diag(hij %*% hkl)), hhh))
}