doubleCenterSS <- function(x) {
  r <- apply(x, 1, mean)
  m <- mean(x)
  return(x - outer(r, r, "+") + m)
}

x <- matrix(rnorm(10), 5, 2)
d <- as.matrix(dist(x))
z <- matrix(rnorm(10), 5, 2)
tilde <- as.matrix(dist(z))

sigma <- function(delta) {
  h <- doubleCenterSS(delta - d)
  return(sum(h^2) / 4)
}

eta <- function(delta, tilde) {
  h <- doubleCenterSS(tilde - d)
  e1 <- sigma(tilde) + .5 * sum(h * (delta - tilde)) + .25 * sum((delta - tilde)^2)
  e2 <- sum((delta - hat(tilde))^2) / 4
  return(c(e1, e2))
}

hat <- function(tilde) {
  h <- doubleCenterSS(tilde - d)
  return(tilde - h)
}

delta <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
