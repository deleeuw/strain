doubleCenterSS <- function(x) {
  r <- apply(x, 1, mean)
  m <- mean(x)
  return(x - outer(r, r, "+") + m)

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
  return((sigma(tilde) + 2 * sum(h * (delta - tilde)) + sum((delta - tilde)^2)) / 4)
}

hat <- function(tilde) {
  h <- doubleCenterSS(tilde - d)
  return(tilde + h)
}

 
print(sigma(tilde))

delta <- as.matrix(dist(matrix(rnorm(10), 5, 2)))
