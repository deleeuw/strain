data(morse, package = "smacof")
morseMatrix <- as.matrix(morse)
morseData <- makeMDSData(morse)
morseLabels <- as.character(attr(morse, "Labels"))

misData <- function(data, howmany) {
  if (howmany == 0) {
    return(data)
  }
  ndat <- length(data[[1]])
  tmis <- sample(1:ndat, howmany)
  return(list(iind = data[[1]][-tmis], 
              jind = data[[2]][-tmis], 
              delta = data[[3]][-tmis], 
              blocks = data[[4]][-tmis], 
              weights = data[[5]][-tmis], 
              nobj = data[[6]], 
              ndat = data[[7]] - howmany))
}
