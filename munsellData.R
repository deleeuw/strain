source("smacofDataUtilities.R")

munsell <- c(-2.37, -.12, -.62, .23, 1.56, 1.09, 2.02, 2.23, 
             -1.01, -1.93, -0.90, .80, -.47, 1.05, 0.78, 0.70, -1.32, -.67, 1.07, .70, 2.62, -.78, 1.25, -1.75, 0.28, -.72, -1.02, -1.23,
             -1.65, 0.49, 0.57, -0.67, 1.88, -1.18, -1.30, 0.42)
munsellMatrix <- matrix(0, 9, 9)
munsellMatrix[outer(1:9, 1:9, ">")] <- munsell
munsellMatrix <- munsellMatrix + t(munsellMatrix)
munsellDist <- as.dist(munsellMatrix)
munsellData <- makeMDSData(munsellDist)