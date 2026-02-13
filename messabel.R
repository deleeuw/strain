messabel<-matrix(c(-1,0,1,-1,0,1,-1,0,1,1,1,1,0,0,0,-1,-1,-1), 9, 2)
messabelDist <- dist(messabel)
messabelData <- makeMDSData(messabelDist - 2)