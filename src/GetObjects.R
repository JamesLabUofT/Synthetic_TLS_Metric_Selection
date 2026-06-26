objectProperties <- read.csv("data/pointCloudProperties.csv")

recentreScan <- function(las){
  meanX <- mean(las$X)
  meanY <- mean(las$Y)
  las$X <- las$X - meanX
  las$Y <- las$Y - meanY
  return(las)
}


lasList <- list()
for (i in 1:nrow(objectProperties)){
  objectLas <- recentreScan(lidR::readLAS(objectProperties[i,]$FilePath))
  lasList <- append(lasList, objectLas)
}

objectProperties$Las <- lasList

Trees <- subset(objectProperties, Type == "Tree")
Trees$bin <- cut(Trees$Height, 0:10*1.8+3.5, labels = 1:10, right = F)

Understory <- subset(objectProperties, Type == "Understory")

Snag <- subset(objectProperties, Type == "Snag")

