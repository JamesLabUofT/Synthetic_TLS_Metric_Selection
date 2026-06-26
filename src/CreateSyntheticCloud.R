source("src/GetObjects.R")

recentreScan <- function(las){
  meanX <- mean(las$X)
  meanY <- mean(las$Y)
  las$X <- las$X - meanX
  las$Y <- las$Y - meanY
  return(las)
}

anglePointCloud <- function(las, angle){
  r <- sqrt(las$X**2 + las$Y**2)
  theta <- atan2(las$Y, las$X)
  
  theta <- theta + angle
  
  las$X <- r * cos(theta)
  las$Y <- r * sin(theta)
  
  return(las)
}

positionPointCloud <- function(las, X, Y){
  las$X <- las$X + X
  las$Y <- las$Y + Y
  return(las)
}

createCloudFromScenario <- function(scenario){
  scenarioScan <- NULL
  #change the following to lapply or mlapply
  for (tree in (1:scenario$numTrees[[1]])){
    treeHeight <- scenario$TreeHeights[[1]][tree]
    
    deadTree <- scenario$Aliveness[[1]][tree] == 0
    
    treeLas <- Trees[Trees$Height == treeHeight,]$Las[[1]]
    
    if (deadTree){
      treeLas <- (Snag %>% filter(Height <= treeHeight) %>% arrange(desc(Height)) %>% slice_head())$Las[[1]]
    }
    
    treeAngle <- scenario$TreeAngles[[1]][tree]
    treeX <- scenario$TreePositions[[1]][tree,][1]
    treeY <- scenario$TreePositions[[1]][tree,][2]
    
    treeLas <- treeLas |> anglePointCloud(treeAngle) |> 
      positionPointCloud(treeX, treeY)
    
    if (is.null(scenarioScan)){
      scenarioScan <- treeLas
    } else {
      scenarioScan <- rbind(scenarioScan, treeLas)
    }
  }
  if (scenario$amountSubcanopy == 0) return(scenarioScan)
  for (understory in (1:scenario$amountSubcanopy)) {
    understoryLas = Understory[scenario$UnderstoryID[[1]][understory],]$Las[[1]]
    
    understoryAngle <- scenario$UnderstoryAngle[[1]][understory]
    understoryX <- scenario$UnderstoryPositions[[1]][understory,][1]
    understoryY <- scenario$UnderstoryPositions[[1]][understory,][2]
    
    understoryLas <- understoryLas |> anglePointCloud(understoryAngle) |>
      positionPointCloud(understoryX, understoryY)
    
    scenarioScan <- rbind(scenarioScan, understoryLas)
  }
  
  return(scenarioScan)
}

#Put all building objects in a line for displaying as a figure in the manuscript
syntheticObjectsPointCloud <- function(){
  finalCloud <- rbind(Understory[1,]$Las[[1]] |> positionPointCloud(0, -3),
                      binOneSnag$Las[[1]] |> positionPointCloud(0, 0),
                      binOneTree$Las[[1]] |> positionPointCloud(0, 3),
                      binTwoSnag$Las[[1]] |> positionPointCloud(0, 6),
                      binTwoTree$Las[[1]] |> positionPointCloud(0, 9),
                      binThreeSnag$Las[[1]] |> positionPointCloud(0, 12),
                      binThreeTree$Las[[1]] |> positionPointCloud(0, 15),
                      binFourSnag$Las[[1]] |> positionPointCloud(0, 18),
                      binFourTree$Las[[1]] |> positionPointCloud(0, 21),
                      binFiveSnag$Las[[1]] |> positionPointCloud(0, 24),
                      binFiveTree$Las[[1]] |> positionPointCloud(0, 27)
  )
  lidR::writeLAS(finalCloud, "tmp/syntheticObjects.las")
}