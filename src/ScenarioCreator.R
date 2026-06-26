library(dplyr)


##ADAPTED FROM the PRACMA package
#Randomly chooses points on a disk such that there are n points, the disk has a radius r
#and it attempts m points at once
poisson_disk_sampling <- function (n, r = 1, m = 10, info = TRUE) 
{
  if (r <= 0) 
    stop("Radius r must be positive reals.")
  if (floor(n) != ceiling(n) || n < 1 || floor(m) != ceiling(m) || 
      m < 1) 
    stop("n and m must be integer numbers.")
  ab <- rep(c(2*r, 2*r), times = c(m, m))
  A <- matrix(0, n, 2)
  
  A[1, ] <- c(2*r, 2*r) * runif(2) #Gets a random point in a box of size 2*r x 2*r
  while (sqrt((A[1,1] - r)**2 + (A[1,2] - r)**2) > r){ #ensures the first point is within the desired radius 
    A[1, ] <- c(2*r, 2*r) * runif(2) #if not, get another one
  }
  i <- 2
  iters <- 0
  while (i <= n) {
    
    if (iters > 100 * n){
      stop("Number of attempts indicates the constraints might be impossible to fulfill")
    }
    iters <- iters + 1
    
    B <- matrix(ab * runif(2 * m), nrow = m, ncol = 2) #Get m random points
    
    g <- which(sqrt((B[,1] - r)**2 + (B[,2] - r)**2) > r) #Remove the points that are outside the disk
    B <- B[-g,]
    
    if (nrow(B) == 0 || is.null(nrow(B))) { #if B is empty, all the random points were outside the disk
      next #Try again
    }
    
    C <- pracma::distmat(B, A[1:(i - 1), ]) #Matrix with the euclidean distances between the random points in the disk and the previously selected points
    k <- which.max(apply(C, 1, min)) #Get the index of the point whose closest point (apply min) is the furthest (which.max)
    A[i, ] <- B[k, ] #Add that point to our final list
    i <- i + 1
  }
  if (info) {
    AA <- pracma::distmat(A, A)
    diag(AA) <- max(AA)
    d <- sqrt(2 * r * r/n)
    cat("Minimal Distance between points: ", min(AA), "\n")
  }
  return(A - c(r,r))
} 

poisson_disk_sampling_V <- Vectorize(poisson_disk_sampling, vectorize.args = "n")

##further ADAPTED FROM the PRACMA package
#Randomly chooses points on a disk such that there are n points, the disk has a radius r
#it attempts m points at once AND all points are at least OtherDist distance away from any given point
#in the matrix Otherpoints.
poisson_disk_sampling_extra_constraints <- function (n, r = 1, m = 10, Otherpoints, OtherDist, info = TRUE) 
{
  if (r <= 0) 
    stop("Radius r must be positive reals.")
  if (floor(n) != ceiling(n) || n < 1 || floor(m) != ceiling(m) || 
      m < 1) {
    if (n == 0){
      return(NA)
    }
    stop("n and m must be integer numbers.")
  }
  ab <- rep(c(2*r, 2*r), times = c(m, m))
  A <- matrix(0, n, 2)
  A[1, ] <- c(2*r, 2*r) * runif(2)
  while (sqrt((A[1,1] - r)**2 + (A[1,2] - r)**2) > r || #ensures the first point is within the desired radius
         length(which(pracma::distmat(A[1,], Otherpoints) < OtherDist )) != 0){ #ensures the first point is far enough from the other points 
    A[1, ] <- c(2*r, 2*r) * runif(2) #if not, get another one
  }
  iters <- 0
  i <- 2
  while (i <= n) {
    
    if (iters > 100 * n){
      stop("Number of attempts indicates the constraints might be impossible to fulfill")
    }
    iters <- iters + 1
    
    B <- matrix(ab * runif(2 * m), nrow = m, ncol = 2) #Get m random points
    
    g <- which(sqrt((B[,1] - r)**2 + (B[,2] - r)**2) > r) #Remove the points that are outside the disk
    B <- B[-g,]
    
    D <- pracma::distmat(B, Otherpoints) #Get the distance between the m random points and the points in Otherpoints
    h <- unique(which(D<OtherDist, arr.ind = T)[, 1]) #Get the indices of the points that are too close to the points in Otherpoints
    B <- B[-h,] #Remove them
    
    if (nrow(B) == 0 || is.null(nrow(B))) { #if B is empty, all the random points were outside the disk
      next #Try again
    }
    
    C <- pracma::distmat(B, A[1:(i - 1), ]) #Matrix with the euclidean distances between the random points in the disk and the previously selected points
    k <- which.max(apply(C, 1, min)) #Get the index of the point whose closest point (apply min) is the furthest (which.max)
    A[i, ] <- B[k, ] #Add that point to our final list
    i <- i + 1
  }
  if (info) {
    AA <- pracma::distmat(A, A)
    diag(AA) <- max(AA)
    d <- sqrt(2 * r * r/n)
    cat("Minimal Distance between points: ", min(AA), "\n")
    
    AOther <- pracma::distmat(A, Otherpoints)
    cat("Minimal Distance with other points: ",min(AOther), "\n")
  }
  return(A - c(r,r))
}

poisson_disk_sampling_extra_constraints_V <- Vectorize(poisson_disk_sampling_extra_constraints, vectorize.args = c("n", "Otherpoints"))

sample_tree_heights <- function(mean_index, height_options, n, sd = 4.4) {
  weights <- dnorm(height_options, mean = mean_index, sd)
  indices <- 1:length(height_options)
  
  probability_distribution <- dnorm(height_options, mean = mean_index, sd = sd)
  probability_distribution <- probability_distribution/sum(probability_distribution)
  
  sampled_indices <- sample(height_options, size = n, replace = T, prob = probability_distribution)
  return(sampled_indices)
}

sample_tree_heights_v <- Vectorize(sample_tree_heights, vectorize.args = c("mean_index", "n"))

sample_dead_stems <- function(stem_probability, n) {
  sample(0:1, size = n, replace = T, prob = c(stem_probability, 1-stem_probability))
}

sample_dead_stems_v <- Vectorize(sample_dead_stems)

sample_angles <- function(number_trees) {
  angle <- sample(36000,number_trees, replace = T)/100
  return(angle)
}

sample_angles_v <- Vectorize(sample_angles)

sample_understory_IDs <- function(number_understory) {
  return(sample(c(1:3), number_understory, replace = T))
}

sample_understory_IDs_v <- Vectorize(sample_understory_IDs)


createScenarioParameters <- function(ParameterSpaceOptions){
  PlotRadius <- 27
  NumberMeterSquared <- pi * PlotRadius**2
  
  
  #SCENARIO PARAMETERS
  TreeHeightDist <- sort(unique(Trees$Height))
  TreeHeightMean <- ParameterSpaceOptions$MeanHeight
  Mortality <- ParameterSpaceOptions$Mortality
  NumberOfTrees <- round(ParameterSpaceOptions$NumTrees * NumberMeterSquared)
  SubcanopyAmount <- round(ParameterSpaceOptions$NumSub * NumberMeterSquared)
  
  ScenarioParameters <- data.frame(heightMean = TreeHeightMean, 
                                    mortalityProportion = Mortality, numTrees = NumberOfTrees, 
                                    amountSubcanopy = SubcanopyAmount)
  
  ScenarioParameters$ID <- 1:nrow(ScenarioParameters)
  
  set.seed(4018)
  
  
  ScenarioParameters <- ScenarioParameters %>% mutate(TreePositions = poisson_disk_sampling_V(numTrees, PlotRadius, 10)) %>% 
    mutate(UnderstoryID = sample_understory_IDs_v(amountSubcanopy)) %>% 
    mutate(UnderstoryPositions = poisson_disk_sampling_extra_constraints_V(amountSubcanopy, PlotRadius, 500, TreePositions, 1)) %>% 
    mutate(TreeHeights = sample_tree_heights_v(heightMean, TreeHeightDist, numTrees)) %>% 
    mutate(Aliveness = sample_dead_stems_v(mortalityProportion, numTrees)) %>% 
    mutate(TreeAngles = sample_angles_v(numTrees)) %>% 
    mutate(UnderstoryAngle = sample_angles_v(amountSubcanopy))
  
  saveRDS(ScenarioParameters, file = "data/Scenarios/SyntheticScenarios.rds") 
}

createScenarioParameters(readRDS("data/ParamaterSpace.rds"))
