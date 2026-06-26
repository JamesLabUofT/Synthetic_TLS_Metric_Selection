library(dplyr)
Lacunarity <- function(las){
  
  voxSize <- 0.5
  
  las <- lidR::voxelize_points(las, voxSize)
  
  
  
  minX <- las@header$`Min X`
  minY <- las@header$`Min Y`
  minZ <- las@header$`Min Z`
  
  las$X <- (las$X - minX) * 1/voxSize + 1
  las$Y <- (las$Y - minY) * 1/voxSize + 1 
  las$Z <- (las$Z - minZ) * 1/voxSize + 1
  
  newArr = array(dim = c(max(las$X), max(las$Y), max(las$Z)), data = 0)
  
  newArr[cbind(las$X, las$Y, las$Z)] <- 1
  
  lacunarity <- computeLacunarity(newArr)
  
  return(tibble(HeterogeneityIndex = lacunarity[[1]], #Scott et al., 2022
                 TotalLacunarity = lacunarity[[2]])) #Cooper 2022
}

#Ehbrecht et al., 2016
ENL <- function(las){
  #Voxelize by 20cm
  las <- lidR::voxelize_points(las, 0.2)
  
  computeENL <- function(Z){
    #Count number of voxels in 1m slices, divided by total number of voxels (an array p)
    p <- table(ceiling(Z))/length(Z)
    
    #Zero_D_ENL <- length(p) This is not an interesting metric, but including it for the sake of completeness
    One_D_ENL <- exp(-sum(p*log(p)))
    Two_D_ENL <- 1/sum(p**2)
    
    return(tibble(OneDENL = One_D_ENL, TwoDENL = Two_D_ENL))
  }
  
  df <- select(las@data, X, Y, Z)
  df <- round(df, 3)
  
  df$X <- ceiling(df$X) #Summarize into 1mx1m columns (with 0.2m horizontal slices)
  df$Y <- ceiling(df$Y)
  
  
  out <- df %>% group_by(X, Y) %>% 
    reframe(computeENL(Z)) %>% 
    ungroup() %>% summarise(OneDENL = mean(OneDENL), TwoDENLCV = sd(TwoDENL)/mean(TwoDENL), TwoDENL = mean(TwoDENL))
  ScanWideENL <- computeENL(df$Z)
  out$ScanOneDENL <- ScanWideENL$OneDENL
  out$ScanTwoDENL <- ScanWideENL$TwoDENL
  return(out)
}

#Aalto et al. 2023
Canopy_ratio <- function(las){
  las <- lidR::voxelize_points(las, 1)
  
  df <- select(las@data, X, Y, Z)
  df <- round(df, 3)
  
  
  
  getRelativeHeights <- function(Z){
    
    RelativeHeights <- quantile(Z, c(.25, .98))
    return(
      tibble(RH25 = RelativeHeights[1][[1]],
             RH98 = RelativeHeights[2][[1]])
    )
  }
  
  df <- df %>% group_by(X, Y) %>% 
    reframe(getRelativeHeights(Z)) %>% 
    ungroup() %>% summarise(RH25 = mean(RH25), RH98 = mean(RH98), 
                            CR = (mean(RH98) - mean(RH25))/mean(RH98))
  
  return(df)
}


#Wilson et al., 2021
WilsonEtAl2021_FSG_CH_CBH_UH <- function(las){
  voxelSize <- 0.25
  las <- lidR::voxelize_points(las, voxelSize)
  
  height_col <- sort(unique(las@data$Z))
  
  df <- las@data
  
  df <- round(df, 3)
  
  fuelStrataGap <- function(Z){
    rle <- rle(as.numeric(height_col %in% Z))
    if (rle$values[length(rle$values)] == 0){
      rle$lengths <- rle$lengths[1:length(rle$lengths) - 1]
      rle$values <- rle$values[1:length(rle$values) - 1]
    }
    rle$cumsum <- cumsum(rle$lengths)
    
    orderedLength <- order(rle$lengths)
    
    largest_empty_index <- orderedLength[max(which(orderedLength %in% which(rle$values == 0)))]
    fsg <- rle$lengths[largest_empty_index]* voxelSize
    ch <- rle$cumsum[length(rle$cumsum)] * voxelSize
    cbh <- rle$cumsum[largest_empty_index] * voxelSize
    uh <- rle$cumsum[largest_empty_index - 1] * voxelSize
    return(
      tibble(fsg = fsg,
             ch = ch,
             cbh = cbh,
             uh = uh)
    )
  }
  
  df <- df %>% group_by(X, Y) %>% 
    reframe(fuelStrataGap(Z)) %>% na.omit %>% 
    ungroup() %>% summarise(fsg = median(fsg), ch = median(ch), 
                            cbh = median(cbh), uh = median(uh))
  
  return(df)
}
