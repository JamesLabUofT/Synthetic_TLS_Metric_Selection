Subcanopy <- read.csv("/Users/leo/Forestry/JPBW_Data/FieldData/AllVegSamplingSubcanopy.csv")
Subcanopy_count <- Subcanopy %>% filter(PROJECT == "JPBW") %>% group_by(SITEPLOT) %>% summarise (n = n()) %>% ungroup()
Empty_plots <- Subcanopy %>% filter(State == "")
Subcanopy_count$n[Subcanopy_count$SITEPLOT %in% Empty_plots$SITEPLOT] <- 0

CanopyData <- read.csv("/Users/leo/Forestry/JPBW_Data/FieldData/AllVegSamplingCanopyFix.csv")
Canopy <- CanopyData %>% filter(PROJECT == "JPBW") %>% group_by(SITEPLOT) %>% summarise(numTrees = n(), mortality = 1-sum(State == "A") / n(), meanHeight = mean(height)) 

PlotsDetails <- cbind(as.numeric(Canopy$numTrees), 
                      as.numeric(Canopy$mortality), as.numeric(Canopy$meanHeight), 
                      as.numeric(Subcanopy_count$n))
colnames(PlotsDetails) <- c("NumTrees", "Mortality", "MeanHeight", "NumSub")
PlotsDetails <- as.data.frame(PlotsDetails)

### Uniform sampling

#Get all combinations of points in the 4D space
AllNumTrees <- round(scales::rescale(1:15, to = c(min(PlotsDetails$NumTrees), max(PlotsDetails$NumTrees))))
AllMortality <- round(scales::rescale(1:15, to = c(min(PlotsDetails$Mortality), max(PlotsDetails$Mortality))), digits = 2)
AllMeanHeight <- round(scales::rescale(1:15, to = c(min(PlotsDetails$MeanHeight), max(PlotsDetails$MeanHeight))), digits = 2)
AllNumSub <- round(scales::rescale(1:15, to = c(min(PlotsDetails$NumSub) ,max(PlotsDetails$NumSub))))

AllOptions <- expand.grid(AllNumTrees, AllMortality, AllMeanHeight, AllNumSub)
colnames(AllOptions) <- c("NumTrees", "Mortality", "MeanHeight", "NumSub")
ConvHull <- geometry::convhulln(PlotsDetails)

### Edge cases are not well represented. I need to expand the data that makes the convex hull to enlarge it just a bit


AddBuffer <- function(PlotsDetails, NumTreeBuffer, MortalityBuffer, MeanHeightBuffer, NumSubBuffer){
  PlotsDetailsWBuffer <- PlotsDetails
  for (plot in 1:nrow(PlotsDetails)){
    plotDetails <- PlotsDetails[plot,] 
    for (a in c(-1,1)){
      for (b in c(-1, 1)){
        for (c in c(-1,1)){
          for (d in c(-1,1))
            bufferPlot <- data.frame(NumTrees = plotDetails$NumTrees + a * NumTreeBuffer,
                                     Mortality = plotDetails$Mortality + b * MortalityBuffer,
                                     MeanHeight = plotDetails$MeanHeight + c * MeanHeightBuffer,
                                     NumSub = plotDetails$NumSub + d * NumSubBuffer)
            PlotsDetailsWBuffer <- rbind(PlotsDetailsWBuffer, bufferPlot)
          }
        }
      }
  }
  return(PlotsDetailsWBuffer)
  }

BufferPoints <- PlotsDetails[sort(unique(append(append(append(ConvHull[,4], ConvHull[,3]), ConvHull[,2]), ConvHull[,1]))),]
NewBuffer <- AddBuffer(BufferPoints,
                       (AllNumTrees[2] - AllNumTrees[1])/2,
                       (AllMortality[2] - AllMortality[1])/3,
                       (AllMeanHeight[2] - AllMeanHeight[1])/2,
                       (AllNumSub[2] - AllNumSub[1])/4)
NewBufferConvHull <- geometry::convhulln(NewBuffer)
NewerValidOptions <- AllOptions[-which(!geometry::inhulln(NewBufferConvHull, as.matrix(AllOptions))),]

NewerParameters <- anti_join(NewerValidOptions, NewValidOptions, by = NULL)
NewerParameters$NumSub <- NewerParameters$NumSub / (pi * 3.99*3.99)
NewerParameters$NumTrees <- NewerParameters$NumTrees / (pi*11.28**2)
saveRDS(NewerParameters, "data/ParameterSpace.rds")


nll_function <- function(sd_value, data) {
  data %>%
    group_by(SITEPLOT) %>%
    summarize(
      nll = sum(-dnorm(height, mean = mean(height), sd = sd_value, log = TRUE))
    ) %>%
    pull(nll) %>%
    sum()
}

# Optimize SD
optimal_sd <- optimize(
  nll_function,
  interval = c(0.5, 10), # Search range for SD
  data = EachCanopy
)$minimum

# Print optimal SD
cat("Optimal SD:", optimal_sd, "\n")

