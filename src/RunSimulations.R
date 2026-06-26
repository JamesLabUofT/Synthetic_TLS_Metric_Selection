source("src/CreateSyntheticCloud.R")
source("src/FuelMetrics.R")
reticulate::source_python("src/python/computeLacunarity.py")

getScenarioString <- function(scenario){
  return(paste0("ID", scenario$ID))
}

getPointCloudsForScenarios <- function(ScenarioRDSFile){
  ScenarioSetName <- tools::file_path_sans_ext(basename(ScenarioRDSFile))
  
  ScenarioParameters <- readRDS(ScenarioRDSFile)
  
  Scenarios <- split(ScenarioParameters, seq(nrow(ScenarioParameters)))
  las<-createCloudFromScenario(Scenarios[[1]])

  do.call(rbind, parallel::mclapply(Scenarios, mc.preschedule = F, function(scenario){
    scenarioName <- paste0(ScenarioSetName, getScenarioString(scenario))
    fileName <- paste0("temp/pointCloud/", scenarioName, ".laz")
    if (!file.exists(fileName)) {
      las<-createCloudFromScenario(scenario)
      lidR::filter_poi(las, X**2 + Y**2 <= 25**2)
      lidR::writeLAS(las, fileName)
      
      gc()
    }
  }, mc.cores = 20))
}
getPointCloudsForScenarios("data/Scenarios/SyntheticScenarios.rds")

getMetricsForPointClouds <- function(ScenarioRDSFile){
  ScenarioSetName <- tools::file_path_sans_ext(basename(ScenarioRDSFile))
  
  ScenarioParameters <- readRDS(ScenarioRDSFile)
  
  Scenarios <- split(ScenarioParameters, seq(nrow(ScenarioParameters)))
  
  reticulate::source_python("src/python/computeLacunarity.py")
  
  gc()
  WilsonMetricsOutputs <- do.call(rbind, parallel::mclapply(Scenarios, mc.preschedule = F, function(scenario){
    outFile <- paste0("temp/simulationOutputs/",ScenarioSetName, scenario$ID, ".csv")
    if (!file.exists(outFile)) {
      scenarioName <- paste0(ScenarioSetName, getScenarioString(scenario))
      fileName <- paste0("temp/pointCloud/", scenarioName, ".laz")
      las<-lidR::readLAS(fileName)
      las <- lidR::filter_poi(las, X**2 + Y**2 <= 25**2)
      
      ENLOuts <- ENL(las)
      
      lacunarityOuts <- Lacunarity(las)
      
      canopyRatioOuts <- Canopy_ratio(las)
      
      wilsonOuts <- WilsonEtAl2021_FSG_CH_CBH_UH(las)
      
      metricOutputs <- cbind(ENLOuts, lacunarityOuts, canopyRatioOuts, wilsonOuts)
      
      metricOutputs$ID <- scenario$ID
      
      
      write.csv(metricOutputs, outFile)
    }
    gc()
  }, mc.cores = 35))
  
  outputFiles <- list.files("temp/simulationOutputs", pattern = paste0(ScenarioSetName, "*"), full.names = T)
  simulationOutputs <- read.csv(outputFiles[1])
  
  for (f in 2:length(outputFiles)){
    file <- outputFiles[f]
    output <- read.csv(file)
    simulationOutputs<-rbind(simulationOutputs, output)
  }
  
  
  write.csv(simulationOutputs, paste0("data/SimulationOutputsRaw/", ScenarioSetName, ".csv"))
  simulationOutputs <- select(simulationOutputs, -c("X"))
  
  simulationOutputs <- merge(simulationOutputs, ScenarioParameters)
  
  
  summarisedOutputs <- simulationOutputs %>% group_by(numTrees, amountSubcanopy, mortalityProportion, heightMean) %>% 
    summarise(fsg = mean(fsg), ch = mean(ch), cbh = mean(cbh), uh = mean(uh),
              totalLacunarity = mean(TotalLacunarity),heterogeneityIndex= mean(HeterogeneityIndex), oneDENL = mean(OneDENL), twoDENL = mean(TwoDENL),
              twoDENLCV = mean(TwoDENLCV), scanOneDENL = mean(ScanOneDENL), scanTwoDENL = mean(ScanTwoDENL), RH25 = mean(RH25), RH98 = mean(RH98),
              CR = mean(CR)) %>% ungroup()
  
  write.csv(summarisedOutputs, paste0("data/SimulationOutputsSummarised/", ScenarioSetName, ".csv"))
}

getMetricsForPointClouds("data/Scenarios/SyntheticScenarios.rds")