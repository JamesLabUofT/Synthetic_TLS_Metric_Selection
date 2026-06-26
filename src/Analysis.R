library(dplyr)
library(ggplot2)

CORRELATION_THRESHOLD <- 0.8

dataRaw <- read.csv("data/SimulationOutputsRaw/SyntheticScenarios.csv")
dataRaw <- dplyr::select(dataRaw, -c("X", "X.1"))


ScenarioParameters <- readRDS("data/Scenarios/SyntheticScenarios.rds")
ScenarioParameters  <- dplyr::select(ScenarioParameters, c("heightMean", "mortalityProportion", "numTrees", "amountSubcanopy", "ID"))

data <- merge(dataRaw, ScenarioParameters, by = "ID")

stopifnot(nrow(data) == 7587)


explanatory <- data[, c(2:8, 11:15)]
MRTexplanatory <- data[, c(2:8, 11:15)]
colnames(explanatory) <- c("1D ENL", "2D ENL CV", "2D ENL", "Scan-wide 1D ENL", "Scan-wide 2D ENL",
                           "Heterogeneity Index", "Total Lacunarity", "CR", "FSG", "CH", "CBH", "UH") #More friendly names

explanatory <- explanatory[, c("CH", "CBH", "FSG", "UH", "CR", "1D ENL", "2D ENL", "2D ENL CV", "Scan-wide 1D ENL", "Scan-wide 2D ENL",
                               "Total Lacunarity", "Heterogeneity Index")] #Reorder to stay consistent with manuscript

response <- data[, 16:19]
MRTresponse <- data[, 16:19]
colnames(response) <- c("Mean Tree Height", "Mortality Proportion", "Stem Density", "Understory Density")
response <- response[, c("Stem Density", "Understory Density", "Mean Tree Height", "Mortality Proportion")]

# Unsupervised - Correlation

get_correlation_heatmap <- function(variables){
  correlation_matrix <- round(cor(variables), 2)
  melted_cormat <- reshape2::melt(correlation_matrix)
  colnames(melted_cormat) <- c("Var1", "Var2", "Pearson correlation coefficient")
  
  var_names <- colnames(variables)
  
  label_map <- setNames(var_names, var_names)
  label_map["FSG"] <- "bold('FSG')"
  label_map["UH"] <- "bold('UH')"
  label_map["CR"] <- "bold('CR')"
  label_map["1D ENL"] <- "phantom()^{1}*'D ENL'"
  label_map["2D ENL"] <- "phantom()^{2}*'D ENL'"
  label_map["2D ENL CV"] <- "bold(phantom()^{2}*'D ENL CV')"
  label_map["Scan-wide 1D ENL"] <- "bold('Scan-wide '^{1}*'D ENL')"
  label_map["Scan-wide 2D ENL"] <- "'Scan-wide '^{2}*'D ENL'"
  label_map["Total Lacunarity"] <- "bold('Total Lacunarity')"
  label_map["Heterogeneity Index"] <- "'Heterogeneity Index'"
  
  return(ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=`Pearson correlation coefficient`)) + 
    geom_tile(color = "white") +
    scale_x_discrete(labels = function(x) parse(text = label_map[x], )) +
    scale_y_discrete(labels = function(x) parse(text = label_map[x])) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab") +
    theme_minimal() +
    theme(axis.text.y = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.text.x = element_text(angle = 45,vjust = 1, 
                                     size = 12, hjust = 1), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())+
    coord_fixed() +
    geom_text(aes(Var1, Var2, label = `Pearson correlation coefficient`), color = "black", size = 4))
}

Correlation_heatmap <- get_correlation_heatmap(explanatory)

noCorExp <- explanatory[,-c(sort(caret::findCorrelation(cor(explanatory), CORRELATION_THRESHOLD)))]

NoCorExp_heatmap <- get_correlation_heatmap(noCorExp)

### threshold sensitivity
thresholds <- seq(0, 1, by = 0.001)
numMetrics <- sapply(thresholds, function(x){length(caret::findCorrelation(cor(explanatory), x))}, simplify = "array")

df <- data.frame(Threshold = thresholds, Retained_Metrics = 12 - numMetrics)

ggplot(df, aes(Threshold, Retained_Metrics)) +
  geom_point() +
  labs(x = "Pearson correlation threshold", y = "Number of retained metrics") +
  geom_vline(xintercept = 0.8, col = "red") + scale_y_continuous(breaks=seq(0,12,2)) +
  scale_x_reverse(breaks=seq(0,1,0.2) ) + theme_minimal() 
  
# Supervised - Multivariate Regression Tree
noCorMRTExp <- MRTexplanatory[,-c(sort(caret::findCorrelation(cor(MRTexplanatory), CORRELATION_THRESHOLD)))]

scaledMRTresponse <- MRTresponse


scaledMRTresponse$numTrees <-  as.numeric(as.factor(MRTresponse$numTrees)) - 1
scaledMRTresponse$amountSubcanopy <- as.numeric(as.factor(MRTresponse$amountSubcanopy)) - 1
scaledMRTresponse$mortalityProportion <- as.numeric(as.factor(MRTresponse$mortalityProportion)) - 1
scaledMRTresponse$heightMean <- as.numeric(as.factor(MRTresponse$heightMean)) - 1

mrt <- mvpart::mvpart(as.matrix(scaledMRTresponse) ~ ., data = noCorMRTExp,
                      xval = nrow(scaledMRTresponse),
                      xvmult = 100, xv = "pick",
                      which = 4, legend = T,
                      prn = FALSE)

mrtGroups<- as.numeric(as.factor(mrt$where))



rdaNoCor <- vegan::rda(response ~ FSG +  `Scan-wide 1D ENL`+ CR + `2D ENL CV` + UH +`Total Lacunarity`, data = explanatory)
rdaMRT <- vegan::rda(response ~ FSG + CR + `2D ENL CV` + `Scan-wide 1D ENL`, data = explanatory)
rdaFullModel <- vegan::rda(response ~ ., data = explanatory)
 
RsquareAdj(rdaNoCor)$adj.r.squared
RsquareAdj(rdaMRT)$adj.r.squared
RsquareAdj(rdaFullModel)$adj.r.squared

noCorVarImportance <- anova.cca(rdaNoCor, by = "term") 
anova.cca(rdaMRT, permutations = 999)

anova.cca(rdaFullModel, rdaNoCor) 
anova.cca(rdaFullModel, rdaMRT)
anova.cca(rdaNoCor, rdaMRT)

# Variable importance
variableNames <- c("FSG",  "Scan-wide 1D ENL", "CR", "2D ENL CV", "UH", "Total Lacunarity", "Residual")
anova_results <- data.frame(Variable = variableNames, F_value = noCorVarImportance$F)
anova_results$label <- c("*","*", "*","*","","","")

# Plot the F-values
ggplot(anova_results[1:nrow(anova_results)-1,], aes(x = reorder(Variable, F_value), y = F_value, fill = F_value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Metric", y = "F-value") +
  scale_fill_gradient(low = "pink", high = "red") +
  theme(legend.position = "none") +
  geom_text(aes(label = label), hjust = -0.5, size = 6)



# Representative parametrs for each cluster 
clusterParameterValues <- as.data.frame(as.data.frame(mrt$frame %>% filter(var == "<leaf>") %>% select(c(yval2)))$yval2)
colnames(clusterParameterValues)<- c("MeanHeight", "Mortality", "NumTrees", "NumSub")

clusterParameterValues[12,] <- c(0,0,0,0)
clusterParameterValues[13,] <- c(14,14,14,14)

clusterParameterValues$MeanHeight <- scales::rescale(clusterParameterValues$MeanHeight, to = c(min(MRTresponse$heightMean), max(MRTresponse$heightMean)))
clusterParameterValues$Mortality <- scales::rescale(clusterParameterValues$Mortality, to = c(min(MRTresponse$mortalityProportion), max(MRTresponse$mortalityProportion)))
clusterParameterValues$NumTrees <- scales::rescale(clusterParameterValues$NumTrees, to = c(min(MRTresponse$numTrees), max(MRTresponse$numTrees)))
clusterParameterValues$NumSub <- scales::rescale(clusterParameterValues$NumSub, to = c(min(MRTresponse$amountSubcanopy), max(MRTresponse$amountSubcanopy)))


clusterParameterValues$NumSub <- clusterParameterValues$NumSub / (pi * 27**2)
clusterParameterValues$NumTrees <- clusterParameterValues$NumTrees / (pi*27**2)


### Get boxplots of cluster parameters
clusterParameterValues$ID <- c(1:13)
clusterParameterValues$Letter <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M")
clusterParameterValues$Mortality <- clusterParameterValues$Mortality * 100
clusterParameterValues$Letter <- factor(clusterParameterValues$Letter, levels = rev(clusterParameterValues$Letter)) 

  
clusterNumTrees <- ggplot(clusterParameterValues[1:11,], aes(x=Letter, y = NumTrees, color = Letter)) +
  geom_segment(aes (x=Letter, xend=Letter, y=clusterParameterValues[12,]$NumTrees, yend=NumTrees)) +
  geom_point(size = 4, alpha = 0.6) + coord_flip() + scale_color_viridis_d() +
  ylim(clusterParameterValues[12,]$NumTrees, clusterParameterValues[13,]$NumTrees) +
  labs(y = expression("Stem density (#/m"^2*")"), x = "Cluster") + theme_minimal() + theme(legend.position = "")
  
clusterMeanHeight <- ggplot(clusterParameterValues[1:11,], aes(x=Letter, y = MeanHeight, color = Letter)) +
  geom_segment(aes (x=Letter, xend=Letter, y=clusterParameterValues[12,]$MeanHeight, yend=MeanHeight)) +
  geom_point(size = 4, alpha = 0.6) + coord_flip() + scale_color_viridis_d() +
  ylim(clusterParameterValues[12,]$MeanHeight, clusterParameterValues[13,]$MeanHeight) +
  labs(y = "Mean tree height (m)", x = "Cluster") + theme_minimal() + theme(legend.position = "")

clusterNumSub <- ggplot(clusterParameterValues[1:11,], aes(x=Letter, y = NumSub, color = Letter)) +
  geom_segment(aes (x=Letter, xend=Letter, y=clusterParameterValues[12,]$NumSub, yend=NumSub)) +
  geom_point(size = 4, alpha = 0.6) + coord_flip() + scale_color_viridis_d() +
  ylim(clusterParameterValues[12,]$NumSub, clusterParameterValues[13,]$NumSub) +
  labs(y = expression("Understory density (#/m"^2*")"), x = "Cluster") + theme_minimal() + theme(legend.position = "")

clusterMortality <- ggplot(clusterParameterValues[1:11,], aes(x=Letter, y = Mortality, color = Letter)) +
  geom_segment(aes (x=Letter, xend=Letter, y=clusterParameterValues[12,]$Mortality, yend=Mortality)) +
  geom_point(size = 4, alpha = 0.6) + coord_flip() + scale_color_viridis_d() +
  ylim(clusterParameterValues[12,]$Mortality, clusterParameterValues[13,]$Mortality) +
  labs(y = "Stem mortality (%)", x = "Cluster") + theme_minimal() + theme(legend.position = "")

allClusterParams <- ggpubr::ggarrange(clusterNumTrees, clusterNumSub, clusterMeanHeight, clusterMortality, nrow = 2, ncol = 2)

