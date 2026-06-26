library(ggplot2)
library(dplyr)
RealAndSimulatedEquivPlots <- readRDS("data/RealAndSimulatedEquivPlots.rds")
getzscore <- function(metric){
  real = metric[length(metric)]
  sim = metric[1:length(metric)-1]
  zscore <- (real - mean(sim))/sd(sim)
  return(zscore)
}
zscores_df <- RealAndSimulatedEquivPlots %>% group_by(ID) %>%
  summarise(fsg = getzscore(fsg),
  CR = getzscore(CR),
  TotalLacunarity = getzscore(TotalLacunarity),
  TwoDENLCV = getzscore(TwoDENLCV),
  ScanOneDENL = getzscore(ScanOneDENL),
  UH = getzscore(uh))

var_names <- c( "fsg", "UH", "TotalLacunarity", "TwoDENLCV","ScanOneDENL", "CR")

# Create a named viridis color palette
viridis_colors <- viridis::viridis(7)

# Generate boxplots dynamically
plots <- lapply(seq_along(var_names), function(i) {
  if (var_names[i] == "fsg"){
    metric_label <- expression("FSG z-scores")
    label = "*"
    lim <- 8.5
  } else if (var_names[i] == "CR") {
    metric_label <- expression("CR z-scores")
    label = "*"
    lim <- 2.5
  } else if (var_names[i] == "TwoDENLCV") {
    metric_label <- expression(""^2*"D ENL CV z-scores")
    label = "*"
    lim <- 2.5
  } else if (var_names[i] == "ScanOneDENL") {
    metric_label <- expression("Scan-wide "^1*"D ENL z-scores")
    label = "*"
    lim <- 2.5
  } else if (var_names[i] == "TotalLacunarity") {
    metric_label <- expression("Total Lacunarity z-scores")
    label = ""
    lim <- 8.5
  } else {
    metric_label <- expression("UH z-scores")
    label = ""
    lim <- 8.5}
  p <- ggplot(zscores_df, aes(x = var_names[i], y = .data[[var_names[i]]])) +
    geom_boxplot(fill = viridis_colors[i+1]) +
    labs(y = metric_label, x = "") +
    theme_minimal() +
    ylim(-lim, lim) +
    theme(axis.text.x = element_blank())
  p + geom_text(data=(data.frame(y=lim, label=label)),
                aes (y=y, label=label),
                size = 6, fontface= "bold", hjust = 6, vjust=0.7)
})

ggpubr::ggarrange(plotlist = plots, ncol = 3, nrow = 2)
