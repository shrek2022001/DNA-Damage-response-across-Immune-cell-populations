install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c(
  "CATALYST",
  "SingleCellExperiment",
  "ComplexHeatmap"
))
install.packages("uwot")
BiocManager::install("FlowSOM")
BiocManager::install("CATALYST", dependencies = TRUE)

library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(uwot)
library(FlowSOM)

