#installing packages
#------------------------------------------#

# this is just for loading examlpe date
# BiocManager::install("flowWorkspaceData")
#BiocManager::install("openCyto")
#flowWorkspaceData

#install packages from BioConductor
install.packages("BiocManager")
install.packages("devtools")

#packages for data cleaning
BiocManager::install("flowCore")
BiocManager::install("ggcyto")
BiocManager::install("flowWorkspace")
BiocManager::install("flowAI")
BiocManager::install("flowClean")
devtools::install_github("jmeskas/flowCut")

# Loading packages
#------------------------------------------#
library(ggcyto)
library(flowCore)
library(flowViz)
library(flowWorkspace)
library(flowAI)
library(flowClean)
library(flowCut)

