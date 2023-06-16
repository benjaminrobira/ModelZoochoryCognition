##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the effect of space competition between trees to counter balance moving rule effect (i.e. shrinkage when moving to targets)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script tests runs illustrative simulations to have the map of the counter effect of spacing tree on resource shrinkage
# 
# Setting up environment --------------------------------------------------

rm(list = ls())

## Parameters --------------------------------------------------------------

source("Scripts/R/0_Parameters.R")

## Libraries ---------------------------------------------------------------

library(Rcpp)
library(readr)
library(ggplot2)
library(ggpubr)
library(adehabitatHR)

## Functions ---------------------------------------------------------------

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")

plotMap <- function(
    df,
    valueSpacing){
  plot <- ggplot(df, aes(x = x, y = y)) +
    geom_circle(aes(x0 = x, y0 = y, 
                    r = valueSpacing * sqrt(mapSize * mapSize / (numberTrees * pi))), 
                linetype = "dotted") + 
    geom_point(size = 0.5) +
    labs(x = "", y = "") +
    theme_void() +
    theme(    
      strip.text.x = element_text(face = "bold", colour = "black"),
      strip.text.y = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(colour = NA, fill = "white"),
      panel.grid.minor = element_line(colour = "white"),
      panel.grid.major = element_line(colour = "white"),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    xlim(c(0, mapSize)) +
    ylim(c(0, mapSize)) +
    coord_fixed()
  return(plot)
}

estimateShrinkage <- function(outputMapFinal){
  
  boundary <- SpatialLines(
    list(
      cbind(c(0,0,1000,1000,0), c(0,1000,1000,0,0)) %>% 
        Line() %>% 
        Lines(ID = "boundary")
    )
  )
  
  #Get the necessary extend to have 95%
  outputTry = "KEEP GOING"
  extentUse = -0.05
  while(outputTry == "KEEP GOING"){
    print(extentUse)
    extentUse = extentUse + 0.05
    outputTry <- tryCatch(
      {
        UD <- kernelUD(xy = SpatialPoints(outputMapFinal[, 1:2]), h = 50,
                       boundary = boundary, extent = extentUse)
        shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/(mapSize*mapSize)
        "STOP"
      }, error = function(e){
        return("KEEP GOING")
      }
    )
  }
  
  UD <- kernelUD(xy = SpatialPoints(outputMapFinal[, 1:2]), h = 50,
                 boundary = boundary, extent = extentUse)
  shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/(mapSize*mapSize)
  return(shrinkage)
}

# Run simulations ---------------------------------------------------------

set.seed(1)
dateStart <- runif(numberTrees, 0, cycleLength)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))

# Low tree space competition

runSimulation(
  cycleLimitNumber = cycleLimitNumber,
  repetitionNumber = 1,
  timeDelayForDispersal = timeDelayForDispersal,
  torporTime = torporTime,
  saveTreeMap = saveTreeMap,
  samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
  nameInit = "Output/SensitivitySpaceTree/MovingTargetTreeSpaceLow",
  mapSize = mapSize,
  quadratSize = quadratSize,
  numberTrees = numberTrees,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + fruitingLength),
  homogeneousDistribution = homogeneousDistribution,
  treeClusterNumber = treeClusterNumber,
  treeClusterSpread = treeClusterSpread,
  maximumFoodToYield_v = maximumFoodToYield_v,
  cycleLength = cycleLength,
  fruitingLength = fruitingLength,
  noReturnTime = noReturnTime,
  whatValueUnknownTemporal = whatValueUnknownTemporal,
  whatRule = whatRule,
  exponentialRate = exponentialRate,
  perceptualRange = perceptualRange,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1,
  speed = speed,
  DispersalProbability = DispersalProbability,
  useProvidedMap = TRUE,
  moveOnlyToFruitingTrees = TRUE,
  moveOnlyToTarget = FALSE,
  intensityCompetitionForSpace = 0.85
)

# Intermediate tree space competition

runSimulation(
  cycleLimitNumber = cycleLimitNumber,
  repetitionNumber = 1,
  timeDelayForDispersal = timeDelayForDispersal,
  torporTime = torporTime,
  saveTreeMap = saveTreeMap,
  samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
  nameInit = "Output/SensitivitySpaceTree/MovingTargetTreeSpaceIntermediate",
  mapSize = mapSize,
  quadratSize = quadratSize,
  numberTrees = numberTrees,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + fruitingLength),
  homogeneousDistribution = homogeneousDistribution,
  treeClusterNumber = treeClusterNumber,
  treeClusterSpread = treeClusterSpread,
  maximumFoodToYield_v = maximumFoodToYield_v,
  cycleLength = cycleLength,
  fruitingLength = fruitingLength,
  noReturnTime = noReturnTime,
  whatValueUnknownTemporal = whatValueUnknownTemporal,
  whatRule = whatRule,
  exponentialRate = exponentialRate,
  perceptualRange = perceptualRange,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1,
  speed = speed,
  DispersalProbability = DispersalProbability,
  useProvidedMap = TRUE,
  moveOnlyToFruitingTrees = TRUE,
  moveOnlyToTarget = FALSE,
  intensityCompetitionForSpace = 0.85
)

# High Tree space competition

runSimulation(
  cycleLimitNumber = cycleLimitNumber,
  repetitionNumber = 1,
  timeDelayForDispersal = timeDelayForDispersal,
  torporTime = torporTime,
  saveTreeMap = saveTreeMap,
  samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
  nameInit = "Output/SensitivitySpaceTree/MovingTargetTreeSpaceHigh",
  mapSize = mapSize,
  quadratSize = quadratSize,
  numberTrees = numberTrees,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + fruitingLength),
  homogeneousDistribution = homogeneousDistribution,
  treeClusterNumber = treeClusterNumber,
  treeClusterSpread = treeClusterSpread,
  maximumFoodToYield_v = maximumFoodToYield_v,
  cycleLength = cycleLength,
  fruitingLength = fruitingLength,
  noReturnTime = noReturnTime,
  whatValueUnknownTemporal = whatValueUnknownTemporal,
  whatRule = whatRule,
  exponentialRate = exponentialRate,
  perceptualRange = perceptualRange,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1,
  speed = speed,
  DispersalProbability = DispersalProbability,
  useProvidedMap = TRUE,
  moveOnlyToFruitingTrees = TRUE,
  moveOnlyToTarget = FALSE,
  intensityCompetitionForSpace = 0.85
)


# Plot the maps -----------------------------------------------------------

outputMapFinalLow <- read_table2("Output/SensitivitySpaceTree/MovingTargetTreeSpaceLow_p15.811388_s1.000000_t1.000000_r1_Map_36500.txt") %>% as.data.frame() 
outputMapFinalIntermediate <- read_table2("Output/SensitivitySpaceTree/MovingTargetTreeSpaceIntermediate_p15.811388_s1.000000_t1.000000_r1_Map_36500.txt") %>% as.data.frame() 
outputMapFinalHigh <- read_table2("Output/SensitivitySpaceTree/MovingTargetTreeSpaceHigh_p15.811388_s1.000000_t1.000000_r1_Map_36500.txt") %>% as.data.frame() 

# 
# outputMapFinal <- read_table2("Output/SensitivityMovingRule/ReTestMovingTargetTrees_p15.811388_s1.000000_t1.000000_r1_Map_36500.txt") %>% as.data.frame() 
# plotMap(outputMapFinal, 0.45)
# estimateShrinkage(outputMapFinal)

#Low competition
mapLow <- plotMap(outputMapFinalLow, 0.05)
shrinkageLow <- estimateShrinkage(outputMapFinalLow)
mapLow <- mapLow + ggtitle(
  paste0("Spacing tree intensity: 0.05, Shrinkage: ", round(shrinkageLow, digits = 3))
)

#Intermediate competition
mapIntermediate <- plotMap(outputMapFinalIntermediate, 0.05)
shrinkageIntermediate <- estimateShrinkage(outputMapFinalIntermediate)
mapIntermediate <- mapIntermediate + ggtitle(
  paste0("Spacing tree intensity: 0.05, Shrinkage: ", round(shrinkageIntermediate, digits = 3))
)

#High competition
mapHigh <- plotMap(outputMapFinalHigh, 0.85)
shrinkageHigh <- estimateShrinkage(outputMapFinalHigh)
mapHigh <- mapHigh + ggtitle(
  paste0("Spacing tree intensity: 0.85, Shrinkage: ", round(shrinkageHigh, digits = 3))
)

# Save the maps -----------------------------------------------------------

allMaps <- ggarrange(
  mapLow,
  mapHigh,
  ncol = 2
)
saveRDS(allMaps, "Renvironment/Plots/mapsMovingTargetSpaceTreeLoworHigh.rds")
