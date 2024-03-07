##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the effect of moving rule
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script runs the simulations modifying the agents' learning locations when new seeds replace known trees.

# Setting up environment --------------------------------------------------

## Parameters --------------------------------------------------------------

source("Scripts/R/0_Parameters.R")
numberRepetitionsSensitivity = numberRepetitions

## Libraries ---------------------------------------------------------------

library(Rcpp)
library(readr)

## Functions ---------------------------------------------------------------

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")

#Note: For the following simulations I am only saving the map at the last time. I also rerun with the same map for testing the moving rule.

# Moving rule difference: moving to all trees vs fruiting trees -----------

for(r in 118:numberRepetitionsSensitivity){
  print(r)
  set.seed(r)
  dateStart <- runif(numberTrees, 0, cycleLength)
  locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
  
  runSimulationMemoryRevised(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityLearning/TestLearning",
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
    spatialKnowledgeRate = 0.25,
    temporalKnowledgeRate = 0.25,
    speed = speed/2,
    DispersalProbability = DispersalProbability,
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = FALSE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  
  set.seed(r)# restart everything to be sure the seed use remains the same
  dateStart <- runif(numberTrees, 0, cycleLength)
  locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
  
  runSimulationMemoryRevised(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityLearning/TestLearning",
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
    spatialKnowledgeRate = 0.5,
    temporalKnowledgeRate = 0.5,
    speed = speed*2,
    DispersalProbability = DispersalProbability,
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = FALSE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  
  set.seed(r)# restart everything to be sure the seed use remains the same
  dateStart <- runif(numberTrees, 0, cycleLength)
  locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
  
  runSimulationMemoryRevised(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityLearning/TestLearning",
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
    spatialKnowledgeRate = 0.75,
    temporalKnowledgeRate = 0.75,
    speed = speed*2,
    DispersalProbability = DispersalProbability,
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = FALSE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
}
