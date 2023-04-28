##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the effect of moving rule
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script tests the different Rcpp function used for running the simulation in order to see whether those were doing what intended.

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

for(r in 1:numberRepetitionsSensitivity){
  print(r)
  set.seed(r)
  dateStart <- runif(numberTrees, 0, cycleLength)
  locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))

  runSimulation(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityMovingRule/TestMovingFruitTrees",
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
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  
  runSimulation(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityMovingRule/TestMovingNotFruitTrees",
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
    moveOnlyToFruitingTrees = FALSE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  
  runSimulation(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
    nameInit = "Output/SensitivityMovingRule/TestMovingTargetTrees",
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
    moveOnlyToFruitingTrees = FALSE,
    moveOnlyToTarget = TRUE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
}
