##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the effect of space competition between trees
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

# Space tree effect -------------------------------------------------------

competitionIntensity_v <- c(0.05, 0.45, 0.85)#seq(0.05 0.85, 0.05)
for(r in 1:numberRepetitionsSensitivity){

  for(c in 1:length(competitionIntensity_v)){
    
    set.seed(r)
    dateStart <- runif(numberTrees, 0, cycleLength)
    locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
    
    runSimulation(
      cycleLimitNumber = cycleLimitNumber,
      repetitionNumber = (r-1)*3 + c,
      timeDelayForDispersal = timeDelayForDispersal,
      torporTime = torporTime,
      saveTreeMap = saveTreeMap,
      samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
      nameInit = "Output/SensitivitySpaceTree/TestSpaceTree",
      mapSize = mapSize,
      quadratSize = quadratSize,
      numberTrees = numberTrees,
      #Wont be used
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
      intensityCompetitionForSpace = competitionIntensity_v[c]
    )
  }
}
