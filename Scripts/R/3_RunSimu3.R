#~~~~~~~~
#Running simulation
#~~~~~~~~

# What does this script do?
# This script runs the simulation for a given spatio-temporal knowledge rate (parameter s chooses which level). This is because it relies on Rcpp which I failed to parallelise (which seems possible only for packages at time of coding).

# Setting up environment --------------------------------------------------

##Libraries ---------------------------------------------------------------

library(Rcpp)
library(readr)
#Parallelising #https://www.r-bloggers.com/2017/09/a-guide-to-parallelism-in-r/
library(parallel)
library(doParallel)


# Loading parameters ------------------------------------------------------
load("Scripts/R/Parameterisation.RData")
spatialKnowledge <- c(0, 0.25, 0.5, 0.75, 1)
temporalKnowledge <- c(0, 0.25, 0.5, 0.75, 1)
s = 3

# Functions ---------------------------------------------------------------
Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")

# Run analysis: Main foraging ---------------------------------------------

#Trick to parallel loop despite using C++ which otherwise causes problem: task 1 failed - "valeur NULL passée comme adresse symbolique"
#https://stackoverflow.com/questions/18245193/doparallel-issue-with-inline-function-on-windows-7-works-on-linux
#https://stackoverflow.com/questions/25062383/cant-run-rcpp-function-in-foreach-null-value-passed-as-symbol-address
#You should create a package for it to work. Not time for now to learn to do that.

#Compare change in distribution
# for(s in 1:length(spatialKnowledge)){
print(s)

#Parallelising
cores=detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)
registerDoSEQ() #if you opt not to run parallelising -> C++ can't be exported with parallelisation apparently. I tried 
environmentPath <- getwd()

foreach(r=1:numberRepetitions, .packages=c('Rcpp'), .inorder = TRUE) %dopar% { #for(r in 1:30){
  #print(r)
  set.seed(r) #Unique random seed for reproducibility
  dateStart <- runif(numberTrees, 0, 365)
  locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
  runSimulation(
    cycleLimitNumber = cycleLimitNumber,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = samplingMapTime_v,
    nameInit = "Output/Main/Main",
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
    spatialKnowledgeRate = spatialKnowledge[s],
    temporalKnowledgeRate = temporalKnowledge[s],
    speed = speed,
    DispersalProbability = DispersalProbability, 
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = TRUE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
}

#Run measure of efficiency on final distribution without dispersion
stopCluster(cl)
gc()

# Run analysis: Efficiency test ---------------------------------------------

#Parallelising
cores=detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)
registerDoSEQ() #if you opt not to run parallelising -> C++ can't be exported with parallelisation apparently. I tried 
environmentPath <- getwd()

foreach(r=1:numberRepetitions, .packages=c('Rcpp'), .inorder = TRUE) %dopar% { 
  ##Compare agent efficiency
  path <- paste0(
    "Output/Main/Main",
    "_p15.811388",
    "_s",
    format(spatialKnowledge[s], nsmall=2),
    "0000_t",
    format(temporalKnowledge[s], nsmall=2),
    "0000_r",
    r
  )
  locTree <- matrix(unlist(matrix(read_table2(paste0(path,"_Map_",cycleLimitNumber*cycleLength,".txt")))), nrow = numberTrees, ncol = 4)
  set.seed(r) #Unique random seed for reproducibility
  #Null agent
  runSimulation(
    cycleLimitNumber = cycleLimitNumber/5,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime ,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0),
    nameInit = paste0("Output/Main/TestEfficiencyNULL",s),
    mapSize = mapSize,
    quadratSize = quadratSize,
    numberTrees = numberTrees,
    treeLocInit_m = locTree[,1:2],
    fruitingTimesInit_m = locTree[,3:4],
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
    spatialKnowledgeRate = 0,
    temporalKnowledgeRate = 0,
    speed = speed,
    DispersalProbability = 0, 
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = TRUE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  set.seed(r) #Unique random seed for reproducibility
  #Intermediate agent
  runSimulation(
    cycleLimitNumber = cycleLimitNumber/5,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime ,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0),
    nameInit = paste0("Output/Main/TestEfficiencyINTERMEDIATE",s),
    mapSize = mapSize,
    quadratSize = quadratSize,
    numberTrees = numberTrees,
    treeLocInit_m = locTree[,1:2],
    fruitingTimesInit_m = locTree[,3:4],
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
    speed = speed, 
    DispersalProbability = 0, 
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = TRUE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  set.seed(r) #Unique random seed for reproducibility
  #Omniscient/prescient agent
  runSimulation(
    cycleLimitNumber = cycleLimitNumber/5,
    repetitionNumber = r,
    timeDelayForDispersal = timeDelayForDispersal,
    torporTime = torporTime ,
    saveTreeMap = saveTreeMap,
    samplingMapTime_v = c(0),
    nameInit = paste0("Output/Main/TestEfficiencyOMNISCIENT",s),
    mapSize = mapSize,
    quadratSize = quadratSize,
    numberTrees = numberTrees,
    treeLocInit_m = locTree[,1:2],
    fruitingTimesInit_m = locTree[,3:4],
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
    DispersalProbability = 0, 
    useProvidedMap = TRUE,
    moveOnlyToFruitingTrees = TRUE,
    moveOnlyToTarget = FALSE,
    intensityCompetitionForSpace = intensityCompetitionForSpace
  )
  
  gc()
  memory.size()
}
# Stop cluster
stopCluster(cl)
gc()
memory.size()
# }  

