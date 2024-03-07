##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the different Rcpp functions to see if it is working well
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script tests the different Rcpp function used for running the simulation in order to see whether those were doing what intended.

# Setting up environment --------------------------------------------------

## Parameters --------------------------------------------------------------

source("Scripts/R/0_Parameters.R")

## Libraries ---------------------------------------------------------------

library(Rcpp)
library(readr)

## Functions ---------------------------------------------------------------

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")

# Spatial distribution tree -----------------------------------------------

#Homogeneous
matrixTree <- distributionTree(
  numberTrees = 100,
  lowerBorder = 0,
  upperBorder = 100
)
plot(matrixTree[,1], matrixTree[,2])

#Heterogeneous
matrixTree <- distributionTree(
  numberTrees = 100,
  lowerBorder = 0,
  upperBorder = 100,
  homogeneousDistribution = FALSE,
  treeClusterNumber = 9,
  treeClusterSpread = 5
)
plot(matrixTree[,1], matrixTree[,2])

# Temporal distribution tree ----------------------------------------------

fruitingDates <- dateTree(
  numberTrees = 100,
  cycleLength = 365,
  fruitingLength = 30
)
fruitingDates

# Move agent --------------------------------------------------------------

#Obvious case
agentCoord = c(0,0)
treeLoc = matrix(c(1,1,2,2,3,3), ncol = 2, byrow = TRUE)
treeFood = c(1,1,1)

a <- moveAgentKnowledgeBased(
  currentCoordinatesAgent = agentCoord,
  treeLoc = treeLoc,
  treeFood = treeFood,
  whatRule = "closest",
  IDcurrentTree = NA
)
a

#Equality: select closest?
treeFood = c(sqrt(2)/sqrt(18),0,1)
b <- moveAgentKnowledgeBased(
  currentCoordinatesAgent = agentCoord,
  treeLoc = treeLoc,
  treeFood = treeFood,
  whatRule = "closest",
  IDcurrentTree = NA
)
b

treeFood = c(sqrt(2)/sqrt(18),0,1)
#Equality: select farthest?
c <- moveAgentKnowledgeBased(
  currentCoordinatesAgent = agentCoord,
  treeLoc = treeLoc,
  treeFood = treeFood,
  whatRule = "farthest",
  IDcurrentTree = NA
)
c

#Agent at one resource
d <- moveAgentKnowledgeBased(
  currentCoordinatesAgent = c,
  treeLoc = treeLoc,
  treeFood = treeFood,
  whatRule = "farthest",
  IDcurrentTree = NA
)
d

# Food update -------------------------------------------------------------
food <- foodTree(
  currentTime = 30,
  startTime = fruitingDates[,1],
  endTime = fruitingDates[,2],
  lengthOfCycle = 365,
  foodConsumed = rep(0, times=nrow(fruitingDates)),
  maximumFoodToYield = rep(1, times=nrow(fruitingDates)),
  linear = TRUE
)
food

food <- foodTree(
  currentTime = 30,
  startTime = fruitingDates[,1],
  endTime = fruitingDates[,2],
  lengthOfCycle = 365,
  foodConsumed = rep(0, times=nrow(fruitingDates)),
  maximumFoodToYield = rep(1, times=nrow(fruitingDates)),
  linear = FALSE
)
food

# Changing coordinate system ----------------------------------------------

newCoords <- changeReferentialCoordinatesAlongTravel(
  departureLocation = c(0,0),
  arrivalLocation = c(3,3),
  treeLoc = cbind(c(1,2,3), c(1,2,3))
)
newCoords

# Is tree on the way to target or at start/end perceptual range--------------------------------------------

isSeen <- visitedTrees(
  previousCoordinatesAgent = c(0,0),
  currentCoordinatesAgent = c(3,3),
  treeLoc = cbind(c(1,2,3), c(1,2,3)),
  perceptualRange = 0.5
)
isSeen

# Which trees were visited ------------------------------------------------

visited_v <- rep(0, times = nrow(matrixTree))
visited_v[sample(1:length(visited_v), 5)] <- 1#Random assignment, not linked to previous visit correctly here, just to see if ordering correct

whichVisited <- which(visited_v == 1)
distWhichVisited <- apply(matrixTree[whichVisited,], 1, function(x){
  (x[1] - matrixTree[1,1])**2 + (x[2] - matrixTree[1,2])**2
})
whichVisited
distWhichVisited

orderedSequence <- visitedTreesSequence(
  matrixTree,
  visited_v,
  matrixTree[1,]
)
orderedSequence

# New coordinates if dispersed --------------------------------------------
# Rcpp::sourceCpp("Script/Rcpp/FunctionsRcpp.cpp")
#
# newCoordinatesDispersal <- newCoordinatesAfterDispersal(
#   treeLoc = cbind(runif(1000, 0, 3), runif(1000, 0, 3)),
#   previousCoordinatesAgent = c(0,0),
#   currentCoordinatesAgent = c(3,3),
#   perceptualRange = 0,
#   mapSize = 1000
# )
# newCoordinatesDispersal

newCoordinatesDispersal <- newCoordinatesAfterDispersal(
  previousCoordinatesAgent = c(0,0),
  currentCoordinatesAgent = c(3,3),
  rangeDispersal = 0,
  mapSize = 1000
)
newCoordinatesDispersal

newCoordinatesDispersal <- lapply(1:1000, function(x){newCoordinatesAfterDispersal(
  previousCoordinatesAgent = c(0,0),
  currentCoordinatesAgent = c(1000,1000),
  rangeDispersal = 30,
  mapSize = 1000
)})

newCoordinatesDispersal <- lapply(1:1000, function(x){newCoordinatesAfterDispersal(
  previousCoordinatesAgent = c(0,0),
  currentCoordinatesAgent = c(1000,1000),
  rangeDispersal = 200,
  mapSize = 1000
)})
newCoordinatesDispersal <- do.call("rbind", newCoordinatesDispersal)

plot(newCoordinatesDispersal[,1], newCoordinatesDispersal[,2])
abline(a=200*sin(45*pi/180)*2, b=1)
abline(a=-200*sin(45*pi/180)*2, b=1)

# Random selection for spatio-temporal knowledge --------------------------

whichTreeKnowledge <- whichDistanceAndFoodKnown(
  numberTrees = 100,
  spatialKnowledgeRate = 0,
  temporalKnowledgeRate = 0
)
whichTreeKnowledge

whichTreeKnowledge <- whichDistanceAndFoodKnown(
  numberTrees = 100,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1
)
whichTreeKnowledge

whichTreeKnowledge <- whichDistanceAndFoodKnown(
  numberTrees = 100,
  spatialKnowledgeRate = 0.5,
  temporalKnowledgeRate = 0.7
)
whichTreeKnowledge

whichTreeKnowledge <- whichDistanceAndFoodKnown(
  numberTrees = 100,
  spatialKnowledgeRate = 0.5,
  temporalKnowledgeRate = 0.4
)
whichTreeKnowledge
whichTreeKnowledge[whichTreeKnowledge[,1]==0&whichTreeKnowledge[,2]==1,] #empty, ok
length(whichTreeKnowledge[whichTreeKnowledge[,1]==1,1])/nrow(whichTreeKnowledge) #ok
length(whichTreeKnowledge[whichTreeKnowledge[,2]==1,2])/nrow(whichTreeKnowledge) #ok

# Checking which trees are within perceptual range ------------------------

whichPerceivedTrees <- whichPerceived(
  treeLoc = cbind(c(1,2,3), c(1,2,3)),
  currentLoc = c(0,0),
  perceptualRange = 0.5
)
whichPerceivedTrees

whichPerceivedTrees <- whichPerceived(
  treeLoc = cbind(c(1,2,3), c(1,2,3)),
  currentLoc = c(0,0),
  perceptualRange = sqrt(2)
)
whichPerceivedTrees

whichPerceivedTrees <- whichPerceived(
  treeLoc = cbind(c(1,2,3), c(1,2,3)),
  currentLoc = c(0,0),
  perceptualRange = 1.5
)
whichPerceivedTrees

# Which trees are fruiting ------------------------------------------------

whichFruitingTrees <- whichInFruit(
  currentTime = 15,
  startTimeVector = c(0,50,10),
  endTimeVector = c(20, 40, 50),
  lengthOfCycle = 365
)
whichFruitingTrees

whichFruitingTrees <- whichInFruit(
  currentTime = 15,
  startTimeVector = c(0,50,10),
  endTimeVector = c(20, 60, 50),
  lengthOfCycle = 365
)
whichFruitingTrees

# Random movement ---------------------------------------------------------

newCoords = matrix(NA, ncol=2, nrow=2000)

for(i in 1:2000){
  newCoords[i,] <- moveAgentRandom(
    currentLocation = c(0,0),
    mapSize = 1000,
    exponentialRate = 0.01
  )
}
newCoords
hist(sqrt(newCoords[,1]**2 + newCoords[,2]**2))
newCoords[newCoords[,1]<0|newCoords[,1]>1000,] #empty ok
newCoords[newCoords[,2]<0|newCoords[,2]>1000,] #empty ok

# Recreate location/food available matrix based on knowledge --------
treeLocTest <- cbind(runif(100, 0, 1000), runif(100, 0, 1000))
foodQuantityTest <- runif(100, 0, 1)
whichVisited <- floor(runif(100, 0, 365))

toCheck <- which(
  whichTreeKnowledge[,1]==1 &
  whichTreeKnowledge[,2]==0 &
  whichVisited >= 365 - 15
)
toCheck

newTreeMatrix <- filterKnowledge(
  treeLoc = treeLocTest,
  foodTree = foodQuantityTest,
  whichDistanceAndFoodKnown = whichTreeKnowledge,
  timeLastVisit_v = whichVisited,
  365,
  noReturnTime = 15,
  whatValueUnknownTemporal = 0.001,
  mapSize = 1000
)
newTreeMatrix
newTreeMatrix[toCheck,]

# Update last time visit and food consumed --------------------------------

matrixUpdateVisits <- updateVisitTimesAndConsumptionCycle(
                                    currentTime = 70,
                                    startTimeVector = c(0, 20, 10, 5),
                                    endTimeVector = c(0, 20, 10, 5) + 15,
                                    lengthOfCycle = 30,
                                    foodConsumedVector = c(0, 0.2, 0, 0),
                                    visitedTrees_v = c(1, 0, 1, 0),
                                    lastTimeVisitTrees_v = c(0,0,0,0),
                                    lastTimeVisitTreesWithFruit_v = c(0,0,0,0),
                                    foodQuantityAtTree_v = c(0.5, 0, 0, 0.25)
)
matrixUpdateVisits


# Nearest neighbour -------------------------------------------------------

NN <- nearestNeighbourDistance(
  locMatrix = cbind(runif(300, 0, 1000), runif(300, 0, 1000)),
  maxDistancePossible = 1000 + 1
)
NN


# Lloyd -------------------------------------------------------------------

distrib = cbind(runif(300, 0, 1000), runif(300, 0, 1000))
lloyd <- lloydIndex(
  treeLoc = distrib,
  mapSize = 1000,
  quadratSize = 50
)
lloyd

#Test for which mean closest from 1 and sd reduced
patchSizeToTest <- c(10, 25, 50, 100, 200)
resultsMeanRandom = rep(0, times=length(patchSizeToTest))
resultsSdRandom = rep(0, times=length(patchSizeToTest))

resultsMeanClustered = rep(0, times=length(patchSizeToTest))
resultsSdClustered = rep(0, times=length(patchSizeToTest))
numberPatchinessRepetition = 500

for(i in 1:length(patchSizeToTest)){

  transitoryoutput_v_Random = rep(0, times=numberPatchinessRepetition)
  transitoryoutput_v_Clustered = rep(0, times=numberPatchinessRepetition)

  for(j in 1:numberPatchinessRepetition){
  #Random
  distrib = cbind(runif(1000, 0, 1000), runif(1000, 0, 1000))
  lloyd <- lloydIndex(
    treeLoc = distrib,
    mapSize = 1000,
    quadratSize = patchSizeToTest[i]
  )
  transitoryoutput_v_Random[j] = lloyd[2]

  #Clustered
  distrib = distributionTree(
    numberTrees = 1000,
    lowerBorder = 0,
    upperBorder = 1000,
    homogeneousDistribution = FALSE,
    treeClusterNumber = 1000/10,
    treeClusterSpread = 25
  )
  lloyd <- lloydIndex(
    treeLoc = distrib,
    mapSize = 1000,
    quadratSize = patchSizeToTest[i]
  )
  transitoryoutput_v_Clustered[j] = lloyd[2]
  }
  resultsMeanRandom[i] = mean(transitoryoutput_v_Random)
  resultsSdRandom[i] = sd(transitoryoutput_v_Random)

  resultsMeanClustered[i] = mean(transitoryoutput_v_Clustered)
  resultsSdClustered[i] = sd(transitoryoutput_v_Clustered)
}

resultsMeanRandom
resultsSdRandom

resultsMeanClustered
resultsSdClustered

#50 is correct then

# Run simulation test ----------------------------------------------------------

#Check if output are working well
Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")
set.seed(10)
dateStart <- runif(numberTrees, 0, 30)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))

runSimulation(
  cycleLimitNumber = 1,
  repetitionNumber = 1,
  timeDelayForDispersal = 0.5,
  torporTime = 1,
  saveTreeMap = TRUE,
  samplingMapTime_v = seq(from = 0, to = 2*30, length.out = 3),
  nameInit = "Output/TestForAssessingModel/WithDispersal_checkingOutput_TEST",
  mapSize = 100,
  quadratSize = 5,
  numberTrees = 1000,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + 3),
  homogeneousDistribution = TRUE,
  treeClusterNumber = 0,
  treeClusterSpread = 0,
  maximumFoodToYield_v = rep(1, times = 10),
  cycleLength = 30,
  fruitingLength = 3,
  noReturnTime = 2,
  whatValueUnknownTemporal = 0.000,
  whatRule = "closest",
  exponentialRate = 0.01,
  perceptualRange = 15,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1,
  speed = 1000,
  DispersalProbability = 0.01, 
  useProvidedMap = FALSE
)

map <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_checkingOutput_TEST_p15.000000_s1.000000_t1.000000_r1_Map_0.txt"))
map
efficiency <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_checkingOutput_TEST_p15.000000_s1.000000_t1.000000_r1_Efficiency.txt"))
efficiency
continuous <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_checkingOutput_TEST_p15.000000_s1.000000_t1.000000_r1_EfficiencyContinuous.txt"))
continuous
routine <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_checkingOutput_TEST_p15.000000_s1.000000_t1.000000_r1_Routine.txt"))
routine #NOTE: doublons are due to random movements, checked

#Check if this works and if shrinkage resources
dateStart <- runif(numberTrees, 0, 365)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))

runSimulation(
        cycleLimitNumber = 100,
        repetitionNumber = 1,
        timeDelayForDispersal = 0.5,
        torporTime = 1,
        saveTreeMap = TRUE,
        samplingMapTime_v = seq(from=0, to=100*365, length.out = 5),
        nameInit = "Output/TestForAssessingModel/WithDispersalTEST",
        mapSize = 1000,
        quadratSize = 50,
        numberTrees = 1000,
        treeLocInit_m = locTree,
        fruitingTimesInit_m = cbind(dateStart, dateStart + 30),
        homogeneousDistribution = TRUE,
        treeClusterNumber = 0,
        treeClusterSpread = 0,
        maximumFoodToYield_v = rep(1, times=1000),
        cycleLength = 365,
        fruitingLength = 30,
        noReturnTime = 2,
        whatValueUnknownTemporal = 0.000,
        whatRule = "closest",
        exponentialRate = 0.01,
        perceptualRange = 15,
        spatialKnowledgeRate = 1,
        temporalKnowledgeRate = 1,
        speed = 1000,
        DispersalProbability = 0.01, 
        useProvidedMap = FALSE,
        intensityCompetitionForSpace = 0.1
      )

#Initial map
map <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersalTEST_p15.000000_s1.000000_t1.000000_r1_Map_0.txt"))
plot(map[,1], map[,2], pch=19, col="black")
library(dplyr)
alignmentPoints_f(map %>% dplyr::select(x,y))#Note: need to be imported from script A.
patchiness_f(
  map %>% dplyr::select(x,y),
  0,
  0,
  mapSize,
  mapSize,
  50
)
#Final map
map2 <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersalTEST_p15.000000_s1.000000_t1.000000_r1_Map_3650.txt"))
#plot(map2[,1], map2[,2], pch=19, col="black")#if plotted without previous
points(map2[,1], map2[,2], pch=1, col="red")#if added
alignmentPoints_f(map2 %>% dplyr::select(x,y))
patchiness_f(
  map2 %>% dplyr::select(x,y),
  0,
  0,
  mapSize,
  mapSize,
  50
)

#Check if distance between tree is correct
sapply(1:nrow(map2), function(x){
  min(
    sqrt(
      (map2[-x,1] - map2[x,1])**2+
        (map2[-x,2] - map2[x,2])**2
    )
  )
})


#Check what is happening with reduced resource trees
dateStart <- runif(numberTrees, 0, 30)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 10), runif(numberTrees, 0, 10)))

runSimulation(
  cycleLimitNumber = 60,
  repetitionNumber = 1,
  timeDelayForDispersal = 0.5,
  torporTime = 1,
  saveTreeMap = TRUE,
  samplingMapTime_v = seq(from = 0, to = 60*30, length.out = 5),
  nameInit = "Output/TestForAssessingModel/WithDispersal_limitedTrees_TEST",
  mapSize = 100,
  quadratSize = 5,
  numberTrees = 10,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + 3),
  homogeneousDistribution = TRUE,
  treeClusterNumber = 0,
  treeClusterSpread = 0,
  maximumFoodToYield_v = rep(1, times = 10),
  cycleLength = 30,
  fruitingLength = 3,
  noReturnTime = 1,
  whatValueUnknownTemporal = 0.000,
  whatRule = "closest",
  exponentialRate = 0.01,
  perceptualRange = 15,
  spatialKnowledgeRate = 1,
  temporalKnowledgeRate = 1,
  speed = 1000,
  DispersalProbability = 1, 
  useProvidedMap = FALSE
)

#Initial map
map <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_limitedTrees_TEST_p15.000000_s1.000000_t1.000000_r1_Map_0.txt"))
plot(map[,1], map[,2], pch=19, col="black")
text(map[,1], map[,2]-2, label = 1:nrow(map))
library(dplyr)
alignmentPoints_f(map %>%  dplyr::select(x,y))

#Final map
map2 <- as.data.frame(read_table2("Output/TestForAssessingModel/WithDispersal_limitedTrees_TEST_p15.000000_s1.000000_t1.000000_r1_Map_1800.txt"))
#plot(map2[,1], map2[,2], pch=19, col="black")#if plotted without previous
points(map2[,1], map2[,2], pch=1, col="red")#if added
text(map2[,1], map2[,2]-2, label = 1:nrow(map))
alignmentPoints_f(map2 %>%  dplyr::select(x,y))

#Have a visual representation of the evolution with time of the tree distribution

# runSimulation(
#   cycleLimitNumber = 2,
#   repetitionNumber = 1,
#   timeDelayForDispersal = 0.5,
#   torporTime = 1,
#   saveTreeMap = TRUE,
#   samplingMapTime_v = seq(from=60, to=2*60, length.out = 60),
#   nameInit = "Output/TestSameFruitingDate/TestContinuousMapping/Test",
#   mapSize = 1000,
#   quadratSize = 50,
#   numberTrees = 1000,
#   treeLocInit_m = locTree,
#   fruitingTimesInit_m = cbind(dateStart, dateStart + 30),
#   homogeneousDistribution = TRUE,
#   treeClusterNumber = 0,
#   treeClusterSpread = 0,
#   maximumFoodToYield_v = rep(1, times=1000),
#   cycleLength = 60,
#   fruitingLength = 30,
#   noReturnTime = 5,
#   whatValueUnknownTemporal = 0.000,
#   whatRule = "closest",
#   exponentialRate = 0.01,
#   perceptualRange = 10,
#   spatialKnowledgeRate = 0.5,
#   temporalKnowledgeRate = 0.5,
#   speed = 1000,
#   DispersalProbability = 10000, 
#   useProvidedMap = FALSE,
#   linear = FALSE
# )

runSimulation(
  cycleLimitNumber = 2,
  repetitionNumber = 1,
  timeDelayForDispersal = 0.5,
  torporTime = 1,
  saveTreeMap = TRUE,
  samplingMapTime_v = seq(from=60, to=2*60, length.out = 60),
  nameInit = "Output/TestSameFruitingDate/TestContinuousMapping/Test",
  mapSize = 1000,
  quadratSize = 50,
  numberTrees = 1000,
  treeLocInit_m = locTree,
  fruitingTimesInit_m = cbind(dateStart, dateStart + 30),
  homogeneousDistribution = TRUE,
  treeClusterNumber = 0,
  treeClusterSpread = 0,
  maximumFoodToYield_v = rep(1, times=1000),
  cycleLength = 60,
  fruitingLength = 30,
  noReturnTime = 5,
  whatValueUnknownTemporal = 0.000,
  whatRule = "closest",
  exponentialRate = 0.01,
  perceptualRange = 10,
  spatialKnowledgeRate = 0.5,
  temporalKnowledgeRate = 0.5,
  speed = 1000,
  DispersalProbability = 10000, 
  useProvidedMap = FALSE,
  linear = FALSE
)

repetitionNumber = 1
path <- paste0(
  "Output/TestSameFruitingDate/TestContinuousMapping/Test",
  "_p10.000000",
  "_s",
  0.50,
  "00000_t",
  0.50,
  "00000_r",
  repetitionNumber
)

for(i in 60:120){
  map <- as.data.frame(read_table2(paste0(path,"_Map_", i, ".txt")))
  plot(map[,1], map[,2], pch=19, col="black")
  if(i > 60){
    map <- as.data.frame(read_table2(paste0(path,"_Map_", 60, ".txt")))
    plot(map[,1], map[,2], pch=19, col="red")
    map <- as.data.frame(read_table2(paste0(path,"_Map_", i - 1, ".txt")))
    points(map[,1], map[,2], pch=19, col="blue")
  }
  map <- as.data.frame(read_table2(paste0(path,"_Map_", i, ".txt")))
  points(map[,1], map[,2], pch=19, col="black")
  
  Sys.sleep(1)
}
  
map <- as.data.frame(read_table2(paste0(path,"_Map_", 60, ".txt")))
plot(map[,1], map[,2], pch=19, col="green")

map <- as.data.frame(read_table2(paste0(path,"_Map_", 85, ".txt")))
points(map[,1], map[,2], pch=19, col="red")

map <- as.data.frame(read_table2(paste0(path,"_Map_", 120, ".txt")))
points(map[,1], map[,2], pch=1, col="blue")

# Run simulation test: new replacement memory due to seed dispersal -------

dateStart <- runif(numberTrees, 0, 365)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
runSimulationMemoryRevised(
  cycleLimitNumber = cycleLimitNumber,
  repetitionNumber = 1,
  timeDelayForDispersal = timeDelayForDispersal,
  torporTime = torporTime,
  saveTreeMap = saveTreeMap,
  samplingMapTime_v = samplingMapTime_v,
  nameInit = "Output/TestSimulationNewMemorySeedDispersal",
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
  spatialKnowledgeRate = 0.5,
  temporalKnowledgeRate = 0.5,
  speed = speed,
  DispersalProbability = DispersalProbability, 
  useProvidedMap = TRUE,
  moveOnlyToFruitingTrees = FALSE,
  moveOnlyToTarget = FALSE,
  intensityCompetitionForSpace = intensityCompetitionForSpace
)
