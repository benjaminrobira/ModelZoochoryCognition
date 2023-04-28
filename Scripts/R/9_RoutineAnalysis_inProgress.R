###~~~~~~~~~~~~~~
## Analysing routine output
###~~~~~~~~~~~~~~

# Setting up environments -------------------------------------------------

rm(list = ls())

## Parameters -------------------------------------------------------------

#Load parameters
load("Scripts/R/Parameterisation.RData")

## Packages ----------------------------------------------------------------

library(readr)
library(progress)
library(dplyr)
library(parallel)
library(doParallel)

## Functions ---------------------------------------------------------------

source("Scripts/R/RiotteLambert_2017_routine_f/entropy_mdf.R")

# Analysis ----------------------------------------------------------------

## Load files --------------------------------------------------------------

#What are the files containing the series of trees visited?
listFiles_v <- list.files("Output/SensitivityMovingRule")
whichRoutine_v <- grep("Routine", listFiles_v)
filesRoutine_v <- listFiles_v[whichRoutine_v]

listFiles_v <- list.files("Output/SensitivityMovingRule")
whichMapFinal_v <- grep("Map_7300", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]

outputRoutine <- mclapply(1:length(filesRoutine_v), mc.cores = 5, function(whichFile){#test for 5, replace with length(filesRoutine_v)
  seqVisits <- read_csv(paste0("Output/SensitivityMovingRule/",  filesRoutine_v[whichFile]), 
                        col_names = FALSE)
  seqVisits <- separate(seqVisits, col = "X1", into = c("seq", "time"), sep = " timer ")
  seqVisits_1 <- sapply(seqVisits[,1], function(x){substr(x, 1, 1)})
  whichRandom_v <- seqVisits[,1] == "r-1000"
  seqVisits_2 <- as.numeric(sapply(seqVisits[,1], function(x){gsub("[^0-9]", "", x)}))
  #Reput if -1000
  seqVisits_2[whichRandom_v] <- "-1000"
  seqVisits <- cbind(seqVisits_1, seqVisits_2) %>% as.data.frame
  seqVisits[,2] <- as.numeric(seqVisits[,2])
  colnames(seqVisits) <- c("moveType", "targetID")
  seqVisits <- seqVisits[seqVisits$targetID != -1000,]
  
  if(grepl("MovingNot", filesRoutine_v[whichFile])){
    whatRule <- "MovingAllTrees"
  }else{
    whatRule <- "MovingFruitTrees"
  }
  
  outputMapFinal <- read_table2(paste0("Output/SensitivityMovingRule/",  filesMapFinal_v[whichFile])) %>% as.data.frame() 
  
  return(list(whatRule, seqVisits, outputMapFinal))
})

length(outputRoutine)

## Routine index estimation ------------------------------------------------

routineIndex_v <- mclapply(1:length(outputRoutine), mc.cores = 5, function(whichApply){
  entropyValue <- entropy_O1(outputRoutine[[whichApply]][[2]][,2])
  return(1 - entropyValue)
})
routineIndex_v

outputRoutine[[whichApply]][[2]][1:5,]


seq <- as.numeric(outputRoutine[[whichApply]][[2]][,2]) + 1
seq <- seq[seq != lag(seq)]
coordsTravels = outputRoutine[[whichApply]][[3]][seq,]


plot(outputRoutine[[whichApply]][[3]][,1], outputRoutine[[whichApply]][[3]][,2], pch = 19)
lines(outputRoutine[[whichApply]][[3]][seq,][1:100,1], outputRoutine[[whichApply]][[3]][seq,][1:100,2], col = "red")
points(outputRoutine[[whichApply]][[3]][seq,][1:100,1], outputRoutine[[whichApply]][[3]][seq,][1:100,2], pch = 19, col = "red")





data <- seqVisits


outputMapFinal

library(ggplot2)
library(scales)
ggplot(outputMapFinal, aes(x = x, y = y)) +
  geom_point(pch = 21, colour = "black", fill = "white", size = 2.5) +
  geom_point(data = outputMapFinal[seqVisits$targetID, ], aes(x = x, y = y), pch = 19, colour = "darkgoldenrod") +
  geom_path(data = outputMapFinal[seqVisits$targetID[1:100], ], aes(x = x, y = y), pch = 19, colour = "darkgoldenrod") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4)) + 
  scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))


Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")
dateStart <- runif(numberTrees, 0, 365)
locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))

runSimulation(
  cycleLimitNumber = 1,
  repetitionNumber = 1,
  timeDelayForDispersal = 0.5,
  torporTime = 1,
  saveTreeMap = TRUE,
  samplingMapTime_v = seq(from=0, to=100*365, length.out = 5),
  nameInit = "Output/TestForAssessingModel/TEST",
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
  spatialKnowledgeRate = 0,
  temporalKnowledgeRate = 0,
  speed = 1000,
  DispersalProbability = 0, 
  useProvidedMap = TRUE,
  intensityCompetitionForSpace = 0.1
)


seqVisits <-  read_csv("Output/TestForAssessingModel/TEST_p15.000000_s0.000000_t0.000000_r1_Routine.txt", 
                       col_names = FALSE)

#Transform the recorded series into a series of only the visited tree ID (no timer, no info if true target or targeted en route)
seqVisits <- separate(seqVisits, col = "X1", into = c("seq", "time"), sep = " timer ")
seqVisits_1 <- sapply(seqVisits[,1], function(x){substr(x, 1, 1)})
whichRandom_v <- seqVisits[,1] == "r-1000"
seqVisits_2 <- as.numeric(sapply(seqVisits[,1], function(x){gsub("[^0-9]", "", x)}))
#Reput if -1000
seqVisits_2[whichRandom_v] <- "-1000"
  
seqVisits <- cbind(seqVisits_1, seqVisits_2) %>% as.data.frame
seqVisits[,2] <- as.numeric(seqVisits[,2])
colnames(seqVisits) <- c("moveType", "targetID")
seqVisits <- seqVisits[seqVisits$targetID != -1000,]




plot(locTree[,1], locTree[,2], pch = 19)
points(locTree[as.numeric(seqVisits[,2][1:6]) + 1,1], locTree[as.numeric(seqVisits[,2][1:6]) + 1,2], pch = 19, cex = 0.5, col = "red")
lines(locTree[as.numeric(seqVisits[,2][1:6]) + 1,1], locTree[as.numeric(seqVisits[,2][1:6]) + 1,2], col = "red")

isSeen <- visitedTrees(
  previousCoordinatesAgent = locTree[as.numeric(seqVisits[,2][3]),],
  currentCoordinatesAgent = locTree[as.numeric(seqVisits[,2][4]),],
  treeLoc = locTree,
  perceptualRange = 15
)

orderedSequence <- visitedTreesSequence(
  locTree,
  isSeen,
  locTree[as.numeric(seqVisits[,2][3]),]
)
orderedSequence

points(locTree[as.numeric(seqVisits[,2][1]),1], locTree[as.numeric(seqVisits[,2][1]),2], pch = 19, col = "green")
points(locTree[as.numeric(seqVisits[,2][2]),1], locTree[as.numeric(seqVisits[,2][2]),2], pch = 19, col = "green")

points(locTree[isSeen == 1,1], locTree[isSeen == 1,2], pch = 19, cex = 1.5, col = "blue")
points(locTree[orderedSequence + 1,1], locTree[orderedSequence + 1,2], pch = 19, cex = 0.5, col = "orange")
