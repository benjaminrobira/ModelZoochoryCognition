## ~~~~~~~~~~~~~~~~
## List of parameters for running simulation
## ~~~~~~~~~~~~~~~~

# What does this script do?
# This script sets up all parameters needed for running the simulations.

#rm(list=ls())

#Parameters

numberRepetitions = 200#Repetition of same conditions for agent (but note environment - spatio-temporal distribution is resampled)

# Environment related -----------------------------------------------------

## General var -------------------------------------------------------------

mapSize = 1000
quadratSize = 50#To calculate Lloyd index within the simulation too
cycleLimitNumber = 100
cycleLength = 365

## Resource ----------------------------------------------------------------

numberTrees = 1000
#treeLocInit_m = locTree
#fruitingTimesInit_m = cbind(dateStart, dateStart + 30)
homogeneousDistribution = TRUE
treeClusterNumber = 0
treeClusterSpread = 0
maxFoodToYield = 1
maximumFoodToYield_v = rep(maxFoodToYield, times = numberTrees)
fruitingLength = 30
intensityCompetitionForSpace = 0.45

# Agent -------------------------------------------------------------------

## Abilities ---------------------------------------------------------------

torporTime = 1
perceptualRange = 1/(2*sqrt(numberTrees/mapSize/mapSize))
spatialKnowledge <- c(0, 0.25, 0.5, 0.75, 1)
temporalKnowledge <- c(0, 0.25, 0.5, 0.75, 1)
speed = 1000
noReturnTime = 2
whatValueUnknownTemporal = 0.000
exponentialRate = 0.01#Exponential for random movement (step length)
whatRule = "closest" #What to do if known locations have equal interest: choose closest or farthest

## Seed Dispersal ----------------------------------------------------------

timeDelayForDispersal = 0.5
DispersalProbability = 0.01

# Saving output options ---------------------------------------------------

saveTreeMap = TRUE
samplingMapTime_v = floor(seq(from = 0, to = cycleLimitNumber*cycleLength, length.out = 5))

# Create table for parameters ---------------------------------------------

parametersValue <- matrix(NA, ncol = 5, nrow = 30)

parametersValue[1,] <- c("Environment","Map size","Length of a side of the square environmental map", mapSize, "su")
parametersValue[2,] <- c("Environment","Quadrat size","Length of a side of a square quadrat to calculate Lloyd index of patchiness", quadratSize, "su")
parametersValue[3,] <- c("Environment","Period length","Length of a period before a given plant starts producing again", cycleLength, "tu")
parametersValue[4,] <- c("Environment","Number of seasons","Number of seasons (with seed dispersal plus without seed dispersal) before the simulation is ended", paste(cycleLimitNumber,5, sep = " + "), "")
parametersValue[5,] <- c("Environment","Number of plants","Number of plants hosted by the environment", numberTrees, "")
parametersValue[6,] <- c("Environment","Fruiting length","Time duration of the fruiting period of each plant", fruitingLength, "tu")
parametersValue[7,] <- c("Environment","Maximum food yielded at a plant","Food quantity that a plant might yield at best (peak of the triangular-shaped food distribution)", maxFoodToYield,"fu")
parametersValue[8,] <- c("Environment","Spacing intensity", "Relative length of a square map whose area would correspond to the area of exclusive spaces of all plants without overlapping", paste(c(0.05, 0.45, 0.85)*100, collapse = ", "), "-")
parametersValue[9,] <- c("Agent","Speed","Speed at which the forager moves", paste0("(", paste(c(speed/2, speed, speed*2), collapse = ","), ")"), "su/tu")
parametersValue[10,] <- c("Agent","Torpor time","Time duration for which the forager stops foraging in case no food is available in the environment", torporTime, "tu")
parametersValue[11,] <- c("Agent","Perceptual range","Distance at which the forager is aware of the environment", round(perceptualRange, digits = 2), "su")
parametersValue[12,] <- c("Agent","Knowledge rate","Proportion of plants of the environment for which the forager knows the location and prodution timing", paste0("(", paste(spatialKnowledge, collapse = ", "), ")"),"")
parametersValue[13,] <- c("Agent","No-return time","Time delay before a forager mentally decides to target a previously visited plant", noReturnTime, "tu")
parametersValue[14,] <- c("Agent","Dispersal time","Time duration during which seeds from a previously ingested fruit can be dispersed", timeDelayForDispersal, "tu")
parametersValue[15,] <- c("Agent","Probability of dispersal","Probability (per tu) that the seeds is actually dispersed", DispersalProbability/timeDelayForDispersal,"1/tu")
parametersValue[16,] <- c("Agent", 
                          as.character(expression(lambda["step length"])),
                          "Average step length for random movements used to parameterise the exponential distribution", exponentialRate,"su")
# parametersValue[16,] <- c("","","","")
# parametersValue[17,] <- c("","","","")
# parametersValue[18,] <- c("","","","")#c("","","Random distribution to obtain the step length for random movement","")
# parametersValue[19,] <- c("","","","")#c("","","Random distribution to obtain the start of fruiting of each plant","")
# parametersValue[20,] <- c("","","","")#c("","",as.character(expression(italic(x)~and~italic(y)~"coordinates of a plant"~italic(i))),"")

rownames(parametersValue) <- NULL
colnames(parametersValue) <- c("Modelling entity", "Parameter", "Definition", "Value", "Unit")

parametersValue <- parametersValue[!is.na(parametersValue[,1]),]

# Saving image ------------------------------------------------------------

#save.image("Scripts/R/Parameterisation.RData")