##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the effect of sensory range
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script runs the simulations modifying the agents' speed.

# Setting up environment --------------------------------------------------

## Parameters --------------------------------------------------------------

rm(list = ls())
source("Scripts/R/0_Parameters.R")
numberRepetitionsSensitivity = numberRepetitions

## Libraries ---------------------------------------------------------------

library(Rcpp)
library(readr)

## Functions ---------------------------------------------------------------

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")

# Running the simulations -----------

# for(r in 1:numberRepetitionsSensitivity){
#   print(r)
#   set.seed(r)
#   dateStart <- runif(numberTrees, 0, cycleLength)
#   locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
# 
#   runSimulation(
#     cycleLimitNumber = cycleLimitNumber,
#     repetitionNumber = r,
#     timeDelayForDispersal = timeDelayForDispersal,
#     torporTime = torporTime,
#     saveTreeMap = saveTreeMap,
#     samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
#     nameInit = "Output/SensitivitySensoryRange/SensoryRangeIntermediate",
#     mapSize = mapSize,
#     quadratSize = quadratSize,
#     numberTrees = numberTrees,
#     treeLocInit_m = locTree,
#     fruitingTimesInit_m = cbind(dateStart, dateStart + fruitingLength),
#     homogeneousDistribution = homogeneousDistribution,
#     treeClusterNumber = treeClusterNumber,
#     treeClusterSpread = treeClusterSpread,
#     maximumFoodToYield_v = maximumFoodToYield_v,
#     cycleLength = cycleLength,
#     fruitingLength = fruitingLength,
#     noReturnTime = noReturnTime,
#     whatValueUnknownTemporal = whatValueUnknownTemporal,
#     whatRule = whatRule,
#     exponentialRate = exponentialRate,
#     perceptualRange = perceptualRange*2,
#     spatialKnowledgeRate = 0,
#     temporalKnowledgeRate = 0,
#     speed = speed,
#     DispersalProbability = DispersalProbability,
#     useProvidedMap = TRUE,
#     moveOnlyToFruitingTrees = FALSE,
#     moveOnlyToTarget = FALSE,
#     intensityCompetitionForSpace = intensityCompetitionForSpace
#   )
# 
#   set.seed(r)# restart everything to be sure the seed use remains the same
#   dateStart <- runif(numberTrees, 0, cycleLength)
#   locTree <- as.matrix(cbind(runif(numberTrees, 0, 1000), runif(numberTrees, 0, 1000)))
#   
#   runSimulation(
#     cycleLimitNumber = cycleLimitNumber,
#     repetitionNumber = r,
#     timeDelayForDispersal = timeDelayForDispersal,
#     torporTime = torporTime,
#     saveTreeMap = saveTreeMap,
#     samplingMapTime_v = c(0, cycleLimitNumber * cycleLength),
#     nameInit = "Output/SensitivitySensoryRange/SensoryRangeHigh",
#     mapSize = mapSize,
#     quadratSize = quadratSize,
#     numberTrees = numberTrees,
#     treeLocInit_m = locTree,
#     fruitingTimesInit_m = cbind(dateStart, dateStart + fruitingLength),
#     homogeneousDistribution = homogeneousDistribution,
#     treeClusterNumber = treeClusterNumber,
#     treeClusterSpread = treeClusterSpread,
#     maximumFoodToYield_v = maximumFoodToYield_v,
#     cycleLength = cycleLength,
#     fruitingLength = fruitingLength,
#     noReturnTime = noReturnTime,
#     whatValueUnknownTemporal = whatValueUnknownTemporal,
#     whatRule = whatRule,
#     exponentialRate = exponentialRate,
#     perceptualRange = perceptualRange*4,
#     spatialKnowledgeRate = 0,
#     temporalKnowledgeRate = 0,
#     speed = speed,
#     DispersalProbability = DispersalProbability,
#     useProvidedMap = TRUE,
#     moveOnlyToFruitingTrees = FALSE,
#     moveOnlyToTarget = FALSE,
#     intensityCompetitionForSpace = intensityCompetitionForSpace
#   )
# }

# Extracting results and plots --------------------------------------------

### Load complementary data -------------------------------------------------

load("Renvironment/outputTestCognition.RData")

dataLowSensoryRangeNaive <- indicesMain_df %>% 
  filter(knowledgeRate == 0 & time == 36500) %>% 
  dplyr::select(-time) %>% 
  mutate(knowledgeRate = "Low") %>% 
  rename(SensoryRange = "knowledgeRate")


tableIndexSpatial <- read_delim("Renvironment/tableIndexSpatial.txt", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

#Patchiness
patchinessRouteValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Patchiness" &
                                                  tableIndexSpatial$Density == "High" &
                                                  tableIndexSpatial$Distribution == "Route"]
patchinessHeterogeneousValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Patchiness" &
                                                          tableIndexSpatial$Density == "High" &
                                                          tableIndexSpatial$Distribution == "Heterogeneous"]
patchinessHomogeneousValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Patchiness" &
                                                        tableIndexSpatial$Density == "High" &
                                                        tableIndexSpatial$Distribution == "Homogeneous"]

alignmentRouteValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Alignment" &
                                                 tableIndexSpatial$Density == "High" &
                                                 tableIndexSpatial$Distribution == "Route"]
alignmentHeterogeneousValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Alignment" &
                                                         tableIndexSpatial$Density == "High" &
                                                         tableIndexSpatial$Distribution == "Heterogeneous"]
alignmentHomogeneousValue <- tableIndexSpatial$Value[tableIndexSpatial$Index == "Alignment" &
                                                       tableIndexSpatial$Density == "High" &
                                                       tableIndexSpatial$Distribution == "Homogeneous"]


#Spat autocorr

load("Renvironment/TestSpatAutocorr.RData")

spAutoCorrIndex_dflong <- pivot_longer(
  spAutoCorrIndex_df,
  !corr,
  names_to = "Distribution", 
  values_to = "Moran"
) %>% 
  mutate(
    facet = ifelse(corr == 1, "atop", "bottom")
  )

spAutoCorrIndexMean_dflong <- spAutoCorrIndex_dflong %>% 
  group_by(corr, Distribution) %>% 
  summarise(
    Moran = mean(Moran)
  )

moranHomoLow <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHomogeneous" & spAutoCorrIndexMean_dflong$corr == 0]
moranHomoIntermediate <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHomogeneous" & spAutoCorrIndexMean_dflong$corr == 0.5]
moranHomoHigh <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHomogeneous" & spAutoCorrIndexMean_dflong$corr == 1]

moranHeteroLow <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHeterogeneous" & spAutoCorrIndexMean_dflong$corr == 0]
moranHeteroIntermediate <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHeterogeneous" & spAutoCorrIndexMean_dflong$corr == 0.5]
moranHeteroHigh <- spAutoCorrIndexMean_dflong$Moran[spAutoCorrIndexMean_dflong$Distribution == "spAutoCorrHeterogeneous" & spAutoCorrIndexMean_dflong$corr == 1]

## Packages ----------------------------------------------------------------

library(readr)
library(progress)
library(Rcpp)
library(adehabitatLT)
library(adehabitatHR)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)

## Functions ---------------------------------------------------------------

source("Scripts/R/RiotteLambert_2017_routine_f/entropy_mdf.R")

### Spatial autocorrelation (Moran, Geary) ----------------------------------

autocorrSp <-
  function(x,
           weight,
           index = c("Moran", "Geary"),
           circular = FALSE,
           na.rm = TRUE) {
    #Adapted from the ape package to have it circularized, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870747/, Moran.I function
    if (dim(weight)[1] != dim(weight)[2])
      stop("'weight' must be a square matrix")
    n <- length(x)
    if (dim(weight)[1] != n)
      stop("'weight' must have as many rows as observations in 'x'")
    ei <- -1 / (n - 1)
    nas <- is.na(x)
    if (any(nas)) {
      if (na.rm) {
        x <- x[!nas]
        n <- length(x)
        weight <- weight[!nas,!nas]
      }
      else {
        warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
        return(list(
          observed = NA,
          expected = ei,
          sd = NA,
          p.value = NA
        ))
      }
    }
    ROWSUM <- rowSums(weight)
    ROWSUM[ROWSUM == 0] <- 1
    weight <- weight / ROWSUM
    s <- sum(weight)
    y <- 0
    obs <- NA
    if (index == "Moran") {
      if (circular) {
        #Mean values
        cosmean = mean(cos(x))
        sinmean = mean(sin(x))
        valuemean = atan2(sinmean, cosmean)
        y <-
          sapply(x, function(value) {
            atan2(sin(value - valuemean), cos(value - valuemean))
          })
      } else{
        m <- mean(x)
        y <- x - m
      }
      v <- sum(y ^ 2)
      cv <- sum(weight * y %o% y)
      obs <- (n / s) * (cv / v)
    } else if (index == "Geary") {
      if (circular) {
        #Mean values
        cosmean = mean(cos(x))
        sinmean = mean(sin(x))
        valuemean = atan2(sinmean, cosmean)
        y <-
          sapply(x, function(value) {
            atan2(sin(value - valuemean), cos(value - valuemean))
          })
        z <-
          sapply(x, function(value) {
            atan2(sin(value - x), cos(value - x))
          })
      } else{
        m <- mean(x)
        y <- x - m
        z <- sapply(x, function(value) {
          (value - x) ** 2
        })
      }
      v <- sum(y ^ 2)
      cv <- sum(weight * z ** 2)
      obs <- (n - 1) / (2 * s) * (cv / v)
    }
    return(obs)
  }

### Reticulation analysis function; following description of index f --------
#https://repositories.lib.utexas.edu/bitstream/handle/2152/68127/ADDIS-DISSERTATION-2017.pdf?sequence=1

# install.packages("cccd")
library(cccd)
library(igraph)

# locations_df <- as.matrix(cbind(runif(500,0,500), runif(500,0,500)))
# distanceType = "Euclidean"

reticulationIndex_f <-
  function(locations_df, distanceType = "Euclidean") {
    #Compute the relative neighbouring graph
    relativeNeighbouringGraph <-
      rng(
        x = locations_df,
        r = 1,
        method = distanceType,
        usedeldir = TRUE,
        open = FALSE,
        algorithm = 'cover_tree'
      )
    
    #Compute the matrix of adjacency = which locations are linked by an edge
    adjacencyMatrix <- as_adjacency_matrix(relativeNeighbouringGraph)
    adjacencyMatrix <- as.matrix(adjacencyMatrix)
    
    #Get the distance matrix between locations
    distanceMatrix <-
      dist(locations_df,
           method = tolower(distanceType),
           upper = TRUE)
    
    #Get the adjacency matrix that now indicates the distance between the points ~ weighted graph
    distanceMatrix_m <- as.matrix(distanceMatrix)
    relativeNeighbouringMatrix <-
      adjacencyMatrix * distanceMatrix_m#Element by element multiplication!
    #relativeNeighbouringMatrix[relativeNeighbouringMatrix==0] <- NA#For removing "links" that do not exist
    
    #Recreate graph of adjacency
    g <-
      graph.adjacency(relativeNeighbouringMatrix,
                      weighted = TRUE,
                      mode = "undirected")
    
    #Calculate the shortest path for each pair
    s.paths <- shortest.paths(g, algorithm = "dijkstra")
    s.paths.lower <- s.paths
    s.paths.lower[lower.tri(s.paths)] <- 0
    s.paths.lower[is.infinite(s.paths.lower)] <- 0
    
    #Calculate the reticulation index:
    
    #Diameter of the graph
    whichMax <-
      as.data.frame(which(s.paths.lower == max(s.paths.lower, na.rm = TRUE), arr.ind = TRUE))
    valueDiameter <-
      unique(s.paths.lower[whichMax[, 1], whichMax[, 2]])[1]
    
    #Euclidean distance between the points forming diameter of the graph
    vectorDistanceBetweenLocDiameter <-
      distanceMatrix_m[whichMax[, 1], whichMax[, 2]]
    
    #Max distance between edges
    maxEdgeDist <- max(distanceMatrix_m)
    
    reticulationIndex <-
      mean((valueDiameter - vectorDistanceBetweenLocDiameter) / maxEdgeDist)
    
    return(reticulationIndex)
  }

###  Reticulation based on dirichlet tesselation ----------------------------
## mean/var of areas of polygons, divided by number of polygons to account for sampling size. All bordering polygons were removed.
reticulationIndexDirichlet_f <- function(locations_df) {
  locations_df <- as.data.frame(locations_df)
  colnames(locations_df) <- c("x", "y")
  locations_df <- unique(locations_df)
  
  library(deldir)
  library(ggvoronoi)
  library(sp)
  library(rgeos)
  polygons <- voronoi_polygon(locations_df, x = "x", y = "y")
  
  #Identify polygons at the border of square map, and remove
  xlim <- range(locations_df[, 1])
  ylim <- range(locations_df[, 2])
  ##Create spatial lines of squared frame
  coordsBoundaries.df <-
    cbind(rep(xlim, each = 2), c(ylim, rev(ylim)))
  coordsBoundaries.df <-
    rbind(coordsBoundaries.df, coordsBoundaries.df[1, ])
  
  boundries_line <- sp::Line(coordsBoundaries.df)
  boundries_lines <-
    sp::Lines(list(boundries_line), ID = "noImportance")
  boundries_spatialLines <- sp::SpatialLines(list(boundries_lines))
  
  ##Intersect tesselation and squared frame
  whichToRemove <-
    sp::over(polygons, boundries_spatialLines, byid = TRUE)
  
  ##Remove polygons in intersection
  polygons <- polygons[which(is.na(whichToRemove)), ]
  plot(polygons)
  points(distribution.df[, 1], distribution.df[, 2], pch = 19)
  
  #Calculate variance in areas of remaining polygons, and total area
  ##Calculate areas of singular polygon and total
  areaByPolygon.v <- gArea(polygons, byid = TRUE)
  totArea <- sum(areaByPolygon.v)
  
  ##Calculate variance/mean in all polygons
  varArea <- var(areaByPolygon.v)
  meanArea <- mean(areaByPolygon.v)
  
  return(varArea / meanArea / length(areaByPolygon.v))
}

### Pseudo fractal dimension ------------------------------------------------

fractalDimension_f <-
  function(locations_df,
           radius.gwfa,
           bandwith.gwfa,
           sample_size.gwfa,
           cell_size.gwfa)
  {
    library(gwfa) #Check for help: Applying two fractal methods to characterise the local and global deviations from scale invariance of built patterns throughout mainland France
    
    locations_df <- as.data.frame(locations_df)
    colnames(locations_df) <- c("x", "y")
    #estimate the fractal dimension on the different radius: see vignette package gwfa
    test = gwfa(
      points = locations_df,
      q = 0,
      radius = radius.gwfa,
      bandwith = bandwith.gwfa,
      sample_size = sample_size.gwfa,
      cell_size = cell_size.gwfa
    )
    
    dfFractal <- (do.call(cbind, test[, 4:(4 + length(radius.gwfa) - 1)]))
    dfFractal[which(is.infinite(dfFractal))] <- NA
    
    X = cbind(rep(1, length(test@radius)), log2(test@radius))
    fit_frac_dim = dfFractal %*% t(solve(t(X) %*% X) %*% t(X))
    test$dimfrac = fit_frac_dim[, 2]
    # plot(rep(X[,2], each=nrow(do.call(cbind,test[,4:(4+length(radius.gwfa)-1)]))),  as.vector(do.call(cbind,test[,4:(4+length(radius.gwfa)-1)])))
    
    #Return average fractal dimension
    return(mean(test$dimfrac, na.rm = TRUE))
  }

### Lloyd index of patchiness -----------------------------------------------

patchiness_f <-
  function(locations_df,
           minx,
           miny,
           maxx,
           maxy,
           cellsize) {
    #Square cells
    locations_df <- as.data.frame(locations_df)
    locations_df$xcoord_rasterCount <-
      floor((locations_df[, 1] - minx) / cellsize)
    locations_df$ycoord_rasterCount <-
      floor((locations_df[, 2] - miny) / cellsize)
    
    library(tidyr)
    library(dplyr)
    countNumber <-
      locations_df %>% group_by(xcoord_rasterCount, ycoord_rasterCount) %>% count
    
    countNumber <- as.data.frame(countNumber)
    #Readd those with 0
    presentQuadrat <- paste(countNumber[, 1], countNumber[, 2], sep = "_")
    allPossibleQuadrat.df <-
      as.data.frame(tidyr::crossing(
        seq(
          from = minx,
          to = maxx - cellsize,
          by = cellsize
        ) / cellsize,
        seq(
          from = miny,
          to = maxy - cellsize,
          by = cellsize
        ) / cellsize
      ))
    allPossibleQuadrat.v <-
      paste(allPossibleQuadrat.df[, 1], allPossibleQuadrat.df[, 2], sep = "_")
    absentQuadrat <-
      which(!(allPossibleQuadrat.v %in% presentQuadrat))
    
    toAdd <-
      cbind(allPossibleQuadrat.df[absentQuadrat, ], rep(0, times = length(absentQuadrat)))
    colnames(toAdd) <- colnames(countNumber)
    countNumber <- rbind(countNumber,
                         toAdd)
    
    values <- as.numeric(as.character(countNumber$n))
    patchinessIndex <-
      sum(values * (values - 1), na.rm = TRUE) / sum(values, na.rm = TRUE)
    patchinessIndex <- patchinessIndex / mean(values, na.rm = TRUE)
    return(patchinessIndex)
  }

### Point alignement (route index) --------------------------------------------------------

alignmentPoints_f <- function(locations_df, buffer = 0, Npoints = NA){
  locations_df <- as.matrix(locations_df)
  alignment_v <- sapply(1:nrow(locations_df), function(row) {
    #Find the Npoints = 2 closest points
    locations_df_rdc <- locations_df[-row, ]
    diffLoc <- lapply(1:nrow(locations_df_rdc), function(rw) {
      locations_df_rdc[rw, ] - locations_df[row, ]
    })
    diffLoc <- do.call("rbind", diffLoc)
    distance_df <- apply(diffLoc ** 2, 1, sum)
    
    if(!is.na(Npoints)){
      pointsIDincluded <- rep(NA, times = Npoints)
      for (pts in 1:Npoints) {
        if (pts == 1) {
          pointsID <- which(distance_df == min(distance_df))[1]
          pointsIDincluded[pts] <- pointsID
        } else{
          removePoints <- pointsIDincluded[which(!is.na(pointsIDincluded))]
          pointsID <-
            which(distance_df == min(distance_df[-removePoints]))[1]
          pointsID <- pointsID[which(!(pointsID %in% removePoints))][1]
          pointsIDincluded[pts] <- pointsID
        }
      }
    }else{
      whichToBeIncluded_v <- which(distance_df < buffer * buffer)
      if(length(whichToBeIncluded_v) < 2){
        Npoints2 = 2
        pointsIDincluded <- rep(NA, times = Npoints2)
        for (pts in 1:Npoints2) {
          if (pts == 1) {
            pointsID <- which(distance_df == min(distance_df))[1]
            pointsIDincluded[pts] <- pointsID
          } else{
            removePoints <- pointsIDincluded[which(!is.na(pointsIDincluded))]
            pointsID <-
              which(distance_df == min(distance_df[-removePoints]))
            pointsID <- pointsID[which(!(pointsID %in% removePoints))][1]
            
            pointsIDincluded[pts] <- pointsID
          }
        }
      }else{
        pointsIDincluded <- whichToBeIncluded_v
      }
    }
    whichCombination_m <- combn(pointsIDincluded, 2)
    
    alignmentPossible_v <- apply(whichCombination_m, 2, function(pointsIDincluded_rdc){
      compareDf <- rbind(locations_df[row, ],
                         locations_df_rdc[pointsIDincluded_rdc, ]) %>%  as.data.frame()
      colnames(compareDf) <- c("x", "y")
      
      #Calculate the three angles of the triangle formed by those points
      # Law of cosine https://en.wikipedia.org/wiki/Law_of_cosines
      
      #Calculate distances
      x <- as.numeric(compareDf$x)
      y <- as.numeric(compareDf$y)
      dist12 = sqrt((y[2] - y[1])**2 + (x[2] - x[1])**2)
      dist13 = sqrt((y[3] - y[1])**2 + (x[3] - x[1])**2)
      dist32 = sqrt((y[2] - y[3])**2 + (x[2] - x[3])**2)
      
      angle1 <- acos(min(max(-1, (dist12**2 + dist13**2 - dist32**2)/(2*dist12*dist13)), 1))#Using max and min because errors due to rounding otherwise
      angle2 <- acos(min(max(-1, (dist32**2 + dist13**2 - dist12**2)/(2*dist32*dist13)), 1))#Using max and min because errors due to rounding otherwise
      angle3 <- acos(min(max(-1, (dist12**2 + dist32**2 - dist13**2)/(2*dist12*dist32)), 1))#Using max and min because errors due to rounding otherwise
      
      if(round(angle1 + angle2 + angle3, digits = 2) != round(pi, digits = 2)){
        print("Error in the calculus")
        print(row)
      }
      # #Rank them along x and calculate angle diff
      # angleDiffBasedOnx <- compareDf %>%
      #   arrange(x) %>%
      #   summarise(diffAngle = abs(atan((y[2] - y[1]) / (x[2] - x[1])) - atan((y[2] - y[3]) /
      #                                                                          (x[2] - x[3]))),
      #             distance1 = sqrt((y[2] - y[1])**2 + (x[2] - x[1])**2),
      #             distance2 = sqrt((y[2] - y[3])**2 + (x[2] - x[3])**2)
      #   )
      # 
      # # angleDiffBasedOnx$diffAngle[1] <-
      # #   ifelse(angleDiffBasedOnx$diffAngle[1] == pi / 2,
      # #          pi / 2,
      # #          angleDiffBasedOnx$diffAngle[1] %% pi / 2)
      # #Rank them along y and calculate angle diff
      # angleDiffBasedOny <- compareDf %>%
      #   arrange(y) %>%
      #   summarise(diffAngle = abs(atan((y[2] - y[1]) / (x[2] - x[1])) - atan((y[2] - y[3]) /
      #                                                                          (x[2] - x[3]))),
      #             distance1 = sqrt((y[2] - y[1])**2 + (x[2] - x[1])**2),
      #             distance2 = sqrt((y[2] - y[3])**2 + (x[2] - x[3])**2)
      #   )
      # # angleDiffBasedOny$diffAngle[1] <-
      # #   ifelse(angleDiffBasedOny$diffAngle[1] == pi / 2,
      # #          pi / 2,
      # #          angleDiffBasedOny$diffAngle[1] %% pi / 2)
      # 
      # 
      # alignmentBasedOnx <- sin(angleDiffBasedOnx$diffAngle[1])#angleDiffBasedOnx$alignment[1]
      # alignmentBasedOny <- sin(angleDiffBasedOny$diffAngle[1])#angleDiffBasedOny$alignment[1]
      # minValue <- min(alignmentBasedOnx, alignmentBasedOny)
      # 
      # #Calculate weight: ratio between length of min/max edge
      # library(sf)
      # 
      # distance_v <- sapply(1:nrow(compareDf), function(pt){
      #   st_length(st_linestring(compareDf[-pt,] %>% as.matrix()))
      # })
      
      #print(minValue)
      return(min(sin(angle1), sin(angle2), sin(angle3))[1])# * min(distance_v)/max(distance_v))
    })
    
    alignmentCorrected <- 1 - min(alignmentPossible_v)
    alignmentCorrected <- (alignmentCorrected - (1 - sin(pi/3))) / (1 - (1 - sin(pi/3)))
    return(alignmentCorrected)
    
  })
  # alignmentValue = c(mean(alignment_v, na.rm = TRUE),
  #                    median(alignment_v, na.rm = TRUE),
  #                    sd(alignment_v, na.rm = TRUE)
  # )
  # 
  #Interestingly, I noticed that the distribution is differently shaped, and this might be perhaps more informative to look at that.
  #An index could be the difference between the 3rd bins and 5th bins:
  # alignmentValue = (length(alignment_v[alignment_v >= 0.8 & alignment_v <= 1]) - length(alignment_v[alignment_v >= 0.4 & alignment_v <= 0.6]))/(length(alignment_v[alignment_v >= 0.8 & alignment_v <= 1]))
  # Or more straightforwardly, the skewness
  library(modeest)
  alignmentValue = -skewness(alignment_v)
  
  return(alignmentValue)
}

#Test a "grid" pattern"
coords <- t(combn(seq(from = 0, to = 1000, length.out = 32), 2))
coords <- rbind(coords , cbind(coords[,2], coords[,1]))
coords <- rbind(coords, cbind(seq(from = 0, to = 1000, length.out = 32), seq(from = 0, to = 1000, length.out = 32)))
coords <- unique(coords)
nrow(coords)
plot(coords[,1], coords[,2], asp = 1)
alignmentPoints_f(coords)

### Result plot -------------------------------------------------------------

plotResults <- function(
    yAxisName,
    xAxisName,
    xVar,
    yVar,
    df,
    boxplotOnly = TRUE,
    categoricalX = FALSE,
    levelsOldNameX = NA,
    levelsNewNameX = NA,
    groupVar = NA,
    colourGroup_v = NA,
    nameGroup = NA,
    levelOrderGroup = NA,
    widthBox = NA,
    differentMeanShapePoints = FALSE,
    valueSegments = NA,
    lineTypeSegments = NA,
    colourSegments = NA,
    nameSegments = NA
){
  #For benchmark lines (added during revision)
  dfSegments = data.frame(
    valueSegments = valueSegments,
    lineTypeSegments = lineTypeSegments,
    colourSegments = colourSegments,
    nameSegments = nameSegments
  )
  
  #Main plotting
  df <- df %>% mutate(
    x = .data[[xVar]],
    y = .data[[yVar]] 
  )
  if(categoricalX){
    for(i in 1:length(levelsOldNameX)){
      df <- df %>% mutate(
        x = ifelse(x == levelsOldNameX[i], levelsNewNameX[i], x)
      )  
    }
    df <- df %>% 
      mutate(
        x = factor(x, levels = levelsNewNameX)
      )
  }
  
  if(!is.na(groupVar)){
    df <- df %>% mutate(
      group = .data[[groupVar]],
      group = factor(group, levels = levelOrderGroup)
    )
    dodge <- position_dodge(width = 0.5)
  }
  
  howManyShapes <- length(unique(df$group))
  if(howManyShapes > 5){
    print("Cannot deal with a groupin variable > 5 levels and different mean point shapes")
    break
  }
  if(boxplotOnly){
    if(is.na(widthBox)){
      widthBox = 0.25
    }
    if(!is.na(groupVar)){
      if(differentMeanShapePoints){
        plot <- ggplot(df, aes(x = x, y = y)) +
          geom_segment(dfSegments, mapping = aes(x = -Inf, xend = Inf, y = valueSegments, yend = valueSegments, linetype = lineTypeSegments, color = colourSegments)) +
          geom_text(dfSegments, mapping = aes(x = Inf, y = valueSegments, label = nameSegments, colour = colourSegments), vjust = -0.5, hjust = 1, fontface = 3) +
          geom_boxplot(aes(x = x, y = y, group = interaction(x, group), fill = group), width = widthBox, notch = TRUE, position = dodge) +
          scale_colour_manual(values = colourSegments) +
          scale_linetype_identity() +
          guides(linetype = "none", colour = "none") +
          stat_summary(
            #fun.data = give.n,
            geom = "text",
            fun.y = "mean",
            size = 4,
            hjust = -0.55,
            position = dodge,
            aes(label = round(after_stat(y), 3), group = group)
          ) +
          stat_summary(
            geom = "point",
            fun.y = "mean",
            size = 2,
            fill = "white",
            aes(group = group, shape = group),
            position = dodge,
            colour = "black"
          ) +
          scale_shape_manual(name = nameGroup,values = c(24, 21, 25, 22, 23)[1:howManyShapes]) +
          labs(x = xAxisName, y = yAxisName) +
          scale_fill_manual(name = nameGroup, values = colourGroup_v) + 
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
          scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
      }else{
        plot <- ggplot(df, aes(x = x, y = y)) +
          geom_segment(dfSegments, mapping = aes(x = -Inf, xend = Inf, y = valueSegments, yend = valueSegments, linetype = lineTypeSegments, color = colourSegments)) +
          geom_text(dfSegments, mapping = aes(x = Inf, y = valueSegments, label = nameSegments, colour = colourSegments), vjust = -0.5, hjust = 1, fontface = 3) +
          scale_colour_manual(values = colourSegments) +
          scale_linetype_identity() +
          guides(linetype = "none", colour = "none") +
          geom_boxplot(aes(x = x, y = y, group = interaction(x, group), fill = group), width = widthBox, notch = TRUE, position = dodge) +
          stat_summary(
            #fun.data = give.n,
            geom = "text",
            fun.y = "mean",
            size = 4,
            hjust = -0.55,
            position = dodge,
            aes(label = round(after_stat(y), 3), group = group)
          ) +
          stat_summary(
            geom = "point",
            fun.y = "mean",
            size = 2,
            fill = "white",
            shape = 21,
            aes(group = group),
            position = dodge,
            colour = "black"
          ) +
          labs(x = xAxisName, y = yAxisName) +
          scale_fill_manual(name = nameGroup, values = colourGroup_v) + 
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
          scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
      }
    }else{
      plot <- ggplot(df, aes(x = x, y = y)) +
        geom_segment(dfSegments, mapping = aes(x = -Inf, xend = Inf, y = valueSegments, yend = valueSegments, linetype = lineTypeSegments, color = colourSegments)) +
        geom_text(dfSegments, mapping = aes(x = Inf, y = valueSegments, label = nameSegments, colour = colourSegments), vjust = -0.5, hjust = 1, fontface = 3) +
        scale_colour_manual(values = colourSegments) +
        scale_linetype_identity() +
        guides(linetype = "none", colour = "none") +
        geom_boxplot(aes(x = x, y = y, group = x), width = 0.25, notch = TRUE, fill = "grey90") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.55,
          aes(label = round(after_stat(y), 3))
        ) +
        stat_summary(
          geom = "point",
          fun.y = "mean",
          size = 2,
          shape = 21,
          fill = "white",
          colour = "black"
        ) +
        labs(x = xAxisName, y = yAxisName) +
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
        scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
    }
  }else{
    if(!is.na(groupVar)){
      plot <- ggplot(df, aes(x = x, y = y)) +
        geom_segment(dfSegments, mapping = aes(x = -Inf, xend = Inf, y = valueSegments, yend = valueSegments, linetype = lineTypeSegments, color = colourSegments)) +
        geom_text(dfSegments, mapping = aes(x = Inf, y = valueSegments, label = nameSegments, colour = colourSegments), vjust = -0.5, hjust = 1, fontface = 3) +
        geom_violin(aes(x = x, y = y, fill = group), position = dodge) +
        scale_colour_manual(values = colourSegments) +
        scale_linetype_identity() +
        guides(linetype = "none", colour = "none") +
        geom_boxplot(aes(x = x, y = y, group =  interaction(x, group)), position = dodge, width = 0.05, fill = "black") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.55,
          position = dodge,
          aes(label = round(after_stat(y), 3), group = group)
        ) +
        stat_summary(
          geom = "point",
          fun.y = "mean",
          size = 2,
          shape = 21,
          fill = "white",
          aes(group = group),
          position = dodge,
          colour = "black"
        ) +
        labs(x = xAxisName, y = yAxisName) +
        scale_fill_manual(name = nameGroup, values = colourGroup_v) + 
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
        scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
    }else{
      plot <- ggplot(df, aes(x = x, y = y)) +
        geom_segment(dfSegments, mapping = aes(x = -Inf, xend = Inf, y = valueSegments, yend = valueSegments, linetype = lineTypeSegments, color = colourSegments)) +
        geom_text(dfSegments, mapping = aes(x = Inf, y = valueSegments, label = nameSegments, colour = colourSegments), vjust = -0.5, hjust = 1, fontface = 3) +
        geom_violin(aes(x = x, y = y, group = x), fill = "grey90", adjust = 1) +
        scale_colour_manual(values = colourSegments) +
        scale_linetype_identity() +
        guides(linetype = "none", colour = "none") +
        geom_boxplot(aes(x = x, y = y, group = x), width = 0.01, fill = "black") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.55,
          aes(label = round(after_stat(y), 3))
        ) +
        stat_summary(
          geom = "point",
          fun.y = "mean",
          size = 2,
          shape = 21,
          fill = "white",
          colour = "black"
        ) +
        labs(x = xAxisName, y = yAxisName) +
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
        scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
    }
  }
  
  if(!categoricalX){
    plot <- plot +
      scale_x_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4))
  }
  
  return(plot)
}

library(ggforce)

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
    ) +
    coord_fixed()
  return(plot)
}

## Extracting output simulations ---------------------------------------------------------

listFiles_v <- list.files("Output/SensitivitySensoryRange")
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichRoutine_v <- grep("Routine", listFiles_v)
filesRoutine_v <- listFiles_v[whichRoutine_v]

length(filesMapFinal_v)

indices_l <- lapply(
  1:length(filesMapFinal_v),
  function(whatFile){
    #print(whatFile)
    system(paste("echo '", as.character(whatFile), "'"))
    fileOfInterest <- filesMapFinal_v[whatFile]
    if(grepl("High", fileOfInterest)){
      whatRule <- "High"
    }else if(grepl("Intermediate", fileOfInterest)){
      whatRule <- "Intermediate"
    }
    
    outputMapFinal <- read_table2(paste0("Output/SensitivitySensoryRange/", fileOfInterest)) %>% as.data.frame() 
    #plot(outputMapFinal[,1], outputMapFinal[,2])
    
    ### Patchiness ---------------------------------------------------------------
    
    patchiness = patchiness_f(outputMapFinal,
                              0,
                              0,
                              mapSize,
                              mapSize,
                              quadratSize)
    
    ### Alignment ---------------------------------------------------------------
    
    alignment = alignmentPoints_f(
      outputMapFinal
    )
    
    ### Spatial autocorr --------------------------------------------------------
    
    #Create distance matrix
    treeDistance_m <- as.matrix(dist(outputMapFinal[, 1:2]))
    #Use the inverse of distance as weight
    treeDistanceInverse_m <- 1 / (treeDistance_m)
    #Remove pb with self distance
    diag(treeDistanceInverse_m) <- 0
    
    #Moran/Geary Index considering circular data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870747/
    startFruitingInRadian_v <-
      as.data.frame(outputMapFinal)[, 3] * 2 * pi / 365
    
    #Transform to -pi pi interval
    startFruitingInRadian_v <-
      ifelse(
        startFruitingInRadian_v > pi,
        startFruitingInRadian_v - 2 * pi,
        startFruitingInRadian_v
      )
    
    #Calcul spatial autocorr
    spatialAutocorr = autocorrSp(
      startFruitingInRadian_v,
      treeDistanceInverse_m,
      circular = TRUE,
      index = "Moran"
    )
    ### Shrinkage resource ------------------------------------------------------
    
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
          shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/((mapSize + extentUse * mapSize)**2)
          "STOP"
        }, error = function(e){
          return("KEEP GOING")
        }
      )
    }
    
    UD <- kernelUD(xy = SpatialPoints(outputMapFinal[, 1:2]), h = 50,
                   boundary = boundary, extent = extentUse)
    shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/((mapSize + extentUse * mapSize)**2)
    
    ### Routine ---------------------------------------------------------------
    
    seqVisits <- read_csv(paste0("Output/SensitivitySensoryRange/",  filesRoutine_v[whatFile]), 
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
    
    routine <- 1 - entropy_O1(seqVisits[,2])
    
    return(c(whatRule, patchiness, alignment, spatialAutocorr, routine, shrinkage))
  }
)

save.image("Renvironment/resultsSensoryRange.RData")

indicesSensoryRange_df <- do.call("rbind", indices_l) %>% as.data.frame()
colnames(indicesSensoryRange_df) <- c("SensoryRange", "patchiness", "alignment", "spatialAutocorr", "routine", "shrinkage")
indicesSensoryRange_df$patchiness <- as.numeric(indicesSensoryRange_df$patchiness)
indicesSensoryRange_df$shrinkage <- as.numeric(indicesSensoryRange_df$shrinkage)
indicesSensoryRange_df$alignment <- as.numeric(indicesSensoryRange_df$alignment)
indicesSensoryRange_df$spatialAutocorr <- as.numeric(indicesSensoryRange_df$spatialAutocorr)
indicesSensoryRange_df$routine <- as.numeric(indicesSensoryRange_df$routine)

#Add main data when sensory range is low
indicesSensoryRange_df_rdc <- indicesSensoryRange_df
indicesSensoryRange_df <- rbind(indicesSensoryRange_df, dataLowSensoryRangeNaive)

## Plotting ----------------------------------------------------------------

library(ggplot2)
library(scales)

plotPatchinessMoving <- plotResults(
  yAxisName = "Patchiness",
  xAxisName = "Sensory Range",
  xVar = "SensoryRange",
  yVar = "patchiness",
  df = indicesSensoryRange_df %>% mutate(patchiness = patchiness*(1-shrinkage)),#Normalise patchiness by shrinkage
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSensoryRange_df$SensoryRange)[c(3,2,1)],
  levelsNewNameX = c("Large", "Intermediate", "Small")[c(3,2,1)],
  valueSegments = c(patchinessHomogeneousValue, patchinessHeterogeneousValue, patchinessRouteValue),
  lineTypeSegments = c("dotted", "dashed", "solid"),
  colourSegments =  c("black", "black", "black"),
  nameSegments = c("Homogeneous", "Heterogeneous", "Route")
)

plotAlignmentMoving <- plotResults(
  yAxisName = "Alignment",
  xAxisName = "Sensory Range",
  xVar = "SensoryRange",
  yVar = "alignment",
  df = indicesSensoryRange_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSensoryRange_df$SensoryRange)[c(3,2,1)],
  levelsNewNameX = c("Large", "Intermediate", "Small")[c(3,2,1)],
  valueSegments = c(alignmentHomogeneousValue, alignmentRouteValue),
  lineTypeSegments = c("dotted", "solid"),
  colourSegments =  c("black","black"),
  nameSegments = c("Homogeneous",  "Route")
)

plotRoutineMoving <- plotResults(
  yAxisName = "Routine",
  xAxisName = "Sensory Range",
  xVar = "SensoryRange",
  yVar = "routine",
  df = indicesSensoryRange_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSensoryRange_df$SensoryRange)[c(3,2,1)],
  levelsNewNameX = c("Large", "Intermediate", "Small")[c(3,2,1)],
  valueSegments = c(0.5),
  lineTypeSegments = c(NA),
  colourSegments =  c(NA),
  nameSegments = c("")
)

plotSpatAutocorrMoving <- plotResults(
  yAxisName = "Spatial autocorrelation",
  xAxisName = "Sensory Range",
  xVar = "SensoryRange",
  yVar = "spatialAutocorr",
  df = indicesSensoryRange_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSensoryRange_df$SensoryRange)[c(3,2,1)],
  levelsNewNameX = c("Large", "Intermediate", "Small")[c(3,2,1)],
  valueSegments = c(0, moranHomoIntermediate, moranHomoHigh, moranHeteroHigh),
  lineTypeSegments = c("longdash", "dotted", "solid", "dashed"),
  colourSegments =  c("black", "black", "black", "black"),
  nameSegments = c("", "Int. synchro. homo", "High. synchro. homo.", 
                   "High. synchro. hetero.")
)

plotShrinkageSensoryRange <- plotResults(
  yAxisName = "Shrinkage",
  xAxisName = "Sensory Range",
  xVar = "SensoryRange",
  yVar = "shrinkage",
  df = indicesSensoryRange_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSensoryRange_df$SensoryRange)[c(3,2,1)],
  levelsNewNameX = c("Large", "Intermediate", "Small")[c(3,2,1)],
  valueSegments = c(0.5),
  lineTypeSegments = c(NA),
  colourSegments =  c(NA),
  nameSegments = c("")
)

library(ggpubr)
mergedPlot <- ggarrange(
  plotPatchinessMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.8,1.6)),
  plotAlignmentMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.1,0.8)),
  plotSpatAutocorrMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.05,0.05)),
  plotRoutineMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.5,1)),
  nrow = 4, 
  ncol = 1
)
mergedPlotSensoryRange <- annotate_figure(mergedPlot,
                                        top = text_grob("Sensory range", face = "bold", size = 16))
mergedPlotSensoryRange

saveRDS(mergedPlotSensoryRange, "Renvironment/Plots/mergedPlotSensoryRange.rds")

