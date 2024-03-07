# What does this script do?
# This script simulates distributions with different spatial autocorrelation in fruiting dates to estimate reference values for comparisons to model output.

rm(list = ls())

# Environment set up ------------------------------------------------------

## Libraries ---------------------------------------------------------------

library(gstat)
library(dplyr)
library(sf)
library(circular)

## Functions ---------------------------------------------------------------

### Source Rcpp functions ---------------------------------------------------

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp") #(contains distributionTree)

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

### Calculate matrix of inverse distance between all points -----------------

matrixDistance_f <- function(data){
  matrixDistances <- lapply(1:nrow(data), function(row){
    distance_v <- sqrt((data[,1] - data[row, 1])**2 + (data[,2] - data[row, 2])**2)
    return(distance_v)
  })
  matrixDistances <- do.call("cbind",  matrixDistances)
  matrixDistances <- 1/matrixDistances
  diag(matrixDistances) <- 0#Because no correlation with itself
  return(matrixDistances)
}

### Create spatial distribution -----------------------------------------------------

createDistribution <- function(nTree, mapSize, clusterNumber = nTree / 10, spreading = 50, type = c("Homogeneous", "Heterogeneous")){
  if(type == "Homogeneous"){
    xy <- cbind(
      runif(nTrees, 0, mapSize),
      runif(nTrees, 0, mapSize)
    ) %>%
      as.data.frame() %>% 
      rename(
        x = "V1", 
        y = "V2"
      )
  }else{
    xy <- distributionTree(
      numberTrees = nTree,
      lowerBorder = 0,
      upperBorder = mapSize,
      homogeneousDistribution = FALSE,
      treeClusterNumber = clusterNumber,
      treeClusterSpread = spreading
    )
    colnames(xy) <- c("x", "y")
  }
  return(xy)
}

### Create phenology --------------------------------------------------------

simulateCircularAutocorr <- function(data, #xy data, with colnames x and y
                                     corr,
                                     seed#coefficient of corr between 0 and 1
){
  library(svMisc)
  print("Initialising, it will start shortly.")
  set.seed(seed)
  #Matrix of distance to neighbours, to speed up the calculus
  matrixDistances <- matrixDistance_f(data)
  
  #Simulate the circular variable: 1) Initialise step by step with the closest neighbours already simulated 2) Repeat 1 time to stabilise the distribution a bit (i.e. that the first simulated values also are correlated to their neighbours as desired)
  data$circularVar <- NA
  data$circularVar[1:floor(0.05*nrow(data))] <- rwrappednormal(floor(0.05*nrow(data)), mu = runif(1, 0, 2*pi), rho = 0.5)
  
  #First go to the closest neighbour of the point of interest successively to sample a value
  pointOfInterest <- 1
  hasBeenUsed_v <- rep(FALSE, times = nrow(data))
  hasBeenUsed_v[1:floor(0.05*nrow(data))] <- TRUE
  for (i in floor(0.05*nrow(data)):nrow(data)){
    #Simulate the circular var value from the closest point
    pointOfInterest <- which.max(matrixDistances[,pointOfInterest]*ifelse(hasBeenUsed_v, -1000, 1))#Take the closest point not used
    pointClosest <- which.max(matrixDistances[,pointOfInterest]*ifelse(hasBeenUsed_v, 1, -1000))#Take the closest point used
    data$circularVar[pointOfInterest] <- rwrappednormal(1, mu = data$circularVar[pointClosest], rho = corr)
    hasBeenUsed_v[pointOfInterest] <- TRUE
    # print(pointOfInterest)
    # print(pointClosest)
  }
  
  #In matrix distance, reput itself as "closest" value
  diag(matrixDistances) <- 1 
  
  #Then do one homogeneisation based on neighbouring points (5% of all points) initial value
  newValues <- sapply(1:nrow(data), function(row){
    progress(row/nrow(data)*100)
    pointsOfInterest_v <- which(matrixDistances[,row] >= quantile(matrixDistances[,row], 0.95))#1:nrow(data)#
    newValues <- sapply(pointsOfInterest_v, function(row2){
      rwrappednormal(1, mu = data$circularVar[row2], rho = corr)
    })
    circularValue <- weighted.mean.circular(x = newValues, w = matrixDistances[pointsOfInterest_v,row]/sum(as.numeric(matrixDistances[pointsOfInterest_v,row])), na.rm = TRUE)
    return(as.numeric(circularValue))
  })
  
  data$circularVar <- newValues
  return(data)
}

### Plotting map for spatial autocorrelation visualisation ------------------

plotMapAutocorr <- function(data, mapSize, gradientColours, pointSize){
  library(ggplot2)
  ggplot() +
    geom_path(aes(x = c(0,0,mapSize,mapSize,0), y = c(0, mapSize, mapSize, 0, 0))) +
    geom_point(data, mapping = aes(x = x, y = y, fill = circularVar), pch = 21, col = "black", size = pointSize) + 
    scale_fill_gradient2(low = gradientColours[1], mid = gradientColours[2], high = gradientColours[3], midpoint = pi) +
    guides(fill = "none") +
    theme_void() +
    theme(    
      strip.text.x = element_text(face = "bold", colour = "black"),
      strip.text.y = element_text(face = "bold", colour = "black"),
      strip.background = element_rect(colour = NA, fill = "white"),
      panel.grid.minor = element_line(colour = "white"),
      panel.grid.major = element_line(colour = "white"),
    ) +
    coord_fixed()
}

## Parameters --------------------------------------------------------------

mapSize <- 1000
nTrees <- 1000
nSimulations <- 200
corrValues_v <- c(0, 0.5, 1)

# Simulation of index for different conditions ----------------------------

library(parallel)
simulatedDateSpAutocorr <- mclapply(1:length(corrValues_v),
                                    mc.cores = 3,
                                    function(i){
  
  dataOutput <- data.frame(
    spAutoCorrHomogeneous = rep(NA, times = nSimulations),
    spAutoCorrHeterogeneous = rep(NA, times = nSimulations)
  )
  
  set.seed(42)
  for(j in 1:nSimulations){
    # Create trees
    xyHomogeneous <- createDistribution(nTree = nTrees, mapSize = mapSize, clusterNumber = nTrees/10, spreading = 50, type = "Homogeneous")
    xyHeterogeneous <- createDistribution(nTree = nTrees, mapSize = mapSize, clusterNumber = nTrees/10, spreading = 50, type = "Heterogeneous")
    
    # Simulate the circular variable
    # The measure of correlation will be related to the value of the variance of the circular Gaussian used
    
    #Homogeneous
    dataOfInterest <- simulateCircularAutocorr(xyHomogeneous, corrValues_v[i], seed = j)
    autoCorrValueHomogeneous <- autocorrSp(
      x = dataOfInterest$circularVar,
      weight = matrixDistance_f(dataOfInterest),
      index = "Moran",
      circular = TRUE,
      na.rm = TRUE
    )
    
    #Heterogeneous
    dataOfInterest <- simulateCircularAutocorr(xyHeterogeneous %>%  as.data.frame(), corrValues_v[i], seed = j)
    autoCorrValueHeterogeneous <- autocorrSp(
      x = dataOfInterest$circularVar,
      weight = matrixDistance_f(dataOfInterest),
      index = "Moran",
      circular = TRUE,
      na.rm = TRUE
    )
    
    dataOutput[j,] <- c(autoCorrValueHomogeneous, autoCorrValueHeterogeneous)
  }
  dataOutput$corr <- corrValues_v[i]
  return(dataOutput)
})

spAutoCorrIndex_df <- do.call("rbind", simulatedDateSpAutocorr)

# Visual illustration of environments -------------------------------------

plotExamplesSpatCorr_l <- lapply(1:length(corrValues_v), function(i){

  set.seed(42)
  # Create trees
  xyHomogeneous <- createDistribution(nTree = nTrees, mapSize = mapSize, clusterNumber = nTrees/10, spreading = 50, type = "Homogeneous")
  xyHeterogeneous <- createDistribution(nTree = nTrees, mapSize = mapSize, clusterNumber = nTrees/10, spreading = 50, type = "Heterogeneous")
  
  # Simulate the circular variable
  # The measure of correlation will be related to the value of the variance of the circular Gaussian used
  
  #Homogeneous
  dataOfInterest <- simulateCircularAutocorr(xyHomogeneous, corrValues_v[i], seed = 42)
  plotHomogeneous <- plotMapAutocorr(dataOfInterest, mapSize, gradientColours = c("white", "black", "white"), pointSize = 4)
  autoCorrValueHomogeneous <- autocorrSp(
    x = dataOfInterest$circularVar,
    weight = matrixDistance_f(dataOfInterest),
    index = "Moran",
    circular = TRUE,
    na.rm = TRUE
  )
  plotHomogeneous <- plotHomogeneous +
    ggtitle(paste0("Moran's I = ", round(autoCorrValueHomogeneous, digits = 3))) +
    theme(plot.title = element_text(face = 3, hjust = 0.05, size = 22))
  
  #Heterogeneous
  dataOfInterest <- simulateCircularAutocorr(xyHeterogeneous %>%  as.data.frame(), corrValues_v[i], seed = 42)
  plotHeterogeneous <- plotMapAutocorr(dataOfInterest, mapSize, gradientColours = c("white", "black", "white"), pointSize = 4)
  autoCorrValueHeterogeneous <- autocorrSp(
    x = dataOfInterest$circularVar,
    weight = matrixDistance_f(dataOfInterest),
    index = "Moran",
    circular = TRUE,
    na.rm = TRUE
  )
  plotHeterogeneous <- plotHeterogeneous +
    ggtitle(paste0("Moran's I = ", round(autoCorrValueHeterogeneous, digits = 3))) +
    theme(plot.title = element_text(face = 3, hjust = 0.05, size = 22))
          
  return(list(plotHomogeneous, plotHeterogeneous))
    
})

save.image("Renvironment/TestSpatAutocorr.RData")
