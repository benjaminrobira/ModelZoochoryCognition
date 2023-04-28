##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing the different spatial indices to describe resource distribution pattern
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# What does this script do?
# This script tests different spatial parameters to see whether it is possible to differentiate homogeneous, heterogeneous and route organisation of points (from a point pattern).

rm(list = ls())

# Setting up environment --------------------------------------------------

## Parameters -------------------------------------------------------------

#Load parameters
load("Scripts/R/Parameterisation.RData")
nTreeLow = 100
nTreeHigh = 1000
numberOfRepetitionsForMetrics = 200
distribution.df_tot <- c()

## Packages ----------------------------------------------------------------

library(readr)
library(dplyr)
library(progress)
library(Rcpp)
Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")
library(adehabitatLT)

## Functions ---------------------------------------------------------------

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
      as.data.frame(crossing(
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


### Angle variations --------------------------------------------------------

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
#Test a "grid" pattern"
coords <- t(combn(seq(from = 0, to = 1000, length.out = 32), 2))
coords <- rbind(coords , cbind(coords[,2], coords[,1]))
coords <- rbind(coords, cbind(seq(from = 0, to = 1000, length.out = 32), seq(from = 0, to = 1000, length.out = 32)))
coords <- unique(coords)
nrow(coords)
plot(coords[,1], coords[,2], asp = 1)

alignmentPoints_f(coords)

# Test distance to use ----------------------------------------------------

# resultDistanceReticulation <- matrix(NA, nrow=numberOfRepetitionsForMetrics*3*4, ncol=3)
# pb <- progress_bar$new(format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
#                        total = numberOfRepetitionsForMetrics,
#                        complete = "=",   # Completion bar character
#                        incomplete = "-", # Incomplete bar character
#                        current = ">",    # Current bar character
#                        clear = FALSE,    # If TRUE, clears the bar when finish
#                        width = numberOfRepetitionsForMetrics)      # Width of the progress bar
#
# for(i in 1:floor(numberOfRepetitionsForMetrics)){
#
#   # Random distribution
#   toAddRow <- which(is.na(resultDistanceReticulation[,1]))[1]
#   distribution.df <- cbind(runif(nTreeHigh, 0, 1000), runif(nTreeLow, 0, 1000))
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),1] <- rep("Homo.", times=4)
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),2] <- c("Euclidean", "Manhattan", "Chebychev", "Canberra")
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),3] <- c(reticulationIndex_f(distribution.df, distanceType ="Euclidean"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Manhattan"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Maximum"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Canberra"))
#
#   # Clustered distribution
#   toAddRow <- which(is.na(resultDistanceReticulation[,1]))[1]
#   distribution.df <- distributionTree(
#     numberTrees = nTreeHigh,
#     lowerBorder = 0,
#     upperBorder = 1000,
#     homogeneousDistribution = FALSE,
#     treeClusterNumber = nTreeHigh/10,
#     treeClusterSpread = 50
#   )
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),1] <- rep("Hete.", times=4)
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),2] <- c("Euclidean", "Manhattan", "Chebychev", "Canberra")
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),3] <- c(reticulationIndex_f(distribution.df, distanceType = "Euclidean"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Manhattan"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Maximum"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Canberra"))
#
#   # Route distribution
#   toAddRow <- which(is.na(resultDistanceReticulation[,1]))[1]
#   distribution.df <- c()
#   for(r in 1:10){
#     xCoord <- 2000
#     yCoord <- 2000
#     ##Simulate routes
#     while(max(xCoord) > 1000 | max(yCoord) > 1000 | min(xCoord) < 0 | min(yCoord) < 0){
#       test <- simm.crw(date=seq(from=1, to=100, by=1), h = 20, r = 0.9,
#                        x0=c(0,0), id="A1", burst="A1",
#                        typeII=TRUE)
#
#       xCoord <- test[[1]]$x
#       yCoord <- test[[1]]$y
#       xCoord <- xCoord - min(xCoord)
#       yCoord <- yCoord - min(yCoord)
#
#       #To not have all routes starting from same point, changing starting point = flipping the distribution so that starting point can be one of the four corner
#       whatToDo <- ceiling(r/3 + 0.33)#sample(c(1,2,3,4), 1)
#       if(whatToDo == 1){#Top left corner = vertical flip
#         yCoord <- 1000 - yCoord
#       }else if(whatToDo == 2){#Top right corner = horizontal + vertical flip
#         yCoord <- 1000 - yCoord
#         xCoord <- 1000 - xCoord
#       }else if(whatToDo == 3){#Bottom right corner = horizontal flip
#         xCoord <- 1000 - xCoord
#       }else{#Keep bottom left corner
#       }
#     }
#     distribution.df <- rbind(distribution.df, as.data.frame(cbind(xCoord, yCoord)))
#   }
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),1] <- rep("Route", times=4)
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),2] <- c("Euclidean", "Manhattan", "Chebychev", "Canberra")
#   resultDistanceReticulation[toAddRow:(toAddRow+4-1),3] <- c(reticulationIndex_f(distribution.df, distanceType = "Euclidean"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Manhattan"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Maximum"),
#                                                              reticulationIndex_f(distribution.df, distanceType = "Canberra"))
#   pb$tick()
# }
#
# resultDistanceReticulation <- as.data.frame(resultDistanceReticulation)
# colnames(resultDistanceReticulation) <- c("Distribution", "Distance", "Value")
# resultDistanceReticulation$Value <- as.numeric(resultDistanceReticulation$Value)
#
# resultDistanceReticulation <- resultDistanceReticulation[!is.na(resultDistanceReticulation[,1]),]
# library(ggplot2)
# library(viridis)
# # function for mean labels
# mean.n <- function(x){
#   return(c(y = mean(x) + 0.1*mean(x), label = round(mean(x),2)))
#   # experiment with the multiplier to find the perfect position
# }
#
# resultDistanceReticulation$Distribution <- factor(resultDistanceReticulation$Distribution,      # Reordering group factor levels
#                                                   levels = c("Homo.", "Hete.", "Route"))
#
# distanceReticulation.ggplot <- ggplot(resultDistanceReticulation , aes(x=Distribution, Value, fill=Distribution)) +
#   facet_wrap(~ Distance, scale="free")+#, switch="y")+
#   #geom_jitter(aes(color = Distribution),
#   #width=0.15, alpha = 0.6)+
#   geom_violin(position = position_nudge(x = 0, y = 0)) +
#
#   # geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
#   # geom_jitter(aes(color = Distribution),
#   #             width=0.15, alpha = 0.6)+
#
#   stat_summary(
#     fun.data = mean.n,
#     geom = "text",
#     fun.y = "mean",
#     size = 3
#   ) +
#   stat_summary(
#     geom = "point",
#     fun.y = "mean",
#     size = 2,
#     shape = 19
#   ) +
#   labs(x="DISTRIBUTION", y="VALUE")+
#   #scale_y_continuous(sec.axis = sec_axis(~ . , name = "INDEX", breaks = NULL, labels = NULL)) +#, guide=guide_axis(position="left")
#   #scale_x_discrete(sec.axis = sec_axis(~ . , name = "DENSITY", breaks = NULL, labels = NULL)) +
#   scale_fill_viridis(discrete=TRUE)+
#   scale_color_viridis(discrete=TRUE)+
#
#   theme(legend.position="none",
#         #Grid
#         panel.grid.major = element_line(colour = "gray93", size=0.65),
#         panel.grid.minor = element_line(colour = "gray97"),
#         #Facet wrap
#         strip.text.x = element_text(face="bold", colour="white"),
#         strip.text.y = element_text(face="bold", colour="white"),
#         strip.background = element_rect(colour="black", fill="black")
#   )
#
# distanceReticulation.ggplot
#
# save.image("T:/Saved_PhD/Model_dispersalSeed/ModelZoochoryCognition/Renvironment/resultsTestDistance.RData")

# load("T:/Saved_PhD/Model_dispersalSeed/ModelZoochoryCognition/Renvironment/resultsTestDistance.RData")

# Test index -------------------------------------------------------

## Random distribution -----------------------------------------------------

randomReticulationDirichletLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
randomReticulationDirichletHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

randomReticulationLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
randomReticulationHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

randomFractalLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
randomFractalHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

randomPatchinessLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
randomPatchinessHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

randomAlignmentLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
randomAlignmentHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

### Homogeneous distribution low density ------------------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar
pb$tick(0)
for (i in 1:numberOfRepetitionsForMetrics) {
  distribution.df <-
    cbind(runif(nTreeLow, 0, 1000), runif(nTreeLow, 0, 1000))
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Homogeneous", times = nrow(distribution.df)),
                   rep("Low", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  randomReticulationDirichletLowDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  randomReticulationLowDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Manhattan")
  randomFractalLowDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  randomPatchinessLowDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  randomAlignmentLowDensity[i] <- alignmentPoints_f(distribution.df)
  
  pb$tick()
}

### Homogeneous distribution high density ----------------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar

for (i in 1:numberOfRepetitionsForMetrics) {
  distribution.df <-
    cbind(runif(nTreeHigh, 0, 1000), runif(nTreeHigh, 0, 1000))
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Homogeneous", times = nrow(distribution.df)),
                   rep("High", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  randomReticulationDirichletHighDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  randomReticulationHighDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Manhattan")
  randomFractalHighDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  randomPatchinessHighDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  randomAlignmentHighDensity[i] <-
    alignmentPoints_f(distribution.df)
  
  pb$tick()
}

## Clustered distribution --------------------------------------------------

clusteredReticulationDirichletLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
clusteredReticulationDirichletHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

clusteredReticulationLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
clusteredReticulationHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

clusteredFractalLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
clusteredFractalHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

clusteredPatchinessLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
clusteredPatchinessHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

clusteredAlignmentLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
clusteredAlignmentHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

### Clustered distribution low density --------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar

for (i in 1:numberOfRepetitionsForMetrics) {
  distribution.df <- distributionTree(
    numberTrees = nTreeLow,
    lowerBorder = 0,
    upperBorder = 1000,
    homogeneousDistribution = FALSE,
    treeClusterNumber = nTreeLow / 10,
    treeClusterSpread = 50
  )
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Heterogeneous", times = nrow(distribution.df)),
                   rep("Low", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  clusteredReticulationDirichletLowDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  clusteredReticulationLowDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Euclidean")
  clusteredFractalLowDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  clusteredPatchinessLowDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  clusteredAlignmentLowDensity[i] <-
    alignmentPoints_f(distribution.df)
  
  pb$tick()
}

### Clustered distribution high density --------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar

for (i in 1:numberOfRepetitionsForMetrics) {
  distribution.df <- distributionTree(
    numberTrees = nTreeHigh,
    lowerBorder = 0,
    upperBorder = 1000,
    homogeneousDistribution = FALSE,
    treeClusterNumber = nTreeHigh / 10,
    treeClusterSpread = 50
  )
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Heterogeneous", times = nrow(distribution.df)),
                   rep("High", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  clusteredReticulationDirichletHighDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  clusteredReticulationHighDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Manhattan")
  clusteredFractalHighDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  clusteredPatchinessHighDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  clusteredAlignmentHighDensity[i] <-
    alignmentPoints_f(distribution.df)
  
  pb$tick()
}

## Route distribution ------------------------------------------------------

routeReticulationDirichletLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
routeReticulationDirichletHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

routeReticulationLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
routeReticulationHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

routeFractalLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
routeFractalHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

routePatchinessLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
routePatchinessHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

routeAlignmentLowDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)
routeAlignmentHighDensity <-
  rep(NA, times = numberOfRepetitionsForMetrics)

### Route distribution low density ------------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar

for (i in 1:numberOfRepetitionsForMetrics) {
  xCoord <- 2000
  yCoord <- 2000
  ##Simulate routes
  while (max(xCoord) > 1000 | max(yCoord) > 1000) {
    test <- simm.crw(
      date = seq(from = 1, to = 100, by = 1),
      h = 20,
      r = 0.9,
      x0 = c(0, 0),
      id = "A1",
      burst = "A1",
      typeII = TRUE
    )
    xCoord <- test[[1]]$x
    yCoord <- test[[1]]$y
    xCoord <- xCoord - min(xCoord)
    yCoord <- yCoord - min(yCoord)
  }
  distribution.df <- as.data.frame(cbind(xCoord, yCoord))
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Route", times = nrow(distribution.df)),
                   rep("Low", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  routeReticulationDirichletLowDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  routeReticulationLowDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Manhattan")
  routeFractalLowDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  routePatchinessLowDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  routeAlignmentLowDensity[i] <- alignmentPoints_f(distribution.df)
  
  pb$tick()
}

### Route distribution high density ------------------------------------------

pb <-
  progress_bar$new(
    format = "(:spin)[:bar] :percent[Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = numberOfRepetitionsForMetrics,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = numberOfRepetitionsForMetrics
  )      # Width of the progress bar

for (i in 1:numberOfRepetitionsForMetrics) {
  distribution.df <- c()
  for (r in 1:10) {
    xCoord <- 2000
    yCoord <- 2000
    ##Simulate routes
    while (max(xCoord) > 1000 |
           max(yCoord) > 1000 | min(xCoord) < 0 | min(yCoord) < 0) {
      
      proba <- floor(runif(1, 0, 2))
      test <- simm.crw(
        date = seq(
          from = 1,
          to = 100,
          by = 1
        ),
        h = 20,
        r = 0.9,
        x0 = c(ifelse(proba == 0, 0, runif(1, 0, 1000)), ifelse(proba == 0, runif(1, 0, 1000), 0)),
        id = "A1",
        burst = "A1",
        typeII = TRUE
      )
      
      xCoord <- test[[1]]$x
      yCoord <- test[[1]]$y
      xCoord <- xCoord - min(xCoord)
      yCoord <- yCoord - min(yCoord)
      
      
      #To not have all routes starting from same point, changing starting point = flipping the distribution so that starting point can be one of the four corner
      whatToDo <- ceiling(r / 3 + 0.33)#sample(c(1,2,3,4), 1)
      if (whatToDo == 1) {
        #Top left corner = vertical flip
        yCoord <- 1000 - yCoord
      } else if (whatToDo == 2) {
        #Top right corner = horizontal + vertical flip
        yCoord <- 1000 - yCoord
        xCoord <- 1000 - xCoord
      } else if (whatToDo == 3) {
        #Bottom right corner = horizontal flip
        xCoord <- 1000 - xCoord
      } else{
        #Keep bottom left corner
      }
    }
    distribution.df <-
      rbind(distribution.df, as.data.frame(cbind(xCoord, yCoord)))
  }
  if (i == 1) {
    colnames(distribution.df) <- c("x", "y")
    toAdd <- cbind(rep("Route", times = nrow(distribution.df)),
                   rep("High", times = nrow(distribution.df)),
                   distribution.df)
    colnames(toAdd) <- c("Distribution", "Density", "x", "y")
    distribution.df_tot <- rbind(distribution.df_tot, toAdd)
  }
  routeReticulationDirichletHighDensity[i] <-
    reticulationIndexDirichlet_f(distribution.df)
  routeReticulationHighDensity[i] <-
    reticulationIndex_f(distribution.df, distanceType = "Manhattan")
  routeFractalHighDensity[i] <-
    fractalDimension_f(
      distribution.df,
      radius.gwfa = seq(
        from = 20,
        to = 60,
        length.out = 10
      ),
      bandwith.gwfa = 200,
      sample_size.gwfa = 500,
      cell_size.gwfa = 200
    )
  routePatchinessHighDensity[i] <-
    patchiness_f(
      distribution.df,
      minx = 0,
      miny = 0,
      maxx = 1000,
      maxy = 1000,
      cellsize = 50
    )
  routeAlignmentHighDensity[i] <- alignmentPoints_f(distribution.df)
  
  pb$tick()
}


## Plotting the distribution -----------------------------------------------

distribution.df_tot$x <- as.numeric(distribution.df_tot$x)
distribution.df_tot$y <- as.numeric(distribution.df_tot$y)

library(ggplot2)
library(viridis)
theme_set(theme_bw(15))

#Reordering labels for facet_grid

distribution.df_tot$Density <-
  factor(distribution.df_tot$Density,      # Reordering group factor levels
         levels = c("Low", "High"))
distribution.df_tot$Distribution <-
  factor(
    distribution.df_tot$Distribution,
    # Reordering group factor levels
    levels = c("Homogeneous", "Heterogeneous", "Route")
  )
distribution.df_tot$x <- as.numeric(distribution.df_tot$x)
distribution.df_tot$y <- as.numeric(distribution.df_tot$y)

distribution.ggplot <- ggplot(distribution.df_tot, aes(x = x, y = y)) +
  geom_point() +
  facet_grid(Distribution ~ Density, switch = "y") +
  labs(x = "", y = "") +
  scale_y_continuous(sec.axis = sec_axis(
    ~ . ,
    name = "DISTRIBUTION",
    breaks = NULL,
    labels = NULL
  )) + #, guide=guide_axis(position="left")
  scale_x_continuous(sec.axis = sec_axis(
    ~ . ,
    name = "DENSITY",
    breaks = NULL,
    labels = NULL
  )) +
  theme(
    legend.position = "none",
    #Labs
    # axis.title.y = element_blank(),#element_text(face = "italic"),
    # axis.title.x = element_blank(),#element_text(face = "italic"),
    axis.text.x = element_blank(),
    #remove x axis labels
    axis.ticks.x = element_blank(),
    #remove x axis ticks
    axis.text.y = element_blank(),
    #remove y axis labels
    axis.ticks.y = element_blank(),
    #remove y axis ticks
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold", colour = "white"),
    strip.text.y = element_text(face = "bold", colour = "white"),
    strip.background = element_rect(colour = NA, fill = "black")
  )
distribution.ggplot

## Plotting results -------------

##General layout
library(ggforce)
library(ggplot2)
library(viridis)
theme_set(theme_bw(15))

source(
  "https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R"
)

#Patchiness

Patchiness.df <- as.data.frame(rbind(
  #Low density
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomPatchinessLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredPatchinessLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routePatchinessLowDensity
  ),
  
  #High density
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomPatchinessHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredPatchinessHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routePatchinessHighDensity
  )
))

colnames(Patchiness.df) <- c("Density", "Distribution", "Value")
Patchiness.df$Value <- as.numeric(as.character(Patchiness.df$Value))

Patchiness.ggplot <-
  ggplot(Patchiness.df, aes(x = Distribution, Value, fill = Distribution)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_jitter(aes(color = Distribution),
              width = 0.15, alpha = 0.6) +
  labs(x = "Distribution", y = "Patchiness") +
  scale_y_continuous(n.breaks = 10) + #, minor_breaks = seq(from=0, to=3, by=0.0125))+
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap( ~ Density, ncol = 2, scales = "free") +
  theme(
    legend.position = "none",
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold"),
    strip.background = element_rect(colour = NA, fill = "white")
  )

#Fractal

Fractal.df <- as.data.frame(rbind(
  #Low density
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomFractalLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredFractalLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeFractalLowDensity
  ),
  #High density
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomFractalHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredFractalHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeFractalHighDensity
  )
))

colnames(Fractal.df) <- c("Density", "Distribution", "Value")
Fractal.df$Value <- as.numeric(as.character(Fractal.df$Value))

Fractal.ggplot <-
  ggplot(Fractal.df, aes(x = Distribution, Value, fill = Distribution)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_jitter(aes(color = Distribution),
              width = 0.15, alpha = 0.6) +
  labs(x = "Distribution", y = "Fractal") +
  scale_y_continuous(n.breaks = 10) + #, minor_breaks = seq(from=0, to=3, by=0.0125))+
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap( ~ Density, ncol = 2, scales = "free") +
  theme(
    legend.position = "none",
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold"),
    strip.background = element_rect(colour = NA, fill = "white")
  )

#Reticulation

Reticulation.df <- as.data.frame(rbind(
  #Low density
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomReticulationLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredReticulationLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeReticulationLowDensity
  ),
  #High density
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomReticulationHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredReticulationHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeReticulationHighDensity
  )
))

colnames(Reticulation.df) <- c("Density", "Distribution", "Value")
Reticulation.df$Value <-
  as.numeric(as.character(Reticulation.df$Value))

Reticulation.ggplot <-
  ggplot(Reticulation.df,
         aes(x = Distribution, Value, fill = Distribution)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_jitter(aes(color = Distribution),
              width = 0.15, alpha = 0.6) +
  labs(x = "Distribution", y = "Reticulation") +
  scale_y_continuous(n.breaks = 10) + #, minor_breaks = seq(from=0, to=3, by=0.0125))+
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap( ~ Density, ncol = 2, scales = "free") +
  theme(
    legend.position = "none",
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold"),
    strip.background = element_rect(colour = NA, fill = "white")
  )


#ReticulationDirichlet

ReticulationDirichlet.df <- as.data.frame(rbind(
  #Low density
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomReticulationDirichletLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredReticulationDirichletLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeReticulationDirichletLowDensity
  ),
  #High density
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomReticulationDirichletHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredReticulationDirichletHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeReticulationDirichletHighDensity
  )
))

colnames(ReticulationDirichlet.df) <-
  c("Density", "Distribution", "Value")
ReticulationDirichlet.df$Value <-
  as.numeric(as.character(ReticulationDirichlet.df$Value))

ReticulationDirichlet.ggplot <-
  ggplot(ReticulationDirichlet.df,
         aes(x = Distribution, Value, fill = Distribution)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_jitter(aes(color = Distribution),
              width = 0.15, alpha = 0.6) +
  labs(x = "Distribution", y = "ReticulationDirichlet") +
  scale_y_continuous(n.breaks = 10) + #, minor_breaks = seq(from=0, to=3, by=0.0125))+
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap( ~ Density, ncol = 2, scales = "free") +
  theme(
    legend.position = "none",
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold"),
    strip.background = element_rect(colour = NA, fill = "white")
  )

Alignment.df <- as.data.frame(rbind(
  #Low density
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomAlignmentLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredAlignmentLowDensity
  ),
  cbind(
    rep("Low", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeAlignmentLowDensity
  ),
  #High density
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Homo.", times = numberOfRepetitionsForMetrics),
    randomAlignmentHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Hete.", times = numberOfRepetitionsForMetrics),
    clusteredAlignmentHighDensity
  ),
  cbind(
    rep("High", times = numberOfRepetitionsForMetrics),
    rep("Route", times = numberOfRepetitionsForMetrics),
    routeAlignmentHighDensity
  )
))

colnames(Alignment.df) <- c("Density", "Distribution", "Value")
Alignment.df$Value <- as.numeric(as.character(Alignment.df$Value))

Alignment.ggplot <-
  ggplot(Alignment.df, aes(x = Distribution, Value, fill = Distribution)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  geom_jitter(aes(color = Distribution),
              width = 0.15, alpha = 0.6) +
  labs(x = "Distribution", y = "Alignment") +
  scale_y_continuous(n.breaks = 10) + #, minor_breaks = seq(from=0, to=3, by=0.0125))+
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap( ~ Density, ncol = 2, scales = "free") +
  theme(
    legend.position = "none",
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold"),
    strip.background = element_rect(colour = NA, fill = "white")
  )


library(gridExtra)
grid.arrange(
  Patchiness.ggplot,
  Fractal.ggplot,
  Reticulation.ggplot,
  ReticulationDirichlet.ggplot,
  Alignment.ggplot,
  ncol = 1
)
#I don't like that x axes are always present, trying with facet_grid

df1 <-     cbind(rep("Patchiness", times = nrow(Patchiness.df)),
                 Patchiness.df)
df2 <-     cbind(rep("Fractal", times = nrow(Fractal.df)),
                 Fractal.df)
df3 <-     cbind(rep("Reticulation", times = nrow(Reticulation.df)),
                 Reticulation.df)
df4 <-     cbind(rep("ReticulationDirichlet", times = nrow(ReticulationDirichlet.df)),
                 ReticulationDirichlet.df)
df5 <-     cbind(rep("Alignment", times = nrow(Alignment.df)),
                 Alignment.df)

#rbind with different colnames
global.df <- mapply(c, df1, df2, df3, df4, df5)
global.df <- as.data.frame(global.df)
colnames(global.df)[1] <- "Index"
global.df$Value <- as.numeric(as.character(global.df$Value))

global.df$Distribution <-
  factor(global.df$Distribution,      # Reordering group factor levels
         levels = c("Homo.", "Hete.", "Route"))
global.df$Index <-
  factor(
    global.df$Index,
    # Reordering group factor levels
    levels = c(
      "Patchiness",
      "Fractal",
      "Reticulation",
      "ReticulationDirichlet",
      "Alignment"
    )
  )
global.df$Density <-
  factor(global.df$Density,      # Reordering group factor levels
         levels = c("Low", "High"))

meanValue.df <-
  global.df %>% dplyr::group_by(across(c(-Value))) %>% dplyr::summarise(meanValue = mean(Value, na.rm =
                                                                                           TRUE))#Check why NA

#function for number of observations
give.n <- function(x){
  return(c(y = mean(x)*1.10, label = round(mean(x), digit = 2)))
  # experiment with the multiplier to find the perfect position
}

global.df_rdc <- global.df[!is.na(global.df$Value), ]

global.ggplotLow <-
  ggplot(global.df_rdc %>% filter(Density == "Low"), aes(x = Distribution, Value, fill = Distribution)) +
  facet_grid(Index ~ ., scales = "free")+
  #geom_jitter(aes(color = Distribution),
  #width=0.15, alpha = 0.6)+
  geom_violin(position = position_nudge(x = 0, y = 0)) +
  
  # geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  # geom_jitter(aes(color = Distribution),
  #             width=0.15, alpha = 0.6)+
  
  stat_summary(
    #fun.data = give.n,
    geom = "text",
    fun.y = "mean",
    size = 3,
    vjust = -0.5,
    aes(label = round(after_stat(y), 2))
  ) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    size = 2,
    shape = 19
  ) +
  labs(x = "DISTRIBUTION", y = "VALUE") +
  ggtitle("Low") +
  #scale_y_continuous(sec.axis = sec_axis(~ . , name = "INDEX", breaks = NULL, labels = NULL)) +#, guide=guide_axis(position="left")
  #scale_x_discrete(sec.axis = sec_axis(~ . , name = "DENSITY", breaks = NULL, labels = NULL)) +
  # scale_fill_viridis(discrete = TRUE) +
  # scale_color_viridis(discrete = TRUE) +
  scale_fill_manual(values = c("white", "grey90", "grey80")) +
  
  theme(
    legend.position = "none",
    title = element_text(face = "bold", hjust = 0.5),
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold", colour = "white"),
    strip.text.y = element_text(face = "bold", colour = "white"),
    strip.background = element_rect(colour = "white", fill = "white")
  ) 

global.ggplotLow

global.ggplotHigh <-
  ggplot(global.df_rdc %>% filter(Density == "High"), aes(x = Distribution, Value, fill = Distribution)) +
  facet_grid(Index ~ ., scales = "free")+
  #geom_jitter(aes(color = Distribution),
  #width=0.15, alpha = 0.6)+
  geom_violin(position = position_nudge(x = 0, y = 0)) +
  
  # geom_flat_violin(position = position_nudge(x = .2, y = 0)) +
  # geom_jitter(aes(color = Distribution),
  #             width=0.15, alpha = 0.6)+
  
  stat_summary(
    #fun.data = give.n,
    geom = "text",
    fun.y = "mean",
    size = 3,
    vjust = -0.5,
    aes(label = round(after_stat(y), 2))
  ) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    size = 2,
    shape = 19
  ) +
  labs(x = "DISTRIBUTION", y = "") +
  ggtitle("High") +
  #scale_y_continuous(sec.axis = sec_axis(~ . , name = "INDEX", breaks = NULL, labels = NULL)) +#, guide=guide_axis(position="left")
  #scale_x_discrete(sec.axis = sec_axis(~ . , name = "DENSITY", breaks = NULL, labels = NULL)) +
  # scale_fill_viridis(discrete = TRUE) +
  # scale_color_viridis(discrete = TRUE) +
  scale_fill_manual(values = c("white", "grey90", "grey80")) +
  
  theme(
    legend.position = "none",
    title = element_text(face = "bold", hjust = 0.5),
    #Grid
    panel.grid.major = element_line(colour = "gray93", size = 0.65),
    panel.grid.minor = element_line(colour = "gray97"),
    #Facet wrap
    strip.text.x = element_text(face = "bold", colour = "white"),
    strip.text.y = element_text(face = "bold", colour = "black"),
    strip.background = element_rect(colour = "white", fill = "white")
  ) 

global.ggplotHigh
#Problem is that they do all have same y scale: differences are not seen anymore... leaving it out more the moment

library(ggpubr)
global.ggplot <- ggarrange(
  global.ggplotLow,
  global.ggplotHigh,
  ncol = 2
)
global.ggplot
save.image(
  "Renvironment/resultsTestIndex2.RData"
)

#Ok so only patchiness and alignment to use for best identification.

##~~~~~~~~~~~
## END SCRIPT
##~~~~~~~~~~~