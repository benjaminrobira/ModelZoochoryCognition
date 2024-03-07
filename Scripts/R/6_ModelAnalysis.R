###~~~~~~~~~~~~~~
## Analysing simulation outputs
###~~~~~~~~~~~~~~

# Setting up environments -------------------------------------------------

rm(list = ls())

## Parameters -------------------------------------------------------------

#Load parameters
load("Scripts/R/Parameterisation.RData")

percentForDifference = 5

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

Rcpp::sourceCpp("Scripts/Rcpp/FunctionsRcpp.cpp")
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

### Convergence simulation --------------------------------------------------------

#Function to assess whether the simulation converged on agent's efficiency and resource distribution pattern (YES or NO)
convergenceSimulation <-
  function(output,
           alpha = 0.05,
           numberInterval,
           intervalCompare,
           percentToAssumeDifference = 5) {
    minInterval <- 0
    maxInterval <- max(output[, 1])
    intervalSize <- (maxInterval - minInterval) / numberInterval
    
    firstInterval <-
      c(
        minInterval + (intervalCompare[1] - 1) * intervalSize,
        minInterval + (intervalCompare[1]) * intervalSize
      )
    secondInterval <-
      c(
        minInterval + (intervalCompare[2] - 1) * intervalSize,
        minInterval + (intervalCompare[2]) * intervalSize
      )
    
    compare1.m <-
      output[output[, 1] > firstInterval[1] &
               output[, 1] <= firstInterval[2], ]
    compare2.m <-
      output[output[, 1] > secondInterval[1] &
               output[, 1] <= secondInterval[2], ]
    
    testDifferenceMeanEfficiency <-
      t.test(compare1.m[, 2], compare2.m[, 2])
    
    outputVector <- NA
    outputVector <-
      ifelse(
        testDifferenceMeanEfficiency$p.value[1] < alpha &
          abs(
            testDifferenceMeanEfficiency$estimate[1] - testDifferenceMeanEfficiency$estimate[2]
          ) > percentToAssumeDifference / 100 * min(
            testDifferenceMeanEfficiency$estimate[1],
            testDifferenceMeanEfficiency$estimate[2]
          ),
        "NO",
        "YES"
      )
    
    return(outputVector)
  }

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

### Load complementary data -------------------------------------------------

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

#Alignment
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
  
# Data extraction: consequence of cognition on resource ---------------------------------------------------------

listFiles_v <- list.files("Output/Main")
listFiles_v <- listFiles_v[grep("Main", listFiles_v)]

whichMapFinal_v <- grep("Map", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]

#order for being aligned with routine which will only be taken for the last
filesMapFinalStart_v <- filesMapFinal_v[grep("36500", filesMapFinal_v)]
filesMapFinalEnd_v <- filesMapFinal_v[!(filesMapFinal_v %in% filesMapFinalStart_v)]
filesMapFinal_v <- c(filesMapFinalStart_v, filesMapFinalEnd_v)

whichRoutine_v <- grep("Routine", listFiles_v)
filesRoutine_v <- listFiles_v[whichRoutine_v]

length(filesMapFinal_v)
length(filesRoutine_v)

# whichContinuous_v <- grep("Main_Continuous", listFiles_v)
# filesContinuous_v <- listFiles_v[whichContinuous_v]

library(qdapRegex)
library(doParallel)
library(parallel)

indices_l <- mclapply(
  1:length(filesMapFinal_v),
  mc.cores = 5,
  function(whatFile){
    system(as.character(whatFile))
    tryCatch(
    {
    fileOfInterest <- filesMapFinal_v[whatFile]
    
    time <- qdapRegex::ex_between(fileOfInterest, 
                                            "Map_", 
                                            ".txt", 
                                            extract = TRUE)
    time <- as.numeric(time[[1]][1])
    if(time == 36500){
      lastMap = TRUE
    }else{
      lastMap = FALSE
    }
    
    cognitionLevel <- qdapRegex::ex_between(fileOfInterest, 
                                            "_s", 
                                            "_", 
                                            extract = TRUE)
    cognitionLevel <- as.numeric(cognitionLevel[[1]][1])
    
    outputMapFinal <- read_table2(paste0("Output/Main/", fileOfInterest)) %>% as.data.frame() 
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
    shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/(mapSize*mapSize)
    
    ### Routine ---------------------------------------------------------------
    if(lastMap){
      seqVisits <- read_csv(paste0("Output/Main/",  filesRoutine_v[whatFile]), 
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
      
      return(c(cognitionLevel, time, patchiness, alignment, spatialAutocorr, routine, shrinkage))
    }else{
      return(c(cognitionLevel, time, patchiness, alignment, spatialAutocorr, NA, shrinkage))
    }
    },
    error = function(e){return(NA)})
  }
)

indices_l <- indices_l[sapply(indices_l, length) > 1]

indicesMain_df <- do.call("rbind", indices_l) %>% as.data.frame()
colnames(indicesMain_df) <- c("knowledgeRate", "time", "patchiness", "alignment", "spatialAutocorr", "routine", "shrinkage")
indicesMain_df$knowledgeRate <- as.numeric(indicesMain_df$knowledgeRate)
indicesMain_df$time <- as.numeric(indicesMain_df$time)
indicesMain_df$patchiness <- as.numeric(indicesMain_df$patchiness)
indicesMain_df$alignment <- as.numeric(indicesMain_df$alignment)
indicesMain_df$spatialAutocorr <- as.numeric(indicesMain_df$spatialAutocorr)
indicesMain_df$routine <- as.numeric(indicesMain_df$routine)

save.image("Renvironment/outputTestCognition.RData")

# Data extraction: efficiency ---------------------------------------------

listFiles_v <- list.files("Output/Main")
whichOMNISCIENT_v <- grep("TestEfficiencyOMNISCIENT", listFiles_v)
filesOMNISCIENT_v <- listFiles_v[whichOMNISCIENT_v]
routineOMNISCIENT_v <- filesOMNISCIENT_v[grep("Routine", filesOMNISCIENT_v)]
whichOMNISCIENT_v <- grep("EfficiencyContinuous", filesOMNISCIENT_v)
filesOMNISCIENT_v <- filesOMNISCIENT_v[whichOMNISCIENT_v]

whichINTERMEDIATE_v <- grep("TestEfficiencyINTERMEDIATE", listFiles_v)
filesINTERMEDIATE_v <- listFiles_v[whichINTERMEDIATE_v]
routineINTERMEDIATE_v <- filesINTERMEDIATE_v[grep("Routine", filesINTERMEDIATE_v)]
whichINTERMEDIATE_v <- grep("EfficiencyContinuous", filesINTERMEDIATE_v)
filesINTERMEDIATE_v <- filesINTERMEDIATE_v[whichINTERMEDIATE_v]

whichNULL_v <- grep("TestEfficiencyNULL", listFiles_v)
filesNULL_v <- listFiles_v[whichNULL_v]
routineNULL_v <- filesNULL_v[grep("Routine", filesNULL_v)]
whichNULL_v <- grep("EfficiencyContinuous", filesNULL_v)
filesNULL_v <- filesNULL_v[whichNULL_v]

length(filesOMNISCIENT_v)
length(filesINTERMEDIATE_v)
length(filesNULL_v)

indices_l <- mclapply(
  1:length(filesOMNISCIENT_v),
  mc.cores = 5,
  function(whatFile){
    
    #Determining what initial conditions
    initialCondition <- qdapRegex::ex_between(filesNULL_v[whatFile], 
                                            "NULL", 
                                            "_", 
                                            extract = TRUE)
    initialCondition <- initialCondition[[1]][1]
    
    ## Checking efficiency convergence -----------------------------------------

    outputContinuousNULL <- read_table2(paste0("Output/Main/",  filesNULL_v[whatFile]))
    convergenceNULL <-
      convergenceSimulation(
        outputContinuousNULL,
        numberInterval = 5,
        intervalCompare = c(4, 5),
        percentToAssumeDifference = percentForDifference
      )
    
    outputContinuousINTERMEDIATE <- read_table2(paste0("Output/Main/",  filesINTERMEDIATE_v[whatFile]))
    convergenceINTERMEDIATE <-
      convergenceSimulation(
        outputContinuousINTERMEDIATE,
        numberInterval = 5,
        intervalCompare = c(4, 5),
        percentToAssumeDifference = percentForDifference
      )
    
    outputContinuousOMNISCIENT <- read_table2(paste0("Output/Main/",  filesOMNISCIENT_v[whatFile]))
    convergenceOMNISCIENT <-
      convergenceSimulation(
        outputContinuousOMNISCIENT,
        numberInterval = 5,
        intervalCompare = c(4, 5),
        percentToAssumeDifference = percentForDifference
      )    
    
    ## Extracting final efficiency ---------------------------------------------
  
    efficiencyNULL = outputContinuousNULL[nrow(outputContinuousINTERMEDIATE), 2]
    efficiencyINTERMEDIATE = outputContinuousINTERMEDIATE[nrow(outputContinuousINTERMEDIATE), 2]
    efficiencyOMNISCIENT = outputContinuousOMNISCIENT[nrow(outputContinuousINTERMEDIATE), 2]

    ## Routine -----------------------------------------------------------------

    for(r in 1:3){
      
      if(r == 1){
        fileRoutine = routineNULL_v[whatFile]
      }else if(r == 2){
        fileRoutine = routineINTERMEDIATE_v[whatFile]
      }else{
        fileRoutine = routineOMNISCIENT_v[whatFile]
      }
      
      seqVisits <- read_csv(paste0("Output/Main/",  fileRoutine), 
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
      
      if(r == 1){
        routineNULL = routine
      }else if(r == 2){
        routineINTERMEDIATE = routine
      }else{
        routineOMNISCIENT = routine
      }
      
    }

    output <- 
      data.frame(
        initCondition = rep(initialCondition, 3),
        type = c("Null", "Intermediate", "Omniscient"),
        simulationNb = rep(whatFile, 3),
        convergence = c(convergenceNULL, convergenceINTERMEDIATE, convergenceOMNISCIENT),
        efficiency = as.numeric(c(efficiencyNULL, efficiencyINTERMEDIATE, efficiencyOMNISCIENT)),
        routine = c(routineNULL, routineINTERMEDIATE, routineOMNISCIENT)
      )
    
    return(output)
  }
)

indicesEfficiency_df <- do.call("rbind", indices_l) %>% as.data.frame()

length(indicesEfficiency_df$convergence[indicesEfficiency_df$convergence == "YES"])/nrow(indicesEfficiency_df)
#ok
# Plots -------------------------------------------------------------------

## Final conditions -------------------------------------------------------

plotShrinkage <- plotResults(
  yAxisName = "Shrinkage",
  xAxisName = "Spatiotemporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "shrinkage",
  df = indicesMain_df %>% filter(time == 36500),
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain_df$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  valueSegments = c(0.5),
  lineTypeSegments = c(NA),
  colourSegments =  c(NA),
  nameSegments = c("")
)

plotPatchiness <- plotResults(
  yAxisName = "Patchiness",
  xAxisName = "Spatiotemporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "patchiness",
  df = indicesMain_df %>% filter(time == 36500) %>% mutate(patchiness = patchiness*(1-shrinkage)),#Normalise patchiness by shrinkage
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain_df$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  valueSegments = c(patchinessHomogeneousValue, patchinessHeterogeneousValue, patchinessRouteValue),
  lineTypeSegments = c("dotted", "dashed", "solid"),
  colourSegments =  c("black", "black", "black"),
  nameSegments = c("Homogeneous", "Heterogeneous", "Route")
)

plotAlignment <- plotResults(
  yAxisName = "Alignment",
  xAxisName = "Spatiotemporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "alignment",
  df = indicesMain_df %>% filter(time == 36500),
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain_df$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  valueSegments = c(alignmentHomogeneousValue, alignmentRouteValue),
  lineTypeSegments = c("dotted", "solid"),
  colourSegments =  c("black", "black"),
  nameSegments = c(" Homogeneous", "Route")
)

plotRoutine <- plotResults(
  yAxisName = "Routine",
  xAxisName = "Spatiotemporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "routine",
  df = indicesMain_df %>% filter(time == 36500),
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain_df$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  valueSegments = c(0.5),
  lineTypeSegments = c(NA),
  colourSegments =  c(NA),
  nameSegments = c("")
)

plotSpatAutocorr <- plotResults(
  yAxisName = "Spatial autocorrelation",
  xAxisName = "Spatiotemporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "spatialAutocorr",
  df = indicesMain_df %>% filter(time == 36500),
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain_df$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  valueSegments = c(0, moranHomoIntermediate, moranHomoHigh, moranHeteroHigh),
  lineTypeSegments = c("longdash", "dotted", "solid", "dashed"),
  colourSegments =  c("black", "black", "black", "black"),
  nameSegments = c("", "Int. synchro. homo", "High. synchro. homo.", 
                   "High. synchro. hetero.")
)

library(ggpubr)
mergedPlot <- ggarrange(
  plotPatchiness + rremove("xlab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.8,1.6)),
  plotAlignment + rremove("xlab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.1,0.8)),
  plotSpatAutocorr + rremove("xlab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.05,0.05)),
  plotRoutine + rremove("xlab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.5,0.9)),
  nrow = 4, 
  ncol = 1
)
mergedPlotMain <- annotate_figure(mergedPlot,
                                        top = text_grob("Spatiotemporal knowledge rate", face = "bold", size = 16))
mergedPlotMain

## Map along time ---------------------------------------------------------

# plotResults(
#   yAxisName = "Patchiness",
#   xAxisName = "Time",
#   xVar = "time",
#   yVar = "patchiness",
#   df = indicesMain_df,
#   categoricalX = TRUE,
#   levelsOldNameX = unique(sort(indicesMain_df$time)),
#   levelsNewNameX = as.character(unique(sort(indicesMain_df$time))),
#   groupVar = "knowledgeRate",
#   colourGroup_v = c("white", "grey95", "grey80", "grey65", "grey50"),
#   nameGroup = "Knowledge rate",
#   levelOrderGroup = as.character(unique(indicesMain_df$knowledgeRate))
# )
# 
# plotResults(
#   yAxisName = "Alignment",
#   xAxisName = "Time",
#   xVar = "time",
#   yVar = "alignment",
#   df = indicesMain_df,
#   categoricalX = TRUE,
#   levelsOldNameX = unique(sort(indicesMain_df$time)),
#   levelsNewNameX = as.character(unique(sort(indicesMain_df$time))),
#   groupVar = "knowledgeRate",
#   colourGroup_v = c("white", "grey95", "grey80", "grey65", "grey50"),
#   nameGroup = "Knowledge rate",
#   levelOrderGroup = as.character(unique(indicesMain_df$knowledgeRate))
# )

## Efficiency -------------------------------------------------------------

plotEfficiency <- plotResults(
  yAxisName = "Efficiency (x1000)",
  xAxisName = "Initial resource condition",
  xVar = "initCondition",
  yVar = "efficiency",
  df = indicesEfficiency_df %>% filter(initCondition %in% c(1,3,5)) %>% mutate(efficiency = efficiency*1000) %>% mutate(type = ifelse(type == "Null", "Naive", type)),
  categoricalX = TRUE,
  levelsOldNameX = c("1", "3", "5"),
  levelsNewNameX = c("Naive", "Intermediate", "Omniscient"),
  groupVar = "type",
  colourGroup_v = c("white", "grey80", "black"),
  nameGroup = "Forager type",
  levelOrderGroup = c("Naive", "Intermediate", "Omniscient"),
  widthBox = 0.5,
  differentMeanShapePoints = TRUE
)

plotRoutine2 <- plotResults(
  yAxisName = "Routine",
  xAxisName = "Initial resource condition",
  xVar = "initCondition",
  yVar = "routine",
  df = indicesEfficiency_df %>% filter(initCondition %in% c(1,3,5)) %>% mutate(type = ifelse(type == "Null", "Naive", type)),
  categoricalX = TRUE,
  levelsOldNameX = c("1", "3", "5"),
  levelsNewNameX = c("Naive", "Intermediate", "Omniscient"),
  groupVar = "type",
  colourGroup_v = c("white", "grey80", "black"),
  nameGroup = "Forager type",
  levelOrderGroup = c("Naive", "Intermediate", "Omniscient"),
  widthBox = 0.5,
  differentMeanShapePoints = TRUE
)

# Save --------------------------------------------------------------------

save.image("Renvironment/mainResults.RData")
saveRDS(indicesMain_df, "Renvironment/mainResults_df.rds")

ggsave(plot = plotPatchiness + theme(axis.title = element_text(face = "bold", size = 25),
                                             axis.text.x = element_text(size = 20),
                                             axis.text.y = element_text(size = 20),
                                             legend.text = element_text(size = 10),
                                             legend.title = element_text(size = 12, face = "bold"),
                                             panel.grid.minor = element_line(colour = "grey93"),
                                             panel.grid.major = element_line(colour = "grey93"),
                                             strip.background = element_rect(colour = "white",
                                                                             fill = "white"),
                                             strip.text = element_text(face = "bold", size = 14))
       , filename = "FIG/plotPatchiness.pdf", width = 7, height = 7, dpi = 400)
saveRDS(mergedPlotMain, "Renvironment/Plots/mainPlots.rds")
saveRDS(plotShrinkage, "Renvironment/Plots/mainShrinkage.rds")
saveRDS(plotEfficiency, "Renvironment/Plots/efficiencyPlots.rds")
saveRDS(plotRoutine2, "Renvironment/Plots/routine2Plots.rds")

# Garbage -----------------------------------------------------------------

# 
# vectorHowManyConvergedEfficiency <-
#   rep(0, times = length(temporalKnowledge))
# vectorHowManyConvergedPatchiness <-
#   rep(0, times = length(temporalKnowledge))
# matrixPatchinessAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixAlignmentAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixmoranTestAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixmoranTestAtEndNonCircular <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixgearyTestAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixfractalTestAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixreticulationTestAtEnd <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# matrixDifferenceINTERMEDIATENULL <-
#   matrix(0, nrow = numberRepetitions, ncol = length(temporalKnowledge))
# 
# percentForDifference = 10
# patchinessAtStart = 0
# 
# for (knowledgeLoop in 1:length(spatialKnowledge)) {
#   for (repetitionNumber in 1:(numberRepetitions)) {
#     # #Efficiency
#     # plot(outputContinuous[,1], outputContinuous[,2], type="l")
#     # #NN
#     # plot(outputContinuous[,1], outputContinuous[,3], type="l")
#     # plot(outputContinuous[,1], outputContinuous[,4], type="l")
#     # #Lloyd patchiness
#     # plot(outputContinuous[,1], outputContinuous[,6], type="l")
#     
#     # Analysing convergence ---------------------------------------------------
#     
#     #To analyse convergence, I split the whole evolution into 5 parts and I compared the third and the fifth. I compared the mean (and sd for the second) of mean efficiency and distance to nearest neighbour.
#     #Both should not differ to indicate convergence. Comparison are simply done with ttest.
#     
#     ##IMPORT DATA WITH Dispersal
#     path <- paste0(
#       "Output/Main/Main",
#       "_p15.811388",
#       "_s",
#       format(spatialKnowledge[knowledgeLoop], nsmall = 2),
#       "0000_t",
#       format(temporalKnowledge[knowledgeLoop], nsmall = 2),
#       "0000_r",
#       repetitionNumber
#     )
#     
#     outputContinuous <-
#       read_table2(paste0(path, "_EfficiencyContinuous.txt"))
#     outputContinuous <- as.data.frame(outputContinuous)
#     outputContinuous <- outputContinuous[-nrow(outputContinuous),]
#     
#     patchinessAtStart = patchinessAtStart + outputContinuous[1, 6]
#     
#     tryCatch({
#       outputMapStart <- read_table2(paste0(path, "_Map_0.txt"))
#       outputMapFinal <- read_table2(paste0(path, "_Map_5475.txt"))
#     }, error = function(e) {
#       outputMapStart <- NA
#       outputMapFinal <- NA
#     })
#     
#     ##IMPORT DATA WITHOUT Dispersal
#     pathNULL <- paste0(
#       "Output/NoDispersalNULL",
#       knowledgeLoop,
#       "_p15.811388",
#       "_s",
#       "0.00",
#       "0000_t",
#       "0.00",
#       "0000_r",
#       repetitionNumber
#     )
#     
#     pathINTERMEDIATE <- paste0(
#       "Output/NoDispersalINTERMEDIATE",
#       knowledgeLoop,
#       "_p15.811388",
#       "_s",
#       "0.50",
#       "0000_t",
#       "0.50",
#       "0000_r",
#       repetitionNumber
#     )
#     
#     pathOMNISCIENT <- paste0(
#       "Output/NoDispersalOMNISCIENT",
#       knowledgeLoop,
#       "_p15.811388",
#       "_s",
#       "1.00",
#       "0000_t",
#       "1.00",
#       "0000_r",
#       repetitionNumber
#     )
#     
#     tryCatch({
#       outputContinuousNULL <-
#         read_table2(paste0(pathNULL, "_EfficiencyContinuous.txt"))
#       outputContinuousNULL <- as.data.frame(outputContinuousNULL)
#       
#       outputContinuousINTERMEDIATE <-
#         read_table2(paste0(pathINTERMEDIATE, "_EfficiencyContinuous.txt"))
#       outputContinuousINTERMEDIATE <-
#         as.data.frame(outputContinuousINTERMEDIATE)
#       
#       outputContinuousOMNISCIENT <-
#         read_table2(paste0(pathOMNISCIENT, "_EfficiencyContinuous.txt"))
#       outputContinuousOMNISCIENT <-
#         as.data.frame(outputContinuousOMNISCIENT)
#       
#       matrixDifferenceINTERMEDIATENULL[repetitionNumber, knowledgeLoop] = outputContinuousINTERMEDIATE[nrow(outputContinuousINTERMEDIATE), 2] - outputContinuousNULL[nrow(outputContinuousNULL), 2]
#       
#       matrixNULL[repetitionNumber, knowledgeLoop] = outputContinuousNULL[nrow(outputContinuousINTERMEDIATE), 2]
#       matrixINTERMEDIATE[repetitionNumber, knowledgeLoop] = outputContinuousINTERMEDIATE[nrow(outputContinuousINTERMEDIATE), 2]
#       matrixOMNISCIENT[repetitionNumber, knowledgeLoop] = outputContinuousOMNISCIENT[nrow(outputContinuousINTERMEDIATE), 2]
#       
#       # if(outputContinuousINTERMEDIATE[nrow(outputContinuousINTERMEDIATE),2] - outputContinuousNULL[nrow(outputContinuousNULL),2] < 0){
#       #   print(knowledgeLoop)
#       #   print(repetitionNumber)
#       #   break
#       # }
#     }, error = function(e) {
#       matrixDifferenceINTERMEDIATENULL[repetitionNumber, knowledgeLoop] = NA
#     })
#     
#     resultConvergence <-
#       convergenceSimulation(
#         outputContinuous,
#         numberInterval = 5,
#         intervalCompare = c(4, 5),
#         percentToAssumeDifference = percentForDifference
#       )
#     
#     #Save results of convergence
#     vectorHowManyConvergedEfficiency[knowledgeLoop] <-
#       vectorHowManyConvergedEfficiency[knowledgeLoop] + ifelse(resultConvergence[1] == "YES", 1, 0)
#     vectorHowManyConvergedPatchiness[knowledgeLoop] <-
#       vectorHowManyConvergedPatchiness[knowledgeLoop] + ifelse(resultConvergence[5] == "YES", 1, 0)
#     
#     #### Plot result convergence in function of Spatiotemporal knowledge ####
#     
#     # Analysing plant pattern difference --------------------------------------
#     
#     outputMapFinal <- as.data.frame(outputMapFinal)
#     
#     ## Patchiness --------------------------------------------------------------
# 
#     matrixPatchinessAtEnd[repetitionNumber, knowledgeLoop] =   patchiness_f(outputMapFinal,
#                                                                         0,
#                                                                         0,
#                                                                         mapSize,
#                                                                         mapSize,
#                                                                         quadratSize)
#     #outputContinuous[nrow(outputContinuous), 6]
#     
# 
#     ## Alignment ---------------------------------------------------------------
# 
#     matrixAlignmentAtEnd[repetitionNumber, knowledgeLoop] = alignmentPoints_f(
#       outputMapFinal
#     )
#     
#     ## Moran/Geary Index -------------------------------------------------------
# 
#     #I analyse how trees fruiting dates are spatially correlated at the end: I use Moran Index
#     
#     ##### Circular version: see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870747/
#     
#     #Geary Index
#     # outputMapFinal <- as.data.frame(cbind(1:300, 1:300, 1:300)) #Just to test
#     
#     #Create distance matrix
#     treeDistance_m <- as.matrix(dist(outputMapFinal[, 1:2]))
#     #Use the inverse of distance as weight
#     treeDistanceInverse_m <- 1 / (treeDistance_m)
#     #Remove pb with self distance
#     diag(treeDistanceInverse_m) <- 0
#     
#     #Moran Index considering circular data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870747/
#     startFruitingInRadian_v <-
#       as.data.frame(outputMapFinal)[, 3] * 2 * pi / 365
#     
#     #Transform to -pi pi interval
#     startFruitingInRadian_v <-
#       ifelse(
#         startFruitingInRadian_v > pi,
#         startFruitingInRadian_v - 2 * pi,
#         startFruitingInRadian_v
#       )
#     
#     #Sample
#     sampleSize = length(startFruitingInRadian_v)
#     
#     #Distance function (instead of euclidean distance)
#     distanceCircular <- function(x1, x2) {
#       distanceType = atan2(sin(x1 - x2), cos(x1 - x2))
#       distanceLinear = x1 - x2
#       return(distanceLinear)
#     }
#     # autocorrSp(as.data.frame(outputMapFinal)[,3],treeDistanceInverse_m, index="Moran")
#     # Moran.I(as.data.frame(outputMapFinal)[,3],treeDistanceInverse_m)
#     
#     matrixmoranTestAtEnd[repetitionNumber, knowledgeLoop] = autocorrSp(
#       startFruitingInRadian_v,
#       treeDistanceInverse_m,
#       circular = TRUE,
#       index = "Moran"
#     )
#     matrixgearyTestAtEnd[repetitionNumber, knowledgeLoop] = autocorrSp(
#       startFruitingInRadian_v,
#       treeDistanceInverse_m,
#       circular = TRUE,
#       index = "Geary"
#     )
#     
#     # ##### To check: non circular version
#     # library(ape)
#     # #Compute moran index on start date fruiting for last distribution
#     # moranTest <- Moran.I(as.data.frame(outputMapFinal)[,3],treeDistanceInverse_m)
#     # matrixmoranTestAtEndNonCircular[repetitionNumber, knowledgeLoop] = moranTest$observed
#     
#     ## Fractal dimension -------------------------------------------------------
#     
#     # library(gwfa) #Check for help: Applying two fractal methods to characterise the local and global deviations from scale invariance of built patterns throughout mainland France
#     # outputMapFinal <- as.data.frame(outputMapFinal)
#     #
#     # test=gwfa(points=outputMapFinal,q=0,radius=seq(from=5, to=100, length.out=10),
#     #           bandwith=500,sample_size=500,cell_size=500)
#     #
#     # #estimate the fractal dimension on the different radius: see vignette package gwfa
#     # X=cbind(rep(1,length(test@radius)),log2(test@radius))
#     # fit_frac_dim=(do.call(cbind,test[,4:10]))%*%t(solve(t(X)%*%X)%*%t(X))
#     # test$dimfrac=fit_frac_dim[,2]
#     matrixfractalTestAtEnd[repetitionNumber, knowledgeLoop] = fractalDimension_f(
#       outputMapFinal[, c(1, 2)],
#       radius.gwfa = seq(
#         from = 20,
#         to = 60,
#         length.out = 10
#       ),
#       bandwith.gwfa = 200,
#       sample_size.gwfa = 500,
#       cell_size.gwfa = 200
#     )
#     
#     ## Reticulation ------------------------------------------------------------
#     matrixreticulationTestAtEnd[repetitionNumber, knowledgeLoop] = reticulationIndex_f(outputMapFinal[, c(1, 2)], "Euclidean")
#     
#     # Analysing efficiency differences ----------------------------------------
#     
#     #I check how different are the foraging efficiencies of three agent types at start, middle and end of run // TO BE DONE LATER
#     
# 
#     
#   }
# }
# 
# save.image(
#   "T:/Saved_PhD/Model_dispersalSeed/ModelZoochoryCognition/Renvironment/dataSimulation.RData"
# )
# 
# 
# 
# 
# 
# outputMapFinal <-
#   read_table(
#     "Output_old/WithDispersionRandomReplacement_p10.000000_s1.000000_t1.000000_r30_Map_7300.txt"
#   )
# outputMapStart <-
#   read_table(
#     "Output_old/WithDispersionRandomReplacement_p10.000000_s1.000000_t1.000000_r30_Map_0.txt"
#   )
# 
# outputMapFinal <- as.data.frame(outputMapFinal)
# outputMapStart <- as.data.frame(outputMapStart)
# 
# plot(outputMapFinal[, 1],
#      outputMapFinal[, 2],
#      pch = 19,
#      col = "black")
# #col=grey(level = abs(startFruitingInRadian_v)/max(abs(startFruitingInRadian_v)), alpha = rep(0.5, times=length(startFruitingInRadian_v))))
# points(outputMapStart[, 1],
#        outputMapStart[, 2],
#        pch = 1,
#        col = "red")
# #col=grey(level = abs(startFruitingInRadian_v)/max(abs(startFruitingInRadian_v)), alpha = rep(0.5, times=length(startFruitingInRadian_v))))
# 
# 
# points(outputMapFinal[1, 1],
#        outputMapFinal[1, 2],
#        pch = 19,
#        col = "red")
# points(outputMapStart[1, 1],
#        outputMapStart[1, 2],
#        pch = 19,
#        col = "blue")
# 
# vectorHowManyConvergedEfficiency / numberRepetitions
# vectorHowManyConvergedPatchiness / numberRepetitions
# 
# #Plot results
# 
# pdf(file = "Presentation/Graphics/patchinessResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# #Import own function
# source(
#   "T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R",
#   local = knitr::knit_global()
# )
# 
# ###~~~~~~~~~~
# ## PATCHINESS
# ###~~~~~~~~~~
# 
# #transform to 2 cols table for boxplot
# dfPatchiness <- as.data.frame(cbind(
#   as.vector(matrixPatchinessAtEnd),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dfPatchiness) <- c("Value", "Knowledge")
# 
# #Patchiness
# maxy = ceiling(max(dfPatchiness$Value))
# miny = 0
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5), ylim = c(0, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = 0,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dfPatchiness,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(0, 0),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = 0,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dfPatchiness$Loc <- as.numeric(as.factor(dfPatchiness$Knowledge))
# dfPatchiness$Loc <- dfPatchiness$Loc + 0.2
# dfPatchiness$LocJittered <- jitter(dfPatchiness$Loc, factor = 0.5)
# points(
#   dfPatchiness$LocJittered,
#   dfPatchiness$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixPatchinessAtEnd, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixPatchinessAtEnd, 2, mean) - maxy / 30,
#   labels = round(apply(matrixPatchinessAtEnd, 2, mean), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(from = 0, to = maxy, by = maxy / 5),
#   labels = seq(from = 0, to = maxy, by = maxy / 5),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = maxy / 2,
#   text = "Index of patchiness",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# pdf(file = "Presentation/Graphics/patchinessResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5), ylim = c(0, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = 0,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(from = 0, to = maxy, by = maxy / 5),
#   labels = seq(from = 0, to = maxy, by = maxy / 5),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = maxy / 2,
#   text = "Index of patchiness",
#   cex = 2,
#   font = 2
# )
# dev.off()
# 
# 
# ###~~~~~~~~~~
# ## MORAN INDEX
# ###~~~~~~~~~~
# pdf(file = "Presentation/Graphics/moranResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# #Spatial correlation
# 
# #transform to 2 cols table for boxplot
# dfmoranTest <- as.data.frame(cbind(
#   as.vector(matrixmoranTestAtEnd),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dfmoranTest) <- c("Value", "Knowledge")
# 
# #moranTest
# maxy = ceiling(max(dfmoranTest$Value) * 10) / 10
# miny = floor(min(dfmoranTest$Value) * 10) / 10
# 
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dfmoranTest,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(miny, miny),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = miny,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dfmoranTest$Loc <- as.numeric(as.factor(dfmoranTest$Knowledge))
# dfmoranTest$Loc <- dfmoranTest$Loc + 0.2
# dfmoranTest$LocJittered <- jitter(dfmoranTest$Loc, factor = 0.5)
# points(
#   dfmoranTest$LocJittered,
#   dfmoranTest$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixmoranTestAtEnd, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixmoranTestAtEnd, 2, mean) - maxy / 30,
#   labels = round(apply(matrixmoranTestAtEnd, 2, mean), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Moran's Index\n(circular adapted; start date of fruiting)",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# 
# pdf(file = "Presentation/Graphics/moranResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Moran's Index\n(circular adapted; start date of fruiting)",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# 
# ###~~~~~~~~~~
# ## GEARY INDEX
# ###~~~~~~~~~~
# pdf(file = "Presentation/Graphics/gearyResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# 
# #transform to 2 cols table for boxplot
# dfgearyTest <- as.data.frame(cbind(
#   as.vector(matrixgearyTestAtEnd),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dfgearyTest) <- c("Value", "Knowledge")
# 
# #gearyTest
# maxy = ceiling(max(dfgearyTest$Value))
# miny = floor(min(dfgearyTest$Value))
# 
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dfgearyTest,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(miny, miny),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = miny,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dfgearyTest$Loc <- as.numeric(as.factor(dfgearyTest$Knowledge))
# dfgearyTest$Loc <- dfgearyTest$Loc + 0.2
# dfgearyTest$LocJittered <- jitter(dfgearyTest$Loc, factor = 0.5)
# points(
#   dfgearyTest$LocJittered,
#   dfgearyTest$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixgearyTestAtEnd, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixgearyTestAtEnd, 2, mean) - maxy / 30,
#   labels = round(apply(matrixgearyTestAtEnd, 2, mean), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Geary's Index\n(circular adapted; start date of fruiting)",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# pdf(file = "Presentation/Graphics/gearyResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Geary's Index\n(circular adapted; start date of fruiting)",
#   cex = 2,
#   font = 2
# )
# dev.off()
# 
# ###~~~~~~~~~~
# ## FRACTAL
# ###~~~~~~~~~~
# 
# pdf(file = "Presentation/Graphics/fractalResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# 
# #transform to 2 cols table for boxplot
# dffractalTest <- as.data.frame(cbind(
#   as.vector(matrixfractalTestAtEnd),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dffractalTest) <- c("Value", "Knowledge")
# 
# #fractalTest
# maxy = ceiling(max(dffractalTest$Value) * 10) / 10
# miny = floor(min(dffractalTest$Value) * 10) / 10
# 
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dffractalTest,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(miny, miny),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = miny,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dffractalTest$Loc <- as.numeric(as.factor(dffractalTest$Knowledge))
# dffractalTest$Loc <- dffractalTest$Loc + 0.2
# dffractalTest$LocJittered <- jitter(dffractalTest$Loc, factor = 0.5)
# points(
#   dffractalTest$LocJittered,
#   dffractalTest$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixfractalTestAtEnd, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixfractalTestAtEnd, 2, mean) - maxy / 50,
#   labels = round(apply(matrixfractalTestAtEnd, 2, mean), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = " \"Fractal dimension\" ",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# 
# pdf(file = "Presentation/Graphics/fractalResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = " \"Fractal dimension\" ",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# ###~~~~~~~~~~
# ## RETICULATION
# ###~~~~~~~~~~
# 
# pdf(file = "Presentation/Graphics/reticulationResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# 
# #transform to 2 cols table for boxplot
# dfreticulationTest <- as.data.frame(cbind(
#   as.vector(matrixreticulationTestAtEnd),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dfreticulationTest) <- c("Value", "Knowledge")
# 
# #reticulationTest
# maxy = ceiling(max(dfreticulationTest$Value) * 10) / 10
# miny = floor(min(dfreticulationTest$Value) * 10) / 10
# 
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dfreticulationTest,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(miny, miny),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = miny,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dfreticulationTest$Loc <-
#   as.numeric(as.factor(dfreticulationTest$Knowledge))
# dfreticulationTest$Loc <- dfreticulationTest$Loc + 0.2
# dfreticulationTest$LocJittered <-
#   jitter(dfreticulationTest$Loc, factor = 0.5)
# points(
#   dfreticulationTest$LocJittered,
#   dfreticulationTest$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixreticulationTestAtEnd, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixreticulationTestAtEnd, 2, mean) - maxy / 50,
#   labels = round(apply(matrixreticulationTestAtEnd, 2, mean), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Reticulation",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# 
# pdf(file = "Presentation/Graphics/reticulationResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Reticulation",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# ###~~~~~~~~~~-
# #Difference efficiency
# ###~~~~~~~~~~-
# 
# pdf(file = "Presentation/Graphics/efficiencyResult.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# 
# #transform to 2 cols table for boxplot
# dfDifferenceINTERMEDIATENULL <- as.data.frame(cbind(
#   as.vector(matrixDifferenceINTERMEDIATENULL),
#   rep(temporalKnowledge, each = numberRepetitions)
# ))
# colnames(dfDifferenceINTERMEDIATENULL) <- c("Value", "Knowledge")
# 
# #DifferenceINTERMEDIATENULL
# 
# #For readability
# dfDifferenceINTERMEDIATENULLInit <- dfDifferenceINTERMEDIATENULL
# dfDifferenceINTERMEDIATENULL$Value <-
#   dfDifferenceINTERMEDIATENULL$Value * 10000
# matrixDifferenceINTERMEDIATENULL <-
#   matrixDifferenceINTERMEDIATENULL * 10000
# 
# maxy = ceiling(max(dfDifferenceINTERMEDIATENULL$Value))
# miny = floor(min(dfDifferenceINTERMEDIATENULL$Value))
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     Value ~ Knowledge,
#     data = dfDifferenceINTERMEDIATENULL,
#     boxwex = 0.25,
#     xlab = "",
#     ylab = "",
#     cex.lab = 1.2,
#     yaxs = "i",
#     xaxs = "i",
#     las = 1,
#     tcl = -0.25,
#     frame.plot = FALSE,
#     xaxt = "n",
#     yaxt = "n",
#     outline = FALSE,
#     add = TRUE
#   )
# 
# 
# for (i in 1:length(temporalKnowledge)) {
#   #Mask half
#   rect(
#     xleft = c(i),
#     xright = c(i + 0.5),
#     ybottom = c(miny, miny),
#     ytop = c(maxy, maxy),
#     col = "white",
#     border = NA
#   )
#   
#   #Readd grid
#   addGrid(
#     xmin = i,
#     xmax = i + 0.5,
#     xintsmall = 0.1,
#     xintbig = 0.5,
#     ymin = miny,
#     ymax = maxy,
#     yintsmall = (maxy - miny) / 20,
#     yintbig = (maxy - miny) / 5,
#     axisPlot = FALSE
#   )#Background grid
#   
# }
# 
# #Read interquartile range
# segments(
#   x0 = c(1:length(temporalKnowledge)),
#   x1 = c(1:length(temporalKnowledge)),
#   y0 = c(boxPlotSaved$stats[1, 1:length(temporalKnowledge)]),
#   y1 = c(boxPlotSaved$stats[5, 1:length(temporalKnowledge)])
# )
# 
# #Add jitter points + link because paired
# 
# dfDifferenceINTERMEDIATENULL$Loc <-
#   as.numeric(as.factor(dfDifferenceINTERMEDIATENULL$Knowledge))
# dfDifferenceINTERMEDIATENULL$Loc <-
#   dfDifferenceINTERMEDIATENULL$Loc + 0.2
# dfDifferenceINTERMEDIATENULL$LocJittered <-
#   jitter(dfDifferenceINTERMEDIATENULL$Loc, factor = 0.5)
# points(
#   dfDifferenceINTERMEDIATENULL$LocJittered,
#   dfDifferenceINTERMEDIATENULL$Value,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixDifferenceINTERMEDIATENULL, 2, mean),
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.2,
#   y = apply(matrixDifferenceINTERMEDIATENULL, 2, mean) - (maxy - miny) /
#     30,
#   labels = round(apply(
#     matrixDifferenceINTERMEDIATENULL, 2, mean
#   ), digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Difference in efficiency\nIntermediate - Null (x 10000)",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# pdf(file = "Presentation/Graphics/efficiencyResultEmpty.pdf")
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# emptyPlot(xlim = c(0.5, length(temporalKnowledge) + 0.5),
#           ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = length(temporalKnowledge) + 0.5,
#   xintsmall = 0.1,
#   xintbig = 0.5,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# #x-axis
# axis(
#   side = 1,
#   line = 0,
#   at = seq(
#     from = 1,
#     to = length(temporalKnowledge),
#     by = 1
#   ),
#   labels = temporalKnowledge,
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 1,
#   line = 2.5,
#   at = (length(temporalKnowledge) + 0.5 + 0.5) / 2,
#   text = "Spatiotemporal knowledge rate",
#   cex = 2,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   labels = seq(
#     from = miny,
#     to = maxy,
#     by = (maxy - miny) / 5
#   ),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3,
#   at = (maxy - miny) / 2 + miny,
#   text = "Difference in efficiency\nIntermediate - Null (x 10000)",
#   cex = 2,
#   font = 2
# )
# 
# dev.off()
# 
# 
# ###~~~~~~~~~~
# ## Drawing model example
# ###~~~~~~~~~~
# 
# library(Rcpp)
# Rcpp::sourceCpp("Script/Rcpp/FunctionsRcpp.cpp")
# 
# cycleLimitN = 1#50
# cycleLength = 50
# Ntrees = 50
# durationFruiting = 3
# #print(r)
# set.seed(1) #Unique random seed for reproducibility
# dateStart <- runif(Ntrees, 0, 50)
# locTree <-
#   as.matrix(cbind(runif(Ntrees, 0, 1000), runif(Ntrees, 0, 1000)))##NOTE: in fact, since a var of Rcpp is called the same, it modifies also the variable here when it is ran!!
# runSimulationForPlot(
#   cycleLimitNumber = cycleLimitN,
#   repetitionNumber = 1,
#   timeDelayForDispersal = 0.05,
#   torporTime = 1,
#   saveTreeMap = TRUE,
#   samplingMapTime_v = c(0),
#   nameInit = "Output/ContinuousForPlot/With_Dispersal_plot_test",
#   mapSize = 1000,
#   quadratSize = 50,
#   numberTrees = Ntrees,
#   treeLocInit_m = locTree,
#   fruitingTimesInit_m = cbind(dateStart, dateStart + durationFruiting),
#   homogeneousDistribution = TRUE,
#   treeClusterNumber = 0,
#   treeClusterSpread = 0,
#   maximumFoodToYield_v = rep(1, times = Ntrees),
#   cycleLength = cycleLength,
#   fruitingLength = durationFruiting,
#   noReturnTime = durationFruiting / 10,
#   whatValueUnknownTemporal = 0.000,
#   whatRule = "closest",
#   exponentialRate = 0.01,
#   perceptualRange = 10,
#   spatialKnowledgeRate = 1,
#   temporalKnowledgeRate = 1,
#   speed = 1000,
#   DispersalProbability = 0.1
# )
# 
# 
# #Load sequence of movements
# path <- paste0(
#   "Output/ContinuousForPlot/With_Dispersal_plot_test",
#   "_p15.811388",
#   "_s",
#   format(1, nsmall = 2),
#   "0000_t",
#   format(1, nsmall = 2),
#   "0000_r",
#   1
# )
# 
# 
# #Initial map
# mapInit <- read_table2(paste0(path, "_Map_", 0, ".txt"))
# mapInit <- as.data.frame(mapInit)
# 
# outputContinuous <-
#   read_table2(paste0(path, "_EfficiencyContinuous.txt"))
# outputContinuous <- as.data.frame(outputContinuous)
# 
# #Load agent figure
# library(png)
# pacmanIMG <- readPNG("Presentation/Drawing/pacman.png")
# 
# #Plot and save each movement and changes in environment
# for (i in 1:nrow(outputContinuous)) {
#   #SAVE UNDER PDF
#   pdf(file = paste0("Presentation/Graphics/Simulation", i, ".pdf"))
#   
#   #Multiple plots
#   layout(
#     mat = cbind(c(1, 1), c(2, 3)),
#     widths = c(5, 5),
#     heights = c(5, 5)
#   )
#   par(mar = c(1, 1, 0, 0),
#       mgp = c(1.5, 1, 0),
#       xpd = TRUE)
#   
#   ###~~~~~~~~~~
#   ## SPATIAL MAP
#   ###~~~~~~~~~~
#   
#   #Load and plot map
#   map <- read_table2(paste0(path, "_Map_", i - 1, ".txt"))
#   map <- as.data.frame(map)
#   
#   plot(
#     mapInit[, 1],
#     mapInit[, 2],
#     xlim = c(0, 1000),
#     ylim = c(0, 1000),
#     main = "",
#     xlab = "",
#     ylab = "",
#     yaxt = "n",
#     xaxt = "n",
#     pch = 21,
#     cex = 1.2,
#     bg = "black",
#     col = "black",
#     xpd = TRUE,
#     asp = 1,
#     frame.plot = FALSE
#   )
#   
#   rect(
#     xleft = 0,
#     xright = 1000,
#     ybottom = 0,
#     ytop = 1000,
#     lwd = 1.5
#   )
#   
#   #Add previous path
#   if (i > 2) {
#     for (j in 1:(i - 1)) {
#       lines(
#         x = outputContinuous[(j - 1):j, 2],
#         y = outputContinuous[(j - 1):j, 3],
#         lwd = 1.5,
#         lty = 2,
#         col = "lightgray"
#       )
#     }
#   }
#   if (i != 1) {
#     lines(
#       x = outputContinuous[(i - 1):i, 2],
#       y = outputContinuous[(i - 1):i, 3],
#       lwd = 1.5,
#       lty = 2
#     )
#   }
#   
#   #Replot points
#   points(
#     map[, 1],
#     map[, 2],
#     pch = 21,
#     cex = 1.2,
#     bg = "darkolivegreen3",
#     col = "black",
#     xpd = TRUE
#   )
#   
#   #Add agent
#   addImg(pacmanIMG,
#          x = outputContinuous[i, 2] ,
#          y = outputContinuous[i, 3],
#          width = 70)
#   
#   #Time
#   text(
#     x = 100,
#     y = -100,
#     labels = paste0("t = ", round(outputContinuous[i, 1], digit = 2)),
#     font = 2,
#     cex = 1.5,
#     xpd = TRUE
#   )
#   
#   #Legend
#   legend(
#     x = 300,
#     y = 0,
#     pch = c(19, 19),
#     col = c("darkolivegreen3", "black"),
#     cex = 1.5,
#     legend = c("Current patch", "Initial patch"),
#     xpd = TRUE,
#     bty = "n"
#   )
#   
#   ###~~~~~~~~~~-
#   ## TEMPORAL MAP
#   ###~~~~~~~~~~-
#   
#   #Plot the temporal distribution of resources
#   par(mar = c(3, 3, 0, 1),
#       mgp = c(1.5, 1, 0),
#       xpd = TRUE)
#   
#   emptyPlot(
#     xlim = c(outputContinuous[i, 1] - 10, outputContinuous[i, 1] + 10),
#     ylim = c(0, 5)
#   )
#   
#   #Plot triangle of availability for each tree
#   for (m in 1:nrow(map)) {
#     start <- map[m, 3]
#     end <- map[m, 4]
#     #First cycle
#     polygon(
#       x = c(start, start + (end - start) / 2, end, start),
#       y = c(0, 1, 0, 0),
#       col = adjustcolor("darkolivegreen3", alpha.f = 0.2)
#     )
#     #Second cycle
#     start = start + cycleLength
#     end = end + cycleLength
#     polygon(
#       x = c(start, start + (end - start) / 2, end, start),
#       y = c(0, 1, 0, 0),
#       col = adjustcolor("darkolivegreen3", alpha.f = 0.2)
#     )
#   }
#   
#   #Hide on the left triangle part too early
#   rect(
#     xleft = -100,
#     xright = outputContinuous[i, 1] - 10,
#     ybottom = -1,
#     ytop = 1.1,
#     lwd = 2,
#     col = "white",
#     border = NA
#   )
#   #Hide on the left triangle part too late
#   rect(
#     xleft = outputContinuous[i, 1] + 1010,
#     xright = outputContinuous[i, 1] + 5,
#     ybottom = -1,
#     ytop = 1.1,
#     lwd = 2,
#     col = "white",
#     border = NA
#   )
#   
#   
#   #Time
#   points(
#     x = outputContinuous[i, 1],
#     y = 0,
#     pch = 19,
#     cex = 2,
#     xpd = TRUE
#   )
#   lines(
#     x = c(outputContinuous[i, 1], outputContinuous[i, 1]),
#     y = c(0, 2),
#     lwd = 2
#   )
#   
#   #Axis
#   arrows(
#     x0 = outputContinuous[i, 1] - 10,
#     y0 = 0,
#     x1 = outputContinuous[i, 1] - 10,
#     y1 = 2,
#     length = 0.1,
#     angle = 20,
#     xpd = TRUE
#   )
#   arrows(
#     x0 = outputContinuous[i, 1] - 10,
#     y0 = 0,
#     x1 = outputContinuous[i, 1] + 10,
#     y1 = 0,
#     length = 0.1,
#     angle = 20,
#     xpd = TRUE
#   )
#   text(
#     x = outputContinuous[i, 1] - 12,
#     y = 1,
#     labels = "Fruit quantity",
#     xpd = TRUE,
#     srt = 90,
#     font = 2,
#     cex = 1.3
#   )
#   text(
#     x = outputContinuous[i, 1],
#     y = -0.3,
#     labels = "t",
#     xpd = TRUE,
#     font = 2,
#     cex = 1.3
#   )
#   
#   ###~~~~~~~~~~
#   ## EFFICIENCY
#   ###~~~~~~~~~~
#   emptyPlot(x = c(0, max(outputContinuous[, 1])), y = c(-max(outputContinuous[, 4] *
#                                                                1000), max(outputContinuous[, 4] * 1000)))
#   
#   #add grid
#   addGrid(
#     xmin = 0,
#     xmax = max(outputContinuous[, 1]),
#     xintsmall = max(outputContinuous[, 1]) / 20,
#     xintbig = max(outputContinuous[, 1]) / 10,
#     ymin = 0,
#     ymax = max(outputContinuous[, 4] * 1000),
#     yintsmall = max(outputContinuous[, 4] * 1000) / 20,
#     yintbig = max(outputContinuous[, 4] * 1000) / 10,
#     axisPlot = FALSE
#   )#Background grid
#   
#   
#   lines(x = outputContinuous[1:i, 1],
#         y = outputContinuous[1:i, 4] * 1000,
#         lwd = 2)
#   #Axis
#   arrows(
#     x0 = 0,
#     y0 = 0,
#     x1 = 0,
#     y1 = ,
#     max(outputContinuous[, 4] * 1000),
#     length = 0.1,
#     angle = 20,
#     xpd = TRUE
#   )
#   arrows(
#     x0 = 0,
#     y0 = 0,
#     x1 = max(outputContinuous[, 1]),
#     y1 = 0,
#     length = 0.1,
#     angle = 20,
#     xpd = TRUE
#   )
#   text(
#     x = -5,
#     y = max(outputContinuous[, 4] * 1000) / 2,
#     labels = "Efficiency",
#     xpd = TRUE,
#     srt = 90,
#     font = 2,
#     cex = 1.3
#   )
#   text(
#     x = max(outputContinuous[, 1]),
#     y = -max(outputContinuous[, 4] * 1000) / 10,
#     labels = "t",
#     xpd = TRUE,
#     font = 2,
#     cex = 1.3
#   )
#   
#   dev.off()
# }