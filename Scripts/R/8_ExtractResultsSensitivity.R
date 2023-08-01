###~~~~~~~~~~~~~~
## Analysing simulation outputs sensitivity
###~~~~~~~~~~~~~~

# Setting up environments -------------------------------------------------

rm(list = ls())

## Parameters -------------------------------------------------------------

#Load parameters
load("Scripts/R/Parameterisation.RData")
#If needed to replot and not running the lapply
#load("Renvironment/outputSensitivityMovingRuleAndSpaceTree.RData")
## Packages ----------------------------------------------------------------

library(readr)
library(progress)
library(dplyr)
library(parallel)
library(doParallel)
library(adehabitatHR)
library(sf)

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
    levelOrderGroup = NA
){
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
  
  if(boxplotOnly){
    if(!is.na(groupVar)){
      plot <- ggplot(df, aes(x = x, y = y)) +
        geom_boxplot(aes(x = x, y = y, group = interaction(x, group), fill = group), width = 0.25, notch = TRUE, position = dodge) +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.65,
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
        geom_boxplot(aes(x = x, y = y, group = x), width = 0.25, notch = TRUE, fill = "grey90") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.65,
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
        geom_violin(aes(x = x, y = y, fill = group), position = dodge) +
        geom_boxplot(aes(x = x, y = y, group =  interaction(x, group)), position = dodge, width = 0.05, fill = "black") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.65,
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
        geom_violin(aes(x = x, y = y, group = x), fill = "grey90", adjust = 1) +
        geom_boxplot(aes(x = x, y = y, group = x), width = 0.01, fill = "black") +
        stat_summary(
          #fun.data = give.n,
          geom = "text",
          fun.y = "mean",
          size = 4,
          hjust = -0.65,
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

# Results extraction ------------------------------------------------------

## Moving rule -------------------------------------------------------------

listFiles_v <- list.files("Output/SensitivityMovingRule")
whichlistFiles_v <- grep("Test", listFiles_v)
listFiles_v <- listFiles_v[whichlistFiles_v]
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichRoutine_v <- grep("Routine", listFiles_v)
filesRoutine_v <- listFiles_v[whichRoutine_v]

length(filesMapFinal_v)

indices_l <- mclapply(
  1:length(filesMapFinal_v),
  mc.cores = 5,
  function(whatFile){
    #print(whatFile)
    system(as.character(whatFile))
    fileOfInterest <- filesMapFinal_v[whatFile]
    if(grepl("MovingNot", fileOfInterest)){
      whatRule <- "MovingAllTrees"
    }else if(grepl("MovingTarget", fileOfInterest)){
      whatRule <- "MovingTarget"
    }else{
      whatRule <- "MovingFruitTrees"
    }
    
    outputMapFinal <- read_table2(paste0("Output/SensitivityMovingRule/", fileOfInterest)) %>% as.data.frame() 
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
    
    seqVisits <- read_csv(paste0("Output/SensitivityMovingRule/",  filesRoutine_v[whatFile]), 
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

indicesMovingRule_df <- do.call("rbind", indices_l) %>% as.data.frame()
colnames(indicesMovingRule_df) <- c("movingRule", "patchiness", "alignment", "spatialAutocorr", "routine", "shrinkage")
indicesMovingRule_df$patchiness <- as.numeric(indicesMovingRule_df$patchiness)
indicesMovingRule_df$shrinkage <- as.numeric(indicesMovingRule_df$shrinkage)
indicesMovingRule_df$alignment <- as.numeric(indicesMovingRule_df$alignment)
indicesMovingRule_df$spatialAutocorr <- as.numeric(indicesMovingRule_df$spatialAutocorr)
indicesMovingRule_df$routine <- as.numeric(indicesMovingRule_df$routine)

## Quick plot --------------------------------------------------------------------

library(ggplot2)
library(scales)

plotPatchinessMoving <- plotResults(
  yAxisName = "Patchiness",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "patchiness",
  df = indicesMovingRule_df %>% mutate(patchiness = patchiness*(1-shrinkage)),#Normalise patchiness by shrinkage
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMovingRule_df$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)]
)

plotAlignmentMoving <- plotResults(
  yAxisName = "Alignment",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "alignment",
  df = indicesMovingRule_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMovingRule_df$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)]
)

plotRoutineMoving <- plotResults(
  yAxisName = "Routine",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "routine",
  df = indicesMovingRule_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMovingRule_df$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)]
)

plotSpatAutocorrMoving <- plotResults(
  yAxisName = "Spatial autocorrelation",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "spatialAutocorr",
  df = indicesMovingRule_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMovingRule_df$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)]
)

plotShrinkageMovingRule <- plotResults(
  yAxisName = "Shrinkage",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "shrinkage",
  df = indicesMovingRule_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMovingRule_df$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)]
)

library(ggpubr)
mergedPlot <- ggarrange(
  plotPatchinessMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.8,1.3)),
  plotAlignmentMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.1,0.8)),
  plotSpatAutocorrMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.5,0.5)),
  plotRoutineMoving + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.5,0.9)),
  nrow = 4, 
  ncol = 1
)
mergedPlotMovingRule <- annotate_figure(mergedPlot,
                top = text_grob("Moving rule", face = "bold", size = 16))
mergedPlotMovingRule

# Spacing tree ------------------------------------------------------------

listFiles_v <- list.files("Output/SensitivitySpaceTree")
whichlistFiles_v <- grep("Test", listFiles_v)
listFiles_v <- listFiles_v[whichlistFiles_v]
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichRoutine_v <- grep("Routine", listFiles_v)
filesRoutine_v <- listFiles_v[whichRoutine_v]

length(filesMapFinal_v)

spacingValue_v <- rep(c(0.05, 0.45, 0.85), times = numberRepetitions)
simuID_v <- paste0("r", 1:(3*numberRepetitions), "_")

indices_l <- mclapply(
  1:length(filesMapFinal_v),
  mc.cores = 5,
  function(whatFile){
    system(as.character(whatFile))
    fileOfInterest <- filesMapFinal_v[whatFile]

    valueSpacing <- spacingValue_v[
      which(sapply(simuID_v, function(x){grepl(x, fileOfInterest)}))[1]
    ]
    
    outputMapFinal <- read_table2(paste0("Output/SensitivitySpaceTree/", fileOfInterest)) %>% as.data.frame() 
    plot(outputMapFinal[,1], outputMapFinal[,2])
    
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
          shrinkage <- 1 - getverticeshr(UD, percent = 95, unin = "m", unout = "m2")@data$area/(mapSize*mapSize)
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
    
    seqVisits <- read_csv(paste0("Output/SensitivitySpaceTree/",  filesRoutine_v[whatFile]), 
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
    
    return(c(valueSpacing, patchiness, alignment, spatialAutocorr, routine, shrinkage))
  }
)
indicesSpacing_df <- do.call("rbind", indices_l) %>% as.data.frame()
colnames(indicesSpacing_df) <- c("valueSpacing", "patchiness", "alignment", "spatialAutocorr", "routine", "shrinkage")

## Quick plot --------------------------------------------------------------------
indicesSpacing_df$shrinkage <- as.numeric(indicesSpacing_df$shrinkage)

plotPatchinessSpacing <- plotResults(
  yAxisName = "Patchiness",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "patchiness",
  df = indicesSpacing_df %>% mutate(patchiness = patchiness*(1-shrinkage)),#Normalise patchiness by shrinkage
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing_df$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85")
)

plotAlignmentSpacing <- plotResults(
  yAxisName = "Alignment",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "alignment",
  df = indicesSpacing_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing_df$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85")
)

plotRoutineSpacing <- plotResults(
  yAxisName = "Routine",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "routine",
  df = indicesSpacing_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing_df$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85")
)

plotSpatAutocorrSpacing <- plotResults(
  yAxisName = "Spatial autocorrelation",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "spatialAutocorr",
  df = indicesSpacing_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing_df$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85")
)

plotShrinkageSpaceTree <- plotResults(
  yAxisName = "Shrinkage",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "shrinkage",
  df = indicesSpacing_df,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing_df$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85")
)

library(ggpubr)
mergedPlot <- ggarrange(
  plotPatchinessSpacing + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.8,2.6)),
  plotAlignmentSpacing + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.1,0.8)),
  plotSpatAutocorrSpacing + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(-0.5,0.5)),
  plotRoutineSpacing  + rremove("xlab") + rremove("ylab") + scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0.5,0.9)),
  nrow = 4, 
  ncol = 1
)
mergedPlotSpaceTree <- annotate_figure(mergedPlot,
                                        top = text_grob("Tree spacing intensity", face = "bold", size = 16))
mergedPlotSpaceTree

# Save output -------------------------------------------------------------

save.image("Renvironment/outputSensitivityMovingRuleAndSpaceTree.RData")
saveRDS(mergedPlotMovingRule, "Renvironment/Plots/movingRulePlots.rds")
saveRDS(mergedPlotSpaceTree, "Renvironment/Plots/spacingTreePlots.rds")
saveRDS(plotShrinkageMovingRule, "Renvironment/Plots/movingRuleShrinkage.rds")
saveRDS(plotShrinkageSpaceTree, "Renvironment/Plots/spacingTreeShrinkage.rds")

# Garbage -----------------------------------------------------------------

# ## True plot (base) --------------------------------------------------------
# 
# library(png)
# 
# layout(mat = rbind(c(1, 2), c(3,4)),
#        widths = c(5, 2),
#        heights = c(5, 5))
# 
# 
# ### Patchiness --------------------------------------------------------------
# 
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# #Import own function
# source("Scripts/R/toolbox.R")
# miny = 1
# maxy <- max(indicesMovingRule_df$patchiness)
# emptyPlot(xlim = c(0.5, 2.5), ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = 2.5,
#   xintsmall = 0.1,
#   xintbig = 1,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     patchiness ~ movingRule,
#     data = indicesMovingRule_df,
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
# for (i in 1:2) {
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
#   x0 = c(1:2),
#   x1 = c(1:2),
#   y0 = c(boxPlotSaved$stats[1, 1:2]),
#   y1 = c(boxPlotSaved$stats[5, 1:2])
# )
# 
# #Add jitter points + link because paired
# 
# indicesMovingRule_df$Loc <- as.numeric(as.factor(indicesMovingRule_df$movingRule))
# indicesMovingRule_df$Loc <- indicesMovingRule_df$Loc + 0.2
# indicesMovingRule_df$LocJittered <- jitter(indicesMovingRule_df$Loc, factor = 0.5)
# points(
#   indicesMovingRule_df$LocJittered,
#   indicesMovingRule_df$patchiness,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:2 + 0.2,
#   y = (indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(patchiness)))$mean,
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.35,
#   y = (indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(patchiness)))$mean,
#   labels = round((indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(patchiness)))$mean, digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# # axis(
# #   side = 1,
# #   line = 1,
# #   at = seq(
# #     from = 0.5,
# #     to = 2.5,
# #     by = 0.5
# #   ),
# #   labels = c("", "Moving to all\nperceived trees", "", "Moving to all\nperceived trees", ""),
# #   tcl = -0.25,
# #   las = 1,
# #   cex.axis = 1.5
# # )
# mtext(
#   side = 1,
#   line = 1,
#   at = c(1, 2),
#   text = c("Moving to all\nperceived trees", "Moving to all\nperceived trees"),
#   cex = 1,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(from = miny, to = maxy, by = (maxy - miny) / 5),
#   labels = round(seq(from = miny, to = maxy, by = (maxy - miny) / 5), digits = 2),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3.5,
#   at = miny + (maxy - miny) / 2,
#   text = "Index of patchiness",
#   cex = 2,
#   font = 2
# )
# 
# par(mar = c(0, 0, 0, 0), mgp = c(2, 0.5, 0))
# emptyPlot()
# patchiness_legendIMG <- readPNG("FIG/Patchiness_legend.png")
# addImg(patchiness_legendIMG,
#        x = 0.5,
#        y = 0.5,
#        width = 0.5)
# 
# ### Alignment ---------------------------------------------------------------
# 
# par(mar = c(4, 7, 0.5, 0.5),
#     mgp = c(3.5, 1, 0),
#     xpd = TRUE)
# #Import own function
# source("Scripts/R/toolbox.R")
# miny = 0.3
# maxy <- max(indicesMovingRule_df$alignment)
# emptyPlot(xlim = c(0.5, 2.5), ylim = c(miny, maxy))
# 
# #add grid
# addGrid(
#   xmin = 0.5,
#   xmax = 2.5,
#   xintsmall = 0.1,
#   xintbig = 1,
#   ymin = miny,
#   ymax = maxy,
#   yintsmall = (maxy - miny) / 20,
#   yintbig = (maxy - miny) / 5,
#   axisPlot = FALSE
# )#Background grid
# 
# boxPlotSaved <-
#   boxplot(
#     alignment ~ movingRule,
#     data = indicesMovingRule_df,
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
# for (i in 1:2) {
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
#   x0 = c(1:2),
#   x1 = c(1:2),
#   y0 = c(boxPlotSaved$stats[1, 1:2]),
#   y1 = c(boxPlotSaved$stats[5, 1:2])
# )
# 
# #Add jitter points + link because paired
# 
# indicesMovingRule_df$Loc <- as.numeric(as.factor(indicesMovingRule_df$movingRule))
# indicesMovingRule_df$Loc <- indicesMovingRule_df$Loc + 0.2
# indicesMovingRule_df$LocJittered <- jitter(indicesMovingRule_df$Loc, factor = 0.5)
# points(
#   indicesMovingRule_df$LocJittered,
#   indicesMovingRule_df$alignment,
#   cex = 0.5,
#   pch = 19,
#   col = "grey",
#   xpd = TRUE
# )
# 
# #Add mean points
# points(
#   x = 1:2 + 0.2,
#   y = (indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(alignment)))$mean,
#   pch = 19,
#   cex = 1.25
# )
# text(
#   x = 1:length(temporalKnowledge) + 0.35,
#   y = (indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(alignment)))$mean,
#   labels = round((indicesMovingRule_df %>% group_by(movingRule) %>% summarise(mean = mean(alignment)))$mean, digit = 2),
#   cex = 1.25
# )
# 
# #x-axis
# # axis(
# #   side = 1,
# #   line = 1,
# #   at = seq(
# #     from = 0.5,
# #     to = 2.5,
# #     by = 0.5
# #   ),
# #   labels = c("", "Moving to all\nperceived trees", "", "Moving to all\nperceived trees", ""),
# #   tcl = -0.25,
# #   las = 1,
# #   cex.axis = 1.5
# # )
# mtext(
#   side = 1,
#   line = 1,
#   at = c(1, 2),
#   text = c("Moving to all\nperceived trees", "Moving to all\nperceived trees"),
#   cex = 1,
#   font = 2
# )
# 
# #y-axis
# axis(
#   side = 2,
#   line = 0,
#   at = seq(from = miny, to = maxy, by = (maxy - miny) / 5),
#   labels = round(seq(from = miny, to = maxy, by = (maxy - miny) / 5), digits = 2),
#   tcl = -0.25,
#   las = 1,
#   cex.axis = 1.5
# )
# mtext(
#   side = 2,
#   line = 3.5,
#   at = miny + (maxy - miny) / 2,
#   text = "Index of alignment",
#   cex = 2,
#   font = 2
# )
# 
# par(mar = c(0, 0, 0, 0), mgp = c(2, 0.5, 0))
# emptyPlot()
# alignment_legendIMG <- readPNG("FIG/alignment_legend.png")
# addImg(alignment_legendIMG,
#        x = 0.5,
#        y = 0.5,
#        width = 0.5)
