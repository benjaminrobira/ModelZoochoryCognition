# What does this script?
# This script adds a complement: it calculates the variance in fruiting date 

rm(list = ls())

#Data
library(tidyr)
library(readr)
library(dplyr)
#Graph
library(ggplot2)
library(ggpubr)
library(scales)
#Parallelisation
library(parallel)
library(doParallel)
#Character processing
library(qdapRegex)

# Function


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
    differentMeanShapePoints = FALSE
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
        geom_boxplot(aes(x = x, y = y, group = x), width = 0.25, notch = TRUE, fill = "black") +
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
        geom_violin(aes(x = x, y = y, fill = group), position = dodge) +
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
        geom_violin(aes(x = x, y = y, group = x), fill = "black", adjust = 1) +
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

# Moving rule -------------------------------------------------------------

listFiles_v <- list.files("Output/SensitivityMovingRule")
whichlistFiles_v <- grep("Test", listFiles_v)
listFiles_v <- listFiles_v[whichlistFiles_v]
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichMapInit_v <- grep("Map_0", listFiles_v)
filesMapInit_v <- listFiles_v[whichMapInit_v]

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
    
    fileOfInterest <- filesMapFinal_v[whatFile]
    outputMapFinal <- read_table2(paste0("Output/SensitivityMovingRule/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapFinal$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthFinal <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    fileOfInterest <- filesMapInit_v[whatFile]
    outputMapInit <- read_table2(paste0("Output/SensitivityMovingRule/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapInit$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthInit <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    return(c(meanVectorLengthInit, meanVectorLengthFinal, whatRule))
  }
)
indicesMoving <- do.call("rbind", indices_l) %>% 
  as.data.frame() %>% 
  mutate(
    experience = "Moving rule"
  ) %>% 
  rename(
    Initial = "V1",
    Final = "V2"
  ) %>% 
  pivot_longer(cols = Initial:Final, values_to = "varianceDate", names_to = "whatMap")
colnames(indicesMoving)[1] <- c("movingRule")

indicesMoving$varianceDate <- as.numeric(indicesMoving$varianceDate)
plotVariationDateMoving <- plotResults(
  yAxisName = "Variation in fruiting date",
  xAxisName = "Moving rule",
  xVar = "movingRule",
  yVar = "varianceDate",
  df = indicesMoving,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMoving$movingRule)[c(2,1,3)],
  levelsNewNameX = c("Only fruit trees", "All trees", "Only target trees")[c(2,1,3)],
  groupVar = "whatMap",
  colourGroup_v = c("white", "black"),
  nameGroup = "Conditions",
  levelOrderGroup = c("Initial", "Final"),
  widthBox = 0.4,
  differentMeanShapePoints = TRUE
)

# Spacing tree -----------------------------------------------------
numberRepetitions = 200
spacingValue_v <- rep(c(0.05, 0.45, 0.85), times = numberRepetitions)
simuID_v <- paste0("r", 1:(3*numberRepetitions), "_")

listFiles_v <- list.files("Output/SensitivitySpaceTree")
whichlistFiles_v <- grep("Test", listFiles_v)
listFiles_v <- listFiles_v[whichlistFiles_v]
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichMapInit_v <- grep("Map_0", listFiles_v)
filesMapInit_v <- listFiles_v[whichMapInit_v]

indices_l <- mclapply(
  1:length(filesMapFinal_v),
  mc.cores = 5,
  function(whatFile){
    #print(whatFile)
    system(as.character(whatFile))
    fileOfInterest <- filesMapFinal_v[whatFile]
    valueSpacing <- spacingValue_v[
      which(sapply(simuID_v, function(x){grepl(x, fileOfInterest)}))[1]
    ]
    
    fileOfInterest <- filesMapFinal_v[whatFile]
    outputMapFinal <- read_table2(paste0("Output/SensitivitySpaceTree/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapFinal$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthFinal <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    fileOfInterest <- filesMapInit_v[whatFile]
    outputMapInit <- read_table2(paste0("Output/SensitivitySpaceTree/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapInit$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthInit <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    return(c(meanVectorLengthInit, meanVectorLengthFinal, valueSpacing))
  }
)
indicesSpacing <- do.call("rbind", indices_l) %>% 
  as.data.frame() %>% 
  mutate(
    experience = "Spacing tree intensity"
  ) %>% 
  rename(
    Initial = "V1",
    Final = "V2"
  ) %>% 
  pivot_longer(cols = Initial:Final, values_to = "varianceDate", names_to = "whatMap")
colnames(indicesSpacing)[1] <- c("valueSpacing")


plotVariationDateSpacing <- plotResults(
  yAxisName = "Variation in fruiting date",
  xAxisName = "Tree spacing intensity",
  xVar = "valueSpacing",
  yVar = "varianceDate",
  df = indicesSpacing,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesSpacing$valueSpacing),
  levelsNewNameX = c("0.05", "0.45", "0.85"),
  groupVar = "whatMap",
  colourGroup_v = c("white", "black"),
  nameGroup = "Conditions",
  levelOrderGroup = c("Initial", "Final"),
  widthBox = 0.4,
  differentMeanShapePoints = TRUE
)


# Spatio-temporal knowledge -----------------------------------------------
listFiles_v <- list.files("Output/Main")
listFiles_v <- listFiles_v[grep("Main", listFiles_v)]
whichMapFinal_v <- grep("Map_36500", listFiles_v)
filesMapFinal_v <- listFiles_v[whichMapFinal_v]
whichMapInit_v <- grep("Map_0", listFiles_v)
filesMapInit_v <- listFiles_v[whichMapInit_v]

indices_l <- mclapply(
  1:length(filesMapFinal_v),
  mc.cores = 5,
  function(whatFile){
    #print(whatFile)
    system(as.character(whatFile))
    fileOfInterest <- filesMapFinal_v[whatFile]
    cognitionLevel <- qdapRegex::ex_between(fileOfInterest, 
                                            "_s", 
                                            "_", 
                                            extract = TRUE)
    cognitionLevel <- as.numeric(cognitionLevel[[1]][1])
    
    
    fileOfInterest <- filesMapFinal_v[whatFile]
    outputMapFinal <- read_table2(paste0("Output/Main/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapFinal$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthFinal <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    fileOfInterest <- filesMapInit_v[whatFile]
    outputMapInit <- read_table2(paste0("Output/Main/", fileOfInterest)) %>% as.data.frame() 
    dateAngle_v <- (outputMapInit$startFruit * (2*pi) / 365) %% (2*pi)
    #Taken from trajr
    meanVectorLengthInit <- Mod(mean(complex(modulus = 1, argument = dateAngle_v)))
    
    return(c(meanVectorLengthInit, meanVectorLengthFinal, cognitionLevel))
  }
)
indicesMain <- do.call("rbind", indices_l) %>% 
  as.data.frame() %>% 
  mutate(
    experience = "Spatio-temporal knowledge rate"
  ) %>% 
  rename(
    Initial = "V1",
    Final = "V2"
  ) %>% 
  pivot_longer(cols = Initial:Final, values_to = "varianceDate", names_to = "whatMap")
colnames(indicesMain)[1] <- c("knowledgeRate")
  
plotVariationDateMain <- plotResults(
  yAxisName = "Variation in fruiting date",
  xAxisName = "Spatio-temporal knowledge rate",
  xVar = "knowledgeRate",
  yVar = "varianceDate",
  df = indicesMain,
  categoricalX = TRUE,
  levelsOldNameX = unique(indicesMain$knowledgeRate),
  levelsNewNameX = as.character(c(0, 0.25, 0.5, 0.75, 1)),
  groupVar = "whatMap",
  colourGroup_v = c("white", "black"),
  nameGroup = "Conditions",
  levelOrderGroup = c("Initial", "Final"),
  widthBox = 0.4,
  differentMeanShapePoints = TRUE
)

mergedPlot <- ggarrange(
  plotVariationDateMain +  
    theme(legend.text = element_text(size = 18), 
          legend.title = element_text(size = 22)) +
    scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0, 0.3)),
  plotVariationDateMoving + rremove("ylab") + 
    theme(legend.text = element_text(size = 18), 
          legend.title = element_text(size = 22)) +
    scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0, 0.3)),
  plotVariationDateSpacing + rremove("ylab") + 
    theme(legend.text = element_text(size = 18), 
          legend.title = element_text(size = 22)) +
    scale_y_continuous(breaks = extended_breaks(n = 4), minor_breaks = extended_breaks(n = 6*4), limits = c(0, 0.3)),
  nrow = 1, 
  ncol = 3,
  common.legend = TRUE,
  legend = "right"
)
mergedPlot

saveRDS(mergedPlot , "Renvironment/Plots/plotVariationFruitingDates.rds")
