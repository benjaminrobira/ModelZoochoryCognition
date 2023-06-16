#What does this script do?
#This scripts builds the plots for having illustrative legend of y axis for the main result graph.

library(ggplot2)
library(png)
library(grid)

#Homogeneous
image <- readPNG("FIG/homogeneous.png")
rasterhomogeneous <- rasterGrob(image, interpolate = TRUE)

#Heterogeneous
image <- readPNG("FIG/heterogeneous.png")
rasterheterogeneous <- rasterGrob(image, interpolate = TRUE)

#Route
image <- readPNG("FIG/route.png")
rasterroute <- rasterGrob(image, interpolate = TRUE)

#Routine
image <- readPNG("FIG/routine.png")
rasterroutine <- rasterGrob(image, interpolate = TRUE)

#No routine
image <- readPNG("FIG/noroutine.png")
rasternoroutine <- rasterGrob(image, interpolate = TRUE)

#Sp corr
image <- readPNG("FIG/spcorr.png")
rasterspcorr <- rasterGrob(image, interpolate = TRUE)

#No sp corr
image <- readPNG("FIG/nospcorr.png")
rasternospcorr <- rasterGrob(image, interpolate = TRUE)

#Patchiness legend
patchinessLegend <- ggplot() +
  geom_segment(mapping = aes(x = 0, xend = 0, y = 0.1, yend = 1),
               arrow = arrow(length = unit(0.03, "npc"), ends="last"),
               color = "black",
               size = 0.7) +
  annotation_custom(rasterheterogeneous, xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1) + 
  annotation_custom(rasterhomogeneous, xmin = 0, xmax = 0.2, ymin = 0.1, ymax = 0.3) + 
  xlim(c(0, 0.4)) + 
  ylim(c(0, 1)) + 
  theme_void()

#Alignment legend
alignmentLegend <- ggplot() +
  geom_segment(mapping = aes(x = 0, xend = 0, y = 0.1, yend = 1),
               arrow = arrow(length = unit(0.03, "npc"), ends="last"),
               color = "black",
               size = 0.7) +
  annotation_custom(rasterroute, xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1) + 
  annotation_custom(rasterhomogeneous, xmin = 0, xmax = 0.2, ymin = 0.1, ymax = 0.3) + 
  xlim(c(0, 0.4)) + 
  ylim(c(0, 1)) + 
  theme_void()

#Sp corr legend
spcorrLegend <- ggplot() +
  geom_segment(mapping = aes(x = 0, xend = 0, y = 0.1, yend = 1),
               arrow = arrow(length = unit(0.03, "npc"), ends="last"),
               color = "black",
               size = 0.7) +
  annotation_custom(rasterspcorr, xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1) + 
  annotation_custom(rasternospcorr, xmin = 0, xmax = 0.2, ymin = 0.1, ymax = 0.3) + 
  xlim(c(0, 0.4)) + 
  ylim(c(0, 1)) + 
  theme_void()

#Routine legend
routineLegend <- ggplot() +
  geom_segment(mapping = aes(x = 0, xend = 0, y = 0.1, yend = 1),
               arrow = arrow(length = unit(0.03, "npc"), ends="last"),
               color = "black",
               size = 0.7) +
  annotation_custom(rasterroutine, xmin = 0, xmax = 0.2, ymin = 0.8, ymax = 1) + 
  annotation_custom(rasternoroutine, xmin = 0, xmax = 0.2, ymin = 0.1, ymax = 0.3) + 
  xlim(c(0, 0.4)) + 
  ylim(c(0, 1)) + 
  theme_void()

#Merge plots
library(ggpubr)
plotLegend <- ggarrange(
  patchinessLegend, 
  alignmentLegend, 
  spcorrLegend, 
  routineLegend, 
  nrow = 4, 
  ncol = 1
)
saveRDS(plotLegend, "FIG/plotLegend.rds")
