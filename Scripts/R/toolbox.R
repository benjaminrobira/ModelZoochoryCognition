## ~~~~~~~
# TOOLBOX
## ~~~~~~~

## What is it?

#This toolbox gathers many handmade functions (of my own, or burrowed, and potentially modified from the Internet!) that were created because repeatedly used. 
#Note that I have plenty of others linked to movement analyses specifically (in case you need).
#How each function work is described, please enjoy!

## OH NO! :( AN ERROR!

#Found a mistake or an error? I am sad (really, this is my highest fright), and I truly hope this is no big deal. 
#Please contact me at benjamin.robira@normalesup.org so we can solve it together!

## ~~~~~~~
#Function list:

# ~~~~~
###Data processing
# ~~~~~

#firstup
## Capitalizes the first letter of a string.

#as.numcharac
## This function allows you to combine as.numeric and as.character.

#'%nin%'
## This function allows you to select what is not included in a given variable.

#getmode
## This function calculates the mode (the value the most represented).

#cbind.fill
## This function allows you to bind columns of different sizes.

#replaceNAtable
## This function allows you to change all the NA in a df with a defined symbol.

#round_preserve_sum
## This function allows you to round a vector that sums to one, while preserving the sum to one.

#getIDMatrixFromCoords
## This function allows you to extract the vector coordinates of a matrix.

# ~~~~~
###P-value related
# ~~~~~
#pvalueToText
## This function will transform the p-value into the classical symbols "*" "**" "***" "n.s.". Default is 0.05, 0.01, 0.001. > 0.05

#pvalueRound 
## This function will return as character p=value of the pvalue. In case it's below the treshold for very significant, in order to avoid very small pvalue it will be displayed as "p<threshold". 

#textTestOneVar
## This function will allow you to directly write (e.g. for Rmarkdown documents) the statistics, df and p-value of a linear model full vs null comparison.

#textEstOneVar
## This function will allow you to directly write (e.g. for Rmarkdown documents) the statistics and p-value for a predictor test.

#roundIntelligent
## To round in a smart way (with scientific writing or not): the rounding will occur for the value before the power of ten in the scientific writing

# ~~~~~
###Linear models
# ~~~~~

#reportGamSmooth
## Report the result of smoothed terms for GAMs.

#mainResults
## Extract the results of a (generalized) linear model.

#stabResults 
## Extract the stability of a (generalized) linear model.

#diagnostics.plot.dharma
## Plot the different diagnostics for a (generalised) linear model. Does not work for some (e.g. beta models).

# ~~~~~
###Graphic related
# ~~~~~

# The graphic related functions were written to be used, for most, with the base R package. Later, I finally switched to ggplot, and those cannot be used in this case.

#addImg
## This function is used for optimising display of jpeg or png images directly into the graphs. it was taken from the internet 

#pastellize -> seems not to work well... 
## Creates a pastel version of a given colour.

#addGrid
## This function add a background grid.

#addCorner 
## To add segment corners to the map.

#errorBars
## To add the sd/error/CI bars to a given plot.

#addLabel
## To add a circled labels to number the different panels of a multi-panel plot.

#emptyPlot
## To create an empty base plot.

#patternLayer
## Fill a given polygon with a given pattern.

#patternLayerCoordinates
## Fill a given SPATIAL polygon with a given pattern.

#legend2 (from : https://bitbucket.org/dominikfroehlich/legend2/src/master/legend2/legend2.r)
## patched version of R plot legend function to allow for variable spacings in horizontal legend

# scientific_10
## This function will allow you to transform labels in elegant scientific writing: e.g: 5x10^2 instead of 5e10^2

# ~~~~~
### R markdown
# ~~~~~

# R mardown related functions are here to help display results in a R markdown file, or add functionalities (e.g. word count etc.).

# rangePrint
## This function allows you to print range between [,]

# numberToLetter
## This function allows you to transform number <= 10 to letters. 

# RmdWords
## Calculates number of words in a Rmarkdown document. Seems rather inaccurate.

# citeR
## Merge two bib files: one existing (referenced literature), and on created for R packages used. 

# extractRversion
## To extract the current version of R used for analysis. 

## ~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## ~~~~~~~
# Data processing --------------------------------------------------------------------
## ~~~~~~~

## ~~~~~~~
## firstup --------------------------------------------------------------------
## ~~~~~~~

# Capitalizes the first letter of a string.

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# ~~~~~

## ~~~~~~~
## as.numcharac --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to combine as.numeric and as.character.

as.numcharac <- function(x){
	return(as.numeric(as.character(x)))
}

# ~~~~~

## ~~~~~~~
## %nin% --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to select what is not included in a given variable.

'%nin%' <- Negate('%in%')

# ~~~~~

## ~~~~~~~
## getmode --------------------------------------------------------------------
## ~~~~~~~

# This function calculates the mode (the value the most represented).

# Create the function.
getmode <- function(
v #a numerical vector
) {
   uniqv <- unique(v)
   return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

# ~~~~~

## ~~~~~~~
## Concatenate vector of different size (fill with NA) --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to bind columns of different sizes.

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(NA, n-nrow(x), ncol(x))))) 
}

# ~~~~~

## ~~~~~~~
## replaceNAtable --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to change all the NA in a df with a defined symbol.

replaceNAtable <- function(input, replaceWith){
input <- apply(input, 2, function(x){
  toChange <- x
  toChange[is.na(toChange)] <- replaceWith
  toChange <- gsub(" ", "", toChange)
  return(toChange)
  })
 return(input)
}
# ~~~~~

## ~~~~~~~
## Rounding while preserving summing to 1 --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to round a vector that sums to one, while preserving the sum to one.

#https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum

round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

# ~~~~~

## ~~~~~~~
## Extracting coordinates of matrix when matrix is considered as a vector --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to extract the vector coordinates of a matrix.

getIDMatrixFromCoords <- function(mat, coords1, coords2){
  n.row <- nrow(mat)
  n.col <- ncol(mat)
  ID <- apply(cbind(coords1, coords2), 1, function(x)
  {
  (x[1] - 1) * n.col + x[2]
  }
  )
  return(ID)
}
## ~~~~~

## ~~~~~~~
# P-value related --------------------------------------------------------------------
## ~~~~~~~

## ~~~~~~~
## pvalueToText --------------------------------------------------------------------
## ~~~~~~~

# This function will transform the p-value into the classical symbols "*" "**" "***" "n.s.". Default is 0.05, 0.01, 0.001. > 0.05

pvalueToText <- function(
p, #The pvalue as numeric
thresholdVerySign=0.001,
thresholdIntermediateSign=0.01,
thresholdWeaklySign=0.05
){
	if(p<=thresholdVerySign){
		return("***");
	}
	else if(p<=thresholdIntermediateSign){
		return("**");
	}
	else if (p<=thresholdWeaklySign){
		return("*");
	}
	else{
		return("n.s.");
	}
}

# ~~~~~

## ~~~~~~~
## pvalueRound --------------------------------------------------------------------
## ~~~~~~~

# This function will return as character p=value of the pvalue. In case it's below the treshold for very significant, in order to avoid very small pvalue it will be displayed as "p<threshold". 

pvalueRound <- function(
p,
thresholdVerySign=0.001,
text=TRUE
){
	if(p<=thresholdVerySign){
		if(text){
			return(paste("p < ", thresholdVerySign, sep=""))
		}else{
			return(paste("< ", thresholdVerySign, sep=""))
		}
	}
	else{
		if(text){
			return(paste("p = ", round(p, digit=3), sep=""))
		}else{
			return(round(p, digit=3))
		}
	}
}

# ~~~~~

## ~~~~~~~
## textModelFullNull --------------------------------------------------------------------
## ~~~~~~~

# This function will allow you to directly write (e.g. for Rmarkdown documents) the statistics, df and p-value of a linear model full vs null comparison.

# comparison: the dataframe of the lmtest::lrtest(null, full)
# Note: for incompatibility, the chisq term has to be inserted prior the function... I am working on it!

textModelFullNull <- function(comparison, digit=2){
  output <- paste(
    #as.character("$\chi$"),
    " = ",
    round(comparison[2,4], digit=digit),
    ", df = ",
    comparison[2,3],
    ", ",
    pvalueRound(comparison[2,5]),
    sep=""
  )
  return(output)
}

## ~~~~~~~
## textTestOneVar --------------------------------------------------------------------
## ~~~~~~~

# This function will allow you to directly write (e.g. for Rmarkdown documents) the statistics and p-value for a predictor test.

#It necessitates:
# -coeffLine: a vector containing the statistics and the p-value (note: use that of drop1 for multiple levels!), and the df (for df version)
# -statistics: the "statistics used: z, t, chisq, F, deviance etc... that you want to be displayed

textTestOneVar <- function(coeffLine, digit=2, statistics){
  output <- paste(
    paste(statistics, " = ", sep=""),
	ifelse(abs(coeffLine[1])<1/(10**digit),format(coeffLine[1], digits=digit+1, scientific = T),round(coeffLine[1], digit=digit)),
    ", ",
    pvalueRound(coeffLine[2]),
    sep=""
  )
  return(output)
}

textTestOneVarWithdf <- function(coeffLine, digit=2, statistics){
  output <- paste(
	paste(statistics, " = ", sep=""),
	ifelse(abs(coeffLine[1])<1/(10**digit),format(coeffLine[1], digits=digit+1, scientific = T),round(coeffLine[1], digit=digit)),
    ", ",
	"df = ",
	coeffLine[3],
	", ",
    pvalueRound(coeffLine[2]),
    sep=""
  )
  return(output)
}

## ~~~~~~~
## textEstOneVar --------------------------------------------------------------------
## ~~~~~~~

# This function will allow you to directly write (e.g. for Rmarkdown documents) the statistics and p-value for a predictor test.
#It necessitates:
# -confintLine: a vector containing the estimate, the lower border of the confint, the upper border of the confint in that order.

textEstOneVar <- function(confintLine, digit=2){
  output <- paste(
    "est. = ",
    ifelse(abs(confintLine[1])<1/(10**digit),format(confintLine[1], digits=digit+1, scientific = T),round(confintLine[1], digit=digit)),
    ", CI95% = [",
    ifelse(abs(confintLine[2])<1/(10**digit),format(confintLine[2], digits=digit+1, scientific = T),round(confintLine[2], digit=digit)),
    ",",
    ifelse(abs(confintLine[3])<1/(10**digit),format(confintLine[3], digits=digit+1, scientific = T),round(confintLine[3], digit=digit)),
    "]",
    sep=""
  )
  return(output)
}  

## ~~~~~~~
## roundIntelligent --------------------------------------------------------------------
## ~~~~~~~

# To round in a smart way (with scientific writing or not): the rounding will occur for the value before the power of ten in the scientific writing
roundIntelligent <- function(x, digit=2){
	ifelse(abs(x)<1/(10**digit),format(x, digits=digit+1, scientific = T),round(x, digit=digit))
}

## ~~~~~~~

## ~~~~~~~
# Linear models --------------------------------------------------------------------
## ~~~~~~~

## ~~~~~~~
## reportGamSmooth --------------------------------------------------------------------
## ~~~~~~~

# Report the result of smoothed terms for GAMs.

#The input vector should be the line displaying result of the smoothed terms. This can be accessed if fitting the gam with the gam function of the mgcv package
#using summary(model)$s.table

reportGamSmooth <- function(vectorSmoothTermResults, digit=2){
  paste(
    "edf = ",
    ifelse(abs(vectorSmoothTermResults[1])<1/(10**digit),format(vectorSmoothTermResults[1], digits=digit+1, scientific = T),round(vectorSmoothTermResults[1], digit=digit)),
    ", ",
    "Ref.df = ", 
    vectorSmoothTermResults[2],
    ", ",
    "F = ",
    ifelse(abs(vectorSmoothTermResults[3])<1/(10**digit),format(vectorSmoothTermResults[3], digits=digit+1, scientific = T),round(vectorSmoothTermResults[3], digit=digit)),
    ", ",
    pvalueRound(vectorSmoothTermResults[4]),
    sep=""
  )
}

## ~~~~~~~
## mainResults --------------------------------------------------------------------
## ~~~~~~~

# Extract the results of a (generalized) linear model.

mainResults <- function(
model,
modelName
){

	coeffModel <- as.data.frame(summary(model)$coefficients)
	confintModel <- as.data.frame(confint(model))
	vectorVarNb <- nrow(coeffModel)
	obsNumber <- nrow(model$model)
	family <- family(model)

	if(family=="gaussian"){
	drop1Model <- as.data.frame(drop1(model, test="F"))
	}
	else{
	drop1Model <- as.data.frame(drop1(model, test="Chisq"))
	}

	output <- matrix(NA, ncol=10, nrow=vectorVarNb)
	colnames(output) <- c("Model", "Adj. R²","Variables", "Est.", "Sd", "Lower CI (95%)", "Upper CI (95%)", #"RSS", 
	"Stat.", "Df", "p-value")

	output[,1] <- paste(modelName, " (N=", obsNumber,")", sep="")

	library(rsq)
	output[,2] <- rsq(model)
	output[,3] <- rownames(coeffModel)
	output[,4] <- coeffModel[,1]
	output[,5] <- coeffModel[,2]
	output[,6] <- confintModel[,1]
	output[,7] <- confintModel[,2]

	#Determine which variable is factor

	output[1,8] <- drop1Model[1,5]
	output[1,9] <- drop1Model[1,1]
	output[1,10] <- drop1Model[1,6]

	counterFactor=0
	for(i in 2:ncol(model$model)){
		if(!is.factor(model$model[,i])){
			output[i+counterFactor,8] <- drop1Model[i,5]
			output[i+counterFactor,9] <- drop1Model[i,1]
			output[i+counterFactor,10] <- drop1Model[i,6]
		}
		else{
			for(j in 1:(length(levels(model$model[,i]))-1)){
				output[i+counterFactor+j-1,8] <- drop1Model[i,5]
				output[i+counterFactor+j-1,9] <- drop1Model[i,1]
				output[i+counterFactor+j-1,10] <- drop1Model[i,6]
			}
			counterFactor=counterFactor+length(levels(model$model[,i]))-2
		}
	}
	return(output)
}
## ~~~~~~~

## ~~~~~~~
## stabResults --------------------------------------------------------------------
## ~~~~~~~

# Extract the stability of a (generalized) linear model.

stabResults <- function(
model,
modelName
){
	# vectorVarNb <- nrow(coeffModel)
	obsNumber <- nrow(model$model)
	family <- family(model)

	dfBetasmodel <- round(cbind(coefficients(model), coefficients(model)+
				  t(apply(X=dfbeta(model), MARGIN=2, FUN=range))), 5)

	output <- matrix(NA, ncol=10, nrow=vectorVarNb)
	colnames(output) <- c("Model","N", "VIF", "dffits", "lev", "CD", "Disp.","Variable", "dfbeta (min)", "dfbeta (max)")

	output[,1] <- modelName

	library(rsq)
	output[,2] <- obsNumber
	output[,3] <- max(vif(model))
	output[,4] <- max(dffits(model))
	output[,5] <- max(as.vector(influence(model)$hat))
	output[,6] <- max(cooks.distance(model))
	source("diagnostics.R")
	if(family=="poisson"){output[,7] <- overdisp.test(model)[4]}
	
	output[,8] <-  rownames(summary(model)$coefficients)
	output[,9] <- dfBetasmodel[,1]
	output[,10] <- dfBetasmodel[,2]
	
	return(output)
}
	
## ~~~~~~~~~

## ~~~~~~~
## Hist residuals, QQplot, and Residuals vs fitted values (homogeneity) --------------------------------------------------------------------
## ~~~~~~~

# Plot the different diagnostics for a (generalised) linear model. Does not work for some (e.g. beta models).

diagnostics.plot.dharma <-function(mod.res, col=grey(level=0.25, alpha=0.5), breaks.histo=20, quantreg=TRUE){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 3, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="", breaks=breaks.histo)
  mtext(text="Histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  
  library(DHARMa)
  simulationOutput <- simulateResiduals(fittedModel = mod.res, plot = FALSE)
  
  plotQQunif(simulationOutput) # left plot in plot.DHARMa()
  plotResiduals(simulationOutput, quantreg=quantreg)
  #Old way without dharma, from Roger Mundry
  # qqnorm(residuals(mod.res), main="", pch=19)
  # qqline(residuals(mod.res))
  # mtext(text="qq-plot of residuals", side=3, line=0)
  # plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
  # abline(h=0, lty=2)
  # mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

## ~~~~~~~~~

## ~~~~~~~~~
# Graphic related --------------------------------------------------------------------
## ~~~~~~~~~

## ~~~~~~~~~
## addImg --------------------------------------------------------------------
## ~~~~~~~~~

# This function is used for optimising display of jpeg or png images directly into the graphs. it was taken from the internet 
#https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r

addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate, xpd=TRUE)
}

## ~~~~~~~~~

## ~~~~~~~
## pastellizing colors --------------------------------------------------------------------
## ~~~~~~~

# Creates a pastel version of a given colour.

#pastellizeIntensity should vary between 0 (no pastellization) and 1 (full pastellization)
#from: https://datascienceconfidential.github.io/r/graphics/2017/11/23/soothing-pastel-colours-in-r.html

pastellize <- function(x, pastellizeIntensity=0.5){
  
  library(ColorPalette)
  
  # x is a colour
  # p is a number in [0,1]
  # p = 1 - pastellizeIntensity = 1 will give no pastellization
  
  # convert hex or letter names to rgb
  if (is.character(x)) x <- col2rgb(x)/255
  
  # convert vector to rgb
  if (is.numeric(x)) x <- matrix(x, nr=3)
  
  col <- rgb2hsv(x, maxColorValue=1)
  col[2,1] <- col[2,1]*(1-pastellizeIntensity)
  col <- hsv2rgb(h=col[1,1], s=col[2,1], v=col[3,1])
  
  # return in convenient format for plots
  return(col)
}

# ~~~~~

## ~~~~~~~~~
## addGrid --------------------------------------------------------------------
## ~~~~~~~~~

# This function add a background grid.

addGrid <- function(
xmin, #min of abscisse
xmax, #max of abscisse
xintsmall, #interval on abscisse for the thiner lines
xintbig, #interval on abscisse for the thicker lines
ymin, #min or ordinate
ymax, #max or ordinate
yintsmall, #interval on ordinate for the thiner lines
yintbig, #interval on ordinate for the thicker lines
colsmall="gray97", #colour for the thiner lines
colbig="gray93", #colour for the thicker lines
axisPlot=TRUE, #redraw or not the axes based on the sequence given for the ordinate and abscisse
lty=1,
lwdsmall=1,
lwdbig=1,
round=FALSE,#Rounding labels for axis?
digit=c(0,0),#Rounding number for x then y,
cexAxisX=1,
cexAxisY=1,
contour=FALSE,
colContour="black"
){

segments(x0=xmin, x1=xmax, 
         y0=seq(from=ymin, to=ymax, by=yintsmall), 
         y1=seq(from=ymin, to=ymax, by=yintsmall),
         col=colsmall, lty=lty, lwd=lwdsmall)	

segments(x0=seq(from=xmin, to=xmax, by=xintsmall), x1=seq(from=xmin, to=xmax, by=xintsmall), 
         y0=ymin, 
         y1=ymax,
         col=colsmall, lty=lty, lwd=lwdsmall)
		
segments(x0=xmin, x1=xmax, 
         y0=seq(from=ymin,to=ymax, by=yintbig), 
         y1=seq(from=ymin,to=ymax, by=yintbig),
         col=colbig, lty=lty, lwd=lwdbig)

segments(x0=seq(from=xmin, to=xmax, by=xintbig), x1=seq(from=xmin, to=xmax, by=xintbig), 
         y0=ymin, 
         y1=ymax,
         col=colbig, lty=lty, lwd=lwdbig)
		 
if(axisPlot&round==FALSE){
	axis(side=1, at=seq(from=xmin, to=xmax, by=xintbig), labels=seq(from=xmin, to=xmax, by=xintbig), las=1, tcl=-0.25, cex=cexAxisX)
	axis(side=2, at=seq(from=ymin, to=ymax, by=yintbig), labels=seq(from=ymin, to=ymax, by=yintbig), las=1, tcl=-0.25, cex=cexAxisY)	
} else if(axisPlot&round){
	0
	axis(side=2, at=round(seq(from=ymin, to=ymax, by=yintbig), digit=digit[2]), labels=round(seq(from=ymin, to=ymax, by=yintbig), digit=digit[2]), las=1, tcl=-0.25)	
}

if(contour==TRUE){
#Contour
rect(xleft=xmin, xright=xmax, ybottom=ymin, ytop=ymax, col=NA, border=colContour)
}
}

# ~~~~~

## ~~~~~~~~~
## addCorner --------------------------------------------------------------------
## ~~~~~~~~~

# To add segment corners to the map.

#This function add a background grid

addCorner <- function(
xmin, #min of abscisse
xmax, #max of abscisse
ymin, #min or ordinate
ymax, #max or ordinate
lengthSegment=1/10,#proportion of grid that should be covered by the segment
predefinedLength=NA,#will erase proportion, if not square it's to be used
col="black",#colour of the corners
lty=1,#lty of the corners
lwd=1#lwd of the corner
){

if(is.na(predefinedLength)){
#upper left
segments(
x0=xmin,
x1=xmin + (xmax-xmin)*lengthSegment,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymax - (ymax-ymin)*lengthSegment,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#upperright
segments(
x0=xmax,
x1=xmax - (xmax-xmin)*lengthSegment,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymax - (ymax-ymin)*lengthSegment,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#lowerleft
segments(
x0=xmin,
x1=xmin + (xmax-xmin)*lengthSegment,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + (ymax-ymin)*lengthSegment,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

#lowerright
segments(
x0=xmax,
x1=xmax - (xmax-xmin)*lengthSegment,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymin + (ymax-ymin)*lengthSegment,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)
}
else{
#upper left
segments(
x0=xmin,
x1=xmin + predefinedLength,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymax - predefinedLength,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#upperright
segments(
x0=xmax,
x1=xmax - predefinedLength,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymax - predefinedLength,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#lowerleft
segments(
x0=xmin,
x1=xmin + predefinedLength,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + predefinedLength,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

#lowerright
segments(
x0=xmax,
x1=xmax - predefinedLength,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + predefinedLength,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)
}
}
 
## ~~~~~~~~~
## Create the sd bars --------------------------------------------------------------------
## ~~~~~~~-

# To add the sd/error/CI bars to a given plot.

errorBars <- function(location, meanPt, barValue, refUnit, minValue, maxValue, upperBarValue=NA, lowerBarValue=NA, col="black", lty=1, lwd=1, horiz=FALSE, symmetrical=TRUE, xpd=TRUE){
#Location indicates the loc on the x or y axis (if vert or horiz plot)
#bar value is a unique value, and should indicate the sd/se to ADD to the mean.library(tidyr)
#If not symmetrical, barValue should be set to a random value (won't be used) and you should use upperBarValue and lowerBarValue instead.
#ref unit should indicates the unit to have 10% to draw the little border lines of the bars.
#min and max value should indicate whether there are borders to the variable range (or plot) to adjust the bars
	if(symmetrical){
		if(meanPt - barValue < minValue){
			minPt=minValue
		} else{
			minPt=meanPt - barValue
		}
		if(meanPt + barValue > maxValue){
			maxPt=maxValue
		} else{
			maxPt=meanPt + barValue
		}
		if(horiz){
		  segments(y0=location, y1=location, x0=minPt, x1=maxPt, col=col, lty=lty, lwd=lwd, xpd=TRUE)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=maxPt, x1=maxPt, col=col, lwd=lwd, xpd=TRUE)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=minPt, x1=minPt, col=col, lwd=lwd, xpd=TRUE)
		} else{
		  segments(x0=location, x1=location, y0=minPt, y1=maxPt, col=col, lty=lty, lwd=lwd, xpd=TRUE)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=maxPt, y1=maxPt, col=col, lwd=lwd, xpd=TRUE)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=minPt, y1=minPt, col=col, lwd=lwd, xpd=TRUE)
		}
	} else{
		if(length(lowerBarValue < minValue)==0){
			minPt=lowerBarValue 
		} else{
			lowerBarValue[lowerBarValue < minValue] <- minValue
			minPt=lowerBarValue
		}
		if(length(upperBarValue > maxValue)==0){
			maxPt=upperBarValue
		} else{
			upperBarValue[upperBarValue > maxValue]=maxValue
			maxPt=upperBarValue
		}
		if(horiz){
		  segments(y0=location, y1=location, x0=minPt, x1=maxPt, col=col, lty=lty, lwd=lwd, xpd=TRUE)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=maxPt, x1=maxPt, col=col, lwd=lwd, xpd=TRUE)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=minPt, x1=minPt, col=col, lwd=lwd, xpd=TRUE)
		} else{
		  segments(x0=location, x1=location, y0=minPt, y1=maxPt, col=col, lty=lty, lwd=lwd, xpd=TRUE)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=maxPt, y1=maxPt, col=col, lwd=lwd, xpd=TRUE)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=minPt, y1=minPt, col=col, lwd=lwd, xpd=TRUE)
		}
	}
}

# ~~~~~	
 
## ~~~~~~~~~
## Labelling a panel in a regular manner based on par definition --------------------------------------------------------------------
## ~~~~~~~~~

# To add a circled labels to number the different panels of a multi-panel plot.

#modified from: https://seananderson.ca/2013/10/21/panel-letters/
  
#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.

#circle: plot circle or not
#radiuscircle, size of the radius, as in draw.circle, in x units
#circle.border/bg, colour of border/bg of circle if plotted
#font.col #font colour
#font.size #font size
addLabel <- function(xfrac, yfrac, label, circle = FALSE, radiuscircle=NA, circle.bg=NA, circle.border="black", font.col="black", font.size=1, ...) {
  if(circle){
  	  u <- par("usr")
	  x <- u[1] + xfrac * (u[2] - u[1])
	  y <- u[4] - yfrac * (u[4] - u[3])
	  
	  library(plotrix)
	  draw.circle(x=x, y=y, radius=radiuscircle, border=circle.border, col=circle.bg)
	  text(x, y, label, col=font.col, cex=font.size, ...)

  }else{
	  u <- par("usr")
	  x <- u[1] + xfrac * (u[2] - u[1])
	  y <- u[4] - yfrac * (u[4] - u[3])
	  text(x, y, label, col=font.col, cex=font.size, adj=c(0.5, 0.5),...)
  }
}
  
  
# ~~~~~	
 
## ~~~~~~~~~--

## ~~~~~~~~~
## Create empty plot --------------------------------------------------------------------
## ~~~~~~~~~

# To create an empty base plot.

emptyPlot <- function(xlim=c(0,1), ylim=c(0,1), asp1=FALSE){
	if(asp1==FALSE){
		plot(x=0, y=0, main="", type="n",
		pch=19, xlab="", ylab="", 
		frame.plot=FALSE, xlim=xlim, ylim=ylim, las=1, tcl=-0.25, xaxs="i", yaxs="i",
		xpd=TRUE, yaxt="n", xaxt="n")
	 }else{
	 	plot(x=0, y=0, main="", type="n", asp=1,
		pch=19, xlab="", ylab="", 
		frame.plot=FALSE, xlim=xlim, ylim=ylim, las=1, tcl=-0.25, xaxs="i", yaxs="i",
		xpd=TRUE, yaxt="n", xaxt="n")
	 }
} 

## ~~~~~~~~~--
 
## ~~~~~~~~~
## Create filled pattern for polygons/Necessitates the splancs package --------------------------------------------------------------------
## ~~~~~~~~~

# Fill a given polygon with a given pattern.

# Please note:
# the tmap package (tm_fill, option pattern) and cartography package (patternLayer function) for more advanced patterns (hexagonal or wave for instance !!)

patternLayer <- function(x, y, pch=19, lwd=1, cex=1, col="black", xspacing, yspacing){

	#Initial polygon coordinates
	polygon_sp <- cbind(as.numeric(as.character(x)),as.numeric(as.character(y)))
	#Extract border points
	minX <- min(polygon_sp[,1])
	maxX <- max(polygon_sp[,1])
	minY <- min(polygon_sp[,2])
	maxY <- max(polygon_sp[,2])
	#Create a grid of the points
	library(splancs)
	pointsToPlot_init <- gridpts(polygon_sp,xs=xspacing, ys=yspacing)
	pointsToPlot_spshifted <- pointsToPlot_init
	pointsToPlot_spshifted[,1] <- pointsToPlot_spshifted[,1] + xspacing/2
	pointsToPlot_spshifted[,2] <- pointsToPlot_spshifted[,2] + yspacing/2 
	#Remove points out of polygon
	library(sp)
	isInPolygon <- sp::point.in.polygon(pointsToPlot_spshifted[,1], pointsToPlot_spshifted[,2] , polygon_sp[,1], polygon_sp[,2], mode.checked=FALSE)
	pointsToPlot_spshifted <- pointsToPlot_spshifted[isInPolygon==1,]
	pointsToPlot <- rbind(pointsToPlot_init, pointsToPlot_spshifted)
	
	#Plot
	points(x=pointsToPlot[,1], y=pointsToPlot[,2], pch=pch, cex=cex, lwd=lwd, col=col)
}
			
# ~~~~~			

## ~~~~~~~~~
## Create a spatial dataframe of points for a filled pattern (usable with ggplot with geom_points) --------------------------------------------------------------------
## ~~~~~~~~~

# Fill a given SPATIAL polygon with a given pattern.

#argProjection should meet the requirement for the CRS() function to work properly
patternLayerCoordinates <- function(x, y, xspacing, yspacing, argProjection){
	#Initial polygon coordinates
	polygon_sp <- cbind(as.numeric(as.character(x)),as.numeric(as.character(y)))
	#Extract border points
	minX <- min(polygon_sp[,1])
	maxX <- max(polygon_sp[,1])
	minY <- min(polygon_sp[,2])
	maxY <- max(polygon_sp[,2])
	#Create a grid of the points
	library(splancs)
	pointsToPlot_init <- gridpts(polygon_sp,xs=xspacing, ys=yspacing)
	pointsToPlot_spshifted <- pointsToPlot_init
	pointsToPlot_spshifted[,1] <- pointsToPlot_spshifted[,1] + xspacing/2
	pointsToPlot_spshifted[,2] <- pointsToPlot_spshifted[,2] + yspacing/2 
	#Remove points out of polygon
	library(sp)
	isInPolygon <- sp::point.in.polygon(pointsToPlot_spshifted[,1], pointsToPlot_spshifted[,2] , polygon_sp[,1], polygon_sp[,2], mode.checked=FALSE)
	pointsToPlot_spshifted <- 	pointsToPlot_spshifted[isInPolygon==1,]
	pointsToPlot <- rbind(pointsToPlot_init, pointsToPlot_spshifted)
	
	#Get the spatialPoints dataframe
	pointsToPlot_sp<- SpatialPointsDataFrame(pointsToPlot, data=as.data.frame(seq(from=1, to=nrow(pointsToPlot), by=1)))
	proj4string(pointsToPlot_sp) <-  CRS(argProjection)
	return(pointsToPlot_sp)
}
			
# ~~~~~		


### ~~~~~~~
## Legend2 -------------------------------------------------------------------- 
### ~~~~~~~

# patched version of R plot legend function to allow for variable spacings in horizontal legen

# from https://bitbucket.org/dominikfroehlich/legend2/src/master/legend2/legend2.r
#  Modified file based on src/library/graphics/R/legend.R
#
#  Copyright (C) 1995-2019 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

legend2 <-
function(x, y = NULL, legend, fill = NULL, col = par("col"), border="black",
         lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"),
         box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"),
     pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
     xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
     text.width = NULL, text.col = par("col"), text.font = NULL,
     merge = do.lines && has.pch, trace = FALSE,
     plot = TRUE, ncol = 1, horiz = FALSE, title = NULL,
     inset = 0, xpd, title.col = text.col, title.adj = 0.5,
         seg.len = 2)
{
    ## the 2nd arg may really be `legend'
    if(missing(legend) && !missing(y) &&
       (is.character(y) || is.expression(y))) {
    legend <- y
    y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)

    if(!missing(xpd)) {
        op <- par("xpd")
        on.exit(par(xpd=op))
        par(xpd=xpd)
    }
    title <- as.graphicsAnnot(title)
    if(length(title) > 1) stop("invalid 'title'")
    legend <- as.graphicsAnnot(legend)
    n.leg <- if(is.call(legend)) 1 else length(legend)
    if(n.leg == 0) stop("'legend' is of length 0")
    auto <-
    if (is.character(x))
        match.arg(x, c("bottomright", "bottom", "bottomleft", "left",
               "topleft", "top", "topright", "right", "center"))
    else NA

    if (is.na(auto)) {
    xy <- xy.coords(x, y, setLab = FALSE)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) stop("invalid coordinate lengths")
    } else nx <- 0

    xlog <- par("xlog")
    ylog <- par("ylog")

    rect2 <- function(left, top, dx, dy, density = NULL, angle, ...) {
    r <- left + dx; if(xlog) { left <- 10^left; r <- 10^r }
    b <- top  - dy; if(ylog) {  top <- 10^top;  b <- 10^b }
    rect(left, top, r, b, angle = angle, density = density, ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx; if(xlog) { x1 <- 10^x1; x2 <- 10^x2 }
    y2 <- y1 + dy; if(ylog) { y1 <- 10^y1; y2 <- 10^y2 }
    segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
    if(xlog) x <- 10^x
    if(ylog) y <- 10^y
    points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
    ## ~~~~~- need to adjust  adj == c(xadj, yadj) ?? --
    if(xlog) x <- 10^x
    if(ylog) y <- 10^y
    text(x, y, ...)
    }
    if(trace) {
    catn <- function(...)
        do.call(cat, c(lapply(list(...),formatC), "\n"))
        fv <- function(...)
            paste(vapply(lapply(list(...), formatC),
                         paste, collapse=",", ""),
                  collapse=", ")
    }
    Cex <- cex * par("cex")        # = the `effective' cex for text

    ## at this point we want positive width even for reversed x axis.

    xyc <- xyinch(par("cin"), warn.log=FALSE) # [uses par("usr") and "pin"]
    xc <- Cex * xyc[1L]
    yc <- Cex * xyc[2L]

    if(is.null(text.width)){
#    text.width <- max(abs(strwidth(legend, units="user",
#                       cex=cex, font = text.font)))
        if(horiz){
            text.width <- c(0, abs(strwidth(legend[1:length(legend)-1], units="user",
                       cex=cex, font = text.font))) + (x.intersp + 1) * max(xc)
            for(i in 2:length(text.width)){
                text.width[i] <- text.width[i] + text.width[i-1]  # sum of the vector is x position
            }
        }else{
            text.width <- max(abs(strwidth(legend, units="user",
                       cex=cex, font = text.font))) + (x.intersp + 1) * max(xc)
        }
    }else if(!all(is.numeric(text.width)) || all(text.width < 0)){
        stop("'text.width' must be numeric, >= 0")
    }else{
        if(horiz){
            text.width <- c(0, text.width[1:length(text.width)-1])
            for(i in 2:length(text.width)){
                text.width[i] <- text.width[i] + text.width[i-1]  # sum of the vector is x position
            }
        }else{
            text.width <- max(text.width)
        }
    }


    if(horiz & length(text.width) == length(legend)){
        secondvector <- (1:length(legend))-1
        text.width <- (text.width / secondvector)
        text.width[1] <- 0
    }

    if(any(n_ <- xc < 0)) text.width[n_] <- -text.width[n_]

    xchar  <- xc
    xextra <- 0
    yextra <- yc * (y.intersp - 1)
    ## watch out for reversed axis here: heights can be negative
    ymax   <- yc * max(1, strheight(legend, units="user", cex=cex)/yc)
    ychar <- yextra + ymax
    if(trace) catn("  xchar=", fv(xchar), "; (yextra, ychar)=", fv(yextra,ychar))

    if(mfill) {
    ##= sizes of filled boxes.
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- max(xbox) ## + x.intersp*xchar
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 0))
         ) || !missing(lwd)

    ## number of ("rbinded") legends _per_ column:
    n.legpercol <-
    if(horiz) {
        if(ncol != 1)
                warning(gettextf("horizontal specification overrides: Number of columns := %d",
                                 n.leg), domain = NA)
        ncol <- n.leg
        1
    } else ceiling(n.leg / ncol)

    has.pch <- !missing(pch) && length(pch) > 0 # -> default 'merge' is available
    if(do.lines) {
    x.off <- if(merge) -0.7 else 0
    } else if(merge)
    warning("'merge = TRUE' has no effect when no line segments are drawn")

    if(has.pch) {
    if(is.character(pch) && !is.na(pch[1L]) &&
           nchar(pch[1L], type = "c") > 1) {
        if(length(pch) > 1)
        warning("not using pch[2..] since pch[1L] has multiple chars")
        np <- nchar(pch[1L], type = "c")
        pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
        ## this coercion was documented but not done in R < 3.0.0
        if(!is.character(pch)) pch <- as.integer(pch)
    }

    if (is.na(auto)) {
    ##- Adjust (x,y) :
    if (xlog) x <- log10(x)
    if (ylog) y <- log10(y)
    }
    if(nx == 2) {
    ## (x,y) are specifiying OPPOSITE corners of the box
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top  <- y[2L]
    w <- diff(x)# width
    h <- diff(y)# height
    w0 <- w/ncol # column width

    x <- mean(x)
    y <- mean(y)
    if(missing(xjust)) xjust <- 0.5
    if(missing(yjust)) yjust <- 0.5

    }
    else {## nx == 1  or  auto
    ## -- (w,h) := (width,height) of the box to draw -- computed in steps
    h <- (n.legpercol + !is.null(title)) * ychar + yc
    xch1 <- max(xchar)
    w0 <- text.width + (x.intersp + 1) * xch1
    if(mfill)    w0 <- w0 + dx.fill
    if(do.lines)    w0 <- w0 + (seg.len + x.off)*xch1
    if(horiz){
    w <- sum(w0) + .5* xch1
    }else{
        w <- ncol*w0 + .5* xch1
    }
    if (!is.null(title)
        && (abs(tw <- strwidth(title, units="user", cex=cex) + 0.5*xchar)) > abs(w)) {
        xextra <- (tw - w)/2
        w <- tw
    }

    ## ~~~~~ (w,h) are now the final box width/height.

    if (is.na(auto)) {
        left <- x - xjust * w
        top     <- y + (1 - yjust) * h
    } else {
        usr <- par("usr")
        inset <- rep_len(inset, 2)
        insetx <- inset[1L]*(usr[2L] - usr[1L])
        left <- switch(auto, "bottomright" =,
               "topright" =, "right" = usr[2L] - w - insetx,
               "bottomleft" =, "left" =, "topleft" = usr[1L] + insetx,
               "bottom" =, "top" =, "center" = (usr[1L] + usr[2L] - w)/2)
        insety <- inset[2L]*(usr[4L] - usr[3L])
        top <- switch(auto, "bottomright" =,
              "bottom" =, "bottomleft" = usr[3L] + h + insety,
              "topleft" =, "top" =, "topright" = usr[4L] - insety,
              "left" =, "right" =, "center" = (usr[3L] + usr[4L] + h)/2)
    }
    }

    if (plot && bty != "n") { ## The legend box :
    if(trace)
        catn("  rect2(", left, ",", top,", w=", w, ", h=", h, ", ...)",
                 sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
              lwd = box.lwd, lty = box.lty, border = box.col)
    }

    ## (xt[],yt[]) := `current' vectors of (x/y) legend text
    xt <- left + xchar + xextra +
      (w0 * rep.int(0:(ncol-1), rep.int(n.legpercol,ncol)))[1L:n.leg]  # n.leg <- length(leged) or 1
    yt <- top -    0.5 * yextra - ymax -
      (rep.int(1L:n.legpercol,ncol)[1L:n.leg] - 1 + !is.null(title)) * ychar

    if (mfill) {        #- draw filled boxes -------------
    if(plot) {
        if(!is.null(fill)) fill <- rep_len(fill, n.leg)
        rect2(left = xt, top=yt+ybox/2, dx = xbox, dy = ybox,
          col = fill,
          density = density, angle = angle, border = border)
    }
    xt <- xt + dx.fill
    }
    if(plot && (has.pch || do.lines))
    col <- rep_len(col, n.leg)

    ## NULL is not documented but people use it.
    if(missing(lwd) || is.null(lwd))
    lwd <- par("lwd") # = default for pt.lwd
    if (do.lines) {            #- draw lines ---------------------
        ## NULL is not documented
    if(missing(lty) || is.null(lty)) lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) & !is.na(lwd)
    if(trace)
        catn("  segments2(",xt[ok.l] + x.off*xchar, ",", yt[ok.l],
         ", dx=", seg.len*xchar, ", dy=0, ...)")
    if(plot)
        segments2(xt[ok.l] + x.off*xchar, yt[ok.l],
                      dx = seg.len*xchar, dy = 0,
              lty = lty[ok.l], lwd = lwd[ok.l], col = col[ok.l])
    # if (!merge)
    xt <- xt + (seg.len+x.off) * xchar
    }
    if (has.pch) {            #- draw points -------------------
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
        ok <- !is.na(pch)
        if (!is.character(pch)) {
            ## R 2.x.y omitted pch < 0
            ok <- ok & (pch >= 0 | pch <= -32)
        } else {
            ## like points
            ok <- ok & nzchar(pch)
        }
    x1 <- (if(merge && do.lines) xt-(seg.len/2)*xchar else xt)[ok]
    y1 <- yt[ok]
    if(trace)
        catn("  points2(", x1,",", y1,", pch=", pch[ok],", ...)")
    if(plot)
        points2(x1, y1, pch = pch[ok], col = col[ok],
            cex = pt.cex[ok], bg = pt.bg[ok], lwd = pt.lwd[ok])
##D    if (!merge) xt <- xt + dx.pch
    }

    xt <- xt + x.intersp * xchar
    if(plot) {
    if (!is.null(title))
            text2(left + w*title.adj, top - ymax, labels = title,
                  adj = c(title.adj, 0), cex = cex, col = title.col)

    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
    }
    invisible(list(rect = list(w = w, h = h, left = left, top = top),
           text = list(x = xt, y = yt)))
}


## ~~~~~~~


## ~~~~~~~
## Elegant scientific writing -------------------------------------------------------------------- 
## ~~~~~~~

# This function will allow you to transform labels in elegant scientific writing: e.g: 5x10^2 instead of 5e10^2
# From: https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales

scientific_10 <- function(x) {
  parse(text=gsub("e", "%*%10^", scales::scientific_format()(x)))
}

#Use it with for instance with ggplot: scale_y_continuous(label = scientific_10)

## ~~~~~~~


## ~~~~~~~
# R markdown related --------------------------------------------------------------------
## ~~~~~~~


## ~~~~~~~
## rangePrint -------------------------------------------------------------------- 
## ~~~~~~~

# This function allows you to print range between [,]

rangePrint <- function(x, digit=2){
  paste(
    "[",
    round(min(x), digit=digit),
    ",",
    round(max(x), digit=digit),
    "]",
    sep=""
  )
}

## ~~~~~~~
## Transform to letter --------------------------------------------------------------------
## ~~~~~~~

# This function allows you to transform number <= 10 to letters.

numberToLetter <- function(x){
  vector.v <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
  if(x<11){
    return(vector.v[x])
  } else{
    return(x)
  }
}

## ~~~~~~~~~
## To count and display word and character number --------------------------------------------------------------------
## ~~~~~~~~~

# Calculates number of words in a Rmarkdown document. Seems rather inaccurate.

#from: https://stackoverflow.com/questions/46317934/exclude-sections-from-word-count-in-r-markdown

#The file is the RMd file. All parts between
#<!---TC:ignore--->
#
#will not be counted. References are not counted either.
#Note that you rmd file need to be up to date. Therefore you have to knit twice.

RmdWords <- function(file) {

	library(stringr)
	library(tidyverse)

  # Creates a string of text
  file_string <- file %>%
    readLines() %>%
    paste0(collapse = " ") %>%
	#Remove commented text 
	str_replace_all("<!-- .*? -->", "") %>%
    # Remove YAML header
    str_replace_all("^<--- .*?--- ", "") %>%    
    str_replace_all("^--- .*?--- ", "") %>%
    # Remove code
    str_replace_all("```.*?```", "") %>%
    str_replace_all("`.*?`", "") %>%
    # Remove LaTeX
    str_replace_all("[^\\\\]\\$\\$.*?[^\\\\]\\$\\$", "") %>%
    str_replace_all("[^\\\\]\\$.*?[^\\\\]\\$", "") %>%
	# Remove citation
	#str_replace_all(" \@ref{.*?}", "") %>%#To fig/table
	#str_replace_all("[@.*?]", "") %>%#General
    # Delete text between tags
    str_replace_all("TC:ignore.*?TC:endignore", "") %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_replace_all("  ", "") %>%
	# Remove formula
	str_replace_all("$ .*? $", "") 
    #str_replace_all("<", "") %>%
    #str_replace_all(">", "") 

	
  # Save several different results
  word_count <- str_count(file_string, "\\S+")
  char_count <- str_replace_all(string = file_string, " ", "") %>% str_count()

   return(list(num_words = word_count, num_char = char_count, word_list = file_string))
}

# ~~~~~	

## ~~~~~~~~~
## Merging the .bib file with articles with a newly created .bib with used packages --------------------------------------------------------------------
## ~~~~~~~~~

# Merge two bib files: one existing (referenced literature), and on created for R packages used.

#To cite R package directly, from: https://stackoverflow.com/questions/60026015/use-citation-in-r-markdown-to-automatically-generate-a-bibliography-of-r-packa

#bibliographyArticle: complete (absolute) path towards the .bib with all cited articles
#bibliographyArticle: complete (absolute) path towards the .bib with all cited articles + added R packages (i.e. OUTPUT)
#then, the different packages to use, separated by commas: e.g. readr, lubridate, Rcpp and so on. You can then cite them within the Rmarkdown file using @packageName. 
#But be careful, in doing so, some packages might be referenced with multiple citation (i.e. no "Rcpp" found because it creates @Rcpp1, @Rcpp2 etc...


citeR <- function(bibliographyArticle, bibliographyOutput, ...)
{

  packages <- unlist(lapply(as.list(match.call()), deparse))[c(-1, -2, -3)]
  #Withdraw argument for capital letters if inserted
  if(packages[1]=="FALSE"){
    packages <- packages[-1]
  }
  names(packages) <- NULL
  #Rbibs <- ""
  library(bibtex)
  write.bib(packages, file = "packages.bib", append = FALSE, verbose = TRUE)

  bib_init <- readLines(bibliographyArticle)
  
  big_bib <- c(bib_init, "\n", readLines("packages.bib"))

  #Add {{  }} for title so as to force keeping capital letters
  big_bib <- 
		sapply(big_bib, function(x) {
			#if(keepCapitalLetters==FALSE){
			  if(length(grep( "title", x, fixed=TRUE))>0){
				x <- gsub("title = { ", "title = {{ ", x, fixed=TRUE)
				x <- gsub("title ={ ", "title = {{ ", x, fixed=TRUE)
				x <- gsub("title= { ", "title = {{ ", x, fixed=TRUE)
				x <- gsub("title={ ", "title = {{ ", x, fixed=TRUE)
				x <- gsub(" },", " }},", x, fixed=TRUE)
				
				#Italic for species
				x <- gsub("(", "(\\textit{", x, fixed=TRUE)
				x <- gsub(")", "})", x, fixed=TRUE)
			  }
			#}else{
				##Keep capital for species name
				#if(length(grep( "title", x, fixed=TRUE))>0){
					#x <- gsub("(", "(\textit{", x, fixed=TRUE)
					#x <- gsub(")", "})", x, fixed=TRUE)
				#}
			#}
		 return(x)
		}
		)
  writeLines(big_bib, bibliographyOutput) 
}

## ~~~~~~~~~
## Extracting the used R version to print --------------------------------------------------------------------
## ~~~~~~~~~

# To extract the current version of R used for analysis.

extractRversion <- function(...){

info <- R.version
info <- data.frame(do.call(cbind,info))

library(tidyr)
info <- separate(info, col = "version.string", into=c("useless", "useless2", "version", "dateVersion"), sep=" ")
return(paste("v.",info$version[1], sep=""))

}
