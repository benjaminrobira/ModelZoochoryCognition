
par(mar=c(4, 7, 0.5, 0.5), mgp=c(3.5, 1, 0), xpd=TRUE)

#transform to 2 cols table for boxplot
dfNULL <- as.data.frame(
  cbind(
    as.vector(matrixNULL),
    rep(temporalKnowledge, each=numberRepetitions),
    rep("Null", times=length(temporalKnowledge)*numberRepetitions)
  )
)
colnames(dfNULL) <- c("Value", "Knowledge", "Agent")

dfINTERMEDIATE <- as.data.frame(
  cbind(
    as.vector(matrixINTERMEDIATE),
    rep(temporalKnowledge, each=numberRepetitions),
    rep("Intermediate", times=length(temporalKnowledge)*numberRepetitions)
  )
)
colnames(dfINTERMEDIATE) <- c("Value", "Knowledge", "Agent")


dfOMNISCIENT <- as.data.frame(
  cbind(
    as.vector(matrixOMNISCIENT),
    rep(temporalKnowledge, each=numberRepetitions),
    rep("Omniscient", times=length(temporalKnowledge)*numberRepetitions)
  )
)
colnames(dfOMNISCIENT) <- c("Value", "Knowledge", "Agent")

tableEfficiencyAgent_nodispersal <- rbind(dfOMNISCIENT, dfINTERMEDIATE, dfNULL)


tableEfficiencyAgent_nodispersal <- cbind(
  runif(numberRepetitions*length(temporalKnowledge)*3, 0, 1),
  rep(rep(temporalKnowledge, each=numberRepetitions), times=3),
  rep(c("Null", "Intermediate", "Omniscient"), each=numberRepetitions*length(temporalKnowledge))
)
tableEfficiencyAgent_nodispersal <- as.data.frame(tableEfficiencyAgent_nodispersal)
colnames(tableEfficiencyAgent_nodispersal) <- c("Value", "Knowledge", "Agent")

tableEfficiencyAgent_nodispersal$Value <- as.numeric(tableEfficiencyAgent_nodispersal$Value)


tableEfficiencyAgent_nodispersal$Agent <- factor(tableEfficiencyAgent_nodispersal$Agent,      # Reordering group factor levels
       levels = c("Null", "Intermediate", "Omniscient"))

library(ggplot2)
ggplot(tableEfficiencyAgent_nodispersal,aes(x=Knowledge,y=Value))+
  geom_boxplot(aes(fill=Agent))+
  stat_summary(aes(group=Agent), fun = "mean", geom = "point", shape = 19, size = 2, position=position_dodge(0.75))+
  scale_fill_manual(values=c("white", "lightgrey", "darkgrey"))+
  theme_bw()+
  theme(
    #panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
  )+
  xlab("Efficiency") + 
  ylab("Spatio-temporal knowledge rate") +
  labs(fill = "Agent type")



  
  








+
  theme(
#Grid
panel.grid.major = element_line(colour = "gray93", size=0.65),
panel.grid.minor = element_line(colour = "gray97"),




maxy = ceiling(max(dfNULL$Value, dfOMNISCIENT)*10)/10
miny = floor(min(dfreticulationTest$Value)*10)/10

emptyPlot(xlim=c(0.5,length(temporalKnowledge)+0.5), ylim=c(miny,maxy))

#add grid
addGrid(
  xmin = 0.5, xmax = length(temporalKnowledge)+0.5, xintsmall = 0.1, xintbig = 0.5,
  ymin = miny, ymax = maxy, yintsmall = (maxy - miny)/20, yintbig = (maxy - miny)/5,
  axisPlot = FALSE)#Background grid

boxPlotSaved <- boxplot(Value ~ Knowledge, data = dfreticulationTest, boxwex=0.25,
                        xlab="", ylab="", cex.lab=1.2,
                        yaxs="i", xaxs="i", las=1, tcl=-0.25,
                        frame.plot=FALSE, xaxt="n", yaxt="n",
                        outline=FALSE, add=TRUE)


for(i in 1:length(temporalKnowledge)){
  #Mask half
  rect(xleft=c(i),
       xright=c(i+0.5),
       ybottom=c(miny,miny),
       ytop=c(maxy, maxy),
       col="white",
       border = NA)
  
  #Readd grid
  addGrid(
    xmin = i, xmax = i+0.5, xintsmall = 0.1, xintbig = 0.5,
    ymin = miny, ymax = maxy, yintsmall = (maxy - miny)/20, yintbig = (maxy - miny)/5,
    axisPlot = FALSE)#Background grid
  
}

#Read interquartile range
segments(x0=c(1:length(temporalKnowledge)), x1=c(1:length(temporalKnowledge)), y0=c(boxPlotSaved$stats[1,1:length(temporalKnowledge)]), y1=c(boxPlotSaved$stats[5,1:length(temporalKnowledge)]))

#Add jitter points + link because paired

dfreticulationTest$Loc <- as.numeric(as.factor(dfreticulationTest$Knowledge))
dfreticulationTest$Loc <- dfreticulationTest$Loc + 0.2
dfreticulationTest$LocJittered <- jitter(dfreticulationTest$Loc, factor=0.5)
points(dfreticulationTest$LocJittered, dfreticulationTest$Value, cex=0.5, pch=19, col="grey", xpd=TRUE)

#Add mean points
points(x = 1:length(temporalKnowledge) + 0.2, y = apply(matrixreticulationTestAtEnd, 2, mean), pch = 19, cex=1.25)
text(x = 1:length(temporalKnowledge) + 0.2, y = apply(matrixreticulationTestAtEnd, 2, mean) - maxy/50, labels = round(apply(matrixreticulationTestAtEnd, 2, mean), digit = 2), cex=1.25)       

#x-axis
axis(side=1, line=0, at=seq(from=1, to=length(temporalKnowledge), by=1), labels=temporalKnowledge, tcl=-0.25, las=1, cex.axis=1.5)
mtext(side=1, line=2.5, at=(length(temporalKnowledge)+0.5+0.5)/2, text="Spatio-temporal knowledge rate", cex=2, font=2)

#y-axis
axis(side=2, line=0, at=seq(from=miny, to=maxy, by=(maxy-miny)/5), labels=seq(from=miny, to=maxy, by=(maxy-miny)/5), tcl=-0.25, las=1, cex.axis=1.5)
mtext(side=2, line=3, at=(maxy-miny)/2+miny, text="Reticulation", cex=2, font=2)

dev.off()


pdf(file="Presentation/Graphics/reticulationResultEmpty.pdf")
par(mar=c(4, 7, 0.5, 0.5), mgp=c(3.5, 1, 0), xpd=TRUE)
emptyPlot(xlim=c(0.5,length(temporalKnowledge)+0.5), ylim=c(miny,maxy))

#add grid
addGrid(
  xmin = 0.5, xmax = length(temporalKnowledge)+0.5, xintsmall = 0.1, xintbig = 0.5,
  ymin = miny, ymax = maxy, yintsmall = (maxy - miny)/20, yintbig = (maxy - miny)/5,
  axisPlot = FALSE)#Background grid
#x-axis
axis(side=1, line=0, at=seq(from=1, to=length(temporalKnowledge), by=1), labels=temporalKnowledge, tcl=-0.25, las=1, cex.axis=1.5)
mtext(side=1, line=2.5, at=(length(temporalKnowledge)+0.5+0.5)/2, text="Spatio-temporal knowledge rate", cex=2, font=2)

#y-axis
axis(side=2, line=0, at=seq(from=miny, to=maxy, by=(maxy-miny)/5), labels=seq(from=miny, to=maxy, by=(maxy-miny)/5), tcl=-0.25, las=1, cex.axis=1.5)
mtext(side=2, line=3, at=(maxy-miny)/2+miny, text="Reticulation", cex=2, font=2)
