## Header and raw data import ----

rm(list=(ls()))
setwd("C:/cygwin64/home/Andrew/Dropbox/shared_AL/Porto_iTRAQ_neutral-sites/")

gfpMeasurements <- read.csv("GFP_measurements.csv")

## Raw data visutalisation ----
plot(gfpMeasurements)
x<-gfpMeasurements[gfpMeasurements$Condition=="WT",]
plot(x)
plot(x[x$Bio.Replicate==1,])

plot(gfpMeasurements$Condition,gfpMeasurements$Value - mean(gfpMeasurements[gfpMeasurements$Condition=="WT",]$Value))

## Pulling values for level-specific variation ----

for(i in 1:3){
  print(mean(x[x$Bio.Replicate==i,1]))
}

for(i in 1:length(levels(gfpMeasurements$Condition))){
  print(levels(gfpMeasurements$Condition)[i])
}

## Generates mean for all readings of each biological replicate of each condition
a<-c()
for(i in 1:length(levels(gfpMeasurements$Condition))){
  x<-gfpMeasurements[
    gfpMeasurements$Condition==levels(gfpMeasurements$Condition)[i],
    ]
  b <- c()
  for(j in 1:length(levels(gfpMeasurements$Bio.Replicate))){
    b<-c(b,
      mean(
        x[x$Bio.Replicate==j,]$Value
        )
      )
  }
  a <- rbind(a,b)
}



## Heatmap Protein Analysis ----
require(gplots)
proteinData <- read.csv("./exportQuants.csv", row.names = 1)

heatmap.colours <- colorRampPalette(colors = c("green","black","red"))

heatmap.2(as.matrix(proteinData), col = heatmap.colours(512), scale="row")
heatmap.2(as.matrix(proteinData), col = heatmap.colours(512))

heatmap.2(as.matrix(proteinData[1:364,c(1:11,13:15)]), col = heatmap.colours(512))
heatmap.2(as.matrix(proteinData[1:364,c(1:11,13:15)]), col = heatmap.colours(512), scale="row")

proteinData[360:364,]
