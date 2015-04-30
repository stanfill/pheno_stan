#Use k-means to group different genetic variables

library(simplepheno)
library(data.table)
library(rpart)
data(phenology_temperate_2011);data(phenology_heat_2011)
data(phenology_temperate_2012);data(phenology_heat_2012)
data(phenology_temperate_2013);data(phenology_heat_2013)
phenology_temperate_2012 <- na.omit(phenology_temperate_2012)

temperatedat <- rbind(phenology_temperate_2011,phenology_temperate_2012)
temperatedat <- rbind(temperatedat,phenology_temperate_2013)
temperatedat <- temperatedat[,-3]

temperatedat$GID <- as.factor(temperatedat$GID)
temperatedat <- as.data.table(temperatedat)
setkey(temperatedat,GID)
plot(temperatedat)

#Tree analysis
temperatedat$numGID <- temperatedat$GID
levels(temperatedat$numGID) <- 1:length(unique(temperatedat$GID))
tres <- rpart(dth~numGID,data=temperatedat)
printcp(tres)
plotcp(tres)

plot(tres, uniform=TRUE, main="Classification Tree for Days to Heading")
text(tres, use.n=TRUE, all=FALSE, cex=.8)

ptres <- prune(tres,cp=.02)
printcp(ptres)
plotcp(ptres)

plot(ptres, uniform=TRUE, main="Classification Tree for Days to Heading")
text(ptres, use.n=TRUE, all=TRUE, cex=.8)


G1GID <- levels(temperatedat$GID)[c(15,16,17,19,20,21,22,23,28,29)]
G2GID <- levels(temperatedat$GID)[c(1,2,3,8,10,13,14,18,24,25)]
G3GID <- levels(temperatedat$GID)[c(4,5,6,7,9,11,12,26,27,30)]




