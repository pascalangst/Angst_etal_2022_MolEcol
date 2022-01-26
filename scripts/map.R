# Draw a tree connected to points on a map
setwd("~/bioinformatics/GeographicPhylogeny/")

library(phytools)

# load tree and coordinates (both is ordered alphabetically)
tree <- read.nexus("finalconcatenation.tre")
coords <- read.delim("coords.txt", row.names = 1)
# combine them to one object
obj<-phylo.to.map(tree,coords,plot=FALSE)
# this can be plotted
plot(obj, ftype="i")

# colors:
# find the tips left and right of the base
tmp<-tree$edge[which(tree$edge[,1]==(Ntip(tree)+1)),2]
left<-tmp[1]
middle<-tmp[2]
# find the tips right of the second node
right<-tmp[2]+1
# empty color matrix
colors<-matrix(NA,nrow(coords),2,dimnames=list(rownames(coords)))
# fill matrix
blue<-tree$tip.label[getDescendants(tree,middle)]
blue<-blue[!is.na(blue)]
for(i in 1:length(blue))
  colors[blue[i],1:2]<-c("steelblue4","steelblue4")
green<-tree$tip.label[getDescendants(tree,right)]
green<-green[!is.na(green)]
for(i in 1:length(green))
  colors[green[i],1:2]<-c("steelblue1","steelblue1")
red<-tree$tip.label[getDescendants(tree,left)]
red<-red[!is.na(red)]
for(i in 1:length(red))
  colors[red[i],1:2]<-c("orangered2","orangered2")
# this can be plotted
plot(obj, ftype="i", colors=colors)

# zoom in and make some space for the axes, remove the points at the tips, draw full lines and make the points in the map slightly bigger
obj<-phylo.to.map(tree,coords, ftype = "i", ylim=c(30,80), xlim=c(-5,110), mar=c(5.1,5.1,2.1,2.1), colors=colors, pts = FALSE, lty = 1, cex.points=c(0,1.5))

# # Change some things in the develpers code
# trace("plot.phylo.to.map",edit=TRUE)
# 
# # line 122: no fill, change border color and width
# map(map, add = TRUE, fill = FALSE, col = "darkgray", lwd = 2, mar = rep(0,
# # line 156: lines should finish before the points
# coords[i, 2] + 2), col = colors[i, 1], lty = lty,
# # line 159: change the color of the points
# points(coords, pch = pch, col = colors[,1], cex = cex.points[2], bg = colors[, 
#   
# # lines 161 ff: change colors of branches
# for (i in 1:5) lines(rep(x[cw$edge[i, 2]], 2), Y[i, 
#                                                  ], lwd = lwd[1], lend = 2, col = "orangered2")
# for (i in 7:17) lines(rep(x[cw$edge[i, 2]], 2), Y[i, 
#                                                   ], lwd = lwd[1], lend = 2, col = "steelblue1")
# for (i in 18:nrow(Y)) lines(rep(x[cw$edge[i, 2]], 
#                                 2), Y[i, ], lwd = lwd[1], lend = 2, col = "steelblue4")
# for (i in 6) lines(rep(x[cw$edge[i, 2]], 2), Y[i, 
#                                                ], lwd = lwd[1], lend = 2, col = "black")
# for (i in c(1, 4) + n) lines(range(x[cw$edge[which(cw$edge[, 
#                                                            1] == i), 2]]), Y[which(cw$edge[, 1] == i), 1], 
#                              lwd = lwd[1], lend = 2, col = "black")
# for (i in 2:3 + n) lines(range(x[cw$edge[which(cw$edge[, 
#                                                        1] == i), 2]]), Y[which(cw$edge[, 1] == i), 1], 
#                          lwd = lwd[1], lend = 2, col = "orangered2")
# for (i in 5:9 + n) lines(range(x[cw$edge[which(cw$edge[, 
#                                                        1] == i), 2]]), Y[which(cw$edge[, 1] == i), 1], 
#                          lwd = lwd[1], lend = 2, col = "steelblue1")
# for (i in 10:cw$Nnode + n) lines(range(x[cw$edge[which(cw$edge[, 
#                                                                1] == i), 2]]), Y[which(cw$edge[, 1] == i), 1], 
#                                  lwd = lwd[1], lend = 2, col = "steelblue4")

# Saved changes to this file. It can be sourced now
source("plot.phylo.to.map_ham.R")

# plot again
obj<-phylo.to.map(tree,coords, ftype = "i", ylim=c(30,80), xlim=c(-10,110), mar=c(5.1,5.1,2.1,2.1), colors=colors, pts = FALSE, lty = 1, cex.points=c(0,1.5))

# create the axes
long.lines<-seq(-40,180,by=20)
lat.lines<-seq(10,80,by=10)
# plot axes
axis(1,at=long.lines,cex.axis=1.5)
axis(2,at=lat.lines,cex.axis=1.5)
# fake axes to make a frame for the map
axis(3,at=c(-50,170), tck =0 , labels = c("",""), pos = 80)
axis(4,at=c(0,80), tck =0 , labels = c("",""))
# axis titles
title(xlab="Longitude (degrees)", cex.lab = 1.5)
title(ylab="Latitude (degrees)                               ", cex.lab = 1.5)
# vertival line D.magna
segments(85,30,85,75, lty = 2, lwd = 2)
# line of H.t. split
segments(120,40,-10,50, lty = 3, lwd = 3, col = "steelblue3")

# map using projection:

#modify functions: 
#  phylo.to.map():
#  map <- map(database, regions, xlim = xlim, ylim = ylim, plot = FALSE, 
#             fill = TRUE, resolution = 0,proj='orth',orient=c(10, 50, 0))
#
#plot.phylo.to.map():
#  lines(c(x[tip.i], mapproject(list(y=obj$coords[,"lat"], x=obj$coords[,"long"]))$x[i]), c(Y[which(cw$edge[, 
#                                                                                                           2] == tip.i), 2] - if (from.tip) 0 else sh[tip.i], 
#                                                                                           mapproject(list(y=obj$coords[,"lat"]+2, x=obj$coords[,"long"]))$y[i]), col = colors[i, 1], lty = lty,
#        lwd = lwd[2])
#}
#points(mapproject(list(y=obj$coords[,"lat"], x=obj$coords[,"long"])), pch = pch, col = colors[,1], cex = cex.points[2], bg = colors[, 
#                                                                                                                                    2])   
#
#map.grid():
#  text(mapproject(expand.grid(x = tx + xinc + 35 * 0.5, 
#                              y = y + yinc * 0.05)), labels = auto.format(y), cex = cex, 
#       adj = c(0, 0), col = col, font = font, ...)	    

source("plot.phylo.to.map_proj_ham.R") # orthographic projection

obj<-phylo.to.map(tree,coords, ftype = "i", ylim=c(20,80), xlim=c(-10,110), mar=c(5.1,5.1,2.1,2.1), colors=colors, pts = FALSE, lty = 1, cex.points=c(0,1.5))
plot(obj, ftype="i",mar=c(6.1,5.1,2.1,2.1), colors=colors, pts = FALSE, lty = 0, cex.points=c(0,1.5))
plot(obj, ftype="i",mar=c(0,0,2.1,0), colors=colors, pts = FALSE, lty = 1, cex.points=c(0,1.5))             
#points(mapproject(list(y=obj$coords[,"lat"], x=obj$coords[,"long"])), pch = 21, col = colors[,1], cex = , cex.points=c(0,1.5), bg = colors[,2]) 
map.grid(lim = c(-50,140,20,80), col="black", lwd=0.3, nx = 5,ny=5, pretty=T,cex = 0.75)
