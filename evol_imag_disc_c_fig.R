#!/usr/bin/env Rscript

rm(list=ls())
args= commandArgs(trailingOnly = TRUE)
setwd("~/desktop")

library(data.table)

out=read.table("outcomp_effects_3oscs_beta_c_2coupling.txt",header=T)

nn=75

bs=seq(0.01,1.5,length=nn)
cs=seq(0.1,1.2,length=nn)
es=seq(0.001,1.0,length=nn)

pdf('3osc_evol_c_beta_plane_2coupling.pdf')
n.colors = 90
color.fun = colorRampPalette(colors = c("thistle1", "thistle4"), bias = 20)
col.key = data.table("color" = color.fun(n = n.colors),
                     "value" = seq(from = min(out), to = max(out), along.with = 1:n.colors))

# make the plot with image- leaving space for legend using fig
par(fig = c(0.1,.8,0.1,1), mar = c(2,2,2,0))
image(cs,bs,as.matrix(out,byrow=T,ncol=nn),
      useRaster=T,col = col.key$color)
points(0.6,1.0,pch=19)
mtext(text = "neutral feedback strength (c1)", side=1, adj = 0.5, line = 2.5, cex = 1.0)
mtext(text = "contribution to fitness (beta1)", side=2, padj = 0.5, line = 3.5, cex = 1.0)
contour(cs,bs,as.matrix(out,byrow=T,ncol=nn),add=T)

par(fig = c(.9,1,.3,.8), mar = c(0.5,0.1,0.1,2.5),new = T)
plot(x = rep(1,length(col.key$value)), 
     y = col.key$value, xlim = c(0,1), col = col.key$color, 
     type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4,lwd = 0, las = 2, line = -.75)
mtext(side=2, text = "Fitness", adj = 0.5, line = 1, cex = 0.8)
dev.off()

max(out)