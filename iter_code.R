#!/usr/bin/env Rscript

rm(list=ls())
args= commandArgs(trailingOnly = TRUE)
setwd("~/desktop")

library(data.table)

#4th order runge-kutta function

rk4=function(dydx,y,t,n){
  for(jj in	1:n)
    yt[jj]=y[jj]+dydx[jj]*h/2
  dyt=derivs(t+h/2,yt)
  for(jj in 1:n)	
    yt[jj]=y[jj]+h/2*dyt[jj]
  dym=derivs(t+h/2,yt)
  for(jj in 1:n)
    yt[jj]=y[jj]+h*dym[jj]
  dym=dyt+dym
  dyt=derivs(t+h,yt)
  for(jj in 1:n)
    yout[jj]=y[jj]+(h/6)*(dydx[jj]+dyt[jj]+2*dym[jj])
  yout
}


pastvalue=function(dyno,lag){
  
  tinc=0.01
  bfl=3000
  
  maxlag=(bfl-1)*tinc
  p=0
  
  if(lag>maxlag){
    if(p==bfl) midpnt=1 else midpnt=p+1
    past=hisy[midpnt,dyno]}
  
  else
    
  {
    dispnt=floor((lag/tinc)+0.5)
    disp=(dispnt*tinc-lag)/tinc
    
    A=0.5*disp*(disp-1)
    B=(1-disp^2)
    C=0.5*disp*(disp+1)
    
    if (dispnt>p) midpnt=p+dispnt else midpnt=bfl+p-dispnt
    
    if(midpnt==1) leftpnt=bfl else leftpnt=midpnt-1
    
    if(midpnt==bfl) rightpnt=1 else rightpnt=midpnt+1
    
    R=hisy[rightpnt,dyno]
    L=hisy[leftpnt,dyno]
    M=hisy[midpnt,dyno]
    
    past=A*L+B*M+C*R
  }
  
  past}


nvar=15

#iterate to compute LCE for different strengths of competition

nn=10
res_cv=array(0,dim=c(nn,nn))

bs=seq(0.5,1.5,length=nn)
cs=seq(0.1,1.2,length=nn)

for(ki in 1:nn){
for(ij in 1:nn){

y=array(dim=(nvar))

derivs=function(x,y){
  
  
  cd=cs[ki]
  beta=cs[ij]
  alpha=5.0
  
  cd1=0.6
  beta1=1.0
  alpha1=5.0
  
  cv=0.001*y[15]
  
  #yx[1]=y[2]
  #yx[2]=cd*(1-y[1]^2)*y[2]-y[1]
  
  #yx[3]=alpha-y[3]+beta*(yx[1]+yx[2])
  
  D=matrix(c(0,1,0,-1,cd*(1-y[1]^2),0,0,0,-1-cv),byrow=T,ncol=3)
  E=c(y[1],y[2],y[3])
  G=D%*%E
  
  yx[1]=G[1]
  yx[2]=G[2]
  yx[3]=alpha+beta*(G[1]+G[2])+G[3]
  
  A=matrix(c(0,1,0,-1-2*cd*y[1]*y[2],cd*(1-y[1]^2),0,-beta+beta*2*cd*y[1]*y[2],
             beta+beta*cd*(1-y[1]^2),-1-cv),byrow=T,ncol=3)
  B=matrix(c(y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12]),byrow=T,ncol=3)
  C=A%*%B
  
  yx[4]=C[1,1]
  yx[5]=C[1,2]
  yx[6]=C[1,3]
  
  yx[7]=C[2,1]
  yx[8]=C[2,2]
  yx[9]=C[2,3]
  
  yx[10]=C[3,1]
  yx[11]=C[3,2]
  yx[12]=C[3,3]
  
  D1=matrix(c(0,1,0,-1,cd1*(1-y[13]^2),0,0,0,-1),byrow=T,ncol=3)
  E1=c(y[13],y[14],y[15])
  G1=D1%*%E1
  
  yx[13]=G1[1]
  yx[14]=G1[2]
  yx[15]=alpha1+beta1*(G1[1]+G1[2])+G1[3]
  
  yx
}

#set the initial conditions 

y[1]=0.0
y[2]=1.0
y[3]=5.0

for(i in 4:12) y[i]=0.0
y[4]=1.0
y[8]=1.0
y[12]=1.0

y[13]=0.0
y[14]=1.0
y[15]=5.0


#set the number of iterations (obs), integration step (h), buffer length (buflen) and start time (v)

obs=7500
buflen=30
v=0.0
h=0.01

#arrays for the RK4 integrator (see below)

yt=array(1:nvar)
dyt=array(1:nvar)
dym=array(1:nvar)
yout=array(1:nvar)
yx=array(1:nvar)
eqn=array(1:nvar)

#these next two lines set up the history variables - it is necessary to carefully define what this history should be and ensure the state variable make sense

hisy=array(0,dim=c(buflen,nvar))
hisy[1:buflen,1]=60
tmp=array(0,dim=c(buflen,nvar))

#define output arrays
t1=array(1:obs)
x1=array(dim=c(obs,nvar))

#iterate the model 


for(i in 1:obs){
  eqn=derivs(v,y)
  y=rk4(eqn,y,v,nvar)
  
  #the next 7 lines of code update the history buffer
  
  for(ii in 1:buflen){
    for(jj in 1:nvar){
      tmp[ii,jj]=hisy[ii,jj]}}
  
  for(ii in 2:buflen){
    for(jj in 1:nvar){
      hisy[ii,jj]=tmp[ii-1,jj]}}
  
  for(jj in 1:nvar) hisy[1,jj]=y[jj]
  
  v=v+h
  x1[i,]=y
  t1[i]=v
  
}

res_cv[ij,ki]=mean(log(abs(x1[,10]+x1[,11]+x1[,12])/abs(1.0))/max(t1))
}
}


out=res_cv

n.colors = 90
color.fun = colorRampPalette(colors = c("darkgreen", "grey"), bias = 20)
col.key = data.table("color" = color.fun(n = n.colors),
                     "value" = seq(from = min(out), to = max(out), along.with = 1:n.colors))

# make the plot with image- leaving space for legend using fig
par(fig = c(0.1,.8,0.1,1), mar = c(2,2,2,0))
image(cs,bs,as.matrix(out,byrow=T,ncol=nn),
      useRaster=T,col = col.key$color)
mtext(text = "neural feedback strength (c)", side=1, adj = 0.5, line = 2.5, cex = 1.0)
mtext(text = "contribution to fitness (beta)", side=2, padj = 0.5, line = 3.5, cex = 1.0)

par(fig = c(.9,1,.3,.8), mar = c(0.5,0.1,0.1,2.85),new = T)
plot(x = rep(1,length(col.key$value)), 
     y = col.key$value, xlim = c(0,1), col = col.key$color, 
     type = "n", xaxs = "i", yaxs = "i", ann = F, axes = F)
segments(x0 = 0, x1 = 1, y0 = col.key$value, y1 = col.key$value, col = col.key$color, lwd = 5)
axis(side = 4,lwd = 0, las = 2, line = -.75)
mtext(text = "Legend", adj = 0, line = 1, cex = 0.8)