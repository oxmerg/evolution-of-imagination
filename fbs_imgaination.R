#imagination optimal control

obs=10000 
nobs=obs+1
tol=0.001

beta=array(0,dim=c(nobs)) 
beta1=array(0,dim=c(nobs)) 
oldbeta=array(0,dim=c(nobs)) 
x=array(0,dim=c(nobs)) 
oldx=array(0,dim=c(nobs))
y=array(0,dim=c(nobs))
oldy=array(0,dim=c(nobs))
M=array(0,dim=c(nobs))
oldM=array(0,dim=c(nobs)) 

lambda1=array(0,dim=c(nobs)) 
lambda2=array(0,dim=c(nobs))
lambda3=array(0,dim=c(nobs))
oldlambda1=array(0,dim=c(nobs))
oldlambda2=array(0,dim=c(nobs))
oldlambda3=array(0,dim=c(nobs)) 

nbeta=array(0,dim=c(nobs)) 
nbeta1=array(0,dim=c(nobs)) 
noldbeta=array(0,dim=c(nobs)) 
nx=array(0,dim=c(nobs)) 
noldx=array(0,dim=c(nobs))
ny=array(0,dim=c(nobs))
noldy=array(0,dim=c(nobs))
N=array(0,dim=c(nobs))
noldN=array(0,dim=c(nobs)) 

nlambda1=array(0,dim=c(nobs)) 
nlambda2=array(0,dim=c(nobs))
nlambda3=array(0,dim=c(nobs))
noldlambda1=array(0,dim=c(nobs))
noldlambda2=array(0,dim=c(nobs))
noldlambda3=array(0,dim=c(nobs)) 

pow=function(x,n){x^n}

x[1]=1.0
y[1]=1.0
M[1]=20.0

nx[1]=1.0
ny[1]=1.0
N[1]=20.0

alpha=1.75
rho=0.5
c=2.25
gamma=0.1

epsilon1=3.0
epsilon2=2.0

nalpha=1.75
nrho=0.5
nc=2.5
ngamma=0.1

nepsilon1=3.0
nepsilon2=2.0

h=0.0025
ht=h/2
count=0

repeat {
  
  for(i in 1:obs){
    oldbeta[i]=beta[i]
    oldx[i]=x[i]
    oldy[i]=y[i]
    oldM[i]=M[i]
    oldlambda1[i]=lambda1[i]
    oldlambda2[i]=lambda2[i]
    oldlambda3[i]=lambda3[i]
    
    noldbeta[i]=nbeta[i]
    noldx[i]=nx[i]
    noldy[i]=ny[i]
    noldN[i]=N[i]
    noldlambda1[i]=nlambda1[i]
    noldlambda2[i]=nlambda2[i]
    noldlambda3[i]=nlambda3[i]
    }
  
  for(j in 1:obs){
    
    k11=y[j]
    k12=c*(1-pow(x[j],2))*y[j]-x[j]
    k13=alpha-rho*M[j]-gamma*M[j]*N[j]+beta[j]*(y[j]+c*(1-pow(x[j],2))*y[j]-x[j])
    
    nk11=ny[j]
    nk12=nc*(1-pow(nx[j],2))*ny[j]-nx[j]
    nk13=nalpha-nrho*N[j]-ngamma*M[j]*N[j]+nbeta[j]*(ny[j]+nc*(1-pow(nx[j],2))*ny[j]-nx[j])
    
    #------------------
    
    k21=y[j]+ht*k12
    k22=c*(1-pow((x[j]+ht*k11),2))*(y[j]+ht*k12)-(x[j]+ht*k11)
    k23=alpha-rho*(M[j]+ht*k13)-gamma*(M[j]+ht*k13)*(N[j]+ht*nk13)+0.5*(beta[j]+beta[j+1])*((y[j]+ht*k12)+c*(1-pow((x[j]+ht*k11),2))*(y[j]+ht*k12)-(x[j]+ht*k11))
    
    nk21=ny[j]+ht*nk12
    nk22=nc*(1-pow((nx[j]+ht*nk11),2))*(ny[j]+ht*nk12)-(nx[j]+ht*nk11)
    nk23=nalpha-nrho*(N[j]+ht*nk13)-ngamma*(M[j]+ht*k13)*(N[j]+ht*nk13)+0.5*(nbeta[j]+nbeta[j+1])*((ny[j]+ht*nk12)+nc*(1-pow((nx[j]+ht*nk11),2))*(ny[j]+ht*nk12)-(nx[j]+ht*nk11))
    
    #------------------
    
    k31=y[j]+ht*k22
    k32=c*(1-pow((x[j]+ht*k21),2))*(y[j]+ht*k22)-(x[j]+ht*k21)
    k33=alpha-rho*(M[j]+ht*k23)-gamma*(M[j]+ht*k23)*(N[j]+ht*nk23)+0.5*(beta[j]+beta[j+1])*((y[j]+ht*k22)+c*(1-pow((x[j]+ht*k21),2))*(y[j]+ht*k22)-(x[j]+ht*k21))
    
    nk31=ny[j]+ht*nk22
    nk32=nc*(1-pow((nx[j]+ht*nk21),2))*(ny[j]+ht*nk22)-(nx[j]+ht*nk21)
    nk33=nalpha-nrho*(N[j]+ht*nk23)-ngamma*(M[j]+ht*k23)*(N[j]+ht*nk23)+0.5*(nbeta[j]+nbeta[j+1])*((ny[j]+ht*nk22)+nc*(1-pow((nx[j]+ht*nk21),2))*(ny[j]+ht*nk22)-(nx[j]+ht*nk21))
    
    #------------------
    
    k41=y[j]+h*k32
    k42=c*(1-pow((x[j]+h*k31),2))*(y[j]+h*k32)-(x[j]+h*k31)
    k43=alpha-rho*(M[j]+h*k33)-gamma*(M[j]+h*k33)*(N[j]+h*nk33)+beta[j+1]*((y[j]+h*k32)+c*(1-pow((x[j]+h*k31),2))*(y[j]+h*k32)-(x[j]+h*k31))
    
    nk41=ny[j]+h*nk32
    nk42=nc*(1-pow((nx[j]+h*nk31),2))*(ny[j]+h*nk32)-(nx[j]+h*nk31)
    nk43=nalpha-nrho*(N[j]+h*nk33)-ngamma*(M[j]+h*k33)*(N[j]+h*nk33)+nbeta[j+1]*((ny[j]+h*nk32)+nc*(1-pow((nx[j]+h*nk31),2))*(ny[j]+h*nk32)-(nx[j]+h*nk31))
    
    #------------------
    
    x[j+1]=x[j]+(h/6.0)*(k11+2*k21+2*k31+k41);
    y[j+1]=y[j]+(h/6.0)*(k12+2*k22+2*k32+k42);
    M[j+1]=M[j]+(h/6.0)*(k13+2*k23+2*k33+k43);
    
    nx[j+1]=nx[j]+(h/6.0)*(nk11+2*nk21+2*nk31+nk41);
    ny[j+1]=ny[j]+(h/6.0)*(nk12+2*nk22+2*nk32+nk42);
    N[j+1]=N[j]+(h/6.0)*(nk13+2*nk23+2*nk33+nk43);
    
      
  }
  
  
  for(i in 1:obs){
    j=obs+2-i;
    
    dbeta=0.5*(beta[j]+beta[j-1])
    ndbeta=0.5*(nbeta[j]+nbeta[j-1])
    
    kl11=-(-2*lambda2[j]*c*x[j]*y[j]- lambda2[j] - 2*lambda3[j]*c*beta[j]*x[j]*y[j]-lambda3[j]*beta[j])
    kl12=-(lambda1[j]+lambda2[j]*c*(1-pow(x[j],2))+lambda3[j]*(beta[j]+beta[j]*c*(1-pow(x[j],2))))
    kl13=-(epsilon1 - lambda3[j]*(rho+gamma*N[j]))
    
    nkl11=-(-2*nlambda2[j]*nc*nx[j]*ny[j] - nlambda2[j] - 2*nlambda3[j]*nc*nbeta[j]*nx[j]*ny[j]-nlambda3[j]*nbeta[j])
    nkl12=-(nlambda1[j]+nlambda2[j]*nc*(1-pow(nx[j],2))+nlambda3[j]*(nbeta[j]+nbeta[j]*nc*(1-pow(nx[j],2))))
    nkl13=-(nepsilon1 - nlambda3[j]*nrho-nlambda3[j]*ngamma*M[j])
    
    #------------------ 
    
    kl21=-(-2*(lambda2[j]-ht*kl12)*c*x[j]*y[j]- (lambda2[j]-ht*kl12) - 2*(lambda3[j]-ht*kl13)*c*dbeta*x[j]*y[j]-(lambda3[j]-ht*kl13)*dbeta)
    kl22=-((lambda1[j]-ht*kl11)+(lambda2[j]-ht*kl12)*c*(1-pow(x[j],2))+(lambda3[j]-ht*kl13)*(dbeta+dbeta*c*(1-pow(x[j],2))))
    kl23=-(epsilon1 - (lambda3[j]-ht*kl13)*(rho+gamma*N[j]))
    
    nkl21=-(-2*(nlambda2[j]-ht*nkl12)*nc*nx[j]*ny[j]- (nlambda2[j]-ht*nkl12) - 2*(nlambda3[j]-ht*nkl13)*nc*ndbeta*nx[j]*ny[j]-(nlambda3[j]-ht*nkl13)*ndbeta)
    nkl22=-((nlambda1[j]-ht*nkl11)+(nlambda2[j]-ht*nkl12)*nc*(1-pow(nx[j],2))+(nlambda3[j]-ht*nkl13)*(ndbeta+ndbeta*nc*(1-pow(nx[j],2))))
    nkl23=-(nepsilon1 - (nlambda3[j]-ht*nkl13)*nrho-(nlambda3[j]-ht*nkl13)*ngamma*M[j])
    
    #------------------
    
    kl31=-(-2*(lambda2[j]-ht*kl22)*c*x[j]*y[j]- (lambda2[j]-ht*kl22) - 2*(lambda3[j]-ht*kl23)*c*dbeta*x[j]*y[j]-(lambda3[j]-ht*kl23)*dbeta)
    kl32=-((lambda1[j]-ht*kl21)+(lambda2[j]-ht*kl22)*c*(1-pow(x[j],2))+(lambda3[j]-ht*kl23)*(dbeta+dbeta*c*(1-pow(x[j],2))))
    kl33=-(epsilon1 - (lambda3[j]-ht*kl23)*(rho+gamma*N[j]))
    
    nkl31=-(-2*(nlambda2[j]-ht*nkl22)*nc*nx[j]*ny[j]- (nlambda2[j]-ht*nkl22) - 2*(nlambda3[j]-ht*nkl23)*nc*ndbeta*nx[j]*ny[j]-(nlambda3[j]-ht*nkl23)*ndbeta)
    nkl32=-((nlambda1[j]-ht*nkl21)+(nlambda2[j]-ht*nkl22)*nc*(1-pow(nx[j],2)) +(nlambda3[j]-ht*nkl23)*(ndbeta+ndbeta*nc*(1-pow(nx[j],2))))
    nkl33=-(nepsilon1 - (nlambda3[j]-ht*nkl23)*nrho-(nlambda3[j]-ht*nkl23)*ngamma*M[j])
    
    #------------------
    
    kl41=-(-2*(lambda2[j]-h*kl32)*c*x[j]*y[j]- (lambda2[j]-h*kl32) - 2*(lambda3[j]-h*kl33)*c*beta[j-1]*x[j]*y[j]-(lambda3[j]-h*kl33)*beta[j-1])
    kl42=-((lambda1[j]-h*kl31)+(lambda2[j]-h*kl32)*c*(1-pow(x[j],2))+(lambda3[j]-h*kl33)*(beta[j-1]+beta[j-1]*c*(1-pow(x[j],2))))
    kl43=-(epsilon1 - (lambda3[j]-h*kl33)*(rho+gamma*N[j]))
    
    nkl41=-(-2*(nlambda2[j]-h*nkl32)*nc*nx[j]*ny[j]- (nlambda2[j]-h*nkl32) - 2*(nlambda3[j]-h*nkl33)*nc*nbeta[j-1]*nx[j]*ny[j]-(nlambda3[j]-h*nkl33)*nbeta[j-1])
    nkl42=-((nlambda1[j]-h*nkl31)+(nlambda2[j]-h*nkl32)*nc*(1-pow(nx[j],2)) +(nlambda3[j]-h*nkl33)*(nbeta[j-1]+nbeta[j-1]*nc*(1-pow(nx[j],2))))
    nkl43=-(nepsilon1 - (nlambda3[j]-h*nkl33)*nrho-(nlambda3[j]-h*nkl33)*ngamma*M[j])
    
    #------------------
  
    lambda1[j-1]=lambda1[j]-(h/6.0)*(kl11+2.0*kl21+2.0*kl31+kl41)
    lambda2[j-1]=lambda2[j]-(h/6.0)*(kl12+2.0*kl22+2.0*kl32+kl42)
    lambda3[j-1]=lambda3[j]-(h/6.0)*(kl13+2.0*kl23+2.0*kl33+kl43)
    
    nlambda1[j-1]=nlambda1[j]-(h/6.0)*(nkl11+2.0*nkl21+2.0*nkl31+nkl41)
    nlambda2[j-1]=nlambda2[j]-(h/6.0)*(nkl12+2.0*nkl22+2.0*nkl32+nkl42)
    nlambda3[j-1]=nlambda3[j]-(h/6.0)*(nkl13+2.0*nkl23+2.0*nkl33+nkl43)
  }
  
  for(kv in 1:obs){
    beta1[kv]=lambda3[kv]*(y[kv]+c*(1-pow(x[kv],2))*y[kv]-x[kv])/(2*epsilon2)
    beta[kv]=(beta1[kv]+oldbeta[kv])/2
    
    nbeta1[kv]=nlambda3[kv]*(ny[kv]+nc*(1-pow(nx[kv],2))*ny[kv]-nx[kv])/(2*nepsilon2)
    nbeta[kv]=(nbeta1[kv]+noldbeta[kv])/2}
  
  st=array(0, dim=c(14)) 
  tmp=array(0, dim=c(7)) 
  
  nst=array(0, dim=c(14)) 
  ntmp=array(0, dim=c(7)) 
  
  for(kv in 1:obs){
    
    #player 1
    st[1]=st[1]+abs(beta[kv])
    st[2]=st[2]+abs(oldbeta[kv]-beta[kv])
    
    st[3]=st[3]+abs(x[kv])
    st[4]=st[4]+abs(oldx[kv]-x[kv])
    
    st[5]=st[5]+abs(y[kv])
    st[6]=st[6]+abs(oldy[kv]-y[kv])
    
    st[7]=st[7]+abs(M[kv])
    st[8]=st[8]+abs(oldM[kv]-M[kv])
    
    st[9]=st[9]+abs(lambda1[kv])
    st[10]=st[10]+abs(oldlambda1[kv]-lambda1[kv])
    
    st[11]=st[11]+abs(lambda2[kv])
    st[12]=st[12]+abs(oldlambda2[kv]-lambda2[kv])
    
    st[13]=st[13]+abs(lambda3[kv])
    st[14]=st[14]+abs(oldlambda3[kv]-lambda3[kv])
    
    #player 2
    nst[1]=nst[1]+abs(nbeta[kv])
    nst[2]=nst[2]+abs(noldbeta[kv]-nbeta[kv])
    
    nst[3]=nst[3]+abs(nx[kv])
    nst[4]=nst[4]+abs(noldx[kv]-nx[kv])
    
    nst[5]=nst[5]+abs(ny[kv])
    nst[6]=nst[6]+abs(noldy[kv]-ny[kv])
    
    nst[7]=nst[7]+abs(N[kv])
    nst[8]=nst[8]+abs(noldN[kv]-N[kv])
    
    nst[9]=nst[9]+abs(nlambda1[kv])
    nst[10]=nst[10]+abs(noldlambda1[kv]-nlambda1[kv])
    
    nst[11]=nst[11]+abs(nlambda2[kv])
    nst[12]=nst[12]+abs(noldlambda2[kv]-nlambda2[kv])
    
    nst[13]=nst[13]+abs(nlambda3[kv])
    nst[14]=nst[14]+abs(noldlambda3[kv]-nlambda3[kv])
  }
  
  tmp[1]=tol*st[1]-st[2]
  tmp[2]=tol*st[3]-st[4]
  tmp[3]=tol*st[5]-st[6]
  tmp[4]=tol*st[7]-st[8]
  tmp[5]=tol*st[9]-st[10]
  tmp[6]=tol*st[11]-st[12]
  tmp[7]=tol*st[13]-st[14]
  
  t1=min(tmp[1],tmp[2])
  t2=min(t1,tmp[3])
  t3=min(t2,tmp[4])
  t4=min(t3,tmp[5])
  t5=min(t4,tmp[6])
  v=min(t5,tmp[7])
  
  #player 2 
  
  ntmp[1]=tol*nst[1]-nst[2]
  ntmp[2]=tol*nst[3]-nst[4]
  ntmp[3]=tol*nst[5]-nst[6]
  ntmp[4]=tol*nst[7]-nst[8]
  ntmp[5]=tol*nst[9]-nst[10]
  ntmp[6]=tol*nst[11]-nst[12]
  ntmp[7]=tol*nst[13]-nst[14]
  
  nt1=min(ntmp[1],ntmp[2])
  nt2=min(nt1,ntmp[3])
  nt3=min(nt2,ntmp[4])
  nt4=min(nt3,ntmp[5])
  nt5=min(nt4,ntmp[6])
  nv=min(nt5,ntmp[7])

  count=count+1
  if(v>0 && nv>0) {break}
}

par(mfrow=c(1,2))
t=seq(0,obs*h,length=nobs)
plot(seq(0,obs*h,length=10),seq(-20,20,length=10),type="n",bty="l",xlab="Time",
     ylab="beta (contribution to fitness)")
lines(t,beta)
lines(t,nbeta,col="blue")

plot(seq(0,obs*h,length=10),seq(0,80,length=10),type="n",bty="l",xlab="Time",
     ylab="Imagination dynamics")
lines(t,M)
lines(t,N,col="blue")