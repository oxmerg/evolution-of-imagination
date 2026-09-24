#imagination optimal control

obs=10000 
nobs=obs+1
tol=0.001

#player 1

beta=array(0,dim=c(nobs)) 
beta1=array(0,dim=c(nobs)) 
oldbeta=array(0,dim=c(nobs)) 

x=array(0,dim=c(nobs)) 
oldx=array(0,dim=c(nobs))
y=array(0,dim=c(nobs))
oldy=array(0,dim=c(nobs))

x2=array(0,dim=c(nobs)) 
oldx2=array(0,dim=c(nobs))
y2=array(0,dim=c(nobs))
oldy2=array(0,dim=c(nobs))

x3=array(0,dim=c(nobs)) 
oldx3=array(0,dim=c(nobs))
y3=array(0,dim=c(nobs))
oldy3=array(0,dim=c(nobs))

M=array(0,dim=c(nobs))
oldM=array(0,dim=c(nobs)) 

lambda1=array(0,dim=c(nobs)) 
lambda2=array(0,dim=c(nobs))

oldlambda1=array(0,dim=c(nobs))
oldlambda2=array(0,dim=c(nobs))


lambda1a=array(0,dim=c(nobs)) 
lambda2a=array(0,dim=c(nobs))

oldlambda1a=array(0,dim=c(nobs))
oldlambda2a=array(0,dim=c(nobs))

lambda1b=array(0,dim=c(nobs)) 
lambda2b=array(0,dim=c(nobs))

oldlambda1b=array(0,dim=c(nobs))
oldlambda2b=array(0,dim=c(nobs))

lambda3=array(0,dim=c(nobs))
oldlambda3=array(0,dim=c(nobs)) 
#------------------

#player 2

nbeta=array(0,dim=c(nobs)) 
nbeta1=array(0,dim=c(nobs)) 
noldbeta=array(0,dim=c(nobs)) 

nx=array(0,dim=c(nobs)) 
noldx=array(0,dim=c(nobs))
ny=array(0,dim=c(nobs))
noldy=array(0,dim=c(nobs))

nx2=array(0,dim=c(nobs)) 
noldx2=array(0,dim=c(nobs))
ny2=array(0,dim=c(nobs))
noldy2=array(0,dim=c(nobs))

nx3=array(0,dim=c(nobs)) 
noldx3=array(0,dim=c(nobs))
ny3=array(0,dim=c(nobs))
noldy3=array(0,dim=c(nobs))

N=array(0,dim=c(nobs))
noldN=array(0,dim=c(nobs)) 

nlambda1=array(0,dim=c(nobs)) 
nlambda2=array(0,dim=c(nobs))

noldlambda1=array(0,dim=c(nobs))
noldlambda2=array(0,dim=c(nobs))

nlambda1a=array(0,dim=c(nobs)) 
nlambda2a=array(0,dim=c(nobs))

noldlambda1a=array(0,dim=c(nobs))
noldlambda2a=array(0,dim=c(nobs))

nlambda1b=array(0,dim=c(nobs)) 
nlambda2b=array(0,dim=c(nobs))

noldlambda1b=array(0,dim=c(nobs))
noldlambda2b=array(0,dim=c(nobs))

nlambda3=array(0,dim=c(nobs))
noldlambda3=array(0,dim=c(nobs)) 

pow=function(x,n){x^n}

x[1]=1.0
y[1]=1.0
x2[1]=1.0
y2[1]=1.0
x3[1]=1.0
y3[1]=1.0
M[1]=20.0

nx[1]=1.0
ny[1]=1.0
nx2[1]=1.0
ny2[1]=1.0
nx3[1]=1.0
ny3[1]=1.0
N[1]=20.0

alpha=1.75
rho=0.5
c=1.25
c2=1.25
c3=1.25
gamma=0.1

epsilon1=3.0
epsilon2=2.0

nalpha=1.75
nrho=0.5
nc=1.0
nc2=1.0
nc3=1.0
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
    oldx2[i]=x2[i]
    oldy2[i]=y2[i]
    oldx3[i]=x3[i]
    oldy3[i]=y3[i]
    oldM[i]=M[i]
    oldlambda1[i]=lambda1[i]
    oldlambda2[i]=lambda2[i]
    oldlambda1a[i]=lambda1a[i]
    oldlambda2a[i]=lambda2a[i]
    oldlambda1b[i]=lambda1b[i]
    oldlambda2b[i]=lambda2b[i]
    oldlambda3[i]=lambda3[i]
    
    noldbeta[i]=nbeta[i]
    noldx[i]=nx[i]
    noldy[i]=ny[i]
    noldx2[i]=nx2[i]
    noldy2[i]=ny2[i]
    noldx3[i]=nx3[i]
    noldy3[i]=ny3[i]
    noldN[i]=N[i]
    noldlambda1[i]=nlambda1[i]
    noldlambda2[i]=nlambda2[i]
    noldlambda1a[i]=nlambda1a[i]
    noldlambda2a[i]=nlambda2a[i]
    noldlambda1b[i]=nlambda1b[i]
    noldlambda2b[i]=nlambda2b[i]
    noldlambda3[i]=nlambda3[i]
    }
  
  for(j in 1:obs){
    
    k11=y[j]
    k12=c*(1-pow(x[j],2))*y[j]-x[j]
    
    k11a=y2[j]
    k12a=c2*(1-pow(x2[j],2))*y2[j]-x2[j]
    
    k11b=y3[j]
    k12b=c3*(1-pow(x3[j],2))*y3[j]-x3[j]
    
    k13=alpha-rho*M[j]-gamma*M[j]*N[j]+beta[j]*(y[j]+c*(1-pow(x[j],2))*y[j]-x[j])+beta[j]*(y2[j]+c2*(1-pow(x2[j],2))*y2[j]-x2[j])+beta[j]*(y3[j]+c3*(1-pow(x3[j],2))*y3[j]-x3[j])
    
    #------------------
    
    nk11=ny[j]
    nk12=nc*(1-pow(nx[j],2))*ny[j]-nx[j]
    
    nk11a=ny2[j]
    nk12a=nc2*(1-pow(nx2[j],2))*ny2[j]-nx2[j]
    
    nk11b=ny3[j]
    nk12b=nc3*(1-pow(nx3[j],2))*ny3[j]-nx3[j]
    
    nk13=nalpha-nrho*N[j]-ngamma*M[j]*N[j]+nbeta[j]*(ny[j]+nc*(1-pow(nx[j],2))*ny[j]-nx[j])+nbeta[j]*(ny2[j]+nc2*(1-pow(nx2[j],2))*ny2[j]-nx2[j])+nbeta[j]*(ny3[j]+nc3*(1-pow(nx3[j],2))*ny3[j]-nx3[j])
    #------------------
    
    k21=y[j]+ht*k12
    k22=c*(1-pow((x[j]+ht*k11),2))*(y[j]+ht*k12)-(x[j]+ht*k11)
    
    k21a=y2[j]+ht*k12a
    k22a=c2*(1-pow((x2[j]+ht*k11a),2))*(y2[j]+ht*k12a)-(x2[j]+ht*k11a)
    
    k21b=y3[j]+ht*k12b
    k22b=c3*(1-pow((x3[j]+ht*k11b),2))*(y3[j]+ht*k12b)-(x3[j]+ht*k11b)
  
    k23=alpha-rho*(M[j]+ht*k13)-gamma*(M[j]+ht*k13)*(N[j]+ht*nk13)+0.5*(beta[j]+beta[j+1])*((y[j]+ht*k12)+c*(1-pow((x[j]+ht*k11),2))*(y[j]+ht*k12)-(x[j]+ht*k11))+0.5*(beta[j]+beta[j+1])*((y2[j]+ht*k12a)+c2*(1-pow((x2[j]+ht*k11a),2))*(y2[j]+ht*k12a)-(x2[j]+ht*k11a))+0.5*(beta[j]+beta[j+1])*((y3[j]+ht*k12b)+c3*(1-pow((x3[j]+ht*k11b),2))*(y3[j]+ht*k12b)-(x3[j]+ht*k11b))
    
    #------------------
    
    nk21=ny[j]+ht*nk12
    nk22=nc*(1-pow((nx[j]+ht*nk11),2))*(ny[j]+ht*nk12)-(nx[j]+ht*nk11)
    
    nk21a=ny2[j]+ht*nk12a
    nk22a=nc2*(1-pow((nx2[j]+ht*nk11a),2))*(ny2[j]+ht*nk12a)-(nx2[j]+ht*nk11a)
    
    nk21b=ny3[j]+ht*nk12b
    nk22b=nc3*(1-pow((nx3[j]+ht*nk11b),2))*(ny3[j]+ht*nk12b)-(nx3[j]+ht*nk11b)
    
    nk23=nalpha-nrho*(N[j]+ht*nk13)-ngamma*(M[j]+ht*k13)*(N[j]+ht*nk13)+0.5*(nbeta[j]+nbeta[j+1])*((ny[j]+ht*nk12)+nc*(1-pow((nx[j]+ht*nk11),2))*(ny[j]+ht*nk12)-(nx[j]+ht*nk11))+0.5*(nbeta[j]+nbeta[j+1])*((ny2[j]+ht*nk12a)+nc2*(1-pow((nx2[j]+ht*nk11a),2))*(ny2[j]+ht*nk12a)-(nx2[j]+ht*nk11a))+0.5*(nbeta[j]+nbeta[j+1])*((ny3[j]+ht*nk12b)+nc3*(1-pow((nx3[j]+ht*nk11b),2))*(ny3[j]+ht*nk12b)-(nx3[j]+ht*nk11b))
    
    #------------------
    
    k31=y[j]+ht*k22
    k32=c*(1-pow((x[j]+ht*k21),2))*(y[j]+ht*k22)-(x[j]+ht*k21)
    
    k31a=y2[j]+ht*k22a
    k32a=c2*(1-pow((x2[j]+ht*k21a),2))*(y2[j]+ht*k22a)-(x2[j]+ht*k21a)
    
    k31b=y3[j]+ht*k22b
    k32b=c3*(1-pow((x3[j]+ht*k21b),2))*(y3[j]+ht*k22b)-(x3[j]+ht*k21b)
    
    k33=alpha-rho*(M[j]+ht*k23)-gamma*(M[j]+ht*k23)*(N[j]+ht*nk23)+0.5*(beta[j]+beta[j+1])*((y[j]+ht*k22)+c*(1-pow((x[j]+ht*k21),2))*(y[j]+ht*k22)-(x[j]+ht*k21))+0.5*(beta[j]+beta[j+1])*((y2[j]+ht*k22a)+c2*(1-pow((x2[j]+ht*k21a),2))*(y2[j]+ht*k22a)-(x2[j]+ht*k21a))+0.5*(beta[j]+beta[j+1])*((y3[j]+ht*k22b)+c3*(1-pow((x3[j]+ht*k21b),2))*(y3[j]+ht*k22b)-(x3[j]+ht*k21b))
    
    #------------------
    
    nk31=ny[j]+ht*nk22
    nk32=nc*(1-pow((nx[j]+ht*nk21),2))*(ny[j]+ht*nk22)-(nx[j]+ht*nk21)
    
    nk31a=ny2[j]+ht*nk22a
    nk32a=nc2*(1-pow((nx2[j]+ht*nk21a),2))*(ny2[j]+ht*nk22a)-(nx2[j]+ht*nk21a)
    
    nk31b=ny3[j]+ht*nk22b
    nk32b=nc3*(1-pow((nx3[j]+ht*nk21b),2))*(ny3[j]+ht*nk22b)-(nx3[j]+ht*nk21b)
    
    nk33=nalpha-nrho*(N[j]+ht*nk23)-ngamma*(M[j]+ht*k23)*(N[j]+ht*nk23)+0.5*(nbeta[j]+nbeta[j+1])*((ny[j]+ht*nk22)+nc*(1-pow((nx[j]+ht*nk21),2))*(ny[j]+ht*nk22)-(nx[j]+ht*nk21))+0.5*(nbeta[j]+nbeta[j+1])*((ny2[j]+ht*nk22a)+nc2*(1-pow((nx2[j]+ht*nk21a),2))*(ny2[j]+ht*nk22a)-(nx2[j]+ht*nk21a))+0.5*(nbeta[j]+nbeta[j+1])*((ny3[j]+ht*nk22b)+nc3*(1-pow((nx3[j]+ht*nk21b),2))*(ny3[j]+ht*nk22b)-(nx3[j]+ht*nk21b))
    
    #------------------
    
    k41=y[j]+h*k32
    k42=c*(1-pow((x[j]+h*k31),2))*(y[j]+h*k32)-(x[j]+h*k31)
    
    k41a=y2[j]+h*k32a
    k42a=c2*(1-pow((x2[j]+h*k31a),2))*(y2[j]+h*k32a)-(x2[j]+h*k31a)
    
    k41b=y3[j]+h*k32b
    k42b=c3*(1-pow((x3[j]+h*k31b),2))*(y3[j]+h*k32b)-(x3[j]+h*k31b)
    
    k43=alpha-rho*(M[j]+h*k33)-gamma*(M[j]+h*k33)*(N[j]+h*nk33)+beta[j+1]*((y[j]+h*k32)+c*(1-pow((x[j]+h*k31),2))*(y[j]+h*k32)-(x[j]+h*k31))+ beta[j+1]*((y2[j]+h*k32a)+c2*(1-pow((x2[j]+h*k31a),2))*(y2[j]+h*k32a)-(x2[j]+h*k31a))+ beta[j+1]*((y3[j]+h*k32b)+c3*(1-pow((x3[j]+h*k31b),2))*(y3[j]+h*k32b)-(x3[j]+h*k31b))
    
    #------------------
    
    nk41=ny[j]+h*nk32
    nk42=nc*(1-pow((nx[j]+h*nk31),2))*(ny[j]+h*nk32)-(nx[j]+h*nk31)
    
    nk41a=ny2[j]+h*nk32a
    nk42a=nc2*(1-pow((nx2[j]+h*nk31a),2))*(ny2[j]+h*nk32a)-(nx2[j]+h*nk31a)
    
    nk41b=ny3[j]+h*nk32b
    nk42b=nc3*(1-pow((nx3[j]+h*nk31b),2))*(ny3[j]+h*nk32b)-(nx3[j]+h*nk31b)
    
    nk43=nalpha-nrho*(N[j]+h*nk33)-ngamma*(M[j]+h*k33)*(N[j]+h*nk33)+nbeta[j+1]*((ny[j]+h*nk32)+nc*(1-pow((nx[j]+h*nk31),2))*(ny[j]+h*nk32)-(nx[j]+h*nk31))+nbeta[j+1]*((ny2[j]+h*nk32a)+nc2*(1-pow((nx2[j]+h*nk31a),2))*(ny2[j]+h*nk32a)-(nx2[j]+h*nk31a))+nbeta[j+1]*((ny3[j]+h*nk32b)+nc3*(1-pow((nx3[j]+h*nk31b),2))*(ny3[j]+h*nk32b)-(nx3[j]+h*nk31b))
    
    #------------------
    
    x[j+1]=x[j]+(h/6.0)*(k11+2*k21+2*k31+k41);
    y[j+1]=y[j]+(h/6.0)*(k12+2*k22+2*k32+k42);
    
    x2[j+1]=x2[j]+(h/6.0)*(k11a+2*k21a+2*k31a+k41a);
    y2[j+1]=y2[j]+(h/6.0)*(k12a+2*k22a+2*k32a+k42a);
    
    x3[j+1]=x3[j]+(h/6.0)*(k11b+2*k21b+2*k31b+k41b);
    y3[j+1]=y3[j]+(h/6.0)*(k12b+2*k22b+2*k32b+k42b);
    
    M[j+1]=M[j]+(h/6.0)*(k13+2*k23+2*k33+k43);
    
    
    nx[j+1]=nx[j]+(h/6.0)*(nk11+2*nk21+2*nk31+nk41);
    ny[j+1]=ny[j]+(h/6.0)*(nk12+2*nk22+2*nk32+nk42);
    
    nx2[j+1]=nx2[j]+(h/6.0)*(nk11a+2*nk21a+2*nk31a+nk41a);
    ny2[j+1]=ny2[j]+(h/6.0)*(nk12a+2*nk22a+2*nk32a+nk42a);
    
    nx3[j+1]=nx3[j]+(h/6.0)*(nk11b+2*nk21b+2*nk31b+nk41b);
    ny3[j+1]=ny3[j]+(h/6.0)*(nk12b+2*nk22b+2*nk32b+nk42b);
    
    N[j+1]=N[j]+(h/6.0)*(nk13+2*nk23+2*nk33+nk43);
    
      
  }
  
  
  for(i in 1:obs){
    j=obs+2-i;
    
    dbeta=0.5*(beta[j]+beta[j-1])
    ndbeta=0.5*(nbeta[j]+nbeta[j-1])
    
    kl11=-(-2*lambda2[j]*c*x[j]*y[j]- lambda2[j] - 2*lambda3[j]*c*beta[j]*x[j]*y[j]-lambda3[j]*beta[j])
    kl12=-(lambda1[j]+lambda2[j]*c*(1-pow(x[j],2))+lambda3[j]*(beta[j]+beta[j]*c*(1-pow(x[j],2))))
    
    kl11a=-(-2*lambda2a[j]*c2*x2[j]*y2[j]- lambda2a[j] - 2*lambda3[j]*c2*beta[j]*x2[j]*y2[j]-lambda3[j]*beta[j])
    kl12a=-(lambda1a[j]+lambda2a[j]*c2*(1-pow(x2[j],2))+lambda3[j]*(beta[j]+beta[j]*c2*(1-pow(x2[j],2))))
    
    kl11b=-(-2*lambda2b[j]*c3*x3[j]*y3[j]- lambda2b[j] - 2*lambda3[j]*c3*beta[j]*x3[j]*y3[j]-lambda3[j]*beta[j])
    kl12b=-(lambda1b[j]+lambda2b[j]*c3*(1-pow(x3[j],2))+lambda3[j]*(beta[j]+beta[j]*c3*(1-pow(x3[j],2))))
    
    kl13=-(epsilon1 - lambda3[j]*(rho+gamma*N[j]))
    
    #------------------
    
    nkl11=-(-2*nlambda2[j]*nc*nx[j]*ny[j] - nlambda2[j] - 2*nlambda3[j]*nc*nbeta[j]*nx[j]*ny[j]-nlambda3[j]*nbeta[j])
    nkl12=-(nlambda1[j]+nlambda2[j]*nc*(1-pow(nx[j],2))+nlambda3[j]*(nbeta[j]+nbeta[j]*nc*(1-pow(nx[j],2))))
    
    nkl11a=-(-2*nlambda2a[j]*nc2*nx2[j]*ny2[j] - nlambda2a[j] - 2*nlambda3[j]*nc2*nbeta[j]*nx2[j]*ny2[j]-nlambda3[j]*nbeta[j])
    nkl12a=-(nlambda1a[j]+nlambda2a[j]*nc2*(1-pow(nx2[j],2))+nlambda3[j]*(nbeta[j]+nbeta[j]*nc2*(1-pow(nx2[j],2))))
    
    nkl11b=-(-2*nlambda2b[j]*nc3*nx3[j]*ny3[j] - nlambda2b[j] - 2*nlambda3[j]*nc3*nbeta[j]*nx3[j]*ny3[j]-nlambda3[j]*nbeta[j])
    nkl12b=-(nlambda1b[j]+nlambda2b[j]*nc3*(1-pow(nx3[j],2))+nlambda3[j]*(nbeta[j]+nbeta[j]*nc3*(1-pow(nx3[j],2))))
    
    nkl13=-(nepsilon1 - nlambda3[j]*nrho-nlambda3[j]*ngamma*M[j])
    
    #------------------ 
    
    kl21=-(-2*(lambda2[j]-ht*kl12)*c*x[j]*y[j]- (lambda2[j]-ht*kl12) - 2*(lambda3[j]-ht*kl13)*c*dbeta*x[j]*y[j]-(lambda3[j]-ht*kl13)*dbeta)
    kl22=-((lambda1[j]-ht*kl11)+(lambda2[j]-ht*kl12)*c*(1-pow(x[j],2))+(lambda3[j]-ht*kl13)*(dbeta+dbeta*c*(1-pow(x[j],2))))
    
    kl21a=-(-2*(lambda2a[j]-ht*kl12a)*c2*x2[j]*y2[j]- (lambda2a[j]-ht*kl12a) - 2*(lambda3[j]-ht*kl13)*c2*dbeta*x2[j]*y2[j]-(lambda3[j]-ht*kl13)*dbeta)
    kl22a=-((lambda1a[j]-ht*kl11a)+(lambda2a[j]-ht*kl12a)*c2*(1-pow(x2[j],2))+(lambda3[j]-ht*kl13)*(dbeta+dbeta*c2*(1-pow(x2[j],2))))
    
    kl21b=-(-2*(lambda2b[j]-ht*kl12b)*c3*x3[j]*y3[j]- (lambda2b[j]-ht*kl12b) - 2*(lambda3[j]-ht*kl13)*c3*dbeta*x3[j]*y3[j]-(lambda3[j]-ht*kl13)*dbeta)
    kl22b=-((lambda1b[j]-ht*kl11b)+(lambda2b[j]-ht*kl12b)*c3*(1-pow(x3[j],2))+(lambda3[j]-ht*kl13)*(dbeta+dbeta*c3*(1-pow(x3[j],2))))
    
    kl23=-(epsilon1 - (lambda3[j]-ht*kl13)*(rho+gamma*N[j]))
    
    #------------------
    
    nkl21=-(-2*(nlambda2[j]-ht*nkl12)*nc*nx[j]*ny[j]- (nlambda2[j]-ht*nkl12) - 2*(nlambda3[j]-ht*nkl13)*nc*ndbeta*nx[j]*ny[j]-(nlambda3[j]-ht*nkl13)*ndbeta)
    nkl22=-((nlambda1[j]-ht*nkl11)+(nlambda2[j]-ht*nkl12)*nc*(1-pow(nx[j],2))+(nlambda3[j]-ht*nkl13)*(ndbeta+ndbeta*nc*(1-pow(nx[j],2))))
    
    nkl21a=-(-2*(nlambda2a[j]-ht*nkl12a)*nc2*nx2[j]*ny2[j]- (nlambda2a[j]-ht*nkl12a) - 2*(nlambda3[j]-ht*nkl13)*nc2*ndbeta*nx2[j]*ny2[j]-(nlambda3[j]-ht*nkl13)*ndbeta)
    nkl22a=-((nlambda1a[j]-ht*nkl11a)+(nlambda2a[j]-ht*nkl12a)*nc2*(1-pow(nx2[j],2))+(nlambda3[j]-ht*nkl13)*(ndbeta+ndbeta*nc2*(1-pow(nx2[j],2))))
    
    nkl21b=-(-2*(nlambda2b[j]-ht*nkl12b)*nc3*nx3[j]*ny3[j]- (nlambda2b[j]-ht*nkl12b) - 2*(nlambda3[j]-ht*nkl13)*nc3*ndbeta*nx3[j]*ny3[j]-(nlambda3[j]-ht*nkl13)*ndbeta)
    nkl22b=-((nlambda1b[j]-ht*nkl11b)+(nlambda2b[j]-ht*nkl12b)*nc3*(1-pow(nx3[j],2))+(nlambda3[j]-ht*nkl13)*(ndbeta+ndbeta*nc*(1-pow(nx3[j],2))))
    
    nkl23=-(nepsilon1 - (nlambda3[j]-ht*nkl13)*nrho-(nlambda3[j]-ht*nkl13)*ngamma*M[j])
    
    #------------------
    
    kl31=-(-2*(lambda2[j]-ht*kl22)*c*x[j]*y[j]- (lambda2[j]-ht*kl22) - 2*(lambda3[j]-ht*kl23)*c*dbeta*x[j]*y[j]-(lambda3[j]-ht*kl23)*dbeta)
    kl32=-((lambda1[j]-ht*kl21)+(lambda2[j]-ht*kl22)*c*(1-pow(x[j],2))+(lambda3[j]-ht*kl23)*(dbeta+dbeta*c*(1-pow(x[j],2))))
    
    kl31a=-(-2*(lambda2a[j]-ht*kl22a)*c2*x2[j]*y2[j]- (lambda2a[j]-ht*kl22a) - 2*(lambda3[j]-ht*kl23)*c2*dbeta*x2[j]*y2[j]-(lambda3[j]-ht*kl23)*dbeta)
    kl32a=-((lambda1a[j]-ht*kl21a)+(lambda2a[j]-ht*kl22a)*c2*(1-pow(x2[j],2))+(lambda3[j]-ht*kl23)*(dbeta+dbeta*c2*(1-pow(x2[j],2))))
    
    kl31b=-(-2*(lambda2b[j]-ht*kl22b)*c3*x3[j]*y3[j]- (lambda2b[j]-ht*kl22b) - 2*(lambda3[j]-ht*kl23)*c3*dbeta*x3[j]*y3[j]-(lambda3[j]-ht*kl23)*dbeta)
    kl32b=-((lambda1b[j]-ht*kl21b)+(lambda2b[j]-ht*kl22b)*c3*(1-pow(x3[j],2))+(lambda3[j]-ht*kl23)*(dbeta+dbeta*c3*(1-pow(x3[j],2))))
    
    kl33=-(epsilon1 - (lambda3[j]-ht*kl23)*(rho+gamma*N[j]))
    
    #------------------
    
    nkl31=-(-2*(nlambda2[j]-ht*nkl22)*nc*nx[j]*ny[j]- (nlambda2[j]-ht*nkl22) - 2*(nlambda3[j]-ht*nkl23)*nc*ndbeta*nx[j]*ny[j]-(nlambda3[j]-ht*nkl23)*ndbeta)
    nkl32=-((nlambda1[j]-ht*nkl21)+(nlambda2[j]-ht*nkl22)*nc*(1-pow(nx[j],2)) + (nlambda3[j]-ht*nkl23)*(ndbeta+ndbeta*nc*(1-pow(nx[j],2))))
    
    nkl31a=-(-2*(nlambda2a[j]-ht*nkl22a)*nc2*nx2[j]*ny2[j]- (nlambda2a[j]-ht*nkl22a) - 2*(nlambda3[j]-ht*nkl23)*nc2*ndbeta*nx2[j]*ny2[j]-(nlambda3[j]-ht*nkl23)*ndbeta)
    nkl32a=-((nlambda1a[j]-ht*nkl21a)+(nlambda2a[j]-ht*nkl22a)*nc2*(1-pow(nx2[j],2))+(nlambda3[j]-ht*nkl23)*(ndbeta+ndbeta*nc2*(1-pow(nx2[j],2))))
    
    nkl31b=-(-2*(nlambda2b[j]-ht*nkl22b)*nc3*nx3[j]*ny3[j]- (nlambda2b[j]-ht*nkl22b) - 2*(nlambda3[j]-ht*nkl23)*nc3*ndbeta*nx3[j]*ny3[j]-(nlambda3[j]-ht*nkl23)*ndbeta)
    nkl32b=-((nlambda1b[j]-ht*nkl21b)+(nlambda2b[j]-ht*nkl22b)*nc3*(1-pow(nx3[j],2))+(nlambda3[j]-ht*nkl23)*(ndbeta+ndbeta*nc*(1-pow(nx3[j],2))))
    
    nkl33=-(nepsilon1 - (nlambda3[j]-ht*nkl23)*nrho-(nlambda3[j]-ht*nkl23)*ngamma*M[j])
    
    #------------------
    
    kl41=-(-2*(lambda2[j]-h*kl32)*c*x[j]*y[j]- (lambda2[j]-h*kl32) - 2*(lambda3[j]-h*kl33)*c*beta[j-1]*x[j]*y[j]-(lambda3[j]-h*kl33)*beta[j-1])
    kl42=-((lambda1[j]-h*kl31)+(lambda2[j]-h*kl32)*c*(1-pow(x[j],2))+(lambda3[j]-h*kl33)*(beta[j-1]+beta[j-1]*c*(1-pow(x[j],2))))
    
    kl41a=-(-2*(lambda2a[j]-h*kl32a)*c2*x2[j]*y2[j]- (lambda2a[j]-h*kl32a) - 2*(lambda3[j]-h*kl33)*c2*beta[j-1]*x2[j]*y2[j]-(lambda3[j]-h*kl33)*beta[j-1])
    kl42a=-((lambda1a[j]-h*kl31a)+(lambda2a[j]-h*kl32a)*c2*(1-pow(x2[j],2))+(lambda3[j]-h*kl33)*(beta[j-1]+beta[j-1]*c2*(1-pow(x2[j],2))))
    
    kl41b=-(-2*(lambda2b[j]-h*kl32b)*c3*x3[j]*y3[j]- (lambda2b[j]-ht*kl32b) - 2*(lambda3[j]-h*kl33)*c3*beta[j-1]*x3[j]*y3[j]-(lambda3[j]-h*kl33)*beta[j-1])
    kl42b=-((lambda1b[j]-h*kl31b)+(lambda2b[j]-h*kl32b)*c3*(1-pow(x3[j],2))+(lambda3[j]-h*kl33)*(beta[j-1]+beta[j-1]*c3*(1-pow(x3[j],2))))
    
    kl43=-(epsilon1 - (lambda3[j]-h*kl33)*(rho+gamma*N[j]))
    
    #------------------
    
    nkl41=-(-2*(nlambda2[j]-h*nkl32)*nc*nx[j]*ny[j]- (nlambda2[j]-h*nkl32) - 2*(nlambda3[j]-h*nkl33)*nc*nbeta[j-1]*nx[j]*ny[j]-(nlambda3[j]-h*nkl33)*nbeta[j-1])
    nkl42=-((nlambda1[j]-h*nkl31)+(nlambda2[j]-h*nkl32)*nc*(1-pow(nx[j],2)) +(nlambda3[j]-h*nkl33)*(nbeta[j-1]+nbeta[j-1]*nc*(1-pow(nx[j],2))))
    
    nkl41a=-(-2*(nlambda2a[j]-h*nkl32a)*nc2*nx2[j]*ny2[j]- (nlambda2a[j]-h*nkl32a) - 2*(nlambda3[j]-h*nkl33)*nc2*nbeta[j-1]*nx2[j]*ny2[j]-(nlambda3[j]-h*nkl33)*nbeta[j-1])
    nkl42a=-((nlambda1a[j]-h*nkl31a)+(nlambda2a[j]-h*nkl32a)*nc2*(1-pow(nx2[j],2))+(nlambda3[j]-h*nkl33)*(nbeta[j-1]+nbeta[j-1]*nc2*(1-pow(nx2[j],2))))
    
    nkl41b=-(-2*(nlambda2b[j]-h*nkl32b)*nc3*nx3[j]*ny3[j]- (nlambda2b[j]-h*nkl32b) - 2*(nlambda3[j]-h*nkl33)*nc3*nbeta[j-1]*nx3[j]*ny3[j]-(nlambda3[j]-h*nkl33)*nbeta[j-1])
    nkl42b=-((nlambda1b[j]-h*nkl31b)+(nlambda2b[j]-h*nkl32b)*nc3*(1-pow(nx3[j],2))+(nlambda3[j]-h*nkl33)*(nbeta[j-1]+nbeta[j-1]*nc*(1-pow(nx3[j],2))))
    
    nkl43=-(nepsilon1 - (nlambda3[j]-h*nkl33)*nrho-(nlambda3[j]-h*nkl33)*ngamma*M[j])
    
    #------------------
  
    lambda1[j-1]=lambda1[j]-(h/6.0)*(kl11+2.0*kl21+2.0*kl31+kl41)
    lambda2[j-1]=lambda2[j]-(h/6.0)*(kl12+2.0*kl22+2.0*kl32+kl42)
    
    lambda1a[j-1]=lambda1a[j]-(h/6.0)*(kl11a+2.0*kl21a+2.0*kl31a+kl41a)
    lambda2a[j-1]=lambda2a[j]-(h/6.0)*(kl12a+2.0*kl22a+2.0*kl32a+kl42a)
    
    lambda1b[j-1]=lambda1b[j]-(h/6.0)*(kl11b+2.0*kl21b+2.0*kl31b+kl41b)
    lambda2b[j-1]=lambda2b[j]-(h/6.0)*(kl12b+2.0*kl22b+2.0*kl32b+kl42b)
    
    lambda3[j-1]=lambda3[j]-(h/6.0)*(kl13+2.0*kl23+2.0*kl33+kl43)
    
    
    nlambda1[j-1]=nlambda1[j]-(h/6.0)*(nkl11+2.0*nkl21+2.0*nkl31+nkl41)
    nlambda2[j-1]=nlambda2[j]-(h/6.0)*(nkl12+2.0*nkl22+2.0*nkl32+nkl42)
    
    nlambda1a[j-1]=nlambda1a[j]-(h/6.0)*(nkl11a+2.0*nkl21a+2.0*nkl31a+nkl41a)
    nlambda2a[j-1]=nlambda2a[j]-(h/6.0)*(nkl12a+2.0*nkl22a+2.0*nkl32a+nkl42a)
    
    nlambda1b[j-1]=nlambda1b[j]-(h/6.0)*(nkl11b+2.0*nkl21b+2.0*nkl31b+nkl41b)
    nlambda2b[j-1]=nlambda2b[j]-(h/6.0)*(nkl12b+2.0*nkl22b+2.0*nkl32b+nkl42b)
    
    nlambda3[j-1]=nlambda3[j]-(h/6.0)*(nkl13+2.0*nkl23+2.0*nkl33+nkl43)
  }
  
  for(kv in 1:obs){
    beta1[kv]=lambda3[kv]*(y[kv]+c*(1-pow(x[kv],2))*y[kv]-x[kv]+y2[kv]+c2*(1-pow(x2[kv],2))*y2[kv]-x2[kv]+y3[kv]+c3*(1-pow(x3[kv],2))*y3[kv]-x3[kv])/(2*epsilon2)
    beta[kv]=(beta1[kv]+oldbeta[kv])/2
    
    nbeta1[kv]=nlambda3[kv]*(ny[kv]+nc*(1-pow(nx[kv],2))*ny[kv]-nx[kv]+ny2[kv]+nc2*(1-pow(nx2[kv],2))*ny2[kv]-nx2[kv]+ny3[kv]+nc3*(1-pow(nx3[kv],2))*ny3[kv]-nx3[kv])/(2*nepsilon2)
    nbeta[kv]=(nbeta1[kv]+noldbeta[kv])/2}
  
  st=array(0, dim=c(30)) 
  tmp=array(0, dim=c(15)) 
  
  nst=array(0, dim=c(30)) 
  ntmp=array(0, dim=c(15)) 
  
  for(kv in 1:obs){
    
    #player 1
    st[1]=st[1]+abs(beta[kv])
    st[2]=st[2]+abs(oldbeta[kv]-beta[kv])
    
    st[3]=st[3]+abs(x[kv])
    st[4]=st[4]+abs(oldx[kv]-x[kv])
    
    st[5]=st[5]+abs(y[kv])
    st[6]=st[6]+abs(oldy[kv]-y[kv])
    
    st[7]=st[7]+abs(x2[kv])
    st[8]=st[8]+abs(oldx2[kv]-x2[kv])
    
    st[9]=st[9]+abs(y2[kv])
    st[10]=st[10]+abs(oldy2[kv]-y2[kv])
    
    st[11]=st[11]+abs(x3[kv])
    st[12]=st[12]+abs(oldx3[kv]-x3[kv])
    
    st[13]=st[13]+abs(y3[kv])
    st[14]=st[14]+abs(oldy3[kv]-y3[kv])
    
    st[15]=st[15]+abs(M[kv])
    st[16]=st[16]+abs(oldM[kv]-M[kv])
    
    st[17]=st[17]+abs(lambda1[kv])
    st[18]=st[18]+abs(oldlambda1[kv]-lambda1[kv])
    
    st[19]=st[19]+abs(lambda2[kv])
    st[20]=st[20]+abs(oldlambda2[kv]-lambda2[kv])
    
    st[21]=st[21]+abs(lambda1a[kv])
    st[22]=st[22]+abs(oldlambda1a[kv]-lambda1a[kv])
    
    st[23]=st[23]+abs(lambda2a[kv])
    st[24]=st[24]+abs(oldlambda2a[kv]-lambda2a[kv])
    
    st[25]=st[25]+abs(lambda1b[kv])
    st[26]=st[26]+abs(oldlambda1b[kv]-lambda1b[kv])
    
    st[27]=st[27]+abs(lambda2b[kv])
    st[28]=st[28]+abs(oldlambda2b[kv]-lambda2b[kv])
    
    st[29]=st[29]+abs(lambda3[kv])
    st[30]=st[30]+abs(oldlambda3[kv]-lambda3[kv])
    
    #player 2
    nst[1]=nst[1]+abs(nbeta[kv])
    nst[2]=nst[2]+abs(noldbeta[kv]-nbeta[kv])
    
    nst[3]=nst[3]+abs(nx[kv])
    nst[4]=nst[4]+abs(noldx[kv]-nx[kv])
    
    nst[5]=nst[5]+abs(ny[kv])
    nst[6]=nst[6]+abs(noldy[kv]-ny[kv])
    
    nst[7]=nst[7]+abs(nx2[kv])
    nst[8]=nst[8]+abs(noldx2[kv]-nx2[kv])
    
    nst[9]=nst[9]+abs(ny2[kv])
    nst[10]=nst[10]+abs(noldy2[kv]-ny2[kv])
    
    nst[11]=nst[11]+abs(nx3[kv])
    nst[12]=nst[12]+abs(noldx3[kv]-nx3[kv])
    
    nst[13]=nst[13]+abs(ny3[kv])
    nst[14]=nst[14]+abs(noldy3[kv]-ny3[kv])
    
    nst[15]=nst[15]+abs(N[kv])
    nst[16]=nst[16]+abs(noldN[kv]-N[kv])
    
    nst[17]=nst[17]+abs(nlambda1[kv])
    nst[18]=nst[18]+abs(noldlambda1[kv]-nlambda1[kv])
    
    nst[19]=nst[19]+abs(nlambda2[kv])
    nst[20]=nst[20]+abs(noldlambda2[kv]-nlambda2[kv])
    
    nst[21]=nst[21]+abs(nlambda1a[kv])
    nst[22]=nst[22]+abs(noldlambda1a[kv]-nlambda1a[kv])
    
    nst[23]=nst[23]+abs(nlambda2a[kv])
    nst[24]=nst[24]+abs(noldlambda2a[kv]-nlambda2a[kv])
    
    nst[25]=nst[25]+abs(nlambda1b[kv])
    nst[26]=nst[26]+abs(noldlambda1b[kv]-nlambda1b[kv])
    
    nst[27]=nst[27]+abs(nlambda2b[kv])
    nst[28]=nst[28]+abs(noldlambda2b[kv]-nlambda2b[kv])
    
    nst[29]=nst[29]+abs(nlambda3[kv])
    nst[30]=nst[30]+abs(noldlambda3[kv]-nlambda3[kv])
  }
  
  tmp[1]=tol*st[1]-st[2]
  tmp[2]=tol*st[3]-st[4]
  tmp[3]=tol*st[5]-st[6]
  tmp[4]=tol*st[7]-st[8]
  tmp[5]=tol*st[9]-st[10]
  tmp[6]=tol*st[11]-st[12]
  tmp[7]=tol*st[13]-st[14]
  tmp[8]=tol*st[15]-st[16]
  tmp[9]=tol*st[17]-st[18]
  tmp[10]=tol*st[19]-st[20]
  tmp[11]=tol*st[21]-st[22]
  tmp[12]=tol*st[23]-st[24]
  tmp[13]=tol*st[25]-st[26]
  tmp[14]=tol*st[27]-st[28]
  tmp[15]=tol*st[29]-st[30]
  
  t1=min(tmp[1],tmp[2])
  t2=min(t1,tmp[3])
  t3=min(t2,tmp[4])
  t4=min(t3,tmp[5])
  t5=min(t4,tmp[6])
  
  t6=min(t5,tmp[7])
  t7=min(t6,tmp[8])
  t8=min(t7,tmp[9])
  t9=min(t8,tmp[10])
  
  t10=min(t9,tmp[11])
  t11=min(t10,tmp[12])
  t12=min(t11,tmp[13])
  t13=min(t12,tmp[14])
  
  v=min(t13,tmp[15])
  
  #player 2 
  
  ntmp[1]=tol*nst[1]-nst[2]
  ntmp[2]=tol*nst[3]-nst[4]
  ntmp[3]=tol*nst[5]-nst[6]
  ntmp[4]=tol*nst[7]-nst[8]
  ntmp[5]=tol*nst[9]-nst[10]
  ntmp[6]=tol*nst[11]-nst[12]
  ntmp[7]=tol*nst[13]-nst[14]
  ntmp[8]=tol*nst[15]-nst[16]
  ntmp[9]=tol*nst[17]-nst[18]
  ntmp[10]=tol*nst[19]-nst[20]
  ntmp[11]=tol*nst[21]-nst[22]
  ntmp[12]=tol*nst[23]-nst[24]
  ntmp[13]=tol*nst[25]-nst[26]
  ntmp[14]=tol*nst[27]-nst[28]
  ntmp[15]=tol*nst[29]-nst[30]
  
  nt1=min(ntmp[1],ntmp[2])
  nt2=min(nt1,ntmp[3])
  nt3=min(nt2,ntmp[4])
  nt4=min(nt3,ntmp[5])
  nt5=min(nt4,ntmp[6])
  
  nt6=min(nt5,ntmp[7])
  nt7=min(nt6,ntmp[8])
  nt8=min(nt7,ntmp[9])
  nt9=min(nt8,ntmp[10])
  nt10=min(nt9,ntmp[11])
  
  nt11=min(nt10,ntmp[12])
  nt12=min(nt11,ntmp[13])
  nt13=min(nt12,ntmp[14])
  
  nv=min(nt13,ntmp[15])

  count=count+1
  if(v>0 && nv>0) {break}
}

par(mfrow=c(1,2))
t=seq(0,obs*h,length=nobs)
plot(seq(0,obs*h,length=10),seq(-50,50,length=10),type="n",bty="l",xlab="Time",
     ylab="beta (contribution to fitness)")
lines(t,beta)
lines(t,nbeta,col="blue")

plot(seq(0,obs*h,length=10),seq(0,1,length=10),type="n",bty="l",xlab="Time",
     ylab="Imagination dynamics")
lines(t,M/max(M))
lines(t,N/max(N),col="blue")

max(N)-max(M)