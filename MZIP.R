

######## Define function ##########
get_HFD = function(sim.wide){
  m <- ncol(sim.wide)-2
  HFD <- apply(sim.wide[,paste0("d.",1:m)],1,function(x) sum(x==0))
  return(HFD)
}
# optimization function for MZINB - shared parameter 
lMZINB = function(dat, par){
  y <- dat$HFD
  t <- ifelse(dat$trt=="cntl",0,1)
  loglike=rep(NA, length(y))
  gamma0=par[1]
  gamma1=par[2]
  phi=par[3] #to make the parameter space to be (-inf, inf)
  beta0=par[4]
  beta1=par[5]
  sigma=exp(phi)#rescale parameter space to be (0, inf) as domain for dispersion
  nu=exp(beta0+beta1*t)#(marginal mean) 
  p = 1/(1+exp(-gamma0-gamma1*t)) #proportion of zero
  mu = nu/(1-p) #mean for count part
  alpha = 1/sigma #param for neg binom
  theta = 1/(1+(mu/alpha)) #param for neg binom (alpha/(alpha+mu))
  for (i in 1:length(y)){
    if (y[i]==0){loglike[i] =log(p[i] + (1-p[i])*(theta[i]^alpha))}
    else {loglike[i] = log(1-p[i]) + lgamma(y[i]+alpha) - lgamma(alpha)+y[i]*log(1-theta[i])+alpha*log(theta[i]) - lgamma(y[i]+1)}
  }
  return(-sum(loglike))
}

lMZIP = function(dat, par){
  y <- dat$HFD
  t <- ifelse(dat$trt=="cntl",0,1)
  loglike=rep(NA, length(y))
  gamma0=par[1]
  gamma1=par[2]
  beta0=par[3]
  beta1=par[4]
  nu=exp(beta0+beta1*t)#(marginal mean) 
  p = 1/(1+exp(-gamma0-gamma1*t)) #proportion of zero
  mu = nu/(1-p) #mean for count part
  for (i in 1:length(y)){
    if (y[i]==0){loglike[i] =log(p[i] + (1-p[i])*exp(-mu[i]))}
    else {loglike[i] = log(1-p[i]) -mu[i] + y[i]*log(mu[i]) - lgamma(y[i]+1)}
  }
  return(-sum(loglike))
}

lMZIQP = function(dat, par){
  y <- dat$HFD
  t <- ifelse(dat$trt=="cntl",0,1)
  loglike=rep(NA, length(y))
  gamma0=par[1]
  gamma1=par[2]
  beta0=par[3]
  beta1=par[4]
  phi=par[5]#to make the parameter space to be (-inf, inf)
  sigma=exp(phi)#rescale parameter space to be (0, inf) as domain for dispersion
  nu=exp(beta0+beta1*t)#(marginal mean) 
  p = 1/(1+exp(-gamma0-gamma1*t)) #proportion of zero
  mu = nu/(1-p) #mean for count part
  for (i in 1:length(y)){
    if (y[i]==0){loglike[i] =log(p[i] + (1-p[i])*exp(-mu[i]))}
    else {loglike[i] = log(1-p[i]) + log(mu[i]) + (y[i]-1) *log(mu[i]+sigma*y[i]) - mu[i] -sigma*y[i] - lgamma(y[i]+1)}
  }
  return(-sum(loglike))
}

lMZIQP_wrapper = function(dat,try=20){
  
  fit1.1 <- zeroinfl(HFD ~ trt,data=dat,dist="poisson")
  fit.zero <- summary(fit1.1)$coefficients$zero
  fit.count <- summary(fit1.1)$coefficients$count
  fit.p0 <- 1-exp(fit.zero[1,"Estimate"])/(1+exp(fit.zero[1,"Estimate"]))
  fit.mu0 <- exp(fit.count[1,"Estimate"])
  fit.p1 <- 1-exp(sum(fit.zero[,"Estimate"]))/(1+exp(sum(fit.zero[,"Estimate"])))
  fit.mu1 <- exp(sum(fit.count[1:2,"Estimate"]))
  HFD0 <- (1-fit.p0)*fit.mu0
  HFD1 <- (1-fit.p1)*fit.mu1
  
  r0= fit.zero[1,1]; r1=fit.zero[2,1]
  b0 <- HFD0 ;b1=HFD1-HFD0
  inits3=c(r0,r1,b0,b1,0)
  inits2=c(r0,r1,0,0,0)
  inits1=c(0,0,0,0,0)
  
  fit <- optim(par=inits1, lMZIQP,dat=dat, method='BFGS',
               hessian=T)
  
  control <- regression.control(fit)
  
  if(control!=2){
    fit <- optim(par=inits2, lMZIQP,dat=dat, method='BFGS',
                 hessian=T)
    control <- regression.control(fit)
  }
  
  
  if(control!=2){
    fit <- optim(par=inits3,lMZIQP,dat=dat, method='BFGS',
                 hessian=T)
    control <- regression.control(fit)
  }
  
  if(control!=2){
    r0= fit.zero[1,1]; r1=fit.zero[2,1];b0 <- HFD0 ;b1=HFD1-HFD0
    a=r0/try;b=r1/try;c=b0/try;d=b1/try
    r0=r1=b0=b1=0
    while(control!=2){
      r0=r0+a;r1=r1+b;b0=b0+c;b1=b1+d
      inits=c(r0,r1,b0,b1,0)
      fit <- optim(par=inits,lMZIQP,dat=dat, method='BFGS',
                   hessian=T)
      control <- regression.control(fit)
    }
  }
  return(fit)
}  

######## illustration ##########
### set up 
library(rstudioapi)
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

sim.wide <- read.csv("simulated_data.csv")[,-1]
sim.wide$HFD <- get_HFD(sim.wide)
# 0 - alive and outside hosp
# 1 - alive and in hosp
# 2 - died

