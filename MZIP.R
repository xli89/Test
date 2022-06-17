
rm(list=ls())
######## Define function ##########
get_HFD = function(sim.wide){
  m <- ncol(sim.wide)-2
  HFD <- apply(sim.wide[,paste0("d.",1:m)],1,function(x) sum(x==0))
  return(HFD)
}
# optimization function for MZINB - shared parameter 
regression.control <- function(fit){
  regression.control <- tryCatch(
    {
      var.cov <- sqrt(diag(solve(fit$hessian)))
      regression.control <- 2
    },error=function(cond) {
      #          print("Regression caused an error")
      #          print("Here's the original error message:")
      #          print(cond)
      # Choose a return value in case of error
      return(regression.control=-1)
    },
    warning=function(cond) {
      #          print("Regression caused a warning:")
      #          print("Here's the original warning message:")
      #          print(cond)
      regression.control=1
      return(regression.control=1)
    }
  )
  return(regression.control)
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
    else {loglike[i] = log(1-p[i]) + lgamma(y[i]+alpha) -       lgamma(alpha)+y[i]*log(1-theta[i])+alpha*log(theta[i]) - lgamma(y[i]+1)}
  }
  return(-sum(loglike))
}

lMZIP = function(data, par){
# par: gamma0, gamma1: parameter for logit model on zero-inflated part;
#      b0, b1: parameter on marginal poisson model: V=exp(b0+b1*trt)
  y <- data$HFD
  t <- ifelse(data$trt=="cntl",0,1)
  loglike=rep(NA, length(y))
  gamma0=par[1]
  gamma1=par[2]
  beta0=par[3]
  beta1=par[4]
  nu=exp(beta0+beta1*t) #(marginal mean) 
  p = 1/(1+exp(-gamma0-gamma1*t)) #proportion of zero
  mu = nu/(1-p) #mean for count part
  for (i in 1:length(y)){
    if (y[i]==0){loglike[i] =log(p[i] + (1-p[i])*exp(-mu[i]))}
    else {loglike[i] = log(1-p[i]) -mu[i] + y[i]*log(mu[i]) - lgamma(y[i]+1)}
  }
  return(-sum(loglike))
}

lMZI_wrapper = function(HFD,trt,data,try=20,family="poisson",inits=1){
  cl <- match.call()
  
  if(is.null(data[,HFD])|!is.numeric(data[,HFD])){
    stop("HFD (Hospital free days) is not fund or is not numeric")
  }
  if(is.null(data[,trt])){
    stop("trt (treatment) is not fund or is not numeric")
  }
  
  if(! family %in% c("poisson","negbin")){
    stop("family not recognized")
  }

  if(length(inits)>1){
    if(family=="poisson" & length(inits)!=4){ 
        stop("initial parameter is not specified in appropriate format. \n
             Please specify initial parameter as c(r0,r1,b0,b1), or \n
             1: use zeroinflated model estimate \n
             0: use 0 \n")
      }else if(family=="negbin" & length(inits)!=5){
          stop("initial parameter is not specified in appropriate format. \n
             Please specify initial parameter as c(r0,r1,log(phi),b0,b1), or \n
             1: use zeroinflated model estimate \n
             0: use 0 \n")
      }
  }else if(!(inits==1|inits==0)){
    stop("initial parameter is not specified in appropriate format. \n
             Please specify initial parameter as c(r0,r1,log(phi),b0,b1), or \n
             1: use zeroinflated model estimate \n
             0: use 0 \n")
  }
  

  require(pscl)
  
  ### initialization ###
  if(length(inits)==1){
  if(inits==1){
    fit1.1 <- zeroinfl(HFD ~ trt,data=data,dist=family)
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
    if(family=="negbin"){
      phi=1/exp(fit.count[3,1])
    }
    
    if(family=="negbin"){
      inits <- c(r0,r1,phi,b0,b1)
    }else{
      inits <- c(r0,r1,b0,b1)
    }
  }else if(inits==0){
    if(family=="negbin"){
      inits <- rep(0,5)
    }else{
      inits <- rep(0,4)
    }
  }
  }
  
  lMZI <- ifelse(family=="negbin","lMZINB","lMZIP")
  
  ## initial fit ##
  a <-  paste0("optim(par=c(",paste0(inits,collapse = ","),"),",lMZI,",dat=data, method='BFGS',hessian=T)")
  fit <- eval(parse(text=a))
  
  control <- regression.control(fit)
  
  if(control!=2){
    a=r0/try;b=r1/try;c=b0/try;d=b1/try
    while(control!=2){
      r0=r0-a;r1=r1-b;b0=b0-c;b1=b1-d
      if(family=="negbin"){
        inits <- c(r0,r1,phi,b0,b1)
      }else{
        inits <- c(r0,r1,b0,b1)
      }
      a <-  paste0("optim(par=c(",paste0(inits,collapse = ","),"),",lMZI,",dat=data, method='BFGS',hessian=T)")
      fit <- eval(parse(text=a))
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

fit0 <- lMZI_wrapper(HFD="HFD",trt="trt",data=sim.wide,try=20,family="poisson",inits=1)

fit1 <- lMZI_wrapper(HFD="HFD",trt="trt",data=sim.wide,try=20,family="negbin",inits=1)

