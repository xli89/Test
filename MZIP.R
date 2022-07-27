# model based property, model based SE
# bootstrap CI coverage vs delta method
# different situation





rm(list=ls())
######## Define function ##########
get_HFD = function(sim.wide){
  m <- ncol(sim.wide)-2
  HFD <- apply(sim.wide[,paste0("d.",1:m)],1,function(x) sum(x==0))
  return(HFD)
}



lMZI_wrapper = function(HFD,trt,data,try=20,family=c("poisson","negbin"),
                        inits=1, # 1 - initialize with zero-inflated model; 0 - initialize with 0s. 
                        method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")){
  cl <- match.call()
  
  family <- match.arg(family)
  method <- match.arg(method)
  
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
  
  # optimization function for MZINB - shared parameter 
  lMZINB = function(data, par){
    # par: gamma0, gamma1: parameter for logit model on zero-inflated part;
    #      b0, b1: parameter on marginal poisson model: V=exp(b0+b1*trt)
    y <- data$HFD
    if(is.character(data$trt)|is.factor(data$trt)){
      t <- ifelse(data$trt=="cntl",0,1)
    }else{t=data$trt}
    loglike=rep(NA, length(y))
    gamma0=par[1]
    gamma1=par[2]
    beta0=par[3]
    beta1=par[4]
    phi=par[5] #to make the parameter space to be (-inf, inf)
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
  # optimization function for MZINB - shared parameter   
  lMZIP = function(data, par){
    # par: gamma0, gamma1: parameter for logit model on zero-inflated part;
    #      b0, b1: parameter on marginal poisson model: V=exp(b0+b1*trt)
    y <- data$HFD
    if(is.character(data$trt)|is.factor(data$trt)){
      t <- ifelse(data$trt=="cntl",0,1)
    }else{t=data$trt}
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
  
  require(pscl)
  
  ### initialization ###
  fit1.1 <- zeroinfl(HFD ~ trt,data=data,dist=family)
  fit.zero <- summary(fit1.1)$coefficients$zero
  fit.count <- summary(fit1.1)$coefficients$count
  fit.p0 <- exp(fit.zero[1,"Estimate"])/(1+exp(fit.zero[1,"Estimate"]))
  fit.mu0 <- exp(fit.count[1,"Estimate"])
  fit.p1 <- exp(sum(fit.zero[,"Estimate"]))/(1+exp(sum(fit.zero[,"Estimate"])))
  fit.mu1 <- exp(sum(fit.count[1:2,"Estimate"]))
  HFD0 <- (1-fit.p0)*fit.mu0
  HFD1 <- (1-fit.p1)*fit.mu1
  
  r0= fit.zero[1,1]; r1=fit.zero[2,1]
  b0 <- log(HFD0) ;b1=log(HFD1) - log(HFD0)
  if(family=="negbin"){
    phi=1/exp(fit.count[3,1])
  }
  
  if(length(inits)==1){
    if(inits==1){ 
      parameter <- switch(family,"negbin"=c(r0,r1,b0,b1,phi),"poisson"=c(r0,r1,b0,b1))
      
    }else if(inits==0){
      parameter <- switch(family,"negbin"=rep(0,5),"poisson"=rep(0,4))
    }
  }
  
  lMZI <- switch(family,"negbin"=lMZINB,"poisson"=lMZIP)
  
  ## initial fit ##
  fit <- optim(par=parameter,lMZI,dat=data,method=method,hessian=TRUE)
  
  
  control <- regression.control(fit)
  
  if(control!=2){
    a=r0/try;b=r1/try;c=b0/try;d=b1/try
    while(control!=2){
      r0=switch(paste0(inits),"1"=r0-a,"0"=r0+a)
      r1=switch(paste0(inits),"1"=r1-a,"0"=r1+a)
      b0=switch(paste0(inits),"1"=b0-a,"0"=b0+a)
      b1=switch(paste0(inits),"1"=b1-a,"0"=b1+a)
      
      parameter <- switch(family,"negbin"=c(r0,r1,b0,b1,phi),"poisson"=c(r0,r1,b0,b1))
      
      fit <- optim(par=par,lMZI,dat=data,method=method,hessian=TRUE)
      control <- regression.control(fit)
    }
  }
  
  
  
  est <- fit$par #r0,r1,dispersion,b0,b1: the marginal HFD model:log(HFD)=b0+b1*trt
  var.cov <- solve(fit$hessian)
  HFD.cntl <- exp(est[3])
  HFD.trt <- exp(sum(est[3:4]))
  HFD.cntl.var <- exp(est[3])^2*var.cov[3,3]
  HFD.trt.var <- exp(sum(est[3:4]))^2*(var.cov[3,3]+var.cov[4,4]+2*var.cov[3,4])
  diff.var=HFD.cntl.var+HFD.trt.var
  
  
  fit$HFD <- c("trt"=HFD.trt,"cntl"=HFD.cntl,"trt-cntl"=HFD.trt-HFD.cntl,"diff.var"=diff.var)
  names(fit$par)[1:4] <- c("r0","r1","beta0","beta1")
  return(fit)
}   


two.part <- function(HFD,trt,data,family=c("poisson","negbin")){
  
  if(is.null(data[,HFD])|!is.numeric(data[,HFD])){
    stop("HFD (Hospital free days) is not fund or is not numeric")
  }
  if(is.null(data[,trt])){
    stop("trt (treatment) is not fund or is not numeric")
  }
  
  if(! family %in% c("poisson","negbin")){
    stop("family not recognized")
  }
  
  var_tp<-function(par_count,par_zero,var){
    a<-exp(sum(par_count))
    b<-exp(par_count[1])
    c<-exp(sum(par_zero))
    d<-exp(par_zero[1])
    v<-matrix(c(a/(1+c)-b/(1+d),
                a/(1+c),
                -c/((1+c)^2)*a+d/((1+d)^2)*b, 
                -c/((1+c)^2)*a), ncol=4)
    var<-v%*%var%*%t(v)
    return(var)
  }
  
  var_tp.r<-function(par_count,par_zero,var,m){
    a<-exp(sum(par_count))
    b<-exp(par_count[1])
    c<-exp(sum(par_zero))
    d<-exp(par_zero[1])
    v<-matrix(c(-a/(1+c)+b/(1+d),
                -a/(1+c),
                -c/((1+c)^2)*(m-a)+d/((1+d)^2)*(m-b), 
                -c/((1+c)^2)*(m-a)), ncol=4)
    var<-v%*%var%*%t(v)
    return(var)
  }
  
  cl <- match.call()
  family <- match.arg(family)
  
  data$HFD0 <- ifelse(data$HFD==0,1,0)
  
  require(multcomp)
  fit1 <- glm(HFD0 ~ trt,data=data,family="binomial")
  fit.zero <- summary(fit1)$coefficients
  fit.p.cntl <- exp(fit.zero[1,"Estimate"])/(1+exp(fit.zero[1,"Estimate"]))
  
  fit.logit.trt <- confint(glht(fit1, linfct = c("(Intercept) + trt = 0")))
  fit.p.trt <- exp(fit.logit.trt$confint[[1]])/(1+exp(fit.logit.trt$confint[[1]]))

  # Analysis approach2: Fit two part model- poisson 
  fit2 <- switch(family,"poisson"=glm(HFD~trt,data=data[which(data$HFD0==0),],family="poisson"),
                         "negbin"=glm.nb(HFD~trt,data=data[which(data$HFD0 ==0),]))
  fit.count <- summary(fit2)$coefficients
  fit.mu.cntl <- exp(fit.count[1,"Estimate"]) # X2
  fit.count.trt <- confint(glht(fit2, linfct = c("(Intercept) + trt = 0")))
  fit.mu.trt <- exp(fit.count.trt$confint[[1]]) # X2
  
  
  # calculate HFD
  HFD.cntl <- (1-fit.p.cntl)*fit.mu.cntl
  HFD.trt <- (1-fit.p.trt)*fit.mu.trt
  HFD.diff <- HFD.trt - HFD.cntl
  
  # Calculate var(HFD|trt-HFD|cntl) by delta method
  var <- cbind(rbind(vcov(fit2),matrix(0,ncol=2,nrow=2)),rbind(matrix(0,ncol=2,nrow=2),vcov(fit1)))
  diff.var <- var_tp(fit.count[,"Estimate"],fit.zero[,"Estimate"],var)
  
  out <- list(coefficients=list(count=summary(fit2)$coefficients,zero=summary(fit1)$coefficients),
              HFD=c("trt"=HFD.trt,"cntl"=HFD.cntl,"trt-cntl"=HFD.diff,"var"=diff.var))
  return(out)
}
######## illustration ##########
### set up 
library(rstudioapi)
path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

sim.wide <- read.csv("simulated_data.csv")[,-1]
sim.wide <- read.csv(paste0("SimData_",j,".csv"))[,-1]
sim.wide$HFD <- get_HFD(sim.wide)
# 0 - alive and outside hosp
# 1 - alive and in hosp
# 2 - died
library(dplyr)
sim.wide %>%
  group_by(trt) %>%
  summarise(mean=mean(HFD),
            var=var(HFD))

data %>%
  group_by(trt) %>%
  summarise(mean=mean(HFD),
            var=var(HFD))

fit0 <- lMZI_wrapper(HFD="HFD",trt="trt",data=sim.wide,try=20,family="poisson",method="BFGS",inits=0)

fit1 <- lMZI_wrapper(HFD="HFD",trt="trt",data=sim.wide,try=20,family="poisson",method="Nelder-Mead",inits=1)

fit2 <- lMZI_wrapper(HFD="HFD",trt="trt",data=sim.wide,try=20,family="negbin",method="L-BFGS-B",inits=0)

fit.tp1 <- two.part(HFD="HFD",trt="trt",data=sim.wide,family="negbin")

fit.tp0 <- two.part(HFD="HFD",trt="trt",data=sim.wide,family="poisson")

