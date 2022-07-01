# fit marginalized zero inflated model (poisson or negative binomial)
lMZI_wrapper = function(HFD,trt,data,try=20,family="poisson",inits=1){
  cl <- match.call()
  
  if(is.null(data$HFD)|!is.numeric(data$HFD)){
    stop("HFD (Hospital free days) is not fund or is not numeric")
  }
  if(is.null(data$trt)){
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
