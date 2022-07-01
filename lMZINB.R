# likelihood function for MZINB - shared parameter 
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