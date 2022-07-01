# likelihood function for MZIP - shared parameter 
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