#Survival
S_mix <- function(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2){
  (p*exp(- lambda1 * t^(gamma1)*exp(X*beta1*gamma1)) + (1-p)*exp(- lambda2 * t^(gamma2)*exp(X*beta1*gamma2)))^(exp(beta2*X))
}
#simulate data from a mixture Weibull GH model
sim_mix_weib_gh<- function(n,p,lambda1,lambda2,gamma1,gamma2,beta1,beta2,cens_prob=0.1,tmax=10){
  X <- rbinom(n, 1, 0.5)
  U <- runif(n)
  T_true <- numeric(n)
  for(i in 1:n) {
    target <- function(t) S_mix(t, p,lambda1,lambda2,gamma1,gamma2,X[i],beta1,beta2) - U[i]
    # Search from 0 to a very large number (e.g., 100) to ensure we find the tail
    T_true[i] <- uniroot(target, interval = c(0, 100), extendInt = "yes")$root
  }
  t_cens = runif(n,min=0,max=tmax) 
  T = pmin(T_true,t_cens)
  cens = ifelse(T_true > T, 1, 0)
  return(data.frame(time=T,X=X,event=1-cens,T_true=T_true))
}

# Cumulative Hazard
H_mix <- function(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2){
  -log(S_mix(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2))
}

# hazard
h_mix_aft <- function(t,p,lambda1,lambda2,gamma1,gamma2,X,beta){
  # H_t = H_mix_aft(t,p,lambda1,lambda2,gamma1,gamma2,X,beta)
  S_t = S_mix(t,p,lambda1,lambda2,gamma1,gamma2,X,beta,0)
  return( (p*lambda1*gamma1*t^(gamma1-1)*exp(X*beta*gamma1)*exp(- lambda1 * t^(gamma1)*exp(X*beta*gamma1)) +
             (1-p)*lambda2*gamma2*t^(gamma2-1)*exp(X*beta*gamma2)*exp(- lambda2 * t^(gamma2)*exp(X*beta*gamma2)) ) / S_t )
}
h_mix <- function(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2){
  h_mix_aft(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1) * exp(beta2*X)
}

mixWeib = function(res = "S", t, p, lambda1, lambda2, gamma1, gamma2, X, beta1, beta2) {
  if(res == "S") {
    return(S_mix(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2))
  } else if(res == "H") {
    return(H_mix(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2))
  } else if(res == "h") {
    return(h_mix(t,p,lambda1,lambda2,gamma1,gamma2,X,beta1,beta2))
  } else {
    stop("Invalid result type")
  }
}

mixWeibSz = function(szenario, res = "S", t, X=0, beta1=0, beta2=0) {
  if(szenario == 1) {
    return(mixWeib(res, t, p=0.8, lambda1=0.1, lambda2=0.1, gamma1=3, gamma2=1.6, X=X, beta1=beta1, beta2=beta2))
  } else if(szenario == 2) {
    return(mixWeib(res, t, p=0.5, lambda1=1, lambda2=1, gamma1=1.5, gamma2=0.5, X=X, beta1=beta1, beta2=beta2))
  } else if(szenario == 3) {
    return(mixWeib(res, t, p=0.26, lambda1=0.02, lambda2=0.5, gamma1=3, gamma2=0.7, X=X, beta1=beta1, beta2=beta2))
  } else {
    stop("Invalid szenario number")
  }
}

sim_szenario = function(szenario, beta1, beta2, cens_prob = 0.1, n = 1000, tmax=10) {
  if(szenario == 1) {
    return(sim_mix_weib_gh(n=n,p=0.8,lambda1=0.1,lambda2=0.1,gamma1=3,gamma2=1.6,beta1=beta1,beta2=beta2,cens_prob=cens_prob,tmax=tmax))
  } else if(szenario == 2) {
    return(sim_mix_weib_gh(n=n,p=0.5,lambda1=1,lambda2=1,gamma1=1.5,gamma2=0.5,beta1=beta1,beta2=beta2,cens_prob=cens_prob,tmax=tmax))
  } else if(szenario == 3) {
    return(sim_mix_weib_gh(n=n,p=0.26,lambda1=0.02,lambda2=0.5,gamma1=3,gamma2=0.7,beta1=beta1,beta2=beta2,cens_prob=cens_prob,tmax=tmax))
  } else {
    stop("Invalid szenario number")
  }
}