# Imports needed in the masterfile
library(survPen)
library(dplyr)
library(survival) #surv
library(pracma) #gausslegendre
library(dplyr)
library(splines) #ns)

# C++ spline implementation
library(Rcpp)
library(RcppArmadillo)
sourceCpp(file="crs_cpp.cpp")

# ----

#' Main work horse function, calculates various quantities depending on the res argument, including the negative log likelihood and its gradient 
#' and hessian, as well as the effective degrees of freedom and the LCV criterion for optimizing lambda.
#' Also called to fit a model. The function is quite versatile, but should not be called directly by the user, bur rather through "post" or "fit_genhaz".
#' 
#' @param theta <- vector of parameters, including spline coefficients and covariate effects, in the order: 
#'  intercept, spline coefficients, beta1 (for AFT effects), beta2 (for PH effects). Used as init if result is "fit" or "fit_LCV".
#' @param time <- vector of observed times for fitting or desired times for post estimation
#' @param X <- matrix of covariates, used for both aft and ph effects
#' @param knots <- vector of knots for the spline basis, including the boundary knots
#' @param Z <- projection matrix for the spline basis, required for post estimation
#' @param event <- event indicator vector, required for fitting
#' @param lambda <- smoothing parameter for the spline coefficients
#' @param res <- character string specifying the desired result, one of "log_h","gradient_log_h","hessian_log_h","h",
#'  "gradient_h","hessian_h","H","gradient_H","hessian_H","negll","negll_upen","scores","gradient_negll","hessian_negll",
#'  "hessian_negll_upen","basis","edf","LCV","fit","fit_LCV". 
#'  The last two fit a model, "fit" for a given lambda, "fit_LCV" optimized lambda accordng to modified LCV criterion and returns it's fit.
#' @param model_type <- character vector specifying the model type for each covariate, one of "GH", "PH", "AFT", "AH". 
#'  If length 1, repeated for all covariates. GH is the default and corresponds to the full general hazard model, 
#'  PH corresponds to proportional hazards (beta1=0), AFT corresponds to accelerated failure time (beta2=0) and AH corresponds to additive hazards (beta1=-beta2).
#' @param time2 <- vector of second time variable, required for left truncation (time is event time, time2 is truncation time) and interval censoring (time is t1, time2 is t2)
#' @param cens_type <- type of censoring, one of "rc" for right censoring, "lt_rc" for left truncation and right censoring, "ic" for interval censoring. Default is "rc".
#' @param gaussified <- whether to use Gaussian quadrature for the integration in H and its derivatives, or to use numerical integration. Default is TRUE, as numerical integration is very slow.
#' @param nquad <- number of quadrature points to use if gaussified is TRUE. Default is 25, which should be sufficient for most applications.
#' @param control <- list of control parameters for the optimization, including trace, x.tol, rel.tol and max_iter. Default is list(trace=0, x.tol=1e-12, rel.tol=1e-12, max_iter=100).
#' @param interval <- interval for optimizing lambda in LCV, only used if res is "fit_LCV". Default is c(0,13), which corresponds to lambda between 1 and exp(13) ~ 442413.4, a wide range that should cover most reasonable values of lambda.
#' @param tol_LCV <- tolerance for optimizing lambda in LCV, only used if res is "fit_LCV". Default is 1e-3, which should be sufficient for most applications.
#' @param ... 
#'
#' @returns <- see res. If res is "fit" or "fit_LCV", returns a list with the fitted model parameters, standard errors, z-values, p-values, AIC, effective degrees of freedom, lambda and other relevant information.
#' @export 
#'
#' @examples
genhaz_work = function(theta,time,X,knots,Z,event=NULL, lambda=1,
                  res=c("log_h","gradient_log_h","hessian_log_h",
                        "h","gradient_h","hessian_h",
                        "H","gradient_H","hessian_H",
                        "negll","negll_upen","scores","gradient_negll","hessian_negll","hessian_negll_upen",
                        "basis","edf","LCV",
                        "fit","fit_LCV"),
                  model_type = "GH",time2=NULL,cens_type="rc",
                  gaussified = TRUE, nquad = 25, 
                  control=list(trace=0, x.tol=1e-12, rel.tol=1e-12, max_iter = 100), interval=c(0,13),tol_LCV = 1e-3,
                  ...) {
  
  
  res = match.arg(res)
  #model_type = match.arg(model_type)
  nb0 = length(knots)  #intercept plus spline coefficents (n_knots - 1 + 1)
  nb1 = (length(theta) - nb0) / 2
  nb2 = nb1
  
  #check censoring type
  if(!(cens_type %in% c("rc","ic","lt_rc"))){
    stop("Invalid censoring type, supported are rc,ic,lt_rc")
  }
  #for rt time is time as usual, for lt_rc time is time and time2 is truncation time, for ic time is t1 and time2 is t2
  #check if X is vector and make it matrix if necessary
  if(is.null(dim(X))) {
    X = matrix(X, ncol = 1)
  }
  
  #some checks to insure input validity 
  if (nrow(X) != length(time)) {
    if(nrow(X) == 1){
      X = X[rep(1, length(time)), , drop = FALSE]
    }
    else{
      stop("Number of rows in X and length of time must be the same.")
    }

  }
  
  allowed_models <- c("GH", "PH", "AFT", "AH")
  
  if (length(model_type) == 1) {
    model_type <- rep(model_type, nb1)  # repeat for all covariates
  } else if (length(model_type) != nb1) { #also assumes length(beta1) == length(beta2)
    stop("Length of 'model_type' must be 1 or equal to the number of elements in beta1.")
  }
  
  # Check validity of model type
  if (!all(model_type %in% allowed_models)) {
    stop("All elements of 'model_type' must be one of: 'GH', 'PH', 'AFT', 'AH'")
  }  
  
  is_gh = ifelse(model_type == "GH", 1, 0)
  is_ph = ifelse(model_type == "PH", 1, 0)
  is_aft = ifelse(model_type == "AFT", 1, 0)
  is_ah = ifelse(model_type == "AH", 1, 0)
  
  is_not_gh = ifelse(model_type == "GH",  0, 1)
  is_not_ph = ifelse(model_type == "PH",  0, 1)
  is_not_aft = ifelse(model_type == "AFT",0,1)
  is_not_ah = ifelse(model_type == "AH",  0,1)
  
  X1_filtered = X * matrix(is_not_ph, #ifelse(model_type == "PH",0,1),
                           nrow = nrow(X), ncol = nb1, byrow = TRUE)
  X2_filtered = X * matrix(is_gh + is_ph,  #ifelse(model_type == "GH" | model_type == "PH",1,0),
                           nrow = nrow(X), ncol = nb2, byrow = TRUE)
  X3_filtered = X1_filtered * matrix(is_not_ah, #ifelse(model_type == "AH",0,1), 
                                     nrow = nrow(X), ncol = nb1, byrow = TRUE)
  
  filter = is_not_ph #ifelse(model_type=="PH",0,1)
  filter1 = is_gh + is_aft  #ifelse(model_type == "GH" | model_type == "AFT",1,0)
  filter2 = is_gh + is_ph #ifelse(model_type == "GH" | model_type == "PH",1,0)
  eta1 = function(theta) {
    beta1 = theta[nb0+1:nb1]
    return(drop(X %*% (beta1 * is_not_ph)))
  }
  eta2 = function(theta) {
    beta1 = theta[nb0+1:nb1]
    beta2 = theta[nb0+nb1+1:nb2]
    return(drop(X %*% (beta1 * filter1 + beta2 * filter2)))
  }
  cc = gaussLegendre(nquad, -1, 1)
  
  smf_base = smf_cpp(log(time+0.0000000001), knots=knots, Z=Z, intercept=FALSE) #recheck this
  S_penal_wo_intercept = attr(smf_base,"pen")
  
  #normalize S_penal_wo_intercept to have trace equal to nb0-1 (number of splines)
  S_penal_wo_intercept = S_penal_wo_intercept * (nb0-1) / sum(diag(S_penal_wo_intercept))
  
  #add 0 row and column for intercept
  S_penal = matrix(0, nrow=nrow(S_penal_wo_intercept)+1, ncol=ncol(S_penal_wo_intercept)+1)
  S_penal[-1,-1] <- S_penal_wo_intercept
  
  B = function(time,...) smf_cpp(time, knots=knots, Z=Z, intercept=FALSE, ...) # still time as this is a function of time

  log_h = function(theta,time) {
    beta0 = theta[2:nb0]
    drop(B(log(time+0.0000000001) + eta1(theta)) %*% beta0 + theta[1] + eta2(theta))
    #drop(B(time*exp(eta1(theta))) %*% beta0  + eta2(theta))
    ## drop(B(time*exp(X %*% beta1)) %*% beta0 + X %*% beta1 + X %*% beta2)
  }
  h = function(theta,time) exp(log_h(theta,time))
  Quadrature = function(f,theta,time) {
    ## Gaussian quadrature for \int_0^{time} f(\theta,t) dt
    Integral = 0
    for(i in 1:nquad) {
      wi = cc$w[i]*time/2 # vector
      ti = (cc$x[i]+1)*time/2 # vector
      Integral = Integral + wi*f(theta,ti)
    }
    Integral
  }
  gradient_log_h = function(theta,time) {
    beta0 = theta[2:nb0]
    eta1 = eta1(theta)
    Bmat = B(log(time+0.0000000001) + eta1)  # changed name to Bmat to avoid a name conflict
    B_prime = B(log(time+0.0000000001) + eta1, derivs=1)
    db0 = Bmat
    db1 = drop(B_prime %*% beta0) * X1_filtered + X3_filtered
    db2 = X2_filtered
    ## bind together columnwise
    cbind(1, db0, db1, db2)
  }
  gradient_h = function(theta,time)
    h(theta,time)*gradient_log_h(theta,time)
  hessian_log_h = function(theta,time) {
    beta0 = theta[2:nb0]
    eta1 = eta1(theta)
    ## Bmat = B(time*exp(eta1))
    B_prime = B(log(time+0.0000000001) + eta1, derivs=1)
    B_prime_2 = B(log(time+0.0000000001) + eta1, derivs=2)
    ## db0 = Bmat # function of beta1
    ## db1 = drop(B_prime %*% beta0) * time * exp(eta1) * X1 + X1
    ## db2 = X1
    Hessian = array(0, c(length(time),length(theta),length(theta)))
    ## db00 = db02 = db12 = db20 = db21 = db22 = 0
    index0 = 2:nb0; index1 = nb0+1:nb1; index2 = nb0+nb1+1:nb2
    iA = index0; iB = index1
    for (k in 1:length(iB)) {
      Hessian[,iA,iB[k]] = Hessian[,iB[k],iA] = B_prime * X1_filtered[,k]
    }
    
    iA = index1
    for (k in 1:length(iA)) {
      Hessian[,iA,iA[k]] = drop(B_prime_2 %*% beta0) * X1_filtered*X1_filtered[,k] + 
        (drop(B_prime %*% beta0) * time * exp(eta1) * X1_filtered*X1_filtered[,k])*0
    }
    return(Hessian)
  }
  hessian_h = function(theta,time) {
    h = h(theta,time) # n*1
    dlogh = gradient_log_h(theta,time) # n*m
    d2logh = hessian_log_h(theta,time) # n*m*m
    m = ncol(dlogh)
    ## dlogh2 = t(apply(dlogh, 1, function(x) outer(x,x))) # slow
    dlogh2 = dlogh[, rep(1:m, each=m)] * dlogh[, rep(1:m, times=m)]
    dlogh2 = array(dlogh2, dim(d2logh))
    Hessian = h*(dlogh2 + d2logh)
    ## Hessian = apply(Hessian, 2:3, function(x) x*h)
    return(Hessian)
  }
  H = function(theta,time) Quadrature(h,theta,time)
  gradient_H = function(theta,time) Quadrature(gradient_h,theta,time)
  hessian_H = function(theta,time) Quadrature(hessian_h,theta,time)
  negll = function(theta,time,penalized=TRUE,...) {
    logh = log_h(theta,time)
    ## recursive call to allow for different calculations of H
    # H = genhaz_work(theta,time,X,knots,Z,res="H", gaussified=gaussified, nquad=nquad,
    #            model_type=model_type)
    H1 = H(theta,time)
    penal = 0
    if(penalized){
      penal = lambda/2 * t(theta[1:nb0]) %*% S_penal %*% theta[1:nb0]
    }
    if(cens_type=="rc"){
      return(-sum(logh*event - H1) + penal)
    }
    
    # H2 = genhaz_work(theta,time2,X,knots,res="H", gaussified=gaussified, nquad=nquad,
    #             model_type=model_type)
    H2 = H(theta,time2)
    
    if(cens_type=="lt_rc"){
      return(-sum(logh*event - H1 + H2) + penal)
    }
    if(cens_type=="ic"){
      return(-sum(log(exp(-H1)-exp(-H2))) + penal)
    }
    
  }
  gradient_negll = function(theta,time,penalized=TRUE,...) {
    ## assumes Gaussian quadrature
    dlogh = gradient_log_h(theta,time)
    dH = gradient_H(theta,time)
    penal = 0
    if(penalized){
      penal = c(lambda * S_penal %*% theta[1:nb0],rep(0,2*nb1))
    }
    if(cens_type=="rc"){
      return(-colSums(dlogh * event - dH) + penal)
    }
    
    dH2 = gradient_H(theta,time2)
    if(cens_type=="lt_rc"){
      return(-colSums(dlogh * event - dH +dH2) + penal)
    }
    
    # H = genhaz_work(theta,time,X,knots,Z,res="H",lambda = lambda, gaussified=gaussified, nquad=nquad,
    #            model_type=model_type)
    # H2 = genhaz_work(theta,time2,X,knots,Z,res="H",lambda = lambda, gaussified=gaussified, nquad=nquad,
    #             model_type=model_type)
    H1 = H(theta,time)
    H2 = H(theta,time2)
    if(cens_type=="ic"){
      s = (exp(-H2)*dH2 - exp(-H1)*dH) / (exp(-H1) - exp(-H2))
      return(-colSums(s) + penal)
    }
    
  }
  scores = function(theta,time) {
    ## assumes Gaussian quadrature
    dlogh = gradient_log_h(theta,time)
    dH = gradient_H(theta,time)
    
    if(cens_type=="rc"){
      return(dlogh * event - dH)
    }
    
    dH2 = gradient_H(theta,time2)
    if(cens_type=="lt_rc"){
      return(dlogh * event - dH +dH2)
    }
    
    # H = genhaz_work(theta,time,X,knots,Z,res="H",lambda = lambda, gaussified=gaussified, nquad=nquad,
    #            model_type=model_type)
    # H2 = genhaz_work(theta,time2,X,knots,Z,res="H",lambda = lambda, gaussified=gaussified, nquad=nquad,
    #             model_type=model_type)
    H1 = H(theta,time)
    H2 = H(theta,time2)
    if(cens_type=="ic"){
      s = (exp(-H2)*dH2 - exp(-H1)*dH) / (exp(-H1) - exp(-H2))
      return(s)
    }
  }
  hessian_negll = function(theta,time,penalized = TRUE,...) {
    ## assumes Gaussian quadrature
    d2logh = hessian_log_h(theta,time)
    d2H = hessian_H(theta,time)
    
    penal = 0
    if(penalized){
      penal = rbind(
        cbind(lambda*S_penal, matrix(0, nrow = nb0, ncol = 2*nb1)), #think about deriv by log(lambda) for lcv for optimizing lambda -> fauvernier
        matrix(0, nrow = 2*nb1, ncol = nb0 + 2*nb1)
      )
    }
    if(cens_type=="rc"){
      s = d2logh * event - d2H # careful with array dimensions
      return(-apply(s,2:3,sum) + penal)
    }
    
    d2H2 = hessian_H(theta,time2)
    if(cens_type=="lt_rc"){
      s = d2logh * event - d2H + d2H2 # careful with array dimensions
      return(-apply(s,2:3,sum) + penal)
    }
    
    # H = genhaz_work(theta,time,X,knots,Z,res="H",lambda = lambda, #check if direct H call correct and faster
    #            gaussified=gaussified, nquad=nquad,
    #            model_type=model_type)
    # 
    # H2 = genhaz_work(theta,time2,X,knots,Z,res="H",lambda = lambda,
    #           gaussified=gaussified, nquad=nquad,
    #             model_type=model_type)
    H1 = H(theta,time)
    H2 = H(theta,time2)
    dH = gradient_H(theta,time)
    dH2 = gradient_H(theta,time2)
    S1 = exp(-H1)
    S2 = exp(-H2)
    if(cens_type=="ic"){ #unholy fumbling
      v = (S2*dH2 - S1*dH) # n x p
      a = array(rep(v, times = length(theta)), dim = c(length(time), length(theta), length(theta))) # n x p x p where last dim rep
      b = aperm(a,c(1,3,2)) # n x p x p where last dim rep but transposed
      p1 = a * b / (S1-S2)^2 # outer product of [(.) x p] [(.) x p] and normalize 
      
      v = dH # again n x p, same procedure
      a = array(rep(v, times = length(theta)), dim = c(length(time), length(theta), length(theta)))
      b = aperm(a,c(1,3,2))
      dH_2 = a * b
      
      v = dH2
      a = array(rep(v, times = length(theta)), dim = c(length(time), length(theta), length(theta)))
      b = aperm(a,c(1,3,2))
      dH2_2 = a * b
      
      p2 = (S1 * d2H - S1 * dH_2 - S2 * d2H2 + S2 * dH2_2)/(S1-S2)
      
      s = p2-p1
      
      
      return(-apply(s,2:3,sum))
    }
    
    
    
  }
  ##
  if (res=="H") {
    # vectorised = TRUE
    # if(!vectorised && !gaussified) {
    #   eta1 = eta1(theta)
    #   eta2 = eta2(theta)
    #   beta0 = theta[1:nb0]
    #   return(sapply(1:length(time), function(i) {
    #     h = function(t) {
    #       exp(B(log(t+0.0000000001) + eta1[i]) %*% beta0 + eta2[i]) 
    #     }
    #     integrate(h, 0, time[i])$value
    #   }))
    # }
    # else if (vectorised && !gaussified) {
    #   inner = function(t) h(theta,t)
    #   return(rstpm2:::vintegrate(inner, rep(0,length(time)), time)$value)
    # }
    # else return(H(theta,time))
    if (!gaussified) {
        inner = function(t) h(theta,t)
        return(rstpm2:::vintegrate(inner, rep(0,length(time)), time)$value)
    }
    else return(H(theta,time))
    
  }
  
  
  edf = function(theta_lambda,roh){
    #modvec = c(rep(0,nb0),ifelse(model_type == "PH",1,0),ifelse(model_type == "AFT"| model_type == "AH",1,0)) #submatrix to Id for forced 0 params
    modvec = c(rep(0,nb0),is_ph,is_aft+is_ah) #submatrix to Id for forced 0 params
    lambda = exp(roh)
    #effective degress of freedom
    tryCatch({
      #edf matrix is (t(B)B + lambda * S_penal)^-1 t(B)B
      H_pen <- genhaz_work(theta=theta_lambda, time=time, X=X, knots=knots, Z=Z, event=event, lambda = exp(roh), 
                           cens_type=cens_type,model_type=model_type,time2=time2, res = "hessian_negll", control=control)
      H_upen <- genhaz_work(theta=theta_lambda, time=time, X=X, knots=knots, Z=Z, event=event, lambda = exp(roh), 
                            cens_type=cens_type,model_type=model_type,time2=time2, res = "hessian_negll_upen", control=control)
      edf = Trace(solve(H_pen + diag(modvec)) %*% H_upen) - sum(is_not_gh) #correct degrees of freedom for submodels
      
    }, error = function(e){
      print(paste0("Error in computing edf for roh=", roh, ": ", e))
      #also print hessians
      print("Penalized hessian:")
      print(hessian_negll(theta_lambda,time,penalized = TRUE,lambda=lambda))
      print("theta_lambda:")
      print(theta_lambda)
      edf = Inf
    })
    return(edf)
  }
  LCV_expl = function(roh,theta) {
    lambda = exp(roh)
    
    theta_lambda = theta
    init <<- theta_lambda # update global init for next optimization step, should speed up convergence
    edf = edf(theta_lambda,roh)
    
    penalty = (lambda / 2) * t(theta_lambda[1:nb0]) %*% S_penal %*% theta_lambda[1:nb0]
    unpen_negll = negll(theta,time) - as.numeric(penalty)
    
    #the neg ll needs to be unpenalized here!
    return(unpen_negll + edf * log(length(time))/2)
  }
  if (res=="fit") {
    ## fit = optim(theta, negll, gradient_negll, method="BFGS", ..., time=time)
    lower = rep(-Inf, length(theta))
    upper = rep(Inf, length(theta))
    for (i in 1:nb1) {
      if (model_type[i] == "PH") {
        lower[nb0+i] = upper[nb0+i] = 0
      }
      if (model_type[i] == "AH" || model_type[i] == "AFT") {
        lower[nb0+nb1+i] = upper[nb0+nb1+i] = 0
      }
    }
    if (!exists("init", inherits = FALSE) || is.null(init)) {
      init <- theta
    }
    fit = nlminb(init, negll, gradient_negll, hessian_negll, lower=lower, upper=upper, time=time, penalized = TRUE)

    fit$par = fit$par * c(rep(1,nb0),
                          is_not_ph,
                          is_ph+is_gh) +
      c(rep(0,nb0+nb1), -is_ah*fit$par[(nb0+1):(nb0+nb1)])
    
    fit$parnames = c("intercept",paste0("s",1:(nb0-1)), paste0("beta1_",1:nb1), paste0("beta2_",1:nb2))
    names(fit$par) = fit$parnames
    
    
    fit$hessian = -hessian_negll(fit$par,time) #observe the minus here
    rownames(fit$hessian) = colnames(fit$hessian) = fit$parnames
    
    modvec = c(rep(0,nb0),is_ph,is_aft+is_ah) #submatrix to Id for forced 0 params
    # fit$hessmod = -fit$hessian + diag(modvec) # for debugging
    tryCatch({
      h_pen_inv = solve(-fit$hessian + diag(modvec)) - diag(modvec)
    }, error = function(e){
      # print("Error in inverting hessian:")
      # print(e)
      stop("Stopping due to non inveritible Hessian")
    })
    fit$var = h_pen_inv
    fit$se = sqrt(diag(fit$var))
    fit$z = fit$par/fit$se
    fit$p_values <- 2 * (1 - pnorm(abs(fit$z)))
    fit$df = length(theta) - sum(is_not_gh) # #knots + 1 + 1(for intercept) for the splines plus the number of non zero params in beta1,beta2
    fit$AIC = 2*fit$objective + fit$df
    #fit$basis = basis
    fit$model_type = model_type
    fit$knots = knots
    
    fit$Z = Z
    fit$S = S_penal
    fit$negll = fit$objective
    
    fit$edf = edf(fit$par, log(lambda))
    fit$negll = fit$objective
    fit$lambda = lambda
    fit$penalty = (lambda / 2) * t(fit$par[1:nb0]) %*% S_penal %*% fit$par[1:nb0]
    fit$unpen_negll = fit$negll - as.numeric(fit$penalty)
    
    fit$LCV = fit$unpen_negll + fit$edf * log(length(time))/2
    fit$n = length(time)
    fit$nb0 = nb0
    fit$nb1 = nb1
    #fit$X_orig = X_orig
    fit$X_center = attr(X, "scaled:center")
    
    # names(fit$par) = names
    names(fit$se) = fit$parnames
    # names(fit$p_values) = names
    
    return(fit)
  }

  LCV = function(roh) {
    lambda = exp(roh)
    
    #fit model for given lambda
    fit_lambda = genhaz_work(init,time,X,knots,Z,
                        event=event, lambda=exp(roh),
                        res="fit",
                        model_type = model_type,
                        time2=time2,cens_type=cens_type,
                        vectorised=vectorised, gaussified = gaussified, nquad = nquad, 
                        control=control)
    
    theta_lambda = fit_lambda$par
    init <<- theta_lambda # update global init for next optimization step, should speed up convergence
    edf = edf(theta_lambda,roh)
    
    # unpen_negll = genhaz_work(theta_lambda,time,X,knots,Z,event=event, lambda=lambda,
    #                     res="negll_upen",model_type=model_type,time2=time2,cens_type=cens_type,
    #                     vectorised=vectorised, gaussified = gaussified, nquad = nquad, control=control)

    penalty = (lambda / 2) * t(theta_lambda[1:nb0]) %*% S_penal %*% theta_lambda[1:nb0]
    unpen_negll = fit_lambda$negll - as.numeric(penalty)

    #the neg ll needs to be unpenalized here !!!!!!!!!!!
    return(unpen_negll + edf * log(length(time))/2)
  }
  if( res=="fit_LCV") {
    # assign init <<- theta, so that the inner fit in optimize will start with updated theta values, which should speed up convergence
    init <<- theta
    LCV_opt = optimize(LCV, interval = interval, tol=tol_LCV)
    lambda_opt = exp(LCV_opt$minimum)
    
    fit = genhaz_work(theta,time,X,knots,Z,
                      event=event, lambda=lambda_opt,
                      res="fit",
                      model_type = model_type,
                      time2=time2,cens_type=cens_type,
                      vectorised=vectorised, gaussified = gaussified, nquad = nquad, 
                      control=control)
    fit$lambda = lambda_opt
    fit$roh_opt = LCV_opt$minimum
    fit$LCV = LCV_opt$objective
    fit$edf = edf(fit$par,fit$roh_opt)
    
    return(fit)
    
  }
  
  out = switch(res,
               basis=B(log(time+0.0000000001) + eta1(theta)),
               log_h=log_h(theta,time),
               gradient_log_h=gradient_log_h(theta,time),
               hessian_log_h=hessian_log_h(theta,time),
               h=h(theta,time),
               gradient_h=gradient_h(theta,time),
               hessian_h=hessian_h(theta,time),
               negll=negll(theta,time,penalized=TRUE),
               negll_upen=negll(theta,time,penalized=FALSE),
               gradient_negll=gradient_negll(theta,time,penalized=TRUE),
               scores=scores(theta,time),
               hessian_negll=hessian_negll(theta,time,penalized=TRUE),
               hessian_negll_upen=hessian_negll(theta,time,penalized=FALSE),
               ## H: see above
               gradient_H=gradient_H(theta,time),
               hessian_H=hessian_H(theta,time),
               LCV = LCV(log(lambda)),
               edf = edf(theta,log(lambda)),
               NULL)
  if (!is.null(out)) return(out)
  
  stop("res not matched: ", res)
}


#' Function handling post estimation
#'
#' @param fit <- fitted model object from genhaz_work with res = "fit" or "fit_LCV"
#' @param X <- covariate matrix for prediction, should have same number of columns as covariates in the fitted model. 
#'  Can be a vector (for single prediction or same covariate pattern across time) or a matrix (for multiple predictions). If NULL, will be treated as a vector of zeros.
#' @param time <- vector of time points for prediction. If length(time) > 1 and X is a single row, the single row of X will be repeated for each time point. 
#'  If length(time) == 1 and X has multiple rows, the single time point will be used for all rows of X.
#' @param res <- type of result to return, same options as in genhaz_work except "fit" and "fit_LCV"
#' @param event <- optional vector of event indicators for each time point, needed for certain types of predictions (e.g. scores) and for correct handling of censoring types. Should have same length as time.
#'
#' @returns <- see res
#' @export
#'
#' @examples
post = function(fit,X,time,res,event=NULL) { 
  if(is.null(X)){
    X = matrix(rep(0,length(time)),ncol=fit$nb1,byrow=TRUE)
    return(genhaz_work(fit$par, time, X, knots=fit$knots,Z=fit$Z, lambda = fit$lambda,
                       event=event,model_type=fit$model_type, res = res, control=list(trace=0, x.tol=1e-12, rel.tol=1e-12)))
  } else if(!is.null(X) && !is.matrix(X)) {
    X = matrix(rep(X,length(time)),ncol=length(X),byrow=TRUE)
  } else if(!is.null(X) && is.matrix(X) && length(time)==1) {
    time = rep(time, nrow(X))
    
  } else if(!is.null(X) && is.matrix(X) && nrow(X)==1) {
    #repeat the single row of X for each time point
    X = matrix(rep(X, each=length(time)), ncol=ncol(X)) 
  }
  return(genhaz_work(fit$par, time, X, knots=fit$knots,Z=fit$Z, lambda = fit$lambda,
                     event=event,model_type=fit$model_type, res = res, control=list(trace=0, x.tol=1e-12, rel.tol=1e-12)))
}




#' Calculates the position of the knots based on event times, option to have the last knot at 100% or 95% quantile of event times. 
#' Default is boundary knots at 5% and 95%
#' @param time <- vector of time points
#' @param event <- vector of event indicators (1 for event, 0 for censored), same length as time
#' @param n_knots <- number of knots to use for the spline basis, including boundary knots. 
#' @param limit_to_95 <- logical, if TRUE (default), the last knot will be placed at the 95% quantile of event times, 
#'  if FALSE, the last knot will be placed at the 100% quantile (max) of event times.
#'
#' @returns <- vector of knot positions on the log time scale, including boundary knots
#' @export
#'
#' @examples
knot_pattern <- function(time, event, n_knots, limit_to_95 = TRUE) {
  # Determine the upper probability bound
  upper_prob <- if (limit_to_95) 0.95 else 1.0
  # Create a sequence of probabilities from 0 to the upper bound
  probs <- seq(0.05, upper_prob, length.out = n_knots)
  # Filter for events and calculate quantiles
  data.frame(time, event) %>%
    filter(event == 1) %>%
    mutate(logt = log(time)) %>%
    pull(logt) %>% 
    quantile(probs = probs)
}

# ----
#' Function for fitting a GH model
#'
#' @param surv <- a Surv object containing the survival data, can be right-censored, left-truncated and right-censored, or interval-censored. The function will automatically detect the type of censoring and handle it accordingly.
#' @param formula <- a formula specifying the covariates to include in the model. This is used for calculating the design matrix X
#' @param data <- a data frame containing the variables in the formula and the survival data. The function will add the time, time2, and event variables to this data frame for internal use.
#' @param n_knots <- number of knots to use for the spline basis for the baseline hazard. Default is 8, which includes the boundary knots. 
#' @param lambda <- smoothing parameter for the penalized spline. Default If profile is TRUE, the function will optimize lambda based on LCV.
#' @param model_type <- character vector specifying the model type for each covariate, one of "GH", "PH", "AFT", "AH". 
#'  If length 1, repeated for all covariates. GH is the default and corresponds to the full general hazard model, 
#'  PH corresponds to proportional hazards (beta1=0), AFT corresponds to accelerated failure time (beta2=0) and AH corresponds to additive hazards (beta1=-beta2).
#' @param control <- list of control parameters for the optimization.
#' @param profile <- logical, if TRUE, the function will optimize lambda based on LCV. Default is FALSE, in which case the function will fit the model for the given lambda without optimizing.
#' @param lambda_surv <- logical, if TRUE, the function will use the survPen function to get a smoothing parameter lambda. 
#' @param init <- optional vector of initial values for the optimization. If NULL, will be set to a vector of zeros with length equal to the number of parameters (intercept + spline coefficients(#knots-1) + 2 * covariate coefficients).
#' @param interval <- interval for optimizing lambda if profile is TRUE, default is c(0,10) on the log scale (i.e. lambda will be optimized between exp(0)=1 and exp(10) ~ 22026)
#' @param nquad <- number of quadrature points to use for gaussian quadrature integration, default is 25.
#' @param tol_LCV <- tolerance for optimizing lambda based on LCV, default is 1e-4. This controls the precision of the optimization of lambda when profile is TRUE.
#' @param knots <- optional vector of knot positions on the log time scale, including boundary knots. If NULL, will be calculated based on the event times and n_knots using the knot_pattern function. 
#'  If provided, should be on the log time scale and should include the boundary knots.
#' @param reKnot <- Integer, multiple fits with knots on shifted event times can be realised here. Absolute of reKnot controls amount of times the model is refit, the sign the direction. 
#'  Negative direction is for testig purposes and not recommendend. Default is reKnot = 1.
#' @param limit_to_95 <- logical, if TRUE (default), the last knot will be placed at the 95% quantile of event times, instead of 100%
#' @param timeIt <- logical, if TRUE, the function will print the execution time for fitting the model. Default is FALSE.
#'
#' @returns <- a fitted GH model object containing the parameter estimates, standard errors, p-values, and other relevant information..
#' @export
#'
#' @examples
fit_genhaz = function(surv, formula, data, n_knots=8, lambda = 0, model_type="GH", 
                      control=list(trace=0, x.tol=1e-12, rel.tol=1e-12, max_iter = 100), 
                      profile = FALSE, lambda_surv=FALSE, init = NULL, interval=c(0,10), nquad = 25, 
                      tol_LCV=tol_LCV, knots=NULL, reKnot =1, limit_to_95 = TRUE, timeIt = FALSE) {
  # f <- as.formula(paste0("surv ~", deparse(formula[[2]])))
  
  if(timeIt){
    start_time <- Sys.time()
    on.exit({
      end_time <- Sys.time()
      print(paste("Execution time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds"))
    })
  }
  
  #parse surv:
  type_surv <- attr(surv, "type")
  if(type_surv == "counting"){
    cens_type = "lt_rc"
    time = surv[, "start"]
    time2 = surv[, "start"]
    event = surv[, "status"]
  } else if(type_surv == "interval"){
    cens_type = "ic"
    time = surv[, "time1"]
    time2 = surv[, "time2"]
    event = surv[, "status"]
  } else if(type_surv == "right"){
    #print("rc")
    cens_type = "rc"
    time = surv[, "time"]
    time2 = time
    event = surv[, "status"]
    
  } else{
    stop("Sorry, this configuration of Surv is not yet supported.")
  }
  
  #print(cens_type)
  
  if(length(time) != nrow(data)) {
    stop("Length of surv and number of rows in data do not match.")
  }
  data$time = time
  data$time2 = time2
  data$event = event
  #knots <- filter(data, event==1) %>% with(quantile(time, (0:((n_knots-1))*0.95)/(n_knots-1))) #### ignore knot at 100%
  if(is.null(knots)){
    #knots <- filter(data, event==1) %>% with(quantile(time, (0:((n_knots-1)))/(n_knots-1))) #### knot at 100%
    knots <- knot_pattern(time, event, n_knots,limit_to_95 = limit_to_95)
  }
  #knots <- filter(data, event==1) %>% with(quantile(time, (0:((n_knots-1)))/(n_knots-1))) #### knot at 100%
  assign("knots_of_never_to_be_discovered_name",knots, envir = .GlobalEnv)
  
  #model matrix: THIS IS THE DESIGN MATRIX FOR THE COVARIATES ONLY, NOT THE BASELINE SPLINES!
  X <- model.matrix(formula, data=data)[,-1,drop=FALSE] #remove intercept
  
  if(is.null(init)& !lambda_surv) {
    init <- rep(0, 1 + n_knots - 1 + 2*ncol(X)) #intercept, spline coefs, covariate coefs for ph and aft/ah
  }
  
  
  if(is.null(init) & lambda_surv) {
    baseline <- ~smf(time,knots=knots_of_never_to_be_discovered_name)
    sp_ph <- as.formula(paste0("~", deparse(formula[[2]]),"+", deparse(baseline[[2]])))
    ph_mod <- survPen(sp_ph, data = data, t1 = time, event = event)
    
    #extract inital hazard ratio estimates
    #get number of covariates in formula that is like cov1 + cov2 + cov3 ... length(formula[[2]]) will return the desired number plus number of "+"
    #add -1 here if intercept will be done explicitly
    
    coef <- ph_mod$coefficients
    n_spline_coeff <- n_knots - 1 #without intercept
    n_covariates <- length(coef) - n_spline_coeff - 1 #minus intercept
    init <- c(coef[1],tail(coef,-(n_covariates+1)),rep(0,n_covariates),coef[2:(n_covariates+1)])*0 #intercept,spline coefs, 0 for aft, ph coefs
  }
  
  #Boundary.knots = filter(data, event==1) %>% with(range(time))
  Z = attr(smf_cpp(log(time+0.0000000001), knots=knots, intercept=FALSE, derivs=0L),"Z") #recheck this
  
  
  if(lambda_surv) {
    fit = genhaz_work(theta=init, time=time, X=X, knots=knots, Z=Z, event=event, lambda = ph_mod$lambda, 
                      cens_type=cens_type,model_type=model_type,time2=time2, res = "fit", control=control, interval=interval, nquad = nquad, tol_LCV=tol_LCV)
  } else {
    fit = genhaz_work(theta=init, time=time, X=X, knots=knots, Z=Z, event=event, lambda = lambda, 
                      cens_type=cens_type,model_type=model_type,time2=time2, res = ifelse(profile,"fit_LCV","fit"), control=control,
                      interval=interval, nquad = nquad, tol_LCV=tol_LCV)
  }
  
  #fit$ph_mod = ph_mod
  if(lambda_surv) {
    fit$roh_opt = log(ph_mod$lambda)
  }

  knot_hist = list()
  par_hist = list()
  
  riktning <- ifelse(reKnot>0,1,-1)
  if(reKnot<0) {
    reKnot = abs(reKnot)
  }
  
  for(i in seq_len(reKnot)) {
    knot_hist[[i]] = fit$knots
    par_hist[[i]] = fit$par
    
    #get part of model matrix X that corresponds to beta1, i.e. the covariates that are not part of the spline basis, and adapt knots accordingly
    new_X = t(t(X) * fit$par[(fit$nb0+1):(fit$nb0+fit$nb1)])
    #rowwise product of time and exp(-X %*% beta1) for the events --- PLUS ---
    new_time = time * exp(riktning*rowSums(new_X))
    #new_knots <- filter(data, event==1) %>% with(quantile(new_time, (0:((n_knots-1)))/(n_knots-1)))
    new_knots <- knot_pattern(new_time, event, n_knots, limit_to_95 = limit_to_95)
    
    fit <- genhaz_work(theta=fit$par, time=time, X=X, knots=new_knots, Z=Z, event=event, lambda = lambda, 
                       cens_type=cens_type,model_type=model_type,time2=time2, res = ifelse(profile,"fit_LCV","fit"), control=control,
                       interval=interval, nquad = nquad, tol_LCV=tol_LCV)
  }
  
  if(reKnot>0){
    fit$knot_hist = knot_hist
    fit$par_hist = par_hist
  }
  
  #extract parnames from X and assign 
  covariate_names_expanded <- colnames(X)
  parnames <- c(
    "intercept", 
    paste0("s", 1:(n_knots - 1)), 
    paste0("beta1_", covariate_names_expanded), 
    paste0("beta2_", covariate_names_expanded)
  )
  names(fit$par) <- parnames
  
  return(fit)
}



# ----
# stuff
# ----
LR = function(fit_nested, fit_general) {
  LR = -2 * (-fit_nested$objective + fit_general$objective)
  p_value <- pchisq(LR, df = fit_general$df - fit_nested$df, lower.tail = FALSE)
  res=(c(LR,p_value))
  names(res) <- c("LR-statistic","p_value")
  return(res)
}
waldCI = function(fit, param, alpha = 0.05) {
  if(!(param %in% fit$parnames)) {
    stop("param not param of fit")
  }
  lower = fit$par[param] + qnorm(alpha/2) * fit$se[param]
  names(lower) = "lower"
  upper = fit$par[param] + qnorm(1-alpha/2) * fit$se[param]
  names(upper) = "upper"
  return(c(lower,upper))
}
waldCI_minus = function(fit, param1, param2, alpha = 0.05) {
  covMatrix = fit$var[c(param1,param2),c(param1,param2)]
  
  var = covMatrix[1,1] - 2 * covMatrix[1,2] + covMatrix[2,2]
  
  lower = fit$par[param1] - fit$par[param2] + qnorm(alpha/2) * sqrt(var)
  names(lower) = "lower"
  upper = fit$par[param1] - fit$par[param2] + qnorm(1-alpha/2) * sqrt(var)
  names(upper) = "upper"
  return(c(lower,upper))
}
CI = function(fit,t,covariate,alpha=0.05) {
  X_ = matrix(rep(covariate,length(t)),ncol=length(covariate),byrow=TRUE)
  grad_H = post(fit,X_,t,"gradient_H")
  H = post(fit,X_,t,"H")
  logH = log(H)
  var_theta = fit$var
  # grad_H is (2001 × p) and var_theta is (p × p):
  var_H <- rowSums((grad_H %*% var_theta) * grad_H)
  #var_logH = var_H/(H**2)
  SE_logH = sqrt(var_H)/H
  
  grad_h = post(fit,X_,t,"gradient_h")
  h = post(fit,X_,t,"h")
  logh = log(h)
  var_h <- rowSums((grad_h %*% var_theta) * grad_h)
  SE_logh = sqrt(var_h)/h
  
  
  z <- qnorm(1-alpha/2)
  
  res = data.frame(time = t, 
                   H=H, lower_H = exp(logH + z*SE_logH), upper_H = exp(logH - z*SE_logH),
                   S=exp(-H), lower_S = exp(-exp(logH + z*SE_logH)),upper_S = exp(-exp(logH - z*SE_logH)),
                   h=h, lower_h = exp(logh + z*SE_logh), upper_h = exp(logh - z*SE_logh))
  
  # sandwich: USES OBSERVED SCORE AND HESSIAN; NEEDS THEREFORE TRUE EVENT;TIME AND COVARIATE PATTERNS 
  #incorrect as of now therfore not copied here
  
  return(res)
  
}


plot_hazard <- function(fit,covariate,time, xlim=NULL,ylim=NULL) {
  X <- matrix(rep(covariate,length(time)),ncol=length(covariate),byrow=TRUE)
  haz <- post(fit,X,time,"h")
  plot(time, haz, type='l', xlab='Time', ylab='Hazard', main='Estimated Hazard Function',xlim = xlim, ylim=ylim)
}

translate_time <- function(fit, covariate0, covariate1, time){
  X <- matrix(rep(covariate0,length(time)),ncol=length(covariate0),byrow=TRUE)
  X1 <- matrix(rep(covariate1,length(time)),ncol=length(covariate1),byrow=TRUE)
  tausimp = function(t){
    H0t = post(fit,covariate0,t,"H")
    H1 = function(u) post(fit,covariate1,u,"H")
    Hdiff = function(u) H1(u) - H0t
    return(rstpm2::vuniroot(Hdiff,lower=0,upper=max(time),tol = .Machine$double.eps)$root) #interval=c(0,max(time)))
  }
  tau = Vectorize(tausimp,vectorize.args = "t")
  
  time_int = time
  tau_int = tau(time_int)

  h0 = post(fit,covariate0,time_int,"h")
  h1 = post(fit,covariate1,tau_int,"h")
  phi = h0/h1

  return(phi)
  
  # left = function(x) tau(x)
  psi = volterra_solve(k=function(x,s){ x*s*0+1},f = tau, 
                       a = 0, b = max(time),num=2000L,method="trapezoid")
  colnames(psi) = c("time","tvta")
  return(psi)
  
}


translate_time2 <- function(fit, covariate0, covariate1, time){
  X <- matrix(rep(covariate0,length(time)),ncol=length(covariate0),byrow=TRUE)
  X1 <- matrix(rep(covariate1,length(time)),ncol=length(covariate1),byrow=TRUE)
  tausimp = function(t){
    H0t = post(fit,covariate0,t,"H")
    H1 = function(u) post(fit,covariate1,u,"H")
    Hdiff = function(u) H1(u) - H0t
    return(rstpm2::vuniroot(Hdiff,lower=0,upper=max(time),tol = .Machine$double.eps)$root) #interval=c(0,max(time)))
  }
  tau = Vectorize(tausimp,vectorize.args = "t")
  
  # left = function(x) tau(x)
  psi = volterra_solve(k=function(x,s){ x*s*0+1},f = tau, 
                       a = 0, b = max(time),num=2000L,method="trapezoid")
  colnames(psi) = c("time","tvta")
  return(psi)
  
}


