#### Imports, adjust paths here ####
#imports for Genhaz
library(survPen) #not strictly needed for our inititial values here
library(dplyr)
library(survival) #surv
library(pracma) #gausslegendre

# C++ spline implementation
library(Rcpp)
library(RcppArmadillo)
sourceCpp(file="crs_cpp.cpp") #change to location of crs_cpp.cpp

path_gh <- "GH_logT.R" #change to location of GH_pen.R
path_weib <- "MixtureWeibull.R" #change to location of MixtureWeibull.R
source(path_gh) #Genhaz
source(path_weib) #For drawing from our baselines

library(parallel)
RNGkind("L'Ecuyer-CMRG")

save_file_name = "GenHaz_SimStudy_3103.rds"

mc_cores = 80 # leave one core free for the system, adjust if you want to use more or less cores
#mc_cores = 1

#### DGM and FIT function ####
simulate <- function(beta1, beta2, szenario, n_obs, tmax,reKnot, n_knots = 8, tol_LCV = 0.001) {
  data_sim = sim_szenario(szenario, beta1, beta2, n = n_obs, tmax = tmax)
  fit = tryCatch({
    fit_genhaz(
      Surv(data_sim$time, data_sim$event),
      ~ X,
      data = data_sim,
      model_type = "GH",
      init = rep(0,n_knots+2), #we are not using the survPen init
      profile = TRUE, #do optimize over roh
      n_knots = n_knots, #number of knots
      tol_LCV = tol_LCV, 
      reKnot = reKnot #double fit to get adjusted knots
    )
  }, error = function(e) {
    cat("Error in fitting model:", e$message, "\n")
    return(NULL)
  })
  return(fit)
}

#### Setup simulation parameters and result data frame ####
beta1_values = c(-0.5)
beta2_values = c(-0.5, 0, 0.5)
szenarios = 1:3
tmax = 9
n_sim = 1000
n_obs = 2000
n_knots = 8
reKnot = 1 
#GH,AFT,AH
sims = expand.grid(sim = 1:n_sim, szenario = szenarios, beta1 = beta1_values, beta2 = beta2_values, tmax = tmax, n_obs = n_obs, n_knots=n_knots, reKnot = reKnot)
#PH
sims = sims %>% bind_rows(
       expand.grid(sim = 1:n_sim, szenario = szenarios, beta1 = 0, beta2 = -0.5, tmax = tmax, n_obs = n_obs, n_knots=n_knots, reKnot = reKnot))
#add custom sim with n_knots = 6,10
sims = sims %>% bind_rows(
       expand.grid(sim = 1:n_sim,szenario = 1,       beta1 = -0.5, beta2 = -0.5, tmax = tmax, n_obs = n_obs, n_knots=c(6,10), reKnot = reKnot))
#reknot = 4 for testing the convergence of reKnot
sims = sims %>% bind_rows(
  expand.grid(sim = 1:n_sim, szenario = 1:2    , beta1 = -0.5, beta2 = -0.5, tmax = tmax, n_obs = n_obs, n_knots=n_knots, reKnot = 4))

# #create dummy sims with only one line for testing locally
#sims = data.frame(sim = 1:1, szenario = 2, beta1 = -0.5, beta2 = -0.5, tmax = tmax, n_obs = n_obs, n_knots=n_knots,reKnot=1)

#initialize list with fits params and var-cov matrices for later use
fits_pars = vector("list", nrow(sims))
fits_vars = vector("list", nrow(sims))
fits_Z = vector("list", nrow(sims))
fits_knots = vector("list", nrow(sims))
knot_hist = vector("list", nrow(sims))
par_hist = vector("list", nrow(sims))

#### Setup randomness, seed is here #####
set.seed(3927315) # for reproducibility
n_tasks = nrow(sims) 
## first, get the seeds
seeds = vector("list", n_tasks)
seeds[[1]] = nextRNGStream(.Random.seed)
for (i in 2:n_tasks) {
  seeds[[i]] <- nextRNGStream(seeds[[i-1]])
}


#### Run simulation ####
print(paste0("Running ", n_tasks, " simulations ..."))
all_results <- mclapply(1:n_tasks, function(i) {
  .Random.seed <<- seeds[[i]]
  
  #Now with shiny and secure tryCatch
  fit <- tryCatch({
    simulate(beta1=sims$beta1[i],beta2=sims$beta2[i], szenario=sims$szenario[i], n_obs=sims$n_obs[i], tmax=sims$tmax[i], n_knots = sims$n_knots[i], reKnot=sims$reKnot[i])
  }, error = function(e) {
    cat("Error in simulation or fitting for simulation", i, ":", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(fit)) {
    list(
      beta1_est = fit$par["beta1_1"],
      beta2_est = fit$par["beta2_1"],
      se1 = sqrt(fit$var["beta1_1","beta1_1"]),
      se2 = sqrt(fit$var["beta2_1","beta2_1"]),
      conv = fit$convergence,
      lcv = fit$LCV,
      edf = fit$edf,
      roh = fit$roh_opt,
      par = fit$par,
      var = fit$var,
      Z = fit$Z,
      knots = fit$knots,
      knot_hist = fit$knot_hist,
      par_hist = fit$par_hist
    )
  } else {
    list(conv = -1)
  }
}, mc.cores=mc_cores, mc.set.seed=FALSE)

#### Save results ####
for(i in 1:n_tasks) {
  res <- all_results[[i]]
  if(!is.null(res) & res$conv != -1) {
    sims$beta1_est[i] <- res$beta1_est
    sims$beta2_est[i] <- res$beta2_est
    sims$beta1_se_est[i] <- res$se1
    sims$beta2_se_est[i] <- res$se2
    sims$converged[i] <- res$conv
    sims$LCV[i] <- res$lcv
    sims$edf[i] <- res$edf
    sims$roh_opt[i] <- res$roh
    fits_pars[[i]] <- res$par
    fits_vars[[i]] <- res$var
    fits_Z[[i]] <- res$Z
    fits_knots[[i]] <- res$knots
    knot_hist[[i]] <- res$knot_hist
    par_hist[[i]] <- res$par_hist
  }
}
saveRDS(list(sims=sims, fits_pars=fits_pars, fits_vars=fits_vars, fits_Z=fits_Z, fits_knots = fits_knots, knot_hist=knot_hist, par_hist=par_hist), file=save_file_name)
cat("Simulation results saved to", file.path(getwd(), save_file_name), "\n")
