path_gh <- "GH_logT.R" #change to location of GH_pen.R
path_weib <- "MixtureWeibull.R" #change to location of MixtureWeibull.R

source(path_gh)
source(path_weib)

library(survival)
library(rstpm2)
library(survPen)
library(biostat3)

set.seed(23234)


mel <- biostat3::melanoma %>%
  mutate(X = ifelse(stage == "Localised", 0, 1)) %>%
  mutate(event = ifelse(status == "Dead: cancer", 1, 0)) %>%
  mutate(time = surv_mm) %>%
  mutate(period = ifelse(year8594=="Diagnosed 75-84",0,1))



####################### log

source(path_gh)
fit_mel_unadj <- fit_genhaz(Surv(mel$time, mel$event), ~ X, data = mel, model_type = "GH", profile = TRUE, n_knots = 8, 
                      tol_LCV = 0.001, timeIt = TRUE) #,init = rep(0,10))
fit_mel_adj <- fit_genhaz(Surv(mel$time, mel$event), ~ X + period + agegrp + sex, data = mel, model_type = "GH", profile = TRUE, n_knots = 8, 
                            tol_LCV = 0.001, timeIt = TRUE) #,init = rep(0,10))

#debug(fit_genhaz)
#save
saveRDS(list(fit_mel_unadj=fit_mel_unadj, fit_mel_adj = fit_mel_adj), file="Melanoma_fit.rds")
res <- readRDS("Melanoma_fit.rds")
fit_mel_unadj <- res$fit_mel_unadj
fit_mel_adj <- res$fit_mel_adj
fit_mel_unadj$par
fit_mel_adj$par
fit_mel_adj$se

#get time varying hazard ratio from survPen
fit_pen_unadj <- survPen(Surv(time, event) ~ smf(time) + smf(time,by=X), data = mel, t1 = time, event = event)
#fit_pen_unadj

fit_pen_adj <- survPen(Surv(time, event) ~ smf(time) + smf(time,by=X) +
                       period + agegrp + sex, 
                       data = mel, t1 = time, event = event)
#fit_pen_adj


tlim=max(mel$time)
new.time = seq(0,tlim,by=0.01)
CIs <- CI(fit_mel_unadj,seq(0,tlim,by=0.01),0,alpha=0.05)
CIs_exp <- CI(fit_mel_unadj,seq(0,tlim,by=0.01),1,alpha=0.05)
pred.pen.exp <- predict(fit_pen_unadj,data.frame(time=new.time,X=1))
pred.pen.unexp <- predict(fit_pen_unadj,data.frame(time=new.time,X=0))

png("hazards_plot_mel_GH_unadj.png", width=800, height=600)
#the hazards CIs$h for both the unexp and exp in one plot
plot(CIs$time, CIs$h, type="l", col="blue", ylim=c(0,0.018),xlim=c(0,200),
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9,    # title bigger
     xlab="Time (months)", ylab="Hazard", main="Estimated Hazard from GH Model")
lines(new.time, pred.pen.unexp$haz, col="blue", lty=2)
lines(CIs_exp$time, CIs_exp$h, col="red")
lines(new.time, pred.pen.exp$haz, col="red", lty=2)
legend("topright", legend=c("Unexposed-GH","Unexposed-survPen", "Exposed-GH","Exposed-survPen"), col=c("blue","blue", "red", "red"), 
       lty=c(1,2,1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()

CIs_adj <- CI(fit_mel_adj,seq(0,tlim,by=0.01),c(0,1,0,1,0,0),alpha=0.05)
CIs_exp_adj <- CI(fit_mel_adj,seq(0,tlim,by=0.01),c(1,1,0,1,0,0),alpha=0.05)
pred.pen.exp_adj <- predict(fit_pen_adj,data.frame(time=new.time,X=1,sex="Male",agegrp="60-74",period=1))
pred.pen.unexp_adj <- predict(fit_pen_adj,data.frame(time=new.time,X=0,sex="Male",agegrp="60-74",period=1))

png("hazards_plot_mel_GH_adj.png", width=800, height=600)
#the hazards CIs$h for both the unexp and exp in one plot
plot(CIs_adj$time, CIs_adj$h, type="l", col="blue", ylim=c(0,max(CIs_adj$h,CIs_exp_adj$h)),xlim=c(0,max(mel$time)),
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9,    # title bigger
     xlab="Time (months)", ylab="Hazard", main="Estimated Hazard from GH Model")
lines(new.time, pred.pen.unexp_adj$haz, col="blue", lty=2)
lines(CIs_exp_adj$time, CIs_exp_adj$h, col="red")
lines(new.time, pred.pen.exp_adj$haz, col="red", lty=2)
legend("topright", legend=c("Unexposed-GH","Unexposed-survPen", "Exposed-GH","Exposed-survPen"), col=c("blue","blue", "red", "red"), 
       lty=c(1,2,1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()



#get the time varying hazard ratio from the GH model
hr_gh_unadj <- CIs_exp$h / CIs$h
hr_gh_adj <- CIs_exp_adj$h / CIs_adj$h

hr_pen_unadj <- pred.pen.exp$haz / pred.pen.unexp$haz
hr_pen_adj <- pred.pen.exp_adj$haz / pred.pen.unexp_adj$haz


png("HR_plot_mel_GH_unadj.png", width=800, height=600)
#the hazards CIs$h for both the unexp and exp in one plot
plot(new.time, hr_gh_unadj, type="l", col="purple", xlab="Time (months)", 
     ylab="Hazard Ratio", main="Estimated Time-Varying Hazard Ratio from GH Model",
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9, xlim=c(0,200),ylim=c(0,max(hr_gh_unadj, hr_pen_unadj)))
lines(new.time, hr_pen_unadj, type="l", col="purple", xlab="Time (months)", ylab="Hazard Ratio", 
      main="Estimated Time-Varying Hazard Ratio from survPen", lty=2)
legend("topright", legend=c("HR-GH","HR-survPen"), col=c("purple","purple"), 
       lty=c(1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()

png("HR_plot_mel_GH_adj.png", width=800, height=600)
#the hazards CIs$h for both the unexp and exp in one plot
plot(new.time, hr_gh_adj, type="l", col="purple", xlab="Time (months)", 
     ylab="Hazard Ratio", main="Estimated Time-Varying Hazard Ratio from GH Model",
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9, xlim=c(0,200),ylim=c(0,max(hr_gh_adj, hr_pen_adj)))
lines(new.time, hr_pen_adj, type="l", col="purple", xlab="Time (months)", ylab="Hazard Ratio", 
      main="Estimated Time-Varying Hazard Ratio from survPen", lty=2)
legend("topright", legend=c("HR-GH","HR-survPen"), col=c("purple","purple"), 
       lty=c(1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()


png("survival_plot_mel_GH_unadj.png", width=800, height=600)
plot(CIs$time, CIs$S, type="l", col="blue", ylim=c(0,max(CIs$S,CIs_exp$S,pred.pen.exp$surv,pred.pen.unexp$surv)),xlim=c(0,200),
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9,    # title bigger
     xlab="Time (months)", ylab="Survival", main="Estimated Survival from GH Model")
lines(new.time, pred.pen.unexp$surv, col="blue", lty=2)
lines(new.time, CIs_exp$S, col="red")
lines(new.time, pred.pen.exp$surv, col="red", lty=2)
legend("topright", legend=c("Unexposed-GH","Unexposed-survPen", "Exposed-GH","Exposed-survPen"), col=c("blue","blue", "red", "red"), 
       lty=c(1,2,1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()

png("survival_plot_mel_GH_adj.png", width=800, height=600)
plot(new.time, CIs_adj$S, type="l", col="blue", ylim=c(0,max(CIs_adj$S,CIs_exp_adj$S,pred.pen.exp_adj$surv,pred.pen.unexp_adj$surv)),xlim=c(0,max(mel$time)),
     cex.lab=1.5,     # axis labels bigger
     cex.axis=1.5,    # tick labels bigger
     cex.main=1.9,    # title bigger
     xlab="Time (months)", ylab="Survival", main="Estimated Survival from GH Model")
lines(new.time, pred.pen.unexp_adj$surv, col="blue", lty=2)
lines(new.time, CIs_exp_adj$S, col="red")
lines(new.time, pred.pen.exp_adj$surv, col="red", lty=2)
legend("topright", legend=c("Unexposed-GH","Unexposed-survPen", "Exposed-GH","Exposed-survPen"), col=c("blue","blue", "red", "red"), 
       lty=c(1,2,1,2), cex = 1.7)
#rug(mel$time[mel$X==0], col="blue")
#rug(mel$time[mel$X==1], col="red")
dev.off()

