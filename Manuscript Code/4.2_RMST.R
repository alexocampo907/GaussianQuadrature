# Preamble ----------------------------------------------------------------
set.seed(103)
start.time.overall <- proc.time()

library(survival)
library(spatstat.random)
library(gridExtra)

#### Set Parameter load.sim ####
# load.sim: toggle T or F to either run the simulation (load.sim=F,~13 minutes) or load the data set with simulation results (load.sim=T,~few seconds)
load.sim = TRUE

# Begin -------------------------------------------------------------------

iter <- 1000  # number of simulations
Ns <- N <- c(1e6)  # number of subjects per arm
k <- 3  # expected TTE
A_flag <- T  # include direct effect
M_flag <- T  # include indirect effect
times <- c(k)  # time-points to evaluate survival probs / RMST
mu_1 <- -1  # expectation of M|(A=1)
pi_c <- 0.2  # proportion censored

results <- data.frame()

# Calculate true RMSTs ---------------------------------------------------------
tru_rmsts <- data.frame(TE = rep(NA, length(times)), NDE = NA, NIE = NA)

S_func <- function(M, A, t) {  # survival probability
  return(exp(-t/k * (k/(k+1))^(A*A_flag) * ((k+1)/k)^(M*M_flag)))
}

R_func <- function(M, A, t) {  # RMST
  return((k * (k/(k+1))^(-A*A_flag) * ((k+1)/k)^(-M*M_flag)) * (1-S_func(M, A, t)))
}

start.time <- proc.time()
for(i in 1:length(times)) {
  rmst_trt_m1 <- gauss.hermite(R_func, mu_1, 1, A = 1, t = times[i])
  rmst_trt_m0 <- gauss.hermite(R_func, 0, 1, A = 1, t = times[i]) 
  rmst_untrt_m0 <- gauss.hermite(R_func, 0, 1, A = 0, t = times[i])
  
  tru_rmsts[i,"TE"] <- rmst_trt_m1 - rmst_untrt_m0
  tru_rmsts[i,"NDE"] <- rmst_trt_m0 - rmst_untrt_m0
  tru_rmsts[i,"NIE"] <- rmst_trt_m1 - rmst_trt_m0
}
(proc.time() - start.time)[3]

tru_rmsts

rmst_KM <- function(time, surv, tau) {  # RMST estimator for multiple times
  cutoffs <- findInterval(tau, time)
  cumsum(c(surv[1],surv) * diff(c(0,time,0)))[cutoffs] + 
    surv[cutoffs] * (tau - time[cutoffs])
}


# MCI ---------------------------------------------------------------------
if(!load.sim){
start.time <- proc.time()

true_MCI <- matrix(NA,nrow=iter,ncol = 3)
colnames(true_MCI) <- c("MC_TE","MC_NDE","MC_NIE")
time_MCI <- numeric(iter)

for(sim in 1:iter) { 
  # sim=1
  start.time <- proc.time()
  
  M0 <- rnorm(2*N, 0, 1)
  M1 <- rnorm(2*N, mu_1, 1)
  
  lambda00 <- 1/k * (k/(k+1))^(A_flag*0) * ((k+1)/k)^(M_flag*M0)
  lambda11 <- 1/k * (k/(k+1))^(A_flag*1) * ((k+1)/k)^(M_flag*M1)
  lambda10 <- 1/k * (k/(k+1))^(A_flag*1) * ((k+1)/k)^(M_flag*M0)
  
  # RMST functional
  mu00 <- mean((1/lambda00) * (1-exp(-lambda00*times)))
  mu11 <- mean((1/lambda11) * (1-exp(-lambda11*times)))
  mu10 <- mean((1/lambda10) * (1-exp(-lambda10*times)))
  
  MC_TE <- mu11 - mu00
  MC_NDE <- mu10 - mu00
  MC_NIE <- mu11 - mu10
  # tru_rmsts
  
  true_MCI[sim,] <- c(MC_TE,MC_NDE,MC_NIE)
  time_MCI[sim] <- (proc.time() - start.time)[3]
}
saveRDS(true_MCI,"Results/R Objects/Section 4/results_trueMCI.rds")
saveRDS(time_MCI,"Results/R Objects/Section 4/results_timeMCI.rds")
}else{
true_MCI <- readRDS("./Results/R Objects/Section 4/results_trueMCI.rds")
time_MCI <- readRDS("./Results/R Objects/Section 4/results_timeMCI.rds")
}

# Results -----------------------------------------------------------------

# Wrangle
true_MCI <- data.frame(true_MCI)

# EDA
dim(true_MCI)
colMeans(true_MCI)
head(true_MCI)

mean(time_MCI)
mean(time_MCI)/0.073

tru_rmsts
colMeans(true_MCI)-tru_rmsts
apply(true_MCI,2,sd)*2

# elapsed 
# 0.205
head(true_MCI)

true_MCI <- data.frame(true_MCI)

# Figure -------------------------------------------------------------

# TE
psi_tru = true_psi = NA
psi_quad=tru_rmsts[1]

results_summ <- data.frame(method = c("MC integration","Quadrature"),
                           psi_median = unlist(c(median(true_MCI$MC_TE), tru_rmsts[1])),
                           psi_mean = unlist(c(mean(true_MCI$MC_TE), tru_rmsts[1])),
                           lo = c(mean(true_MCI$MC_TE) - 1.96*sd(true_MCI$MC_TE), NA),
                           hi = c(mean(true_MCI$MC_TE) + 1.96*sd(true_MCI$MC_TE), NA))
results_summ$method <- factor(results_summ$method, levels = c("Quadrature","MC integration"), ordered = T)
results_summ 

plot_TE <- ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi)) +
    # geom_hline(yintercept = 12, linetype = "dashed", color = "blue", size=.3) +
    geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
    ylab(expression(TE(tau))) +
    xlab("") +
    ggtitle("Total Effect") +
    theme_light() +
    coord_flip(ylim = c(min(results_summ$lo), max(results_summ$lo)), clip = 'off') +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(1,7,1,1), "lines"))

# NDE
psi_tru = true_psi = NA
psi_quad=tru_rmsts[2]

results_summ <- data.frame(method = c("MC integration","Quadrature"),
                           psi_median = unlist(c(median(true_MCI$MC_NDE), tru_rmsts[2])),
                           psi_mean = unlist(c(mean(true_MCI$MC_NDE), tru_rmsts[2])),
                           lo = c(mean(true_MCI$MC_NDE) - 1.96*sd(true_MCI$MC_NDE), NA),
                           hi = c(mean(true_MCI$MC_NDE) + 1.96*sd(true_MCI$MC_NDE), NA))
results_summ$method <- factor(results_summ$method, levels = c("Quadrature","MC integration"), ordered = T)
results_summ 

(plot_NDE <- ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi)) +
    geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
    ylab(expression(NDE(tau))) +
    xlab("") +
    ggtitle("Natural Direct Effect") +
    theme_light() +
    coord_flip(ylim = c(min(results_summ$lo), max(results_summ$lo)), clip = 'off') +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(1,7,1,1), "lines")))

(plot_NDE_median <- ggplot(results_summ, aes(x = method, y = psi_median, ymin = lo, ymax = hi)) +
    geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
    ylab(expression(NDE(tau))) +
    xlab("") +
    ggtitle("Natural Direct Effect") +
    theme_light() +
    coord_flip(ylim = c(min(results_summ$lo), max(results_summ$lo)), clip = 'off') +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(1,7,1,1), "lines")))

# NIE
psi_tru = true_psi = NA
psi_quad=tru_rmsts[3]

results_summ <- data.frame(method = c("MC integration","Quadrature"),
                           psi_median = unlist(c(median(true_MCI$MC_NIE), tru_rmsts[3])),
                           psi_mean = unlist(c(mean(true_MCI$MC_NIE), tru_rmsts[3])),
                           lo = c(mean(true_MCI$MC_NIE) - 1.96*sd(true_MCI$MC_NIE), NA),
                           hi = c(mean(true_MCI$MC_NIE) + 1.96*sd(true_MCI$MC_NIE), NA))
results_summ$method <- factor(results_summ$method, levels = c("Quadrature","MC integration"), ordered = T)
results_summ 

plot_NIE <- ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi)) +
    # geom_hline(yintercept = 12, linetype = "dashed", color = "blue", size=.3) +
    geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
    ylab(expression(NIE(tau))) +
    xlab("") +
    ggtitle("Natural Indirect Effect") +
    theme_light() +
    coord_flip(ylim = c(min(results_summ$lo), max(results_summ$lo)), clip = 'off') +
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          plot.margin = unit(c(1,7,1,1), "lines"))

# Figure 11 ---------------------------------------------------------------
p_mediation_RMST <- grid.arrange(plot_TE,
             plot_NDE,
             plot_NIE,nrow=3
)

ggsave("Results/Plots/Section 4/Figure 11.pdf", plot = p_mediation_RMST, width = 8, height = 6, device = "pdf")


(proc.time() - start.time.overall)[3]

