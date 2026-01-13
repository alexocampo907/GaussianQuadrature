# Runs simulation study for confounding scenario -------------------------------
library(tidyverse)
library(pracma)
library(mvQuad)
library(gaussquad)

summarize_simresults <- function(results, psi_quad, true_psi = NA) {
  MC_res <- filter(results, method == "MC")$psi
  SMC_res <- filter(results, method == "2SMC")$psi
  MC_mse <- SMC_mse <- quad_mse <- NA  # Compute MSE
  
  if(!is.na(true_psi)) {
    MC_mse <- mean((MC_res - true_psi)^2)
    SMC_mse <- mean((SMC_res - true_psi)^2)
    quad_mse <- (psi_quad - true_psi)^2
  }
  
  results_summ <- data.frame(method = c("MC integration", "Outcome simulation", "Quadrature"),
                             psi_mean = c(mean(MC_res), mean(SMC_res), psi_quad),
                             lo = c(mean(MC_res) - 1.96*sd(MC_res), mean(SMC_res) - 1.96*sd(SMC_res), NA),
                             hi = c(mean(MC_res) + 1.96*sd(MC_res), mean(SMC_res) + 1.96*sd(SMC_res), NA),
                             mse = c(MC_mse, SMC_mse, quad_mse))
  results_summ$method <- factor(results_summ$method, levels = c("Quadrature","MC integration","Outcome simulation"), ordered = T)
  results_summ$msetext <- paste("MSE =", signif(results_summ$mse, 2))
  results_summ
}

run_MCs <- function(CMAT, PrYC, start_time) {
  probs_0 <- PrYC(CMAT, a=0)
  probs_1 <- PrYC(CMAT, a=1)
  
  mid_time <- Sys.time() - start_time
  
  mu0 <- mean(probs_0)
  mu1 <- mean(probs_1)
  psi <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  time_mc <- Sys.time() - start_time
  
  N <- length(probs_0)
  Y0 <- rbinom(N, 1, prob = probs_0)
  Y1 <- rbinom(N, 1, prob = probs_1)
  mu0 <- mean(Y0)
  mu1 <- mean(Y1)
  psi_2 <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  time2 <- (Sys.time() - start_time) + mid_time - time_mc
  
  data.frame(psi = c(psi, psi_2), 
             method = c("MC", "2SMC"), 
             time = as.numeric(c(time_mc, time2), units = "secs"))
  data.frame(psi = psi, 
             method = "MC", 
             time = as.numeric(time_mc, units = "secs"))
}

expit <- function(x) exp(x) / (1+exp(x))

# Figure 4 ---------------------------------------------------------------------
set.seed(101)
iter <- 10000  

PrYc <- function(c,a) expit(1 + 1*a - c)
results <- data.frame()

for(sim in 1:iter) {
  start_time <- Sys.time()
  
  N <- 1000000
  C <- rnorm(N,0,1)
  
  MCres <- run_MCs(C, PrYc, start_time)
  
  results <- rbind(results, MCres)
  if(sim%%100==0) cat("\r",sim*100/iter,"\b% complete")
}
# saveRDS(results, "results_1DGauss.rds")


# Get mvQuad results
quad_time <- rep(NA, 100)

for(i in 1:100) {
  start_time <- Sys.time()
  
  myGrid <- createNIGrid(dim = 1, type = "GHN", level = 20)
  
  mu1 <- quadrature(PrYc, grid = myGrid, a = 1)
  mu0 <- quadrature(PrYc, grid = myGrid, a = 0)
  psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  quad_time[i] <- Sys.time() - start_time
}

# Plot
results_summ <- summarize_simresults(results, psi_quad)

ggplot(results_summ, aes(x=method, y=psi_mean, ymin=lo, ymax = hi)) +
  geom_pointrange(size = 0.6, pch = '|', linewidth = 0.4) +
  coord_flip() +
  ylab(expression(Causal~odds~ratio~psi[italic(OR)])) +
  xlab("") +
  theme_light()



# Figure 5 ---------------------------------------------------------------------
PrYC <- function(cmat,a) expit(1 + 1*a + 0.1*cmat[,1] + 0.1*cmat[,2])

set.seed(101)
iter <- 10000

results <- data.frame()

for(sim in 1:iter) {
  start_time <- Sys.time()
  
  N <- 1000000
  C1 <- rnorm(N,-5,1)
  C2 <- rnorm(N,C1-5,1)
  CMAT <- cbind(C1,C2)
  
  MCres <- run_MCs(CMAT, PrYC, start_time)
  
  results <- rbind(results, MCres)
  if(sim%%100==0) cat("\r",sim*100/iter,"\b% complete")
}
# saveRDS(results, "results_2DGauss.rds")


# Get mvQuad results
quad_time <- rep(NA, 100)

for(i in 1:100) {
  start_time <- Sys.time()
  
  myGrid <- createNIGrid(dim = 2, type = "GHN", level = 20)
  rescale(myGrid, m = c(-10,-5), C = matrix(c(1,1,1,2), 2), dec.type = 1)
  
  mu1 <- quadrature(PrYC, grid = myGrid, a = 1)
  mu0 <- quadrature(PrYC, grid = myGrid, a = 0)
  psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  quad_time[i] <- Sys.time() - start_time
}


# Plot
results_summ <- summarize_simresults(results, psi_quad)

ggplot(results_summ, aes(x=method, y=psi_mean, ymin=lo, ymax = hi)) +
  geom_pointrange(size = 0.6, pch = '|', linewidth = 0.4) +
  coord_flip() +
  ylab(expression(Causal~odds~ratio~psi[italic(OR)])) +
  xlab("") +
  theme_light()



# Figure 6.1 (uniform confounders) ---------------------------------------------
PrYC <- function(cmat,a) expit(-a + 0.5*cmat[,1] + 0.5*cmat[,2])

set.seed(101)
iter <- 10000 

results <- data.frame()

for(sim in 1:iter) {
  start_time <- Sys.time()
  
  N <- 1000000
  C1 <- runif(N, -2, 2)
  C2 <- runif(N, -4, 0)
  CMAT <- cbind(C1,C2)
  MCres <- run_MCs(CMAT, PrYC, start_time)
  
  results <- rbind(results, MCres)
  if(sim%%100==0) cat("\r",sim*100/iter,"\b% complete")
}
# saveRDS(results, "results_2Duni.rds")


# Get true values
e <- exp(1)
mu0_tru <- -(polylog(-1/e^3, 2)-2*polylog(-1/e, 2)-1.806286070444774256651)/4
mu1_tru <- (-2*polylog(-1/e^4, 2)+4*polylog(-1/e^2, 2)+pi^2/6)/8
psi_tru <- (mu1_tru/(1-mu1_tru)) / (mu0_tru/(1-mu0_tru))


# Get mvQuad results
quad_time <- rep(NA, 100)

for(i in 1:100) {
  start_time <- Sys.time()
  
  myGrid <- createNIGrid(dim = 2, type = "GLe", level = 20)
  rescale(myGrid, domain = matrix(c(-2,-4,2,0), 2))
  
  mu1 <- quadrature(PrYC, grid = myGrid, a = 1)/16
  mu0 <- quadrature(PrYC, grid = myGrid, a = 0)/16
  (psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0)))
  
  quad_time[i] <- Sys.time() - start_time
}


# Plot
results_summ <- summarize_simresults(results, psi_quad, true_psi = psi_tru)

ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi, label = msetext)) +
  geom_hline(yintercept = psi_tru, linetype = "dashed", color = "blue", size = .3) +
  geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
  ylab(expression(Causal~odds~ratio~psi[italic(OR)])) +
  xlab("") +
  ggtitle("Uniform confounders") +
  theme_light() +
  geom_text(aes(x = method, y = psi_mean), hjust = -2, size = 3) +
  coord_flip(ylim = c(min(results_summ$lo), max(results_summ$hi)), clip = "off") +
  theme(plot.title  = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(c(1, 7, 1, 1), "lines"))



# Figure 6.2 (exponential confounders) -----------------------------------------
PrYC <- function(cmat,a) expit(-a + 0.5*cmat[,1] + 0.5*cmat[,2])

set.seed(101)
iter <- 10000

results <- data.frame()

for(sim in 1:iter) {
  start_time <- Sys.time()
  
  N <- 1000000
  C1 <- rexp(N,2)
  C2 <- rexp(N,1)
  CMAT <- cbind(C1,C2)
  MCres <- run_MCs(CMAT, PrYC, start_time)
  
  results <- rbind(results, MCres)
  if(sim%%100==0) cat("\r",sim*100/iter,"\b% complete")
}
# saveRDS(results, "results_2Dexp.rds")


# Get true values:
mu0_tru <- 2/3
mu1_tru <- 4/e*(2/3-(1/e-1/2+log(1+e)-log(1+e)/e^2)/e) 
psi_tru <- (mu1_tru/(1-mu1_tru)) / (mu0_tru/(1-mu0_tru))


# Get mvQuad results; first need to reparametrize the function:
PrYC_rscld <- function(cmat,a) expit(-a + 0.25*cmat[,1] + 0.5*cmat[,2])
quad_time <- rep(NA, 100)

for(i in 1:100) {
  start_time <- Sys.time()
  
  myGrid <- createNIGrid(dim = 2, type = "GLa", level = 20)
  
  mu1 <- quadrature(PrYC_rscld, grid = myGrid, a = 1)
  mu0 <- quadrature(PrYC_rscld, grid = myGrid, a = 0)
  psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  quad_time[i] <- Sys.time() - start_time
}


# Plot
results_summ <- summarize_simresults(results, psi_quad, true_psi = psi_tru)

ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi, label = msetext)) +
  geom_hline(yintercept = psi_tru, linetype = "dashed", color = "blue", size = .3) +
  geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
  ylab(expression(Causal~odds~ratio~psi[italic(OR)])) +
  xlab("") +
  ggtitle("Exponential confounders") +
  theme_light() +
  geom_text(aes(x = method, y = psi_mean), hjust = -2, size = 3) +
  coord_flip(ylim = c(min(results_summ$lo), max(results_summ$hi)), clip = "off") +
  theme(plot.title  = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(c(1, 7, 1, 1), "lines"))



# Figure 6.3 (gamma confounders) -----------------------------------------------
PrYC <- function(cmat,a) expit(-a + 0.5*cmat[,1] + 0.5*cmat[,2])

gLagRule.fun <- function(l) {  # Use fact that (Ga(1,2)+Ga(4,2))/2 = Ga(5,1)
  alpha <- 5
  glag_res <- glaguerre.quadrature.rules(l, alpha-1)
  
  n <- glag_res[[l]]$x
  w <- glag_res[[l]]$w / gamma(alpha)  # normalizing constant
  initial.domain <- matrix(c(0, Inf), ncol = 2)
  
  return(list(
    n = as.matrix(n), 
    w = as.matrix(w), 
    features = list(initial.domain = initial.domain)))
}

set.seed(101)
iter <- 10000
results <- data.frame()

for(sim in 1:iter) {
  start_time <- Sys.time()
  
  N <- 1000000
  C1 <- rgamma(N, shape = 1, scale = 2)
  C2 <- rgamma(N, shape = 4, scale = 2)
  
  CMAT <- cbind(C1,C2)
  MCres <- run_MCs(CMAT, PrYC, start_time)
  
  results <- rbind(results, MCres)
  if(sim%%100==0) cat("\r",sim*100/iter,"\b% complete")
}
# saveRDS(results, "results_2Dgam.rds")


# True values
mutru_0 <- 15*zeta(5)/16
mutru_1 <- (-360*polylog(-1/e,5)+3+10*pi^2+7*pi^4)/(360*e)
psi_tru <- (mutru_1/(1-mutru_1)) / (mutru_0/(1-mutru_0))


# Get mvQuad results after reparametrizing the function
PrYC_rscld <- function(cmat,a) expit(cmat - a)
quad_time <- rep(NA, 100)

for(i in 1:100) {
  start_time <- Sys.time()
  
  myGrid <- createNIGrid(dim = 1, type = "gLagRule.fun", level = 20)
  
  mu1 <- quadrature(PrYC_rscld, grid = myGrid, a = 1)
  mu0 <- quadrature(PrYC_rscld, grid = myGrid, a = 0)
  psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  quad_time[i] <- Sys.time() - start_time
}

# Plot
results_summ <- summarize_simresults(results, psi_quad, true_psi = psi_tru)

ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi, label = msetext)) +
  geom_hline(yintercept = psi_tru, linetype = "dashed", color = "blue", size = .3) +
  geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
  ylab(expression(Causal~odds~ratio~psi[italic(OR)])) +
  xlab("") +
  ggtitle("Gamma confounders") +
  theme_light() +
  geom_text(aes(x = method, y = psi_mean), hjust = -2, size = 3) +
  coord_flip(ylim = c(min(results_summ$lo), max(results_summ$hi)), clip = "off") +
  theme(plot.title  = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(c(1, 7, 1, 1), "lines"))



# Figure 7 ---------------------------------------------------------------------
PrYc <- function(c,a) expit(1 + a - 0.5*c)

# Get true results
mu0_tru <- ((e-2)*e + 2*log(e+1))/e^2
mu1_tru <- 1 - 2/e^2 + 2*log(e^2+1)/e^4
psi_tru <- (mu1_tru/(1-mu1_tru)) / (mu0_tru/(1-mu0_tru))

results <- data.frame()
max_gridpoints <- 50


for(i in 2:max_gridpoints) {
  time_vec <- rep(NA, 100000)
  
  for(j in 1:length(time_vec)) {
    start_time <- Sys.time()
    
    myGrid <- createNIGrid(dim = 1, type = "GLa", level = i)
    mu1 <- quadrature(PrYc, grid = myGrid, a = 1)
    mu0 <- quadrature(PrYc, grid = myGrid, a = 0)
    psi_quad <- (mu1/(1-mu1)) / (mu0/(1-mu0))
    
    time_vec[j] <- Sys.time() - start_time
  }
  results <- rbind(results, 
                   data.frame(bias = psi_quad - psi_tru, 
                              gridpoints = i, 
                              time = mean(as.numeric(time_vec, units = "secs"))))
}
# saveRDS(results, "results_biasexp.rds")


# Add some MC results to compare
set.seed(101)
results_mc <- data.frame()
C <- numeric(0)

for(M in c(1e3,1e4,1e5,1e6,1e7,1e8)) {
  start_time <- Sys.time()
  
  C <- c(C, rexp(M,1))
  probs_0 <- PrYc(C, a=0)
  probs_1 <- PrYc(C, a=1)
  mu0 <- mean(probs_0)
  mu1 <- mean(probs_1)
  psi <- (mu1/(1-mu1)) / (mu0/(1-mu0))
  
  time <- Sys.time() - start_time
  results_mc <- rbind(results_mc, 
                      data.frame(bias = psi - psi_tru, 
                                 gridpoints = M, 
                                 time = mean(as.numeric(time, units = "secs"))))
}

# Plot figure 7.1
ggplot(results) +
  geom_point(aes(gridpoints, abs(bias)), size = 1) +
  scale_y_continuous(trans='log2',
                     breaks = c(1e-15,1e-11,1e-7,1e-3)) +
  geom_hline(yintercept = abs(results_mc$bias[4]), linetype="dashed", 
             color = "blue", size=.3) +
  geom_hline(yintercept = abs(results_mc$bias[6]), linetype="dashed", 
             color = "blue", size=.3) +
  annotate("text", x = 40, y = abs(results_mc$bias[4])+0.003, color = "blue",
           label = "Monte~Carlo*','*~italic(N)==10e+05", parse = T, size = 3) +
  annotate("text", x = 40, y = abs(results_mc$bias[6])+0.00004, color = "blue",
           label = "Monte~Carlo*','*~italic(N)==10e+07", parse = T, size = 3) +
  theme_light() +
  xlab("Quadrature points") + ylab("Absolute value of bias")

# Plot figure 7.2
ggplot(results) +
  geom_point(aes(gridpoints, time), size = 1) +
  scale_y_continuous(labels = function(x) format(x, scientific = T)) +
  theme_light() +
  xlab("Quadrature points") + ylab("Runtime (s)")



# Figure 8 ---------------------------------------------------------------------
library(MASS)
nsims <- 500
results <- data.frame()
set.seed(999)

for(d in 2:14) {
  PrYC <- function(cmat,a) expit(1 + 1*a + 0.1*rowSums(cmat))
  
  # MC
  mc_times <- rep(NA, nsims)
  
  for(i in 1:nsims) {
    start_time_mc <- Sys.time()
    
    N <- 1e+06  
    CMAT <- mvrnorm(N, rep(0,d), diag(d))
    probs_0 <- PrYC(CMAT, a=0)
    probs_1 <- PrYC(CMAT, a=1)
    mu0 <- mean(probs_0)
    mu1 <- mean(probs_1)
    psi_mc <- (mu1/(1-mu1)) / (mu0/(1-mu0))
    
    mc_times[i] <- Sys.time() - start_time_mc
  }
  
  # GH
  gh_times <- rep(NA, nsims)
  
  for(i in 1:nsims) {
    start_time_gh <- Sys.time()
    
    myGrid <- createNIGrid(dim = d, type = "GHN", level = 3)  # chosen to match bias
    mu1 <- quadrature(PrYC, grid = myGrid, a = 1)
    mu0 <- quadrature(PrYC, grid = myGrid, a = 0)
    psi_gh <- (mu1/(1-mu1)) / (mu0/(1-mu0))
    
    gh_times[i] <- Sys.time() - start_time_gh
  }
  
  mc_times <- as.numeric(mc_times, units = "secs")
  gh_times <- as.numeric(gh_times, units = "secs")
  
  results <- rbind(results, 
                   data.frame(method = c("Monte Carlo","G-H quadrature"), 
                              d = d, 
                              time = c(median(mc_times), median(gh_times))))
  cat("\r",d*10,"\b% complete")
}
# saveRDS(results, "results_10Dtime.rds")

# Plot 
ggplot(results, aes(x = d, y = time, group = method, color = method)) + 
  geom_line(size = 0.6) +
  scale_y_continuous(trans = 'log10') +
  #scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_color_manual(values = c("#72B879","#787BF5")) +
  ylab("Time (s)") +
  xlab("Dimension of C") +
  theme_light() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) 
