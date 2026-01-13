# rm(list = ls());gc()
set.seed(907)
start.time.overall <- proc.time()
expit <- function(x) exp(x)/(1+exp(x))
library(tidyverse)
library(mvQuad)
library(dplyr)

#### Set Parameter load.sim ####
# load.sim: toggle T or F to either run the simulation (load.sim=F,~13 minutes) or load the data set with simulation results (load.sim=T,~few seconds)
load.sim = F
if(!load.sim){

# results -----------------------------------------------------------------

iter=1000
res <- data.frame(
  psi_MC = rep(NA, iter),
  time_MC = rep(NA, iter),
  psi_GH = rep(NA, iter),
  time_GH = rep(NA, iter)
)

# Repeated Samples -------------------------------------------------------------

#k=1
for(k in 1:iter){
  
  start.time <- proc.time()
  N=1e6 # print(1e6)
  C <- rnorm(N,-10,1)
  U <- rnorm(N,3,1)
  df1 <- data.frame(a=1,C,U)
  df0 <- data.frame(a=0,C,U)
  
  df1 <- df1 %>% mutate(L = 15 + a + U*.1 + rnorm(N))
  df0 <- df0 %>% mutate(L = 15 + a + U*.1 + rnorm(N))
  df1$M <- df0$M <- 0
  
  beta <- c(100,10,2,2,2,2)
  EY10 <- model.matrix(~ a + M + C + L + U, data = df1) %*% beta
  EY00 <- model.matrix(~ a + M + C + L + U, data = df0) %*% beta
  
  (mu_a1m_MC <- mean(EY10))
  (mu_a0m_MC <- mean(EY00))
  (res$psi_MC[k] <- mean(EY10) - mean(EY00)) # Truth is 12
  
  res$time_MC[k] <- (proc.time() - start.time)[3]
  
  # gauss -------------------------------------------------------------------
  
  ECALMU <- function(X,a,m) as.matrix(cbind(1,a,m,X)) %*% beta
  # mean(ECALMU(X1,a=1,m=0))
  # mean(ECALMU(X0,a=0,m=0))
  
  # GHN (multiplies the weights by normal density)
  start.time <- proc.time()
  myGrid <- createNIGrid(dim=3, type="GHN", level=5)
  # cov(X1)
  Sigma_X <- matrix(c(1.01,0,.1,0,1,0,.1,0,1),nrow=3,byrow=T)
  
  # colMeans(X1)
  mvQuad::rescale(myGrid, m=c(16.3,-10,3),C=Sigma_X,dec.type = 1)
  (mu1 <- quadrature(f = ECALMU,grid = myGrid,a=1,m=0))
  
  # mu0
  myGrid <- createNIGrid(dim=3, type="GHN", level=5)
  
  #cov(X0)
  Sigma_X <- matrix(c(1.01,0,.1,0,1,0,.1,0,1),nrow=3,byrow=T)
  
  #colMeans(X0)
  mvQuad::rescale(myGrid, m=c(15.3,-10,3),C=Sigma_X,dec.type = 1 )
  
  (mu0 <- quadrature(f = ECALMU,grid = myGrid,a=0,m=0))
  res$psi_GH[k] <-mu1-mu0
  res$time_GH[k] <- (proc.time() - start.time)[3]
}
saveRDS(res,"Results/R Objects/Section 4/results_CDEcontinuous.rds")
}else{
res <- readRDS("./Results/R Objects/Section 4/results_CDEcontinuous.rds")
}
colMeans(res,na.rm = T)

# Results summary -----------------------------------------------------------------

dim(res)
res <- res %>% filter(!is.na(psi_MC))
dim(res)
head(res)

# Plot -------------------------------------------------------------

psi_tru <- true_psi <-  12
psi_quad=mean(res$psi_GH)

MC_mse <- mean((res$psi_MC - true_psi)^2)
quad_mse <- (12 - true_psi)^2

results_summ <- data.frame(method = c("MC integration","Quadrature"),
                           psi_mean = c(mean(res$psi_MC), psi_quad),
                           lo = c(mean(res$psi_MC) - 1.96*sd(res$psi_MC), NA),
                           hi = c(mean(res$psi_MC) + 1.96*sd(res$psi_MC), NA),
                           mse = c(MC_mse, quad_mse))
results_summ$msetext <- paste("MSE =", signif(results_summ$mse, 2))
results_summ$method <- factor(results_summ$method, levels = c("Quadrature","MC integration"), ordered = T)
results_summ 

# Figure 10 (Top panel) ---------------------------------------------------

(p_continuous <- ggplot(results_summ, aes(x = method, y = psi_mean, ymin = lo, ymax = hi, label = msetext)) +
  geom_hline(yintercept = 12, linetype = "dashed", color = "blue", size=.3) +
  geom_pointrange(size = 0.7, pch = '|', linewidth = 0.4) +
  ylab(expression(Controlled~Direct~Effect~Delta[italic(CDE)])) +
  xlab("") +
  ggtitle("Continuous Outcome") +
  theme_light() +
  geom_text(y = I(1), hjust = -0.15, size = 3) +
  coord_flip(ylim = c(min(results_summ$lo), max(results_summ$lo)), clip = 'off') +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        plot.margin = unit(c(1,7,1,1), "lines")))

ggsave("Results/Plots/Section 4/Figure 10.1.pdf", plot = p_continuous, width = 8, height = 2, device = "pdf")

(proc.time() - start.time.overall)[3]
