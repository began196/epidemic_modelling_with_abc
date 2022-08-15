library(deSolve)
library(EasyABC)
library(MASS)
source("D:/OneDrive - University of Bristol/FRP_content/function.R")

# SIR model setup
times <- seq(0,10)
initial <- c(S=99,I=1,R=0)
param <- c(2,0.3)
output <- ode(y = initial, times = times, func = sir,parms = param)
out <- as.data.frame(output)
out$time <- NULL

matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(7, 60, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")

prior <- list(list("unif",0.5,2.5),list("unif",0,0.5))

sir_model_ss_ini <- function(para){
  return(sir_model_ss(initial,times,para))
}

# our observed summary statistics
summary_stats <- out$I

set.seed(159)

# want our algorithms to generate 1000 samples from the ABC posterior each
start_rej <- Sys.time() 
rej <- abc_rej(sir_model_ss_ini,prior,summary_stats,10000,0.1) 
end_rej <- Sys.time()
rej_time <- end_rej-start_rej # run time

rej_adjust <- local_linear_reg(rej$parameters,rej$summarystats,summary_stats) # local-linear regression adjustment

start_mcmc <- Sys.time()
mcmc <- abc_mcmc(sir_model_ss_ini,prior,summary_stats,n_sample=1000,n_trials = 1000,proportion = 0.01)
end_mcmc <- Sys.time()
mcmc_time <- end_mcmc-start_mcmc

mcmc_adjust <- local_linear_reg(mcmc$parameters,mcmc$summarystats,summary_stats)

smc_prior <- list(c("unif",0.5,2.5),c("unif",0,0.5))
smc <- ABC_sequential("Delmoral",model=sir_model_ss_ini,smc_prior,nb_simul = 1e3,summary_stat_target=summary_stats,tolerance_target=0.001,alpha=0.1)
smc$epsilon

par(mfrow=(c(3,2)))
plot(density(rej$parameters[,1]),main=expression(paste("Posterior of ",beta, " (Rejection)")),ylim=c(0,4))
lines(density(rej_adjust[,1]),col="red")
abline(v=2,lty=2)
plot(density(rej$parameters[,2]),main=expression(paste("Posterior of ",gamma, " (Rejection)")),ylim=c(0,90))
lines(density(rej_adjust[,2]),col="red")
abline(v=0.3,lty=2)
plot(density(mcmc$parameters[,1]),main=expression(paste("Posterior of ",beta, " (MCMC)")),ylim=c(0,8))
lines(density(mcmc_adjust[,1]),col="red")
abline(v=2,lty=2)
plot(density(mcmc$parameters[,2]),main=expression(paste("Posterior of ",gamma, " (MCMC)")),ylim=c(0,270))
lines(density(mcmc_adjust[,2]),col="red")
abline(v=0.3,lty=2)
plot(density(smc$param[,1]),main=expression(paste("Posterior of ",beta, " (SMC)")))
abline(v=2,lty=2)
plot(density(smc$param[,2]),main=expression(paste("Posterior of ",gamma, " (SMC)")))
abline(v=0.3,lty=2)
par(mfrow=c(1,1))

print(c(10000,mcmc$nsim,smc$nsim))
print(c(rej_time,mcmc_time,smc$computime))

mean(smc$param[,1])

smc$epsilon

par(mfrow=c(1,2))
ts.plot(mcmc$parameters[,1],main="",ylab=expression(paste(beta)))
ts.plot(mcmc$parameters[,2],main="",ylab=expression(paste(gamma)))
par(mfrow=c(1,1))

acf(mcmc$parameters[,1])
acf(mcmc$parameters[,2])
