library(deSolve)
library(EasyABC)
library(latex2exp)

sir <- function(time, state, parameters){
  S <- state[1]
  I <- state[2]
  R <- state[3]
  N <-  S+I+R
  beta <- parameters[1]
  gamma <- parameters[2]  
  dS <- -beta * S * I/N
  dI <-  beta * S * I/N - gamma * I
  dR <- gamma * I
  return(list(c(dS, dI, dR)))
}

initial <- c(S = 99, I =1, R = 0)

parameters <- c(1.5,0.2)

R0 <- 1.5/0.2

times      <- seq(0, 10, by = 1)

out <- ode(y = initial, times = times, func = sir, parms = parameters)

out <- as.data.frame(out)

out$time <- NULL

matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(8, 90, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")

prior=list(c("normal",1.5,0.5),c("normal",0.2,0.1))

I_end<- out$I[11]
R_end <- out$R[11]

infec_ss <- function(parameter){
  beta <- parameter[1]
  gamma <- parameter[2]
  output <- ode(y = initial, times = times, func = sir,parms = c(beta,gamma))
  output <- as.data.frame(output)
  output$time <- NULL
  return(c(output$I,output$R))
}

abc_rej <- ABC_rejection(model= infec_ss, prior = prior, nb_simul= 100, summary_stat_target=c(out$I,out$R),tol = 0.1)

plot(density(abc_rej$param[,1]),xlab = expression(paste("Beta, ", beta)),main =expression(paste("Posterior of Beta, ", beta)))
abline(v = 1.5, col = "darkblue", lty = 2)
plot(density(abc_rej$param[,2]),xlab = expression(paste("Gamma, ", gamma)),main =expression(paste("Posterior of Gamma, ", gamma)))
abline(v = 0.2, col = "darkblue", lty = 2)
plot(density(abc_rej$param[,1]/abc_rej$param[,2]),xlab = TeX("R Number, $R_0$"),xlim = c(5,10), main =TeX("Posterior of R Number, $R_0$"))
abline(v = R0,col="darkblue",lty = 2)

simul_I <- abc_rej$stats[,1:11]
simul_R <- abc_rej$stats[,12:22]

plot(times,simul_I[1,],type = "l",lty=1,col="grey",xlab = "Day", ylab ="Number of People", main = "Simulation vs True Data",ylim = c(0,70))
lines(times, simul_R[1,],lty =1, col= "grey")
for(i in 2:dim(simul_I)[1]){
  lines(times,simul_I[i,], lty = 1, col="grey")
  lines(times,simul_R[i,], lty = 1, col="grey")
}
points(x=times, y = out$I, col = "red",pch = 19)
points(x = times, y = out$R, col = "darkgreen",pch = 19,bg = "darkgreen")
legend(0, 50, c("Infected", "Recovered", "Simulated"), pch = c(19,19,19), col = c("red",'darkgreen','grey'), bty = "n")

abc_smc <- ABC_sequential(method = "Emulation", model = infec_ss, prior = prior, summary_stat_target = c(out$I,out$R),nb_simul = 10)

