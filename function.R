
sir <- function(time, state, parameters){
  S <- state[1]
  I <- state[2]
  R <- state[3]
  N <-  S+I+R
  beta <- parameters[1]
  gamma <- parameters[2]  
  dS <- -beta*S*I/N
  dI <-  beta*S*I/N -gamma*I
  dR <- gamma * I
  return(list(c(dS, dI, dR)))
}

sir_model_ss <- function(initial,t,parameter){
  S <- initial[1]
  I <- initial[2]
  R <- initial[3]
  beta <- parameter[1]
  gamma <- parameter[2]
  output <- ode(y = c(S=S,I=I,R=R), times = t, func = sir,parms = c(beta,gamma))
  return((as.data.frame(output)$I))
}

get_priors <- function(prior){
  priors <-list(rep(NA,length(prior))) 
  for(i in 1:length(prior)){
    rprior <- paste0("r",prior[[i]][1])
    argu <- c(1,prior[[i]][-1])
    priors[i] <- do.call(rprior,argu)
  }
  return(array(as.numeric(priors)))
}

proposal_dist <- function(prior,param,range){
  p_prior <- list(rep(NA,length(prior)))
  for(i in 1:length(prior)){
    rprior <- paste0("r",prior[[i]][1])
    argu <- prior[[i]][-1]
    if(rprior == "runif"){
      if(param[i]-range[i]/2<argu[[1]]){
        p_prior[i] <- runif(1,min=(argu[[1]]),max=(param[i]+range[i]/2))
      }else if(param[i]+range[i]/2>argu[[2]]){
        p_prior[i] <- runif(1,min=(param[i]-range[i]/2),max=(argu[[2]]))
      }else{
        p_prior[i] <- runif(1,min=(param[i]-range[i]/2),max=(param[i]+range[i]/2))
      }
      #reflect at boundary
    }else if(rprior == "rnorm"){
      p_prior[i] <- rnorm(1,mean=param[i],sd=range[i]/2)
    }
  }
 return(array(as.numeric(p_prior)))
}

density_prior <- function(prior,para){
  dens <- list(rep(NA,length(prior))) 
  for(i in 1:length(prior)){
    dprior <- paste0("d",prior[[i]][1])
    argu <- prior[[i]][-1]
    dens[[i]] <- do.call(dprior,c(list(para[i]),argu))
  }
  return(prod(array(as.numeric(dens))))
}

density_proposal <- function(prior,para1,para2,range){
  dens <- list(rep(NA,length(prior)))
  for(i in 1:length(prior)){
    dprior <- paste0("d",prior[[i]][1])
    if(dprior == "dunif"){
      argu <- prior[[i]][-1]
      if(param[i]-range[i]/2<argu[[1]]){
        dens[[i]] <- dunif(para1[i],min=(argu[[1]]),max=(para2[i]+range[i]/2))
      }else if(param[i]+range[i]/2>argu[[2]]){
        dens[[i]] <- dunif(para1[i],min=(para2[i]-range[i]/2),max=(argu[[2]]))
      }else{
        dens[[i]] <- dunif(para1[i],min=(para2[i]-range[i]/2),max=(para2[i]+range[i]/2))
      }
    }else if(dprior == "dnorm"){
      dens[[i]] <- dnorm(para1[i],mean=para2[i],sd=range[i]/2)
    }
  }
  return(prod(array(as.numeric(dens))))
}


abc_rej <- function(model,prior,summary_stats,simul_number,proportion = 0.01){
  # Intialise empty lists to store data
  # parameters
  param <- vector("list",length = simul_number) 
  # summary statistics
  sumstats <- vector("list",length = simul_number) 
  # differnce between observed and simulated summary statistics
  diffss <- vector("list",length = simul_number) 
  
  for(jj in 1:simul_number){
    # obtain the parameters from prior
    param[[jj]] <- get_priors(prior) 
    # calculate summary statistics according to the model
    sumstats[[jj]] <- model(param[[jj]]) 
    # calculate euclidean distance
    diffss[[jj]] <- sum((sumstats[[jj]]-summary_stats)^2) 
  }
  # indicate number of simulations to keep
  N <- floor(proportion*simul_number)  
  # sort and select N smallest distance
  indices <- sort(as.numeric(diffss),decreasing=FALSE,index.return= TRUE)$ix[1:N] 
  # get the parameters according to the distance
  param <- matrix(unlist(param),ncol=length(param[[1]]),byrow=T)[indices,] 
  
  d <- length(sumstats[[1]])
  # convert list to matrix
  sumstats <- matrix(unlist(sumstats),ncol=d,byrow=T)[indices,] 
  diff_ss <- matrix(unlist(diffss),ncol=1)[indices,]
  done <- list(param,sumstats,diff_ss) 
  names(done) <- c("parameters","summarystats","diff_sumstats")
  return(done) # return output
}


abc_mcmc <- function(model,prior,summary_stats,n_sample,n_trials=10000, proportion = 0.01){
  # Initialisation step to get epsilon^2
  ini <- abc_rej(model,prior,summary_stats,n_trials,proportion)
  # Get maximum epsilon^2
  max_distance <- max(ini$diff_sumstats)
  # Randomly sample an index of the set of accepted parameters
  index <- sample(1:as.integer(n_trials*proportion),1) 
  # Compute range of proposal distribution
  range <- rep(0,ncol(ini$parameters))
  for(i in 1:length(range)){
    range[i] <- max(ini$parameters[,i])-min(ini$parameters[,i])
  }
  # Create 3 empty matrices
  theta <- matrix(0,nrow=n_sample,ncol=ncol(ini$parameters))
  sum_stat <- matrix(0,nrow=n_sample,ncol=ncol(ini$summarystats))
  diff_stat <- rep(0,length(ini$diff_sumstats))
  # Set starting point of the chain 
  theta[1,] <- ini$parameters[index,]
  sum_stat[1,] <- ini$summarystats[index,]
  diff_stat[1] <- ini$diff_sumstats[index]
  nsim <- n_trials # Number of times model() is ran
  j=1 # Start from j=1
  while(j < n_sample){
    theta_star <- proposal_dist(prior,theta[j,],range) # propose new parameter
    sum_stat[j+1,] <- model(theta_star) # compute summary statistics
    nsim <- nsim+1
    diff <- sum((sum_stat[j+1,]-summary_stats)^2) # calculate squared euclidean distance
    diff_stat[j+1] <- diff 
    if(diff < max_distance){
      # calculate probability of acceptance
      alpha <- (density_prior(prior,theta_star)*density_proposal(prior,theta[j,],theta_star,range))/
        (density_prior(prior,theta[j,])*density_proposal(prior,theta_star,theta[j,],range))  
      dummy <- runif(1)
      # keep parameters with probability alpha if squared distance below tolerance
      if(dummy<=alpha){ 
        theta[j+1,] <- theta_star
        j <- j+1
      }
    }else{next}
  }
 # Put the results into a list
 done <- list(theta,sum_stat,diff_stat,nsim)
 names(done) <- c("parameters","summarystats","diff_sumstats","nsim")
 return(done)
}

local_linear_reg <- function(parameters,simulated,summary_stats){
  n <- nrow(simulated)
  converted_ss <- matrix(rep(summary_stats,n),nrow=n,byrow=T)
  # calculate difference between simulated and observed data
  raw_diff <- simulated-summary_stats
  # perform local-linear regression with uniform kernel
  linear <- lm(parameters~raw_diff)
  # extract beta
  beta <- as.matrix(linear$coefficients[-1,])
  adjusted <- parameters-raw_diff%*%beta
  return(adjusted)
}

sir2 <- function(time, state, parameters){
  # Initial values
  S <- state[1]
  I <- state[2]
  R <- state[3]
  I_new <- state[4]
  N <-  S+I+R
  # Assign parameter values
  beta <- parameters[1]
  gamma <- parameters[2]  
  # SIR model as described
  dS <- -beta*S*I/N
  dI <-  beta*S*I/N -gamma*I
  dR <- gamma * I
  # additional compartment to capture cumulative I
  dI_new <- beta*S*I/N 
  return(list(c(dS, dI, dR,dI_new)))
}

SIR3_model <- function(initial,t,parameter){
  # Assign priors to parameters of our model
  beta <- parameter[1]
  gamma <- parameter[2]
  I <- parameter[3]
  # Initial conditions
  S <- initial[1]
  N <- initial[2]
  R <- N-S-I
  # Outputs of the ode after running for t unit time
  output <- as.data.frame(ode(y = c(S=S,I=I,R=R,I_new=I+R), 
                              times = t, func = sir2,parms = c(beta,gamma)))
  #return I_new to compare with the summary statistics  
  return(output$I_new) 
}

SIR_S_model <- function(initial,t,parameter){
  tau <- t[2]-t[1]
  kk <- 24-length(t)%%24
  m <- length(t)
  beta <- parameter[1]
  gamma <- parameter[2]
  I <- parameter[3]
  S <- initial[1]
  N <- initial[2]
  p1 <- 1-exp(-beta*tau*I/N);p2 <- 1-exp(-gamma*tau)
  ss <- rep(0,length(t));ii <- rep(0,length(t));rr <- rep(0,length(t));i_new <- rep(0,length(t))
  ss[1] <- floor(S); ii[1] <- floor(I); rr[1] <-N-floor(S)-floor(I); i_new <- N-floor(S)
  for(i in 2:m){
    d1 <- rbinom(1,ss[i-1],p1)
    d2 <- rbinom(1,ii[i-1],p2)
    ss[i] <- ss[i-1]-d1
    ii[i] <- ii[i-1]+d1-d2
    rr[i] <- rr[i-1]+d2
    i_new[i] <- i_new[i-1]+d1
    p1 <- 1-exp(-beta*tau*ii[i]/N)
  }
  ss <- matrix(c(ss,rep(NA,kk)),ncol=24,byrow=T); ii <- matrix(c(ii,rep(NA,kk)),ncol=24,byrow=T)
  rr <- matrix(c(rr,rep(NA,kk)),ncol=24,byrow=T);i_new <- matrix(c(i_new,rep(NA,kk)),ncol=24,byrow=T)
  return(list(S=ss[,1],I=ii[,1],R=rr[,1],I_new=i_new[,1]))
}

SIR_S_model_A <- function(initial,t,parameter){
  # calculate small change in time, \tau
  tau <- t[2]-t[1]
  # will be used later on to transform output to a matrix
  kk <- 24-length(t)%%24
  m <- length(t)
  # priors of the parameters
  beta <- parameter[1]
  gamma <- parameter[2]
  I <- parameter[3]
  # initial condition
  S <- initial[1]
  N <- initial[2]
  # initial probability values (p2 is fixed)
  p1 <- 1-exp(-beta*tau*I/N); p2 <- 1-exp(-gamma*tau)
  # initialise arrays for each compartment
  ss <- rep(0,length(t));ii <- rep(0,length(t));rr <- rep(0,length(t));i_new <- rep(0,length(t))
  # starting values for each compartment
  ss[1] <- floor(S); ii[1] <- floor(I); rr[1] <- N-floor(S)-floor(I); i_new <- N-floor(S)
  for(i in 2:m){
    # random binomial samples
    d1 <- rbinom(1,ss[i-1],p1)
    d2 <- rbinom(1,ii[i-1],p2)
    # update compartments
    ss[i] <- ss[i-1]-d1
    ii[i] <- ii[i-1]+d1-d2
    rr[i] <- rr[i-1]+d2
    i_new[i] <- i_new[i-1]+d1
    # update infection probability
    p1 <- 1-exp(-beta*tau*ii[i]/N)
  }
  # clean up the output into a matrix
  i_new <- matrix(c(i_new,rep(NA,kk)),ncol=24,byrow=T)
  # return the first column which corresponds to unit time of days
  return(array(i_new[,1]))
}

SIR_S_model_B <- function(initial,t,parameter){
  tau <- t[2]-t[1]
  kk <- 24-length(t)%%24
  m <- length(t)
  beta <- parameter[1]
  gamma <- parameter[2]
  I <- parameter[3]
  S <- initial[1]
  N <- initial[2]
  p1 <- 1-exp(-beta*tau*I/N);p2 <- 1-exp(-gamma*tau)
  ss <- rep(0,length(t));ii <- rep(0,length(t));rr <- rep(0,length(t));i_new <- rep(0,length(t))
  ss[1] <- floor(S); ii[1] <- floor(I); rr[1] <-N-floor(S)-floor(I); i_new <- N-floor(S)
  for(i in 2:m){
    d1 <- rbinom(1,ss[i-1],p1)
    d2 <- rbinom(1,ii[i-1],p2)
    ss[i] <- ss[i-1]-d1
    ii[i] <- ii[i-1]+d1-d2
    rr[i] <- rr[i-1]+d2
    i_new[i] <- i_new[i-1]+d1
    p1 <- 1-exp(-beta*tau*ii[i]/N)
  }
  ii <- matrix(c(ii,rep(NA,kk)),ncol=24,byrow=T)
  return(array(ii[,1]))
}
