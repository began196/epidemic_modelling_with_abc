library(deSolve)
library(EasyABC)
library(tidyverse)
library(loo)
library(ggplot2)
source("D:/OneDrive - University of Bristol/FRP_content/function.R")

eng <- read.csv(file = "D:/OneDrive - University of Bristol/FRP_content/Data/england_cases.csv", header = T)

eng_dcase <- rev(eng$newCasesBySpecimenDate)
eng_case <- rev(eng$cumCasesBySpecimenDate)
tail(eng,60)
first_ld <- 49
shops_reopened <- 133
tier_system <- 252
second_ld <- 276
second_ld_end <- 304
third_ld <- 336

N <- 56*10^6

plot(eng_case)
abline(v=c(79,108),lty=2)

# #30-49 28th Feb
# prior_eng <-  list(c("normal",0.5,0.01),c("normal",0.25,0.05),c("unif",10,41))
# initial_eng <- c(N-eng_case[30],N)
# t1 <- seq(30,49)
# 
# SIR_model_1 <- function(para){
#   return(SIR3_model(initial_eng,t1,para))
# }
# eng_abc_seq <- ABC_sequential(method = "Drovandi", model = SIR_model_1, prior = prior_eng, summary_stat_target = eng_case[30:49],
#                              nb_simul = 1e3,tolerance_tab=50,alpha=0.95)
# 
# parm_1 <- apply(eng_abc_seq$param,2,mean)
# plot(eng_case[30:49],col="red")
# points(SIR_model_1(parm_1))
# 
# #  0.5003510  0.2351489 23.8650561

# #49-78
# prior_eng2 <-  list(c("unif",0,0.5),c("unif",0,0.5),c("unif",500,4786))
# initial_eng2 <- c(N-eng_case[49],N)
# t2 <- seq(49,78)
# SIR_model_2 <- function(para){
#   return(SIR3_model(initial_eng2,t2,para))
# }
# 
# eng_abc_seq2 <- ABC_sequential(method = "Drovandi", model = SIR_model_2, prior = prior_eng2, summary_stat_target = eng_case[49:78],
#                               nb_simul = 1e4,tolerance_tab=10,alpha=0.9,c=0.001)
# parm_2 <- apply(eng_abc_seq2$param,2,mean)
# plot(eng_case[49:78],col="red")
# points(SIR_model_2(parm_2))
# # c(0.3575739,0.2829968,2339.6338583)
# 
# # prediction
# plot(eng_case[49:88],col="red",ylim=c(0,3.2e5))
# points(SIR3_model(initial_eng2,seq(49,88),parm_2))
# 
# par(mfrow = c(3,1))
# apply(eng_abc_seq2$param,2,hist,breaks=20)
# par(mfrow=c(1,1))

yy <- seq(0,0.4,length.out = 1e5)
plot(yy,dlnorm(yy,log(0.12),0.3),type="l",xlab=TeX(r"($I_0$)"),ylab=TeX(r"($\pi(I_0)$)"))

#79-108
prior_eng3 <-  list(c("lognormal",log(0.12),0.3),c("lognormal",log(0.15),0.3),c("normal",4e4,5e3))
initial_eng3 <-  c(N-eng_case[79],N)
t3 <- seq(79,108)
SIR_model_3 <- function(para){
  return(SIR3_model(initial_eng3,t3,para))
}

eng_abc_seq3 <- ABC_sequential(method = "Delmoral", model = SIR_model_3, prior = prior_eng3,summary_stat_target = eng_case[79:108],
                              nb_simul = 2e4,tolerance_target=10,alpha=0.1)

parm_3 <- apply(eng_abc_seq3$param,2,mean)

plot(eng_case[79:108],type="l",main="",ylab="Cumulative cases",xlab="Day")
points(SIR_model_3(parm_3),pch=20)
lines(SIR_model_3(parm_3),col="red")

rms_eng_d <- sqrt(sum((SIR_model_3(parm_3)-eng_case[79:108])^2))

rms_pred_eng_d <- sqrt(sum((SIR3_model(initial_eng3,seq(79,118),parm_3)[31:40]-eng_case[109:118])^2))

# prediction
plot(seq(30,40),eng_case[108:118],ylim=c(1.75e5,2.25e5),xlim=c(30,40),type="l",main="",ylab="Cumulative cases",xlab="Day",lwd=2)
lines(seq(30,40),SIR3_model(initial_eng3,seq(79,118),parm_3)[30:40],col="blue",lwd=1)
legend(31,2.15e5,legend=c("Actual","Predicted"),col=c("black","blue"),lty=1)

par(mfrow = c(3,1))
apply(eng_abc_seq3$param,2,hist,breaks=20)
par(mfrow=c(1,1))

plot(density(eng_abc_seq3$param[,1]),main="",xlim=c(0,0.3),col="red")
yy <- seq(0,0.3,length.out = 1e4)
lines(yy,dlnorm(yy,log(0.12),0.3))
abline(v=parm_3[1],lty=2,col="red")

plot(density(eng_abc_seq3$param[,2]),col="red",xlim=c(0,0.35),main="")
yy <- seq(0,0.35,length.out = 1e4)
lines(yy,dlnorm(yy,log(0.15),0.3))
abline(v=parm_3[2],lty=2,col="red")

plot(density(eng_abc_seq3$param[,3]),col="red",ylim=c(0,0.00008),main="")
yy <- seq(2e4,6.5e4,length.out = 1e4)
lines(yy,dnorm(yy,4e4,5e3))
abline(v=parm_3[3],lty=2,col="red")

#c(1.124876e-01,1.367132e-01,4.164897e+04)

#Stochastic
# # 30-49
# prior_eng_s <-  list(c("lognormal",log(0.5),0.25),c("lognormal",log(0.25),0.25),c("normal",26,1))
# initial_eng_s <- c(N-eng_case[30],N)
# t1_s <- seq(30,49,by=1/24)
# SIR_model_S1 <- function(para){
#   return(SIR_S_model_A(initial_eng_s,t1_s,para))
# }
# 
# eng_abc_seq_s <- ABC_sequential(method = "Drovandi", model = SIR_model_S1, prior = prior_eng_s,prior_test ="X2<X1" , summary_stat_target = eng_case[30:49],
#                               nb_simul = 10^3,tolerance_tab=10,alpha=0.9)
# 
# parm_1s <- apply(eng_abc_seq_s$param,2,mean)
# parm_1s
# plot(eng_case[30:49],col="red")
# points(SIR_model_S1(parm_1s))
# 
# par(mfrow=c(3,1))
# apply(eng_abc_seq_s$param,2,hist,breaks=20)
# par(mfrow=c(1,1))
# # 0.4982493  0.2467904 26.0079158

# # 49-78
# prior_eng_s2 <-  list(c("unif",0,0.5),c("unif",0,0.5),c("unif",500,4786))
# initial_eng_s2 <- c(N-eng_case[49],N)
# t2_s <- seq(49,78,by=1/24)
# SIR_model_S2 <- function(para){
#   return(SIR_S_model_A(initial_eng_s2,t2_s,para))
# }
# eng_abc_seq_s2 <- ABC_sequential(method = "Drovandi", model = SIR_model_S2, prior = prior_eng_s2,prior_test ="X2<X1",summary_stat_target = eng_case[49:78],
#                                 nb_simul = 1e4,tolerance_tab=10,alpha=0.8,c=0.001)
# 
# parm_2s <- apply(eng_abc_seq_s2$param,2,mean)
# plot(eng_case[49:78],col="red",ylim=c(0,1.8e5))
# points(SIR_model_S2(parm_2s))
# 
# # prediction
# plot(eng_case[49:88],col="red",ylim=c(0,5.5e5))
# points(SIR_S_model_A(initial_eng_s2,seq(49,88,by=1/24),parm_2s))
# 
# par(mfrow=c(3,1))
# apply(eng_abc_seq_s2$param,2,hist,breaks=20)
# par(mfrow=c(1,1))
# # 0.3336490    0.2358737 2811.1016692

# 79-108
prior_eng_s3 <-   list(c("lognormal",log(0.12),0.3),c("lognormal",log(0.15),0.3),c("normal",4e4,5e3))
initial_eng_s3 <- c(N-eng_case[79],N)
t3_s <- seq(79,108,by=1/24)
SIR_model_S3 <- function(para){
  return(SIR_S_model_A(initial_eng_s3,t3_s,para))
}

# eng_abc_seq_s3 <- ABC_sequential(method = "Drovandi", model = SIR_model_S3, prior = prior_eng_s3, summary_stat_target = eng_case[79:133],
#                                  nb_simul = 1e4,tolerance_tab=10,alpha=0.8,c=0.001)
eng_abc_seq_s3 <- ABC_sequential(method = "Delmoral", model = SIR_model_S3, prior = prior_eng_s3,summary_stat_target = eng_case[79:108],
                               nb_simul = 2e4,tolerance_target=10,alpha=0.1)

parm_3s <- apply(eng_abc_seq_s3$param,2,mean)

set.seed(123)
ggg <- matrix(0,nrow=5e3,ncol=40)
for(i in 1:5e3){
  ggg[i,] <- SIR_S_model_A(initial_eng_s3,seq(1,40,by=1/24),parm_3s)
}
plot(eng_case[79:108],main="",xlab="Day",ylab="Cumulative cases",type="l")
for(i in 1:5e3){
  lines(ggg[i,1:30],col="grey")
}
lines(mean_eng[1:30],col="red",lwd=2)
lines(eng_lower[1:30],col="blue",lty=2)
lines(eng_upper[1:30],col="blue",lty=2)
lines(eng_case[79:108],lwd=2)
legend(x=5,y=1.8e5,legend=c("Actual","Simulated","Mean","95% interval"),col=c("black","grey","red","blue"),lty=c(1,1,1,2))

mean_eng <- apply(ggg,2,mean)
sd_eng <- apply(ggg,2,sd)
eng_lower <- mean_eng - sd_eng*qnorm(0.975)
eng_upper <- mean_eng + sd_eng*qnorm(0.975)

plot(eng_case[79:108],main="",xlab="Day",ylab="Cumulative cases",type="l",lwd=2)
lines(mean_eng[1:30],col="red",lwd=2)
lines(eng_lower[1:30],col="blue",lty=2)
lines(eng_upper[1:30],col="blue",lty=2)
legend(x=1,y=1.8e5,legend=c("Actual","Mean","95% interval"),col=c("black","red","blue"),lty=c(1,1,2))

rms_eng_s <- rep(0,5e3)
for(i in 1:5e3){
  rms_eng_s[i] <- sum((ggg[i,1:30]-eng_case[79:108])^2)
}
rms_s_e <- sqrt(mean(rms_eng_s))

rms_eng_dummy <- rep(0,5e3)
for(i in 1:5e3){
  rms_eng_dummy[i] <- sum((ggg[i,31:40]-eng_case[109:118])^2)
}
rms_pred_eng_s <- sqrt(mean(rms_eng_dummy))

# prediction
plot(seq(30,40),eng_case[108:118],ylim=c(1.74e5,2.25e5),xlim=c(30,40),main="",xlab="Day",ylab="Cumulative cases",lwd=2,type="l")
for(i in 1:5e3){
  lines(seq(30,40),ggg[i,30:40],col="lightblue")
}
lines(seq(30,40),eng_case[108:118],col="black",lwd=2)
legend(31,2.2e5,legend=c("Actual","Predicted"),col=c("black","lightblue"),lty=1)

plot(seq(30,40),eng_case[108:118],ylim=c(1.74e5,2.25e5),xlim=c(30,40),main="",xlab="Day",ylab="Cumulative cases",lwd=2,type="l")
for(i in 1:5e3){
  lines(seq(30,40),ggg[i,30:40],col="lightblue")
}
lines(seq(30,40),mean_eng[30:40],col="red",lwd=2)
lines(seq(30,40),eng_upper[30:40],col="blue",lty=2)
lines(seq(30,40),eng_lower[30:40],col="blue",lty=2)
lines(seq(30,40),eng_case[108:118],col="black",lwd=2)
legend(30,2.25e5,legend=c("Actual","Predicted","Mean","95% interval"),col=c("black","lightblue","red","blue"),lty=c(1,1,1,2))

par(mfrow=c(3,1))
apply(eng_abc_seq_s3$param,2,hist,breaks=20)
par(mfrow=c(1,1))
#c(1.191581e-01,1.582567e-01,4.924302e+04)

plot(density(eng_abc_seq_s3$param[,1]),main="",xlim=c(0,0.3),col="red")
yy <- seq(0,0.3,length.out = 1e4)
lines(yy,dlnorm(yy,log(0.12),0.3))
abline(v=parm_3s[1],lty=2,col="red")

plot(density(eng_abc_seq_s3$param[,2]),col="red",xlim=c(0,0.35),main="")
yy <- seq(0,0.35,length.out = 1e4)
lines(yy,dlnorm(yy,log(0.15),0.3))
abline(v=parm_3s[2],lty=2,col="red")

plot(density(eng_abc_seq3$param[,3]),col="red",ylim=c(0,0.00008),main="")
yy <- seq(2e4,6.5e4,length.out = 1e4)
lines(yy,dnorm(yy,4e4,5e3))
abline(v=parm_3s[3],lty=2,col="red")

# Bristol
# nb = 463,400
# 25th March
Nb <- 463400
bris <- read.csv(file = "D:/OneDrive - University of Bristol/FRP_content/Data/bristol_cases.csv", header = T)
bris_dcase <- rev(bris$newCasesBySpecimenDate)
bris_case <- rev(bris$cumCasesBySpecimenDate)

# #30-59
# prior_bris <-  list(c("lognormal",log(0.12),0.5),c("lognormal",log(0.1),0.5),c("unif",50,167))
# initial_bris <- c(Nb-bris_case[30],Nb)
# t1_b <- seq(30,59)
# 
# SIR_model_1b <- function(para){
#   return(SIR3_model(initial_bris,t1_b,para))
# }
# bris_abc_seq <- ABC_sequential(method = "Drovandi", model = SIR_model_1b, prior = prior_bris, summary_stat_target = bris_case[30:59],
#                               nb_simul = 1e4,tolerance_tab=10,alpha=0.8,c=0.001)
# 
# parm_1b <- apply(bris_abc_seq$param,2,mean)
# plot(bris_case[30:59],col="red")
# points(SIR_model_1b(parm_1b))
# 
# # prediction
# plot(bris_case[30:69],col="red")
# points(SIR3_model(initial_bris,seq(30,69),parm_1b))
# 
# # c(0.1197970,0.1014756,126.4402751)
# par(mfrow=c(3,1))
# apply(bris_abc_seq$param,2,hist,breaks=30)
# par(mfrow=c(1,1))

xx <- seq(250,850,length.out = 1e4)
plot(xx,dunif(xx,300,812),type="l",xlab=TeX(r"($I_0$)"),ylab=TeX(r"($\pi(I_0)$)"))

# 60-89
prior_bris2 <- list(c("unif",0,0.1),c("unif",0,0.1),c("unif",300,812))
initial_bris2 <- c(Nb-bris_case[60],Nb)
t2_b <- seq(60,89)

SIR_model_2b <- function(para){
  return(SIR3_model(initial_bris2,t2_b,para))
}
bris_abc_seq2 <- ABC_sequential(method = "Delmoral", model = SIR_model_2b, prior = prior_bris2, summary_stat_target = bris_case[60:89],
                               nb_simul = 5e3,tolerance_target=10,alpha=0.1)

parm_2b <- apply(bris_abc_seq2$param,2,mean)

plot(bris_case[60:89],type="l",main="",xlab="Day",ylab="Cumulative cases")
points(SIR_model_2b(parm_2b),pch=20)
lines(SIR_model_2b(parm_2b),col="red")

rms_bris_d <- sqrt(sum((SIR_model_2b(parm_2b)-bris_case[60:89])^2))

rms_pred_bris_d <- sqrt(sum((SIR3_model(initial_bris2,seq(60,99),parm_2b)[31:40]-bris_case[90:99])^2))

# prediction
plot(seq(30,40),bris_case[89:99],type="l",main="",ylab = "Cumulative cases",xlab= "Day",xlim=c(30,40),ylim=c(1100,1300),lwd=2)
lines(seq(30,40),SIR3_model(initial_bris2,seq(60,99),parm_2b)[30:40],col="blue",lwd=1)
legend(31,1.275e3,legend=c("Actual","Predicted"),col=c("black","blue"),lty=1)

# c(0.03793094,0.07360420,476.66433196)

par(mfrow=c(3,1))
apply(bris_abc_seq2$param,2,hist,breaks=30)
par(mfrow=c(1,1))

plot(density(bris_abc_seq2$param[,1]),main="",xlim=c(0,0.15),col="red")
yy <- seq(-0.1,0.15,length.out = 1e4)
lines(yy,dunif(yy,0,0.1))
abline(v=parm_2b[1],lty=2,col="red")

plot(density(bris_abc_seq2$param[,2]),main="",xlim=c(-0.05,0.15),col="red")
yy <- seq(-0.2,0.15,length.out = 1e4)
lines(yy,dunif(yy,0,0.1))
abline(v=parm_2b[2],lty=2,col="red")

plot(density(bris_abc_seq2$param[,3]),col="red",xlim=c(200,900),main="")
yy <- seq(200,950,length.out = 1e4)
lines(yy,dunif(yy,300,812))
abline(v=parm_2b[3],lty=2,col="red")

# Stochastic
# #30-59
# prior_bris_s <-  list(c("lognormal",log(0.12),0.5),c("lognormal",log(0.1),0.5),c("unif",50,167))
# initial_bris_s <- c(Nb-bris_case[30],Nb)
# t1_b_s <- seq(30,59,by=1/24)
# 
# SIR_model_1b_s <- function(para){
#   return(SIR_S_model_A(initial_bris_s,t1_b_s,para))
# }
# 
# bris_abc_seq_s <- ABC_sequential(method = "Drovandi", model = SIR_model_1b_s, prior = prior_bris_s, summary_stat_target = bris_case[30:59],
#                                nb_simul = 1e4,tolerance_tab=10,alpha=0.8,c=0.001)
# 
# parm_1b_s <- apply(bris_abc_seq_s$param,2,mean)
# plot(bris_case[30:59],col="red")
# points(SIR_model_1b_s(parm_1b_s))
# 
# # prediction
# plot(bris_case[30:69],col="red")
# points(SIR_S_model_A(initial_bris_s,seq(30,69,by=1/24),parm_1b_s))
# 
# par(mfrow=c(3,1))3
# apply(bris_abc_seq_s$param,2,hist,breaks=30)
# par(mfrow=c(1,1))

# c(0.1250080,0.1024393,120)

# 60-89
prior_bris2_s <- list(c("unif",0,0.1),c("unif",0,0.1),c("unif",300,812))
initial_bris2_s <- c(Nb-bris_case[60],Nb)
t2_b_s <- seq(60,89,by=1/24)

SIR_model_2b_s <- function(para){
  return(SIR_S_model_A(initial_bris2_s,t2_b_s,para))
}
bris_abc_seq2_s <- ABC_sequential(method = "Delmoral", model = SIR_model_2b_s, prior = prior_bris2_s, summary_stat_target = bris_case[60:89],
                                nb_simul = 5e3,tolerance_target=10,alpha=0.1)
parm_2b_s <- apply(bris_abc_seq2_s$param,2,mean)
apply(bris_abc_seq2_s$param,2,quantile,prob=0.975)

set.seed(123)
bbb <- matrix(0,nrow=5e3,ncol=40)
for(i in 1:5e3){
  bbb[i,] <- SIR_S_model_A(initial_bris2_s,seq(1,40,by=1/24),parm_2b_s)
}
plot(bris_case[60:89],type="l",main="",ylab="Cumulative cases",xlab="Day",ylim=c(800,1300))
for(i in 1:5e3){
  lines(bbb[i,1:30],col="grey")
}
lines(mean_bris[1:30],col="red",lwd=2)
lines(bris_lower[1:30],lty=2,col="blue")
lines(bris_upper[1:30],lty=2,col="blue")
lines(bris_case[60:89],lwd=2)
legend(x=5,y=1250,legend=c("Actual","Simulated","Mean","95% interval"),col=c("black","grey","red","blue"),lty=c(1,1,1,2))

mean_bris <- apply(bbb,2,mean)
sd_bris <- apply(bbb,2,sd)
bris_upper <- mean_bris+ qnorm(0.975)*sd_bris
bris_lower <- mean_bris- qnorm(0.975)*sd_bris

plot(bris_case[60:89],type="l",main="",ylab="Cumulative cases",xlab="Day",ylim=c(800,1300),lwd=2)
lines(mean_bris[1:30],col="red",lwd=2)
lines(bris_lower[1:30],lty=2,col="blue")
lines(bris_upper[1:30],lty=2,col="blue")
legend(x=1,y=1300,legend=c("Actual","Mean","95% interval"),col=c("black","red","blue"),lty=c(1,1,2))

rms_bris_s <- rep(0,5e3)
for(i in 1:5e3){
  rms_bris_s[i] <- sum((bbb[i,1:30]-bris_case[60:89])^2)
}

rms_s_b <- sqrt(mean(rms_bris_s))

# zscore_eng <- matrix(0,nrow=5e3,ncol=30)
# for(i in 1:5e3){
#   zscore_eng[i,] <- (ggg[i,1:30]-mean_eng[1:30])/sd_eng[1:30]
# }
# boxplot(zscore_eng[,2:30])
# zscore_bris <- matrix(0,nrow=5e3,ncol=30)
# for(i in 1:5e3){
#   zscore_bris[i,] <- (bbb[i,1:30]-mean_bris[1:30])/sd_bris[1:30]
# }
# boxplot(zscore_bris[,2:30])

plot(seq(2,30),(sd_bris/mean_bris)[2:30],type="l",ylab="Coefficient of Variation",xlab="Day",col="red",ylim=c(0,0.025))
lines(seq(2,30),(sd_eng/mean_eng)[2:30],col="blue")
legend(2,0.0225,legend = c("England","Bristol"),lty=1,col=c("blue","red"))

# prediction
plot(seq(30,40),bris_case[89:99],main="",xlab="Day",ylab="Cumulative cases",lwd =2, xlim=c(30,40),ylim=c(1000,1400),type="l")
for(i in 1:5e3){
  lines(seq(30,40),bbb[i,30:40],col="lightblue")
}
lines(seq(30,40),bris_case[89:99],lwd=2)
legend(30,1350,legend=c("Actual","Predicted"),col=c("black","lightblue"),lty=1)

plot(seq(30,40),bris_case[89:99],main="",xlab="Day",ylab="Cumulative cases",lwd =2, xlim=c(30,40),ylim=c(1000,1400),type="l")
for(i in 1:5e3){
  lines(seq(30,40),bbb[i,30:40],col="lightblue")
}
lines(seq(30,40),mean_bris[30:40],col="red",lwd=2)
lines(seq(30,40),bris_lower[30:40],col="blue",lty=2)
lines(seq(30,40),bris_upper[30:40],col="blue",lty=2)
lines(seq(30,40),bris_case[89:99],lwd=2)
legend(x=30,y=1400,legend=c("Actual","Predicted","Mean","95% interval"),col=c("black","lightblue","red","blue"),lty=c(1,1,1,2))

dummy_rms_b <- rep(0,5e3)
for(i in 1:5e3){
  dummy_rms_b[i] <- sum((bris_case[90:99]-bbb[i,31:40])^2)
}
rms_pred_bris_s <- sqrt(mean(dummy_rms_b))

rms <- c(rms_eng_d,rms_s_e,rms_bris_d,rms_s_b)
rms_pred <- c(rms_pred_eng_d,rms_pred_eng_s,rms_pred_bris_d,rms_pred_bris_s)
rms_pred
# c(0.03917407,0.07887510,477.01244784)
par(mfrow=c(3,1))
apply(bris_abc_seq2_s$param,2,hist,breaks=30)
par(mfrow=c(1,1))

plot(density(bris_abc_seq2_s$param[,1]),main="",xlim=c(0,0.15),col="red")
yy <- seq(-0.1,0.15,length.out = 1e4)
lines(yy,dunif(yy,0,0.1))
abline(v=parm_2b_s[1],lty=2,col="red")

plot(density(bris_abc_seq2_s$param[,2]),main="",xlim=c(-0.05,0.15),col="red")
yy <- seq(-0.2,0.15,length.out = 1e4)
lines(yy,dunif(yy,0,0.1))
abline(v=parm_2b_s[2],lty=2,col="red")

plot(density(bris_abc_seq2_s$param[,3]),col="red",xlim=c(200,900),main="")
yy <- seq(200,950,length.out = 1e4)
lines(yy,dunif(yy,300,812))
abline(v=parm_2b_s[3],lty=2,col="red")

# global_case <- read.csv("D:/OneDrive - University of Bristol/FRP_content/Data/time_series_covid19_confirmed_global.csv",header = T)
# rec_case <- read.csv("D:/OneDrive - University of Bristol/FRP_content/Data/time_series_covid19_recovered_global.csv",header= T)
# death_case <- read.csv("D:/OneDrive - University of Bristol/FRP_content/Data/time_series_covid19_deaths_global.csv",header = T)
# #starts at 22nd Jan
# 
# msia_case <- global_case %>%filter(Country.Region=="Malaysia")
# colnames(msia_case) <- NULL
# msia_i <-as.array(apply(msia_case[1,5:ncol(msia_case)],2,as.numeric))
# 
# msia_rec <- rec_case %>%filter(Country.Region=="Malaysia")
# colnames(msia_rec) <- NULL
# msia_rec <-as.array(apply(msia_rec[1,5:ncol(msia_rec)],2,as.numeric))
# 
# msia_death <- death_case %>%filter(Country.Region=="Malaysia")
# colnames(msia_death) <- NULL
# msia_death <-as.array(apply(msia_death[1,5:ncol(msia_death)],2,as.numeric))
# 
# msia_r <- msia_death+msia_rec
# msia_di <- msia_i-msia_r
# plot(msia_di[313:342])
# sum(((msia_di[313:342])*0.01)^2+((msia_di[313:342])*0.01)^2)
# 
# prior_ms <- list(c("unif",0,0.5),c("unif",0,0.5),c("unif",5000,60485))
# tt1 <- seq(313,342)
# n <- 32.7*10^6
# initial_ms <- c(n-msia_i[313],n)
# SIR2_model_1 <- function(para){
#   SIRB_model(initial_ms,tt1,para)
# }
# 
# ms_smc <- ABC_sequential(method = "Drovandi", model = SIR2_model_1, prior=prior_ms,summary_stat_target = c((msia_di[313:342]),(msia_r[313:342])),
#                               nb_simul = 10^4,tolerance_tab=10,alpha=0.3,prior_test = "X1>X2")
# 
# par_1 <- (apply(ms_smc$param,2,median))
# 
# sum_stats1 <- SIRB_model(initial_ms,tt1,par_1)
# 
# par(mfrow=c(2,1))
# plot(msia_di[313:342])
# points(sum_stats1[1:30],col="red")
# plot(msia_r[313:342])
# points(sum_stats1[31:60],col="red")
# 
# par(mfrow=c(3,1))
# apply(ms_rej$param,2,hist,breaks= 20)
# par(mfrow=c(1,1))

# Suppose transmission process is stochastic
# N = 1e2
set.seed(700)
SIR_stoc1 <- SIR_S_model(c(1e2-1,1e2),seq(1,100,by=1/24),c(0.4,0.1,1))
matplot(x = seq(1,100,by=1), y =cbind(SIR_stoc1$S,SIR_stoc1$I,SIR_stoc1$R), type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "Stochastic SIR",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(65, 90, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")

plot(x=seq(1,100),y=SIR_stoc1$I,type="l",xlab="Day",ylab="Number of infected daily",lwd=2)

# First 20 days of daily infected as summary stats
ss_s1 <- SIR_stoc1$I[1:20]

# Deterministic model
SIR_D1 <- function(para){
  sir_model_ss(c(1e2-para[3],para[3],0),seq(1,20),c(para[1],para[2]))
}

# Stochastic model
# SIR_S1 <- function(par){
#   SIR_S_model_B(c(1e3-par[3],1e3),seq(1,20,by=1/24),par)
# }

# Prior
prior_D1 <- list(c("unif",0,1),c("unif",0,0.5),c("lognormal",log(1),3))
# ABC to get parameters
abc_D1 <- ABC_sequential("Delmoral",model=SIR_D1,prior_D1,nb_simul=2e3,summary_stat_target=ss_s1,tolerance_target=1,alpha=0.1)
# abc_S1 <- ABC_sequential("Delmoral",model=SIR_S1,prior_D1,nb_simul=2e3,summary_stat_target=ss_s1,tolerance_target=5,alpha=0.1)
parm_D1 <- apply(abc_D1$param,2,mean)
# parm_S1 <- apply(abc_S1$param,2,mean)

plot(seq(1,20),SIR_D1(parm_D1),type="l",col="blue",xlim=c(1,40),ylim=c(0,50),lwd=2,main="",ylab="Number of infected daily",xlab="Day")
lines(seq(1,40),SIR_stoc1$I[1:40],lwd=2)
lines(seq(20,40),sir_model_ss(c(1e2-parm_D1[3],parm_D1[3],0),seq(1,40),c(parm_D1[1],parm_D1[2]))[20:40],col="red",lwd=2)
legend(1,45,legend=c("Data","Model","Predicted"),lty=1,col=c("black","blue","red"),lwd=2)

# N = 1e7
set.seed(679)
SIR_stoc3 <- SIR_S_model(c(1e7-1,1e7),seq(1,100,by=1/24),c(0.4,0.1,1))
matplot(x = seq(1,100,by=1), y =cbind(SIR_stoc3$S,SIR_stoc3$I,SIR_stoc3$R), type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "Stochastic SIR",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(60, 1e7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")

# First 30 days of daily infected as summary stats
ss_s3 <- SIR_stoc3$I[35:55]

# Deterministic model
SIR_D3 <- function(para){
  return(sir_model_ss(c(1e7-para[3],para[3],0),seq(1,55),c(para[1],para[2]))[35:55])
}

# Prior
prior_D3 <- list(c("unif",0,1),c("unif",0,0.5),c("lognormal",log(1),3))
# ABC to get parameters
abc_D3 <- ABC_sequential("Delmoral",model=SIR_D3,prior_D3,nb_simul=2e3,summary_stat_target=ss_s3,tolerance_target=0.05,alpha=0.1)
parm_D3<- apply(abc_D3$param,2,mean)

plot(seq(1,50),sir_model_ss(c(1e7-parm_D3[3],parm_D3[3],0),seq(1,50),c(parm_D3[1],parm_D3[2])),
     type="l",col="blue",xlim=c(1,70),ylim=c(0,5e6),lwd=2,main="",ylab="Number of infected daily",xlab="Day")
lines(seq(1,70),SIR_stoc3$I[1:70],lwd=2)
lines(seq(50,70),sir_model_ss(c(1e7-parm_D3[3],parm_D3[3],0),seq(1,70),c(parm_D3[1],parm_D3[2]))[50:70],col="red",lwd=2)
legend(1,4.8e6,legend=c("Data","Model","Predicted"),lty=1,col=c("black","blue","red"),lwd=2)

# with data model
prior_eng_d <-  list(c("unif",0,0.2),c("unif",0,0.2),c("normal",4e4,5e3),c("unif",0,1))
initial_eng_d <-  c(N-eng_case[79],N)
t_d <- seq(79,108)
SIR_model_d <- function(para){
  sim <- (SIR3_model(initial_eng_d,t_d,para[1:3]))
  obs <- rnbinom(length(sim),sim,para[4])
  return(obs)
}

eng_abc_seq_d <- ABC_sequential(method = "Delmoral", model = SIR_model_d, prior = prior_eng_d,summary_stat_target = (eng_case[79:108]),
                               nb_simul = 5e3,tolerance_target=0.01,alpha=0.001)

parm_d <- apply(eng_abc_seq_d$param,2,mean)

plot(density(eng_abc_seq_d$param[,1]))

plot(eng_case[79:108],type="l",ylim=c(9e4,2.1e5))
lines(SIR_model_d(parm_d),col="red")


# mode c(0.0725249,0.08766472,40789,0.4963527)
# c(0.1180783,0.1472157,41869.0029579,0.4998256)

prior_eng_ds <- list(c("unif",0,0.2),c("unif",0,0.2),c("normal",4e4,5e3),c("unif",0,1))
initial_eng_ds <- c(N-eng_case[79],N)
t_ds <- seq(79,108,by=1/24)
SIR_model_ds <- function(para){
  sim <- (SIR_S_model_A(initial_eng_ds,t_ds,para[1:3]))
  obs <- rnbinom(length(sim),size=sim,para[4])
  return(obs)
}

eng_abc_seq_ds <- ABC_sequential(method = "Delmoral", model = SIR_model_ds, prior = prior_eng_ds,summary_stat_target = (eng_case[79:108]),
                                 nb_simul = 5e3,tolerance_target=0.01,alpha=0.001)

parm_ds <- apply(eng_abc_seq_ds$param,2,mean)

plot(eng_case[79:108],type="l",ylim=c(9.5e4,2.1e5))
lines(SIR_model_ds(parm_ds),col="red")



prior_bris_d <- list(c("unif",0,0.2),c("unif",0,0.2),c("unif",300,812))
initial_bris_d <- c(Nb-bris_case[60],Nb)
t_d2 <- seq(60,89)

SIR_model_bris <- function(para){
  sim <- (SIR3_model(initial_bris_d,t_d2,para[1:3]))
  dif <- diff(sim)
  # obs <- rnbinom(length(dif),dif,para[4])
  obs <- rpois(length(dif),dif)
  return(obs)
}


bris_abc_d <- ABC_sequential(method = "Delmoral", model = SIR_model_bris, prior = prior_bris_d, summary_stat_target = diff(bris_case[60:89]),
                                nb_simul = 5e3,tolerance_target=0.001,alpha=0.2)

parm_b <- apply(bris_abc_d$param,2,mean)

plot(density(bris_abc_d$param[,2]))

bris_ud <- cumsum(c(bris_case[60],SIR_model_bris(parm_b)))

plot(bris_case[60:89],type="l")
lines(bris_ud,col="red")

prior_bris_ds <- list(c("unif",0,0.2),c("unif",0,0.2),c("unif",300,812),c("unif",0,1))
initial_bris_ds <- c(Nb-bris_case[60],Nb)
t_ds2 <- seq(60,89,by=1/24)

SIR_model_bris_s <- function(para){
  sim <- (SIR_S_model_A(initial_bris_ds,t_ds2,para[1:3]))
  dif <- diff(sim)
  obs <- rnbinom(length(dif),dif,para[4])
  # obs <- rpois(length(dif),dif*para[4])
  return(obs)
}

bris_abc_ds <- ABC_sequential(method = "Delmoral", model = SIR_model_bris_s, prior = prior_bris_ds, summary_stat_target = diff(bris_case[60:89]),
                             nb_simul = 5e3,tolerance_target=0.001,alpha=0.001)

parm_bs <- apply(bris_abc_ds$param,2,mean)




# fff <- matrix(0,ncol=1,nrow=5e3)
# for(i in 1:5e3){
#   fff[i] <- sqrt(sum((SIR_model_ds(eng_abc_seq_ds$param[i,])-eng_case[79:108])^2))
# }
# 
# sort(fff,index.return=T)$ix[1]




# Different R_0 value
# N = 1e3
# R_0 < 1
# Stochastic
set.seed(123)
S_e1 <- matrix(0,ncol=50,nrow=1e3)
for(i in 1:1e3){
  S_e1[i,] <- SIR_S_model(c(1e3-1,1e3),seq(1,50,by=1/24),c(0.2,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0, 50), ylim = c(0, 10))
for(i in 1:1e3){
  lines(S_e1[i,],col="lightblue")
}
duration <- rep(0,1e3)
when_peak <- rep(0,1e3)
peak_i <- rep(0,1e3)
total_i <- rep(0,1e3)
for(i in 1:1e3){
  peak_i[i] <- max(S_e1[i,])
  total_i[i] <- sum(S_e1[i,])
  duration[i]<- match(0,S_e1[i,])
  when_peak[i] <- match(peak_i[i],S_e1[i,])
}

set.seed(123)
S_e2 <- matrix(0,ncol=50,nrow=1e3)
for(i in 1:1e3){
  S_e2[i,] <- SIR_S_model(c(1e7-1,1e7),seq(1,50,by=1/24),c(0.2,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0, 50), ylim = c(0, 10))
for(i in 1:1e3){
  lines(S_e2[i,],col="lightblue")
}
duration_b <- rep(0,1e3)
when_peak_b <- rep(0,1e3)
peak_ib <- rep(0,1e3)
total_ib <- rep(0,1e3)
for(i in 1:1e3){
  peak_ib[i] <- max(S_e2[i,])
  total_ib[i] <- sum(S_e2[i,])
  duration_b[i]<- match(0,S_e2[i,])
  when_peak_b[i] <- match(peak_ib[i],S_e2[i,])
}

# Deterministic
D_e1 <- as.data.frame(ode(y = c(S=1e3-1,I=1,R=0), times = seq(1,50), func = sir,parms = c(0.2,0.5)))
D_e1$time <- NULL
plot(seq(1,50),D_e1$I,type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,2))
print(c(1,floor(sum(D_e1$I[D_e1$I>0])),3,1))

D_f1 <- as.data.frame(ode(y = c(S=1e7-1,I=1,R=0), times = seq(1,50), func = sir,parms = c(0.2,0.5)))
D_f1$time <- NULL
plot(seq(1,50),D_f1$I,type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,2))
print(c(1,floor(sum(D_f1$I[D_f1$I>0])),3,1))

# R_0 = 1
# Stochastic
set.seed(123)
S_f1 <- matrix(0,ncol=80,nrow=1e3)
for(i in 1:1e3){
  S_f1[i,] <- SIR_S_model(c(1e3-1,1e3),seq(1,80,by=1/24),c(0.5,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0, 80), ylim = c(0, 60))
for(i in 1:1e3){
  lines(S_f1[i,],col="lightblue")
}
duration2 <- rep(0,1e3)
when_peak2 <- rep(0,1e3)
peak_i2 <- rep(0,1e3)
total_i2 <- rep(0,1e3)
for(i in 1:1e3){
  peak_i2[i] <- max(S_f1[i,])
  total_i2[i] <- sum(S_f1[i,])
  duration2[i]<- match(0,S_f1[i,])
  when_peak2[i] <- match(peak_i2[i],S_f1[i,])
}

set.seed(123)
S_f2 <- matrix(0,ncol=800,nrow=1e3)
for(i in 1:1e3){
  S_f2[i,] <- SIR_S_model(c(1e7-1,1e7),seq(1,800,by=1/24),c(0.5,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0,800), ylim = c(0, 4e3))
for(i in 1:1e3){
  lines(S_f2[i,],col="lightblue")
}
duration2b <- rep(0,1e3)
when_peak2b <- rep(0,1e3)
peak_i2b <- rep(0,1e3)
total_i2b <- rep(0,1e3)
for(i in 1:1e3){
  peak_i2b[i] <- max(S_f2[i,])
  total_i2b[i] <- sum(S_f2[i,])
  if(min(S_f2[i,])!=0){
    duration2b[i] <- 800
  }else{duration2b[i]<- match(0,S_f2[i,])}
  when_peak2b[i] <- match(peak_i2b[i],S_f2[i,])
} 

# Deterministic
D_e2 <- as.data.frame(ode(y = c(S=1e3-1,I=1,R=0), times = seq(1,200), func = sir,parms = c(0.5,0.5)))
D_e2$time <- NULL
plot(seq(1,200),D_e2$I,type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,2))
print(c(1,floor(sum(D_e2$I)),86,1))

D_f2 <- as.data.frame(ode(y = c(S=1e7-1,I=1,R=0), times = seq(1,100), func = sir,parms = c(0.5,0.5)))
D_f2$time <- NULL
plot(seq(1,100),D_f2$I,type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,2))
print(c(NA,NA,NA,NA))

# R_0>1
# Stochastic
set.seed(123)
S_g1 <- matrix(0,ncol=60,nrow=1e3)
for(i in 1:1e3){
  S_g1[i,] <- SIR_S_model(c(1e3-1,1e3),seq(1,60,by=1/24),c(1.5,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0, 60), ylim = c(0, 500))
for(i in 1:1e3){
  lines(S_g1[i,],col="lightblue")
}
duration3 <- rep(0,1e3)
when_peak3 <- rep(0,1e3)
peak_i3 <- rep(0,1e3)
total_i3 <- rep(0,1e3)
for(i in 1:1e3){
  peak_i3[i] <- max(S_g1[i,])
  total_i3[i] <- sum(S_g1[i,])
  duration3[i]<- match(0,S_g1[i,])
  when_peak3[i] <- match(peak_i3[i],S_g1[i,])
}

set.seed(123)
S_g2 <- matrix(0,ncol=60,nrow=1e3)
for(i in 1:1e3){
  S_g2[i,] <- SIR_S_model(c(1e7-1,1e7),seq(1,60,by=1/24),c(1.5,0.5,1))$I
}
plot(1, type = "n",main= "",xlab = "Day", ylab = "Number of infected daily",xlim = c(0, 60), ylim = c(0, 1e6))
for(i in 1:1e3){
  lines(S_g2[i,],col="lightblue")
}
duration3a <- rep(0,1e3)
when_peak3a <- rep(0,1e3)
peak_i3a <- rep(0,1e3)
total_i3a <- rep(0,1e3)
for(i in 1:1e3){
  peak_i3a[i] <- max(S_g2[i,])
  total_i3a[i] <- sum(S_g2[i,])
  duration3a[i]<- match(0,S_g2[i,])
  when_peak3a[i] <- match(peak_i3a[i],S_g2[i,])
}

# Deterministic
D_e3 <- as.data.frame(ode(y = c(S=1e3-1,I=1,R=0), times = seq(1,1e2), func = sir,parms = c(1.5,0.5)))
D_e3$time <- NULL
plot(seq(1,50),D_e3$I[1:50],type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,300))
print(c(floor(max(D_e3$I)),floor(sum(D_e3$I)),27,9))

D_f3 <- as.data.frame(ode(y = c(S=1e7-1,I=1,R=0), times = seq(1,1e2), func = sir,parms = c(1.5,0.5)))
D_f3$time <- NULL
plot(seq(1,40),D_f3$I[1:40],type="l",main="",xlab="Day",ylab="Number of infected daily",ylim=c(0,3e6))
print(c(floor(max(D_f3$I)),floor(sum(D_f3$I)),57,18))

par(mfrow=c(2,2))
hist((peak_i2),main="Peak infected",xlab="")
abline(v=1,lty=2,col="red")
hist((total_i2),main="Final size",xlab="")
abline(v=NA,lty=2,col="red")
hist((duration2),main="Epidemic duration",xlab="")
abline(v=NA,lty=2,col="red")
hist(when_peak2,main="Peak day",xlab="")
abline(v=NA,lty=2,col="red")

# Vary initial infected
ini_inf <- seq(1,15,by=1)
det <- matrix(0,ncol=4,nrow=15)
for(i in 1:15){
  asf <-  as.data.frame(ode(y = c(S=1e3-ini_inf[i],I=ini_inf[i],R=0), times = seq(1,1e2), func = sir,parms = c(1.5,0.5)))
  det[i,] <- c(floor(max(asf$I)),floor(sum(asf$I)),length(asf$I[asf$I>=0.5]),match(max(asf$I),asf$I))
}

set.seed(123)
stoc <- matrix(0,ncol=5,nrow=15)
for(i in 1:15){
  kka <- matrix(0,nrow=1e3,ncol=4)
  count <- 0
  for(j in 1:1e3){
    efw <- SIR_S_model(c(1e3-ini_inf[i],1e3),seq(1,30,by=1/24),c(1.5,0.5,ini_inf[i]))$I
    kka[j,] <- c(max(efw),sum(efw),length(efw[efw>0])+1,match(max(efw),efw))
    if(max(efw)==1){
      count <- count+1
    }
  }
  stoc[i,] <- c(apply(kka,2,mean),(1e3-count)/1e3)
}

set.seed(123)
stoc_2 <- matrix(0,ncol=5,nrow=15)
for(i in 1:15){
  kka_2 <- matrix(0,nrow=1e3,ncol=4)
  count <- 0
  for(j in 1:1e3){
    efw_2 <- SIR_S_model(c(1e7-ini_inf[i],1e7),seq(1,70,by=1/24),c(1.5,0.5,ini_inf[i]))$I
    kka_2[j,] <- c(max(efw_2),sum(efw_2),length(efw_2[efw_2>0])+1,match(max(efw_2),efw_2))
    if(max(efw_2)==1){
      count <- count+1
    }
  }
  stoc_2[i,] <- c(apply(kka_2,2,mean),(1e3-count)/1e3)
}

par(mfrow=c(2,2))
plot(stoc[,1],main="Peak infected",type="l",ylab="",xlab="No. of initial infected")
lines(det[,1],col="red")
plot(stoc[,2],main="Final size",type="l",ylab="",xlab="No. of initial infected")
lines(det[,2],col="red")
plot(stoc[,3],main="Duration of epidemic",type="l",ylim=c(17,27),ylab="",xlab="No. of initial infected")
lines(det[,3],col="red")
plot(stoc[,4],main="Peak day",type="l",ylim=c(5,10),ylab="",xlab="No. of initial infected")
lines(det[,4],col="red")
par(mfrow=c(1,1))

plot(stoc[,5])

# Vary recovery time
rec_time <- seq(1,15,by=1)
det_rec <- matrix(0,ncol=4,nrow=15)
for(i in 1:15){
  asfd <-  as.data.frame(ode(y = c(S=1e3-1,I=1,R=0), times = seq(1,1.5e2), func = sir,parms = c(1.5,1/rec_time[i])))
  det_rec[i,] <- c(floor(max(asfd$I)),floor(sum(asfd$I)),length(asfd$I[asfd$I>=0.5]),match(max(asfd$I),asfd$I))
}

set.seed(123)
stoc_rec <- matrix(0,ncol=5,nrow=15)
for(i in 1:15){
  kka3 <- matrix(0,nrow=1e3,ncol=4)
  count <- 0
  for(j in 1:1e3){
    efw3 <- SIR_S_model(c(1e3-1,1e3),seq(1,200,by=1/24),c(1.5,1/rec_time[i],1))$I
    kka3[j,] <- c(max(efw3),sum(efw3),length(efw3[efw3>0])+1,match(max(efw3),efw3))
    if(max(efw3)==1){
      count <- count+1
    }
  }
  stoc_rec[i,] <- c(apply(kka3,2,mean),(1e3-count)/1e3)
}

par(mfrow=c(2,2))
plot(stoc_rec[,1],main="Peak infected",type="l",ylab="",xlab="Recovery time (Days)",ylim=c(0,820))
lines(det_rec[,1],col="red")
plot(stoc_rec[,2],main="Final size",type="l",ylab="",xlab="Recovery time (Days)",ylim=c(0,1.5e4))
lines(det_rec[,2],col="red")
plot(stoc_rec[,3],main="Duration of epidemic",type="l",ylab="",xlab="Recovery time (Days)",ylim=c(0,120))
lines(det_rec[,3],col="red")
plot(stoc_rec[,4],main="Peak day",type="l",ylab="",xlab="Recovery time (Days)",ylim=c(4,13))
lines(det_rec[,4],col="red")
par(mfrow=c(1,1))

plot(stoc_rec[,5],type="p",ylab="Probability",xlab="Recovery time (Days)")
lines(stoc_rec[,5])
