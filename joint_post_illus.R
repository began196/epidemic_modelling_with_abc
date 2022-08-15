library(latex2exp)
library(ggplot2)
library(gridExtra)
set.seed(100)
data <- rnorm(100,5,3)

ss_mean <- mean(data) 
N <- 1e5

theta <- rep(0,N)
simulated <- rep(0,N)
euclid_dist <- rep(0,N)

for(i in 1:N){
  theta[i] <- runif(1,0,15)
  simulated[i] <- mean(rnorm(100,theta[i],3))
  euclid_dist[i] <- abs(simulated[i]-ss_mean)
}

ed <- sort(euclid_dist,index.return=T)
parm_ed <- theta[ed$ix[1:692]]
ss <- simulated[ed$ix[1:692]]
dist <- ss-ss_mean

linear_m <- lm(parm_ed~dist)
bb <- linear_m$coefficients[-1]
aa <- linear_m$coefficients[1]
adj_parm <- parm_ed-bb*dist


scatter <- ggplot()+geom_point(aes(y=theta[ed$ix[1:1e3]],x=simulated[ed$ix[1:1e3]]))
hist_right <- ggplot()+geom_histogram(aes(parm_ed))+coord_flip()+xlab(TeX(r"($\theta$)"))

scat <- scatter+geom_vline(xintercept = c(4.95,5,5.05), color = c("blue","red","blue"))+xlab(TeX(r"(Summary statistics, $S(y)$)"))+ylab(TeX(r"(Model parameters, $\theta$)"))+geom_segment(aes(x = 4.95, y = 3.8, xend = 5, yend = 3.8),arrow = arrow(length = unit(0.5, "cm")))
scat2 <- scat+geom_segment(aes(x = 5.05, y = 3.8, xend = 5, yend = 3.8),arrow = arrow(length = unit(0.5, "cm")))+geom_segment(aes(x = 5, y = 3.8, xend = 5.05, yend = 3.8),arrow = arrow(length = unit(0.5, "cm")))+geom_segment(aes(xend = 4.95, yend = 3.8, x = 5, y = 3.8),arrow = arrow(length = unit(0.5, "cm")))
scat3 <- scat2+annotate(geom="text",x=4.975,y=3.85,label=TeX(r"($-\epsilon$)"),color="black")+annotate(geom="text",x=5.025,y=3.85,label=TeX(r"($+\epsilon$)"),color="black")
grid.arrange(scat3,hist_right, ncol=2, nrow=1, widths=c(4, 1))

plot(y=theta[ed$ix[1:1e3]],x=simulated[ed$ix[1:1e3]],xlab=TeX(r"(Summary statistics, $S(y)$)"),ylab=TeX(r"(Model parameters, $\theta$)"))
abline(v=c(4.95,5.05,5),col=c("blue","blue","red"),lty=2)
arrows(4.95,3.8,5,3.8,code=3,length=0.1)
arrows(5,3.8,5.05,3.8,code=3,length=0.1)
text(4.975,3.85,TeX(r"($-\epsilon$)"))
text(5.025,3.85,TeX(r"($+\epsilon$)"))

plot(density(parm_ed))
lines(density(adj_parm),col="red")
