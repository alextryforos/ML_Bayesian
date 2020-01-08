library(msm)
colors = c("red","blue","green") #colors for intercept, slope, and variance respectively
data = read.csv("Data_FinalProject.csv") #39 observations of Height/Weight Data for men
x = data$Height #Create data vector for height
y = data$Weight #Create data vector for weight
#Intercept currently has meaning of "When weight is 0, average height is -intercept-"
x = scale(x,center=TRUE,scale=FALSE) #center the data to change interpretation of intercept
#Intercept NOW has meaning of "When weight is at the mean, average height is -intercept-"
#Iterative process using monte carlo approximation to find choice of prior
set.seed(1)
#Step 1, Select shape parameter
shape_B0 = 185.93 #prior shape
#Step 2, Select rate parameter that satisfies shape/rate = 175
rate_B0 = 1.0625 #prior rate
#Step 3, Generate Random Data
sample=rgamma(10000,shape=shape_B0,rate=rate_B0) 
#Step 4, calculate probability between 150 & 200, want to be approximately 0.99
mean(sample>150 & sample<200)
#Creating prior function for MCMC
prior_B0 = function(B0,shape_B0,rate_B0){
  dgamma(B0,shape = shape_B0, rate = rate_B0)
}
curve(dgamma(x,shape_B0,rate=rate_B0),xlim = c(0,270),xlab = "Intercept (B0)", ylab = "Density", main ="Gamma Prior", col=colors[1])
#slope prior (normal)
#would like mean of slope to be 8 (increase inch by 1 unit leads 
#to a 8 pound increase, on average, in weight)
mu_B1 = 8 #prior mean
sd_B1 = 5 #prior sd
#Prior function for MCMC
prior_B1 = function(B1,mu_B1,sd_B1){
  dnorm(B1,mu_B1,sd_B1)
}
curve(dnorm(x,mean=mu_B1,sd = sd_B1),xlim = c(mu_B1-2.5*sd_B1,mu_B1+2.5*sd_B1),xlab = "Slope (B1)", ylab = "Density", main = "Normal Prior", col=colors[2])
#variance prior (gamma)
shape_var = .5 #prior mean
rate_var = .01 #prior sd
#Prior function for MCMC
prior_var = function(var,shape_var,rate_var){
  dgamma(var,shape=shape_var,rate=rate_var)
}
curve(dgamma(x,shape_var,rate=rate_var),xlim = c(0,1000),xlab = "Variance", ylab = "Density", main ="Gamma Prior", col=colors[3])
library(png)
library(grid)
img_lik <- readPNG("Liklihood-Function.png")
grid.raster(img_lik)
liklihood = function(y,x,theta_B0 ,theta_B1,theta_var){
  sum(dnorm(y,theta_B0+theta_B1*x,sqrt(theta_var)))
}
#prior predictive check
intercept = c()
slope = c()
plot(y~x,cex=0.0001, ylim=c(100,300))
for (i in 1:15) {
  intercept[i] = rgamma(1,shape=shape_B0,rate=rate_B0)
  slope[i] = rnorm(1,mean=mu_B1,sd=sd_B1)
  abline(a=intercept[i],b=slope[i],col="black")
}
#Housekeeping for MCMC
theta_B0 = 175 #initial value for B0 (intercept), chosen as mean of prior
theta_B1 = 8 #intial value for B1 (slope), chosen as mean of prior
theta_var = 100 # inital value for variance, chosen as mean of prior

n.sim = 15000 #number of iterations during MCMC

n.accept_B0 = 0 #initialize acceptance count for intercept
n.accept_B1 = 0 #initialize acceptance count for slope
n.accept_var = 0 #initialize acceptance count for variance

delta_B0 = 15 #sd of proposal distribution for intercept, 50 got 30%
delta_B1 = 7 #sd of proposal distribution for slope, 15 got 30%
delta_var = 3 #sd of proposal distribution for variance, 25 got 30%

int_mcmc = NULL #vector to store intercept values during MCMC
slope_mcmc = NULL #vector to store slope values during MCMC
var_mcmc = NULL #vector to store variance values during MCMC

#MCMC
for (ite in 1:n.sim) {
  ################ Metropolis-Hastings within GIBBS for intercept  ################### 
  #proposal for intercept
  int.star = rtnorm(1,mean=theta_B0,sd=delta_B0,lower = 0, upper = 1300) #upper bound as 1300 since heaviest person in history is 1300
  #calculate posterior density ratio
  post_B0 = (prior_B0(B0 = int.star,shape_B0,rate_B0)*liklihood(y,x,int.star ,theta_B1,theta_var))/
    (prior_B0(B0 = theta_B0,shape_B0,rate_B0)*liklihood(y,x,theta_B0 ,theta_B1,theta_var))
  #calculate proposal density ratio
  proposal_B0 = dtnorm(theta_B0,mean=int.star,sd=delta_B0,lower = 0, upper = 1300)/
    dtnorm(int.star,mean=theta_B0,sd=delta_B0,lower = 0, upper = 1300)
  #calculate acceptance probability
  alpha_int = min(1,post_B0*proposal_B0)
  #Accept or reject proposal?
  r = runif(1)
  if (r<alpha_int){ #if accepted
    theta_B0 = int.star #update current position of theta_B0 to proposed position
    n.accept_B0 = n.accept_B0 + 1 #updaate number of acceptance of intercept
  }
  int_mcmc = c(int_mcmc,theta_B0) #update int_mcmc vector
  
  ###################### Metropolis within GIBBS for slope  ###################### 
  #proposal for slope
  slope.star = rnorm(1,mean=theta_B1,sd=delta_B1)
  #calculate posterior density ratio
  post_B1 = (prior_B1(B1 = slope.star,mu_B1,sd_B1)*liklihood(y,x,theta_B0 ,theta_B1=slope.star,theta_var))/
    (prior_B1(B1 = theta_B1,mu_B1,sd_B1)*liklihood(y,x,theta_B0 ,theta_B1,theta_var))
  #calculate acceptance probability
  alpha_slope = min(1,post_B1)
  #Accept or reject proposal?
  r = runif(1)
  if (r<alpha_slope){ #if accepted
    theta_B1 = slope.star #update current position of theta_B1 to proposed position
    n.accept_B1 = n.accept_B1 + 1 #updaate number of acceptance of slope
  }
  slope_mcmc = c(slope_mcmc,theta_B1) #update slope_mcmc vector
  
  ################ Metropolis-Hastings within GIBBS for variance  ################### 
  #proposal for variance
  var.star = rtnorm(1,mean=theta_var,sd=delta_var,lower = 0)
  #calculate posterior density ratio
  post_var = (prior_var(var = var.star,shape_var,rate_var)*liklihood(y,x,theta_B0 ,theta_B1,theta_var=var.star))/
    (prior_var(var = theta_var,shape_var,rate_var)*liklihood(y,x,theta_B0 ,theta_B1,theta_var))
  #calculate proposal density ratio
  proposal_var = dtnorm(theta_var,mean=var.star,sd=delta_var,lower = 0)/
    dtnorm(var.star,mean=theta_var,sd=delta_var,lower = 0)
  #calculate acceptance probability
  alpha_var = min(1,post_var*proposal_var)
  #Accept or reject proposal?
  r = runif(1)
  if (r<alpha_var){ #if accepted
    theta_var = var.star #update current position of theta_var to proposed position
    n.accept_var = n.accept_var + 1 #updaate number of acceptance of intercept
  }
  var_mcmc = c(var_mcmc,theta_var) #update var_mcmc vector
}
#Checking acceptance rates & trace plots
theta.mcmc = list(int_mcmc,slope_mcmc,var_mcmc) #MCMC vector
names(theta.mcmc) = c("Intercept","Slope","Variance")
a = c(n.accept_B0,n.accept_B1,n.accept_var)/n.sim #acceptance rate vector
names(a) = c("Intercept","Slope","Variance")
for (i in 1:length(a)) {
  plot(theta.mcmc[[i]],type = "l", main = paste("Trace Plot for", names(theta.mcmc)[i]), xlab = "Iteration", ylab=paste(names(a)[i]), col=colors[i])
  print(paste("The acceptance rate for", names(a)[i], "is", round(a[i],digits = 2)))
}
#burn in period of 3000
burnin = 3000
int_mcmc <- int_mcmc[(burnin+1):n.sim]
slope_mcmc <- slope_mcmc[(burnin+1):n.sim]
var_mcmc <- var_mcmc[(burnin+1):n.sim]
for (i in 1:length(theta.mcmc)) {
  print(paste("----------------- Posterior Information on: The",names(theta.mcmc[i]),"----------------"))
  print(paste("The 95 % credible interval for the",names(theta.mcmc)[i]," is from",round(quantile(theta.mcmc[[i]],c(.025,.975))[1],digits=2), "to", round(quantile(theta.mcmc[[i]],c(.025,.975))[2],digits=2), "."))
  print(paste("The 99 % credible interval for the",names(theta.mcmc)[i]," is from",round(quantile(theta.mcmc[[i]],c(.005,.995))[1],digits=2), "to", round(quantile(theta.mcmc[[i]],c(.005,.995))[2],digits=2), "."))
  print(paste("The mean of the of posterior distribution for the",names(theta.mcmc)[i],"is: ", round(mean(theta.mcmc[[i]]),digits=3), "."))
  print(paste("The median of the of posterior distribution for the",names(theta.mcmc)[i],"is: ", round(median(theta.mcmc[[i]]),digits=3), "."))
  print(paste("The variance of the of posterior distribution for the",names(theta.mcmc)[i]," is: ", round(var(theta.mcmc[[i]]),digits=3), "."))
  hist(theta.mcmc[[i]], main=paste("Histogram of",names(theta.mcmc)[i]), xlab = "", col=colors[i])
}
freqmodel = lm(y~x)
fc = c(freqmodel$coefficients, sqrt(sum(freqmodel$residuals^2)/(length(freqmodel$residuals)-length(freqmodel$coefficients))))
names(fc) = c("Intercept", "Slope", "Mean Squared Error")
bc = c(mean(int_mcmc),mean(slope_mcmc),mean(var_mcmc))
names(bc) = c("Intercept", "Slope", "Variance")
for (i in 1:length(fc)) {
  print(paste("The Frequentist estimate for", names(fc)[i], "is: ", round(fc[i],digits = 2)))
}
print("-------------------------------------------------------")
for (i in 1:length(fc)) {
  print(paste("A Bayesian estimate for", names(bc)[i], "is: ", round(bc[i],digits = 2)))
}