############################################
##                                        ##
##       Residuals and Diagnostics        ##
##   for Multinomial Regression Models    ##
##                                        ##
##  Eric A. E. Gerber and Bruce A. Craig  ##
##                                        ##
############################################

####Packages################################
##  Packages  ##
library(car)          #for confidence envelope in qqPlot
library(doParallel)   #for parallel computing
library(dirmult)      #for estimating dirichlet parameters
library(EnvStats)     #for uniform qqPlot
library(iterators)    #for some parallel options
library(MCMCpack)     #for inverse wishart sampling
library(MGLM)         #for multinomial logistic regression with counts
library(mvnfast)      #for multivariate normal sampling

####Real Data Sets##########################
##  Real Data Sets  ##
require(Lahman)       #for the baseball data
require(MM)           #for the pollen data
require(readxl)       #for reading in bird data

####Functions and Algorithms################
##  Functions and Algorithms from Separate (Messier) File  ##
source("mln_functions.R")

#  Function for Converting Counts to Log-Odds  #
newinitialY <- function(W, base){
  Wmat <- as.matrix(W)
  N <- nrow(Wmat)
  K <- ncol(Wmat)
  M <- apply(Wmat, 1, sum)

  tempZmat <- Wmat/M
  Ymat <- matrix(0, nrow = N, ncol = K-1)

  for(i in 1:N){
    if(length(which(tempZmat[i,]==0)) > 0){
      tempZmat[i,] <- (tempZmat[i,]*(M[i] - 1) + (1/K)) / M[i]
    }
    Ymat[i,] <- log(tempZmat[i,-base]/tempZmat[i,base])
  }

  return(Ymat)

}


####Appendix A##############################
##  Appendix A: Relation to Dunn-Smyth when J = 2  ##

#  Set Seed, Parameters for Simulation  #
set.seed(4)
ntrial <- 50
nsamp <- 100

#  Define Model, Simulate Data  #
x <- rnorm(nsamp)
logit <- -1.0 + .5*x
prob <- exp(logit)/(1+exp(logit))

y <- rbinom(nsamp,ntrial,prob)

#  Fit Data  #
fit1 <- glm(cbind(y,ntrial-y)~x,family=binomial)

#  Generate Various Residuals  #
ndist <- 10000
u <- numeric(length=nsamp)
ua <- numeric(length=nsamp)
u1 <- numeric(length=nsamp)
y1 <- y*(ntrial-1)/ntrial + 0.5
LLu <- numeric(length=nsamp)
ULu <- numeric(length=nsamp)
LLua <- numeric(length=nsamp)
ULua <- numeric(length=nsamp)
LLu1 <- numeric(length=nsamp)
ULu1 <- numeric(length=nsamp)
for(i in 1:nsamp){

  yf <- rbinom(ndist,ntrial,fit1$fitted.values[i])
  #  Adjust for zeros  #
  yf <- yf*(ntrial-1)/ntrial + 0.5
  wf <- log(yf/(ntrial-yf))
  md2 <- (wf - mean(wf))^2/var(wf)
  md2a <- (wf - mean(wf))/sqrt(var(wf))
  md2o <- (log(y1[i]/(ntrial-y1[i])) - mean(wf))^2/var(wf)
  md2oa <- (log(y1[i]/(ntrial-y1[i])) - mean(wf))/sqrt(var(wf))
  if(md2o <= min(md2)) {
    u[i] <- runif(1,0,length(md2[md2 <= min(md2)])/(ndist+1))
    LLu[i] <- 0
    ULu[i] <- length(md2[md2 <= min(md2)])/(ndist+1)
  } else if(md2o >= max(md2)) {
    u[i] <- runif(1,length(md2[md2 <= max(md2)])/(ndist+1),1)
    LLu[i] <- length(md2[md2 <= max(md2)])/(ndist+1)
    ULu[i] <- 1
  } else {
    ind <- min(which(md2 == max(md2[md2 < md2o])))
    u[i] <- runif(1,length(md2[md2 <= md2[ind]])/(ndist+1), length(md2[md2 <= md2o])/(ndist+1))
    LLu[i] <- length(md2[md2 <= md2[ind]])/(ndist+1)
    ULu[i] <- length(md2[md2 <= md2o])/(ndist+1)
  }

  if(md2oa <= min(md2a)) {
    ua[i] <- runif(1,0,length(md2a[md2a <= min(md2a)])/(ndist+1))
    LLua[i] <- 0
    ULua[i] <- length(md2a[md2a <= min(md2a)])/(ndist+1)
  } else if(md2oa >= max(md2a)) {
    ua[i] <- runif(1,length(md2a[md2a <= max(md2a)])/(ndist+1),1)
    LLua[i] <- length(md2a[md2a <= max(md2a)])/(ndist+1)
    ULua[i] <- 1
  } else {
    ind <- min(which(md2a == max(md2a[md2a < md2oa])))
    ua[i] <- runif(1,length(md2a[md2a <= md2a[ind]])/(ndist+1), length(md2a[md2a <= md2oa])/(ndist+1))
    LLua[i] <- length(md2a[md2a <= md2a[ind]])/(ndist+1)
    ULua[i] <- length(md2a[md2a <= md2oa])/(ndist+1)
  }

  #  Dunn and Smyth u value  #
  if(y[i] == 0) {
    u1[i] <- runif(1,0,pbinom(1,ntrial,fit1$fitted.values[i]))
    LLu1[i] <- 0
    ULu1[i] <- pbinom(1,ntrial,fit1$fitted.values[i])
  } else if(y[i] == ntrial) {
    u1[i] <- runif(1,pbinom(ntrial-1,ntrial,fit1$fitted.values[i]),1)
    LLu1[i] <- pbinom(ntrial-1,ntrial,fit1$fitted.values[i])
    ULu1[i] <- 1
  } else {
    u1[i] <- runif(1,pbinom(y[i]-1,ntrial,fit1$fitted.values[i]),pbinom(y[i],ntrial,fit1$fitted.values[i]))
    LLu1[i] <- pbinom(y[i]-1,ntrial,fit1$fitted.values[i])
    ULu1[i] <- pbinom(y[i],ntrial,fit1$fitted.values[i])
  }
}

#  MD2 Residuals  #
zestu <- qnorm(u)
zestuqq <- qqnorm(zestu, col = "lightblue", pch = 19, cex = 2)
abline(a=0,b=1)
text(y ~x, labels=rownames(as.data.frame(zestuqq)),data=zestuqq, cex=0.9, font=2)

#  DS Residuals  #
zestu1 <- qnorm(u1)
zestu1qq <- qqnorm(zestu1, col = "lightblue", pch = 19, cex = 2)
abline(a=0,b=1)
text(y ~x, labels=rownames(as.data.frame(zestu1qq)),data=zestu1qq, cex=0.9, font=2)

#  MD Residuals  #
zestua <- qnorm(ua)
zestuaqq <- qqnorm(zestua, col = "lightblue", pch = 19, cex = 2)
abline(a=0,b=1)
text(y ~x, labels=rownames(as.data.frame(zestuaqq)),data=zestuaqq, cex=0.9, font=2)

#  Agreement between LL + UL MD and DS  #
cor.test(LLu1, LLua, method = "spearman")
png("Figure1a.png", units = "in", width = 5, height = 5, res = 600)
plot(LLu1, LLua, xlab = "Dunn-Smyth Lower Limit", ylab = "MD Lower Limit", main = "Agreement in Lower Limit, DS vs. MD residuals")
abline(0, 1)
dev.off()

cor.test(ULu1, ULua, method = "spearman")
png("Figure1b.png", units = "in", width = 5, height = 5, res = 600)
plot(ULu1, ULua, xlab = "Dunn-Smyth Upper Limit", ylab = "MD Upper Limit", main = "Agreement in Upper Limit, DS vs. MD residuals")
abline(0, 1)
dev.off()

png("Figure1.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
plot(LLu1, LLua, xlab = "Dunn-Smyth Lower Limit", ylab = "MD Lower Limit", main = "Agreement in Lower Limit")
abline(0, 1)
plot(ULu1, ULua, xlab = "Dunn-Smyth Upper Limit", ylab = "MD Upper Limit", main = "Agreement in Upper Limit")
abline(0, 1)
dev.off()


####Section 3###############################
##  Section 3  ##

#  Set Seed, Parameters, Simulate Data  #
set.seed(3)
m <- 1000                  #observations
X <- sort(rnorm(m))        #covariate
Xmat <- cbind(1, X, X^2)   #design matrix
n <- 200                   #exposure
J <- 3                     #number of  categories

# M-Logit with quadratic term #
B <- cbind(c(-1, .5, .25), c(-.75, -1, -.5))
XB <- Xmat%*%B
Pi <- cbind(exp(XB)/(1 + rowSums(exp(XB))), 1/(1 + rowSums(exp(XB))))
Y <- t(apply(Pi, 1, rmultinom, n = 1, size = n))

# Create data frame #
ml <- as.data.frame(cbind(Y, Xmat))

# MLN (overdispersed) with no quadratic #
Xmat2 <- cbind(1, X)
B2 <- cbind(c(-1, .5), c(-.75, -1))
epsilon <- rmvn(n = m, mu = numeric(length = J-1), sigma = diag(J-1))
XB2 <- Xmat2%*%B2 + epsilon
Pi2 <- cbind(exp(XB2)/(1 + rowSums(exp(XB2))), 1/(1 + rowSums(exp(XB2))))
Y2 <- t(apply(Pi2, 1, rmultinom, n = 1, size = n))

# Create data frame #
ml2 <- as.data.frame(cbind(Y2, Xmat2))

# DM (overdispersed/negative correlations) with no quadratic #
get_alpha = dirmult(data = Y2)
master_alpha = get_alpha$gamma
B3 <- rbind(log(master_alpha), c(.25, -.15, .05))
XB3 <- Xmat2%*%B3
alpha <- exp(XB3)
Y3 <- t(apply(alpha, 1, rdirmn, n = 1, size = n))

# Create data frame #
ml3 <- as.data.frame(cbind(Y3, Xmat2))

# MLN (overdispersed) with no quadratic (and positive correlations) #
epsilon2 <- rmvn(n = m, mu = numeric(length = J-1), sigma = matrix(c(1,-.9,-.9,4), byrow=T, ncol=2))
XB4 <- Xmat2%*%B2 + epsilon2
Pi4 <- cbind(exp(XB4)/(1 + rowSums(exp(XB4))), 1/(1 + rowSums(exp(XB4))))
Y4 <- t(apply(Pi4, 1, rmultinom, n = 1, size = n))

# Create data frame #
ml4 <- as.data.frame(cbind(Y4, Xmat2))

# Observed Log-Odds #
ly_true <- newinitialY(Y, base = J)
ly2_true <- newinitialY(Y2, base = J)
ly3_true <- newinitialY(Y3, base = J)
ly4_true <- newinitialY(Y4, base = J)

#  Fitting the Model(s)  #
# MN = Multinomial Logistic Regression, with quadratic term #
correct_mlogit <- MGLMreg(cbind(V1, V2, V3) ~ X + V6, data = ml, dist = "MN")
# without quadratic term #
incorrect_mlogit <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml, dist = "MN")

# Multinomial Logistic Normal Regression #
correct_mln <- mixgibbs(W = Y2, X = Xmat2, Z = diag(1), base = J,  mhscale = 1, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)
# Multinomial Logistic Regression with MLN Data #
incorrect_mln <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml2, dist = "MN")

# Dirichlet-Multinomial Regression #
correct_dm <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml3, dist = "DM")
# MLN Regression with Dirichlet-Multinomial Data #
incorrect_dm <- mixgibbs(W = Y3, X = Xmat2, Z = diag(1), base = J,  mhscale = 1, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)

# Multinomial Logistic Normal Regression with Implied Positive Correlations #
correct_mln_pos <- mixgibbs(W = Y4, X = Xmat2, Z = diag(1), base = J, mhscale = 1, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)
# DM Regression with MLN Data (and implied positive correlations) #
incorrect_mln_dm <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml4, dist = "DM")

#  Constructing Sampling Distributions of Fitted Models  #
# Correct MLOGIT #
sampdists_cmlogit <- list(length = m)
for(i in 1:m){
  sampdists_cmlogit[[i]] <- t(rmultinom(n = 100, size = n, prob = correct_mlogit@fitted[i,]))
}

# Incorrect MLOGIT #
sampdists_imlogit <- list(length = m)
for(i in 1:m){
  sampdists_imlogit[[i]] <- t(rmultinom(n = 100, size = n, prob = incorrect_mlogit@fitted[i,]))
}

# Correct MLN #
burnthin <- seq(401, 1001, 2)

posteriormean.y <- matrix(colMeans(correct_mln$Y[burnthin,]), ncol = J-1, byrow = T)
posteriormean.beta <- matrix(colMeans(correct_mln$Beta[burnthin,]), ncol = J-1, byrow = T)
posteriormean.sigma <- matrix(colMeans(correct_mln$Sigma[burnthin,]), ncol = J-1, byrow = T)

Y_samples <- correct_mln$Y[burnthin,]
Beta_samples <- correct_mln$Beta[burnthin,]
Sigma_samples <- correct_mln$Sigma[burnthin,]

sampdists_cmln <- vector(mode = "list", length = m)
for(i in 1:m){
  sampdists_cmln[[i]] <- posteriorsampling.givenparams(Exposure = sum(Y2[i,]), X = t(as.matrix(Xmat2[i,])), dimW = J, sigma_mean = posteriormean.sigma, beta_mean = posteriormean.beta)
}

# Incorrect MLN #
sampdists_imln <- list(length = m)
for(i in 1:m){
  sampdists_imln[[i]] <- t(rmultinom(n = 10000, size = n, prob = incorrect_mln@fitted[i,]))
}

# Correct DM #
pred_alpha <- exp(Xmat2%*%correct_dm@coefficients)
sampdists_cdm <- list(length = m)
for(i in 1:m){
  sampdists_cdm[[i]] <- rdirmn(n = 10000, size = n, alpha = pred_alpha[i,])
}

# Incorrect DM #
posteriormean.y2 <- matrix(colMeans(incorrect_dm$Y[burnthin,]), ncol = J-1, byrow = T)
posteriormean.beta2 <- matrix(colMeans(incorrect_dm$Beta[burnthin,]), ncol = J-1, byrow = T)
posteriormean.sigma2 <- matrix(colMeans(incorrect_dm$Sigma[burnthin,]), ncol = J-1, byrow = T)

Y_samples2 <- incorrect_dm$Y[burnthin,]
Beta_samples2 <- incorrect_dm$Beta[burnthin,]
Sigma_samples2 <- incorrect_dm$Sigma[burnthin,]

sampdists_idm <- vector(mode = "list", length = m)
for(i in 1:m){
  sampdists_idm[[i]] <- posteriorsampling.givenparams(Exposure = sum(Y3[i,]), X = t(as.matrix(Xmat2[i,])), dimW = J, sigma_mean = posteriormean.sigma2, beta_mean = posteriormean.beta2)
}

# Incorrect MLN (DM) #
pred_alpha_mlndm <- exp(Xmat2%*%incorrect_mln_dm@coefficients)
sampdists_imlndm <- list(length = m)
for(i in 1:m){
  sampdists_imlndm[[i]] <- rdirmn(n = 10000, size = n, alpha = pred_alpha_mlndm[i,])
}


####Note 1##################################
##  Note  ##
# This method will not work if there is no covariance in the sampling distribution
# Luckily, this only happens if some observations have extreme (close to 0 or 1) fitted probabilities
# If this is the case for most of the observations, the method will not work
# However, using the parameter settings above (and most real data situations) it should be rare
# If there are any, we will simply remove these "singular" observations

#  Identifying and Removing Any "Singular" Observations (Should Be Very Few, If Any)  #
# Correct MLOGIT #
singulars_A <- unlist(lapply(sampdists_cmlogit, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_A) > 0){
  ly_true_A <- ly_true[-which(singulars_A == 1),]
  sampdists_A <- sampdists_cmlogit[-which(singulars_A == 1)]
} else{
  ly_true_A <- ly_true
  sampdists_A <- sampdists_cmlogit
}

# Incorrect MLOGIT #
singulars_B <- unlist(lapply(sampdists_imlogit, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_B) > 0){
  ly_true_B <- ly_true[-which(singulars_B == 1),]
  sampdists_B <- sampdists_imlogit[-which(singulars_B == 1)]
} else{
  ly_true_B <- ly_true
  sampdists_B <- sampdists_imlogit
}

# Correct MLN #
singulars_C <- unlist(lapply(sampdists_cmln, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_C) > 0){
  ly_true_C <- ly2_true[-which(singulars_C == 1),]
  sampdists_C <- sampdists_cmln1[-which(singulars_C == 1)]
} else{
  ly_true_C <- ly2_true
  sampdists_C <- sampdists_cmln
}

# Incorrect MLN #
singulars_D <- unlist(lapply(sampdists_imln, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_D) > 0){
  ly_true_D <- ly2_true[-which(singulars_D == 1),]
  sampdists_D <- sampdists_imln[-which(singulars_D == 1)]
} else{
  ly_true_D <- ly2_true
  sampdists_D <- sampdists_imln
}

# Correct DM #
singulars_E <- unlist(lapply(sampdists_cdm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_E) > 0){
  ly_true_E <- ly3_true[-which(singulars_E == 1),]
  sampdists_E <- sampdists_cdm[-which(singulars_E == 1)]
} else{
  ly_true_E <- ly3_true
  sampdists_E <- sampdists_cdm
}

# Incorrect DM #
singulars_F <- unlist(lapply(sampdists_idm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_F) > 0){
  ly_true_F <- ly3_true[-which(singulars_F == 1),]
  sampdists_F <- sampdists_idm[-which(singulars_F == 1)]
} else{
  ly_true_F <- ly3_true
  sampdists_F <- sampdists_idm
}

# Incorrect MLN (DM) #
singulars_G <- unlist(lapply(sampdists_imlndm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))

if(sum(singulars_G) > 0){
  ly_true_G <- ly2_true[-which(singulars_G == 1),]
  sampdists_G <- sampdists_imlndm[-which(singulars_G == 1)]
} else{
  ly_true_G <- ly2_true
  sampdists_G <- sampdists_imlndm
}

##  Calculating the Squared Mahalanobis Distances  ##
##  Note:  ##
# This may take some time depending on the size/dimensionality of the data
# Correct MLOGIT #
mds_A <- vector("list", length(which(singulars_A == 0)))
for(i in 1:length(which(singulars_A == 0))){
  yobs_A <- ly_true_A[i,]
  ypred_A <- colMeans(newinitialY(sampdists_A[[i]], J))
  sigmaypred_A <- cov(newinitialY(sampdists_A[[i]], J))

  obsi_A <- rbind(yobs_A, newinitialY(sampdists_A[[i]], J))

  mds_A[[i]] <- apply(obsi_A, 1, mahalanobis, center = ypred_A, cov = sigmaypred_A)
}

# Incorrect MLOGIT #
mds_B <- vector("list", length(which(singulars_B == 0)))
for(i in 1:length(which(singulars_B == 0))){
  yobs_B <- ly_true_B[i,]
  ypred_B <- colMeans(newinitialY(sampdists_B[[i]], J))
  sigmaypred_B <- cov(newinitialY(sampdists_B[[i]], J))

  obsi_B <- rbind(yobs_B, newinitialY(sampdists_B[[i]], J))

  mds_B[[i]] <- apply(obsi_B, 1, mahalanobis, center = ypred_B, cov = sigmaypred_B)
}

# Correct MLN #
mds_C <- vector("list", length(which(singulars_C == 0)))
for(i in 1:length(which(singulars_C == 0))){
  yobs_C <- ly_true_C[i,]
  ypred_C <- colMeans(newinitialY(sampdists_C[[i]], J))
  sigmaypred_C <- cov(newinitialY(sampdists_C[[i]], J))

  obsi_C <- rbind(yobs_C, newinitialY(sampdists_C[[i]], J))

  mds_C[[i]] <- apply(obsi_C, 1, mahalanobis, center = ypred_C, cov = sigmaypred_C)
}

# Incorrect MLN #
mds_D <- vector("list", length(which(singulars_D == 0)))
for(i in 1:length(which(singulars_D == 0))){
  yobs_D <- ly_true_D[i,]
  ypred_D <- colMeans(newinitialY(sampdists_D[[i]], J))
  sigmaypred_D <- cov(newinitialY(sampdists_D[[i]], J))

  obsi_D <- rbind(yobs_D, newinitialY(sampdists_D[[i]], J))

  mds_D[[i]] <- apply(obsi_D, 1, mahalanobis, center = ypred_D, cov = sigmaypred_D)
}

# Correct DM #
mds_E <- vector("list", length(which(singulars_E == 0)))
for(i in 1:length(which(singulars_E == 0))){
  yobs_E <- ly_true_E[i,]
  ypred_E <- colMeans(newinitialY(sampdists_E[[i]], J))
  sigmaypred_E <- cov(newinitialY(sampdists_E[[i]], J))

  obsi_E <- rbind(yobs_E, newinitialY(sampdists_E[[i]], J))

  mds_E[[i]] <- apply(obsi_E, 1, mahalanobis, center = ypred_E, cov = sigmaypred_E)
}

# Incorrect DM #
mds_F <- vector("list", length(which(singulars_F == 0)))
for(i in 1:length(which(singulars_F == 0))){
  yobs_F <- ly_true_F[i,]
  ypred_F <- colMeans(newinitialY(sampdists_F[[i]], J))
  sigmaypred_F <- cov(newinitialY(sampdists_F[[i]], J))

  obsi_F <- rbind(yobs_F, newinitialY(sampdists_F[[i]], J))

  mds_F[[i]] <- apply(obsi_F, 1, mahalanobis, center = ypred_F, cov = sigmaypred_F)
}

# Incorrect MLN (DM) #
mds_G <- vector("list", length(which(singulars_G == 0)))
for(i in 1:length(which(singulars_G == 0))){
  yobs_G <- ly_true_G[i,]
  ypred_G <- colMeans(newinitialY(sampdists_G[[i]], J))
  sigmaypred_G <- cov(newinitialY(sampdists_G[[i]], J))

  obsi_G <- rbind(yobs_G, newinitialY(sampdists_G[[i]], J))

  mds_G[[i]] <- apply(obsi_G, 1, mahalanobis, center = ypred_G, cov = sigmaypred_G)
}

##  Uniform Residuals (These Will Be Back-Transformed to Standard Normal for the Final Form)  ##
# Correct MLOGIT #
mdps_A <- double(length = length(which(singulars_A == 0)))
for(i in 1:length(which(singulars_A == 0))){
  pctl <- ecdf(mds_A[[i]][-1])

  if(mds_A[[i]][1] <= min(mds_A[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_A[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_A[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_A[[i]][1] > max(mds_A[[i]][-1])){
      minpct <- pctl(max(mds_A[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_A[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_A[[i]][-1])[max(which(sort(mds_A[[i]][-1]) < mds_A[[i]][1]))])
      maxpct <- pctl(mds_A[[i]][1])
      if(minpct == 1){
        L <- length(mds_A[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_A[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_A[i] <- runif(1, minpct, maxpct)
}

# Incorrect MLOGIT #
mdps_B <- double(length = length(which(singulars_B == 0)))
for(i in 1:length(which(singulars_B == 0))){
  pctl <- ecdf(mds_B[[i]][-1])

  if(mds_B[[i]][1] <= min(mds_B[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_B[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_B[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_B[[i]][1] > max(mds_B[[i]][-1])){
      minpct <- pctl(max(mds_B[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_B[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_B[[i]][-1])[max(which(sort(mds_B[[i]][-1]) < mds_B[[i]][1]))])
      maxpct <- pctl(mds_B[[i]][1])
      if(minpct == 1){
        L <- length(mds_B[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_B[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_B[i] <- runif(1, minpct, maxpct)
}

# Correct MLN #
mdps_C <- double(length = length(which(singulars_C == 0)))
for(i in 1:length(which(singulars_C == 0))){
  pctl <- ecdf(mds_C[[i]][-1])

  if(mds_C[[i]][1] <= min(mds_C[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_C[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_C[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_C[[i]][1] > max(mds_C[[i]][-1])){
      minpct <- pctl(max(mds_C[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_C[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_C[[i]][-1])[max(which(sort(mds_C[[i]][-1]) < mds_C[[i]][1]))])
      maxpct <- pctl(mds_C[[i]][1])
      if(minpct == 1){
        L <- length(mds_C[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_C[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_C[i] <- runif(1, minpct, maxpct)
}

# Incorrect MLN #
mdps_D <- double(length = length(which(singulars_D == 0)))
for(i in 1:length(which(singulars_D == 0))){
  pctl <- ecdf(mds_D[[i]][-1])

  if(mds_D[[i]][1] <= min(mds_D[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_D[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_D[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_D[[i]][1] > max(mds_D[[i]][-1])){
      minpct <- pctl(max(mds_D[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_D[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_D[[i]][-1])[max(which(sort(mds_D[[i]][-1]) < mds_D[[i]][1]))])
      maxpct <- pctl(mds_D[[i]][1])
      if(minpct == 1){
        L <- length(mds_D[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_D[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_D[i] <- runif(1, minpct, maxpct)
}

# Correct DM #
mdps_E <- double(length = length(which(singulars_E == 0)))
for(i in 1:length(which(singulars_E == 0))){
  pctl <- ecdf(mds_E[[i]][-1])

  if(mds_E[[i]][1] <= min(mds_E[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_E[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_E[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_E[[i]][1] > max(mds_E[[i]][-1])){
      minpct <- pctl(max(mds_E[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_E[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_E[[i]][-1])[max(which(sort(mds_E[[i]][-1]) < mds_E[[i]][1]))])
      maxpct <- pctl(mds_E[[i]][1])
      if(minpct == 1){
        L <- length(mds_E[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_E[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_E[i] <- runif(1, minpct, maxpct)
}

# Incorrect DM #
mdps_F <- double(length = length(which(singulars_F == 0)))
for(i in 1:length(which(singulars_F == 0))){
  pctl <- ecdf(mds_F[[i]][-1])

  if(mds_F[[i]][1] <= min(mds_F[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_F[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_F[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_F[[i]][1] > max(mds_F[[i]][-1])){
      minpct <- pctl(max(mds_F[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_F[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_F[[i]][-1])[max(which(sort(mds_F[[i]][-1]) < mds_F[[i]][1]))])
      maxpct <- pctl(mds_F[[i]][1])
      if(minpct == 1){
        L <- length(mds_F[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_F[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_F[i] <- runif(1, minpct, maxpct)
}

# Incorrect MLN (DM) #
mdps_G <- double(length = length(which(singulars_G == 0)))
for(i in 1:length(which(singulars_G == 0))){
  pctl <- ecdf(mds_G[[i]][-1])

  if(mds_G[[i]][1] <= min(mds_G[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds_G[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds_G[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds_G[[i]][1] > max(mds_G[[i]][-1])){
      minpct <- pctl(max(mds_G[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds_G[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds_G[[i]][-1])[max(which(sort(mds_G[[i]][-1]) < mds_G[[i]][1]))])
      maxpct <- pctl(mds_G[[i]][1])
      if(minpct == 1){
        L <- length(mds_G[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds_G[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps_G[i] <- runif(1, minpct, maxpct)
}

##  Converting to Standard Normal (Final Form) and Creating Plots  ##
##  Note:  ##
# Some additional plots, color coding observations for use in interpreting lack of fit, are included
# Correct MLOGIT #
mdrs_A <- qnorm(mdps_A)
qqnorm(mdrs_A)
abline(0, 1, col = "red", lwd = 2)
png("Figure2a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_A, col.lines = "red", line = "none", main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure2b.png", units = "in", width = 5, height = 5, res = 600)
plot(X[!singulars_A], mdrs_A, main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure2.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_A, col.lines = "red", line = "none", main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(X[!singulars_A], mdrs_A, main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Incorrect MLOGIT #
mdrs_B <- qnorm(mdps_B)
qqnorm(mdrs_B)
abline(0, 1, col = "red", lwd = 2)
png("Figure3a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_B, col.lines = "red", line = "none", main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure3b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_B, main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure3.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_B, col.lines = "red", line = "none", main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_B, main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Correct MLN #
mdrs_C <- qnorm(mdps_C)
qqnorm(mdrs_C)
abline(0, 1, col = "red", lwd = 2)
png("Figure5a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_C, col.lines = "red", line = "none", main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure5b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_C, main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure5.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_C, col.lines = "red", line = "none", main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_C, main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Incorrect MLN #
mdrs_D <- qnorm(mdps_D)
qqnorm(mdrs_D)
abline(0, 1, col = "red", lwd = 2)
png("Figure4a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_D, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure4b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_D, main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure4.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_D, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_D, main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Correct DM #
mdrs_E <- qnorm(mdps_E)
qqnorm(mdrs_E)
abline(0, 1, col = "red", lwd = 2)
png("Figure6a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_E, col.lines = "red", line = "none", main = "Correctly Modeled DM", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = 'lines'))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure6b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_E, main = "Correctly Modeled DM", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure6.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_E, col.lines = "red", line = "none", main = "Correctly Modeled DM", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = 'lines'))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_E, main = "Correctly Modeled DM", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Incorrect DM #
mdrs_F <- qnorm(mdps_F)
qqnorm(mdrs_F)
abline(0, 1, col = "red", lwd = 2)
png("Figure8a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_F, col.lines = "red", line = "none", main = "Incorrectly Modeled DM", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure8b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_F, main = "Incorrectly Modeled DM", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure8.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_F, col.lines = "red", line = "none", main = "Incorrectly Modeled DM", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_F, main = "Incorrectly Modeled DM", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

# Incorrect MLN (DM) #
mdrs_G <- qnorm(mdps_G)
qqnorm(mdrs_G)
abline(0, 1, col = "red", lwd = 2)
png("Figure7a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_G, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = 'lines'))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure7b.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_G, main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

png("Figure7.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs_G, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = 'lines'))
abline(0, 1, col = "red", lwd = 2)
plot(X, mdrs_G, main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Covariate")
abline(h = 0, col = "red")
dev.off()

#  Investigating Potential Ways of Identifying Lack of Fit  #
colorcodeX <- ifelse(abs(X) > 1.96, "red", ifelse(abs(X) > 1.645, "orange", "blue"))
error_likes <- dmvn(X = epsilon, mu = c(0,0), sigma = diag(2))
colorcodeE <- ifelse(error_likes < .025, "red", ifelse(error_likes < .05, "orange", "blue"))

car::qqPlot(mdrs_A, line = "none", envelope = list(style = "lines"), main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeX[!singulars_A])
plot(X[!singulars_A], mdrs_A, main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeX[!singulars_A])

car::qqPlot(mdrs_B, line = "none", envelope = list(style = "lines"), main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeX)
plot(X, mdrs_B, main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeX)

#png("colorcode_CMLN.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs_C, line = "none", envelope = list(style = "lines"), main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeE)
#dev.off()
#png("colorcode_CMLN2.png", units = "in", width = 5, height = 5, res = 600)
plot(X, mdrs_C, main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeE)
dev.off()

car::qqPlot(mdrs_D, line = "none", envelope = list(style = "lines"), main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeE)
plot(X, mdrs_D, main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeE)

car::qqPlot(mdrs_E, line = "none", envelope = list(style = "lines"), main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeE)
plot(X, mdrs_D, main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeE)

car::qqPlot(mdrs_G, line = "none", envelope = list(style = "lines"), main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", col = colorcodeE)
plot(X, mdrs_G, main = "Incorrectly Modeled MLN (w/ DM)", ylab = "MD2 Residuals", xlab = "Covariate", col = colorcodeE)


####Note 2##################################
##  Note  ##
# Figures 9 and 10 (plots for models fit to data with varying exposure) can be created by:
# Replacing the exposure (n) parameter as one of the below and re-running the model fits (A,B,C,D)
# Adjusting plot names accordingly

# # Varying Exposure #
# n <- round(runif(n = 1000, min = 2, max = 400))   #exposure
# n <- round(runif(n = 1000, min = 2, max = 30))    #exposure
#
# # Correct MLOGIT #
# png("Figure9a.png", units = "in", width = 5, height = 5, res = 600)
# car::qqPlot(mdrs_A, col.lines = "red", line = "none", main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# dev.off()
#
# # Incorrect MLOGIT #
# png("Figure9b.png", units = "in", width = 5, height = 5, res = 600)
# car::qqPlot(mdrs_B, col.lines = "red", line = "none", main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# dev.off()
#
# # Correct MLN #
# png("Figure9c.png", units = "in", width = 5, height = 5, res = 600)
# car::qqPlot(mdrs_C, col.lines = "red", line = "none", main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# dev.off()
#
# # Incorrect MLN #
# png("Figure9d.png", units = "in", width = 5, height = 5, res = 600)
# car::qqPlot(mdrs_D, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# dev.off()
#
# png("Figure9.png", units = "in", width = 15, height = 15, res = 600)
# par(mfrow=c(2,2), cex=1.5)
# car::qqPlot(mdrs_A, col.lines = "red", line = "none", main = "Correctly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# car::qqPlot(mdrs_B, col.lines = "red", line = "none", main = "Incorrectly Modeled MLOGIT", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# car::qqPlot(mdrs_C, col.lines = "red", line = "none", main = "Correctly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# car::qqPlot(mdrs_D, col.lines = "red", line = "none", main = "Incorrectly Modeled MLN", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = list(style = "lines"))
# abline(0, 1, col = "red", lwd = 2)
# dev.off()

####Kolmogorov-Smirnov Tests################
##  Kolmogorov-Smirnov Test Simulation Studies  ##
##  Note  ##
# The Kolmogorov-Smirnov Tests are done in parallel (due to heavy computation time)
# They are provided in a separate (messier) .R file, ks_sim_study.R
# The tables from the manuscript may be replicated using that file

####Real Data Analysis######################
##  Real Data Analysis  ##
##  Note  ##
# Some of the objects within each real data application were named similarly/the same as above
# This was due to laziness, but should not cause a problem as long as the code is run in order
####Pollen Data#############################
##  Pollen Data  ##
data(pollen)
head(pollen)

# Explore a bit #
pollen <- as.data.frame(pollen)
cor(pollen)
colMeans(pollen)

# Number of categories #
J <- ncol(pollen)

# Fit three candidate models #
pollen_mlr <- MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1, data = pollen, dist = "MN")

pollen_dm <- MGLMreg(cbind(Pinus, Abies, Quercus, Alnus) ~ 1, data = pollen, dist = "DM")

pollen_mln <- mixgibbs(W = pollen, X = matrix(1, nrow = nrow(pollen), ncol = 1), Z = diag(1), base = J, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)

# Get Posterior MLN Chains #
burnthin <- seq(401, 1001, 2)

posteriormean.y <- matrix(colMeans(pollen_mln$Y[burnthin,]), ncol = J-1, byrow = T)
posteriormean.beta <- matrix(colMeans(pollen_mln$Beta[burnthin,]), ncol = J-1, byrow = T)
posteriormean.sigma <- matrix(colMeans(pollen_mln$Sigma[burnthin,]), ncol = J-1, byrow = T)

Y_samples <- pollen_mln$Y[burnthin,]
Beta_samples <- pollen_mln$Beta[burnthin,]
Sigma_samples <- pollen_mln$Sigma[burnthin,]

# Sampling/Posterior Distributions of Fitted Models for Each Individual #
sd_mlr <- list(length = nrow(pollen))
for(i in 1:nrow(pollen)){
  sd_mlr[[i]] <- t(rmultinom(n = 10000, size = rowSums(pollen)[i], prob = pollen_mlr@fitted[i,]))
}

pred_alpha <- exp(pollen_dm@coefficients)
sd_dm <- list(length = nrow(pollen))
for(i in 1:nrow(pollen)){
  sd_dm[[i]] <- rdirmn(n = 10000, size = rowSums(pollen)[i], alpha = c(pred_alpha))
}

sd_mln <- vector(mode = "list", length = nrow(pollen))
for(i in 1:nrow(pollen)){
  sd_mln[[i]] <- posteriorsampling.givenparams(Exposure = sum(pollen[i,]), X = t(as.matrix(1)), dimW = J, sigma_mean = posteriormean.sigma, beta_mean = posteriormean.beta)
}

# Observed Log-Odds #
ly_true <- newinitialY(pollen, base = J)
head(ly_true)

# Check for Sampling Distributions without Variability #
singularsA <- unlist(lapply(sd_mlr, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsA) > 0){
  ly_true2A <- ly_true[-which(singularsA == 1),]
  sd_mlr2 <- sd_mlr[-which(singularsA == 1)]
} else{
  ly_true2A <- ly_true
  sd_mlr2 <- sd_mlr
}

singularsB <- unlist(lapply(sd_dm, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsB) > 0){
  ly_true2B <- ly_true[-which(singularsB == 1),]
  sd_dm2 <- sd_dm[-which(singularsB == 1)]
} else{
  ly_true2B <- ly_true
  sd_dm2 <- sd_dm
}

singularsC <- unlist(lapply(sd_mln, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsC) > 0){
  ly_true2C <- ly_true[-which(singularsC == 1),]
  sd_mln2 <- sd_mln[-which(singularsC == 1)]
} else{
  ly_true2C <- ly_true
  sd_mln2 <- sd_mln
}

# Calculate Squared Mahalanobis Distances #
mds1A <- vector("list", length(which(singularsA == 0)))
for(i in 1:length(which(singularsA == 0))){
  yobsA <- ly_true2A[i,]
  ypredA <- colMeans(newinitialY(sd_mlr2[[i]], J))
  sigmaypredA <- cov(newinitialY(sd_mlr2[[i]], J))

  playeriA <- rbind(yobsA, newinitialY(sd_mlr2[[i]], J))

  mds1A[[i]] <- apply(playeriA, 1, mahalanobis, center = ypredA, cov = sigmaypredA)
}

mdps1A <- double(length = length(which(singularsA == 0)))
for(i in 1:length(which(singularsA == 0))){
  pctl <- ecdf(mds1A[[i]][-1])

  if(mds1A[[i]][1] <= min(mds1A[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1A[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1A[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1A[[i]][1] > max(mds1A[[i]][-1])){
      minpct <- pctl(max(mds1A[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1A[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1A[[i]][-1])[max(which(sort(mds1A[[i]][-1]) < mds1A[[i]][1]))])
      maxpct <- pctl(mds1A[[i]][1])
      if(minpct == 1){
        L <- length(mds1A[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1A[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1A[i] <- runif(1, minpct, maxpct)
}

mds1B <- vector("list", length(which(singularsB == 0)))
for(i in 1:length(which(singularsB == 0))){
  yobsB <- ly_true2B[i,]
  ypredB <- colMeans(newinitialY(sd_dm2[[i]], J))
  sigmaypredB <- cov(newinitialY(sd_dm2[[i]], J))

  playeriB <- rbind(yobsB, newinitialY(sd_dm2[[i]], J))

  mds1B[[i]] <- apply(playeriB, 1, mahalanobis, center = ypredB, cov = sigmaypredB)
}

mdps1B <- double(length = length(which(singularsB == 0)))
for(i in 1:length(which(singularsB == 0))){
  pctl <- ecdf(mds1B[[i]][-1])

  if(mds1B[[i]][1] <= min(mds1B[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1B[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1B[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1B[[i]][1] > max(mds1B[[i]][-1])){
      minpct <- pctl(max(mds1B[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1B[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1B[[i]][-1])[max(which(sort(mds1B[[i]][-1]) < mds1B[[i]][1]))])
      maxpct <- pctl(mds1B[[i]][1])
      if(minpct == 1){
        L <- length(mds1B[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1B[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1B[i] <- runif(1, minpct, maxpct)
}

mds1C <- vector("list", length(which(singularsC == 0)))
for(i in 1:length(which(singularsC == 0))){
  yobsC <- ly_true2C[i,]
  ypredC <- colMeans(newinitialY(sd_mln2[[i]], J))
  sigmaypredC <- cov(newinitialY(sd_mln2[[i]], J))

  playeriC <- rbind(yobsC, newinitialY(sd_mln2[[i]], J))

  mds1C[[i]] <- apply(playeriC, 1, mahalanobis, center = ypredC, cov = sigmaypredC)
}

mdps1C <- double(length = length(which(singularsC == 0)))
for(i in 1:length(which(singularsC == 0))){
  pctl <- ecdf(mds1C[[i]][-1])

  if(mds1C[[i]][1] <= min(mds1C[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1C[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1C[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1C[[i]][1] > max(mds1C[[i]][-1])){
      minpct <- pctl(max(mds1C[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1C[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1C[[i]][-1])[max(which(sort(mds1C[[i]][-1]) < mds1C[[i]][1]))])
      maxpct <- pctl(mds1C[[i]][1])
      if(minpct == 1){
        L <- length(mds1C[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1C[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1C[i] <- runif(1, minpct, maxpct)
}

# Convert to Standard Normal and Make Plots #
# MLOGIT Model #
mdrs1A <- qnorm(mdps1A)

png("Figure11a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1A, col.lines = "red", line = "none", main = "MLOGIT Pollen Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure11b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1A, main = "MLOGIT Pollen Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure11.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1A, col.lines = "red", line = "none", main = "MLOGIT Pollen Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1A, main = "MLOGIT Pollen Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# DM Model #
mdrs1B <- qnorm(mdps1B)

png("Figure12a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1B, col.lines = "red", line = "none", main = "DM Pollen Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure12b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1B, main = "DM Pollen Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure12.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1B, col.lines = "red", line = "none", main = "DM Pollen Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1B, main = "DM Pollen Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# MLN Model #
mdrs1C <- qnorm(mdps1C)

png("pollen_mln_QQ.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1C, col.lines = "red", line = "none", main = "MLN Pollen Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("pollen_mln_RES.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1C, main = "MLN Pollen Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# Kolmogorov-Smirnov Tests #
ks.test(mdrs1A, "pnorm", mean = 0, sd = 1)
ks.test(mdrs1B, "pnorm", mean = 0, sd = 1)
ks.test(mdrs1C, "pnorm", mean = 0, sd = 1)

####Baseball Data###########################
##  Baseball Data  ##
# Get a subset of the data #
modernBatting <- Batting[Batting$yearID > 1960,]
modernBattingB <- modernBatting[modernBatting$AB > 199,]
set.seed(42)
modernBattingC <- modernBattingB[order(sample(1:nrow(modernBattingB), 1000, FALSE)),]
rownames(modernBattingC) <- NULL
truethree <- cbind(modernBattingC[,c(12,16:17)], OTHER = modernBattingC$AB - rowSums(modernBattingC[,c(12,16:17)]))
head(truethree)

J <- 4

# Fit the candidate models #
ballmlrA <- MGLMreg(cbind(HR, BB, SO, OTHER) ~ 1, data = truethree, dist = "MN")

ballmlrB <- mixgibbs(W = truethree, X = matrix(1, nrow = nrow(truethree), ncol = 1), Z = diag(1), base = J, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)

ballmlrC <- MGLMreg(cbind(HR, BB, SO, OTHER) ~ 1, data = truethree, dist = "DM")

# MLN Posterior Chains #
burnthin <- seq(401, 1001, 2)

posteriormean.y <- matrix(colMeans(ballmlrB$Y[burnthin,]), ncol = J-1, byrow = T)
posteriormean.beta <- matrix(colMeans(ballmlrB$Beta[burnthin,]), ncol = J-1, byrow = T)
posteriormean.sigma <- matrix(colMeans(ballmlrB$Sigma[burnthin,]), ncol = J-1, byrow = T)

Y_samples <- ballmlrB$Y[burnthin,]
Beta_samples <- ballmlrB$Beta[burnthin,]
Sigma_samples <- ballmlrB$Sigma[burnthin,]

# Construct the sampling distribution of the fitted model for each individual #
sampdistsA <- list(length = nrow(truethree))
for(i in 1:nrow(truethree)){
  sampdistsA[[i]] <- t(rmultinom(n = 10000, size = rowSums(truethree)[i], prob = ballmlrA@fitted[i,]))
}

postsamples <- vector(mode = "list", length = nrow(truethree))
for(i in 1:nrow(truethree)){
  postsamples[[i]] <- posteriorsampling.givenparams(Exposure = sum(truethree[i,]), X = t(as.matrix(1)), dimW = J, sigma_mean = posteriormean.sigma, beta_mean = posteriormean.beta)
}

pred_alpha <- exp(ballmlrC@coefficients)
sampdistsC <- list(length = nrow(truethree))
for(i in 1:nrow(truethree)){
  sampdistsC[[i]] <- rdirmn(n = 10000, size = rowSums(truethree)[i], alpha = c(pred_alpha))
}

# Observed Log-Odds #
ly_true <- newinitialY(truethree, base = 4)

# Check for sampling distributions without variability #
singularsA <- unlist(lapply(sampdistsA, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsA) > 0){
  ly_true2A <- ly_true[-which(singularsA == 1),]
  sampdists2A <- sampdistsA[-which(singularsA == 1)]
} else{
  ly_true2A <- ly_true
  sampdists2A <- sampdistsA
}

singularsB <- unlist(lapply(postsamples, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsB) > 0){
  ly_true2B <- ly_true[-which(singularsB == 1),]
  postsamples2 <- postsamples[-which(singularsB == 1)]
} else{
  ly_true2B <- ly_true
  postsamples2 <- postsamples
}

singularsC <- unlist(lapply(sampdistsC, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singularsC) > 0){
  ly_true2C <- ly_true[-which(singularsC == 1),]
  sampdists2C <- sampdistsC[-which(singularsC == 1)]
} else{
  ly_true2C <- ly_true
  sampdists2C <- sampdistsC
}

# Calculate Mahalanobis distances #
mds1A <- vector("list", length(which(singularsA == 0)))
for(i in 1:length(which(singularsA == 0))){
  yobsA <- ly_true2A[i,]
  ypredA <- colMeans(newinitialY(sampdists2A[[i]], J))
  sigmaypredA <- cov(newinitialY(sampdists2A[[i]], J))

  playeriA <- rbind(yobsA, newinitialY(sampdists2A[[i]], J))

  mds1A[[i]] <- apply(playeriA, 1, mahalanobis, center = ypredA, cov = sigmaypredA)
}

mdps1A <- double(length = length(which(singularsA == 0)))
for(i in 1:length(which(singularsA == 0))){
  pctl <- ecdf(mds1A[[i]][-1])

  if(mds1A[[i]][1] <= min(mds1A[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1A[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1A[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1A[[i]][1] > max(mds1A[[i]][-1])){
      minpct <- pctl(max(mds1A[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1A[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1A[[i]][-1])[max(which(sort(mds1A[[i]][-1]) < mds1A[[i]][1]))])
      maxpct <- pctl(mds1A[[i]][1])
      if(minpct == 1){
        L <- length(mds1A[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1A[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1A[i] <- runif(1, minpct, maxpct)
}

mds1B <- vector("list", length(which(singularsB == 0)))
for(i in 1:length(which(singularsB == 0))){
  yobsB <- ly_true2B[i,]
  ypredB <- colMeans(newinitialY(postsamples2[[i]], J))
  sigmaypredB <- cov(newinitialY(postsamples2[[i]], J))

  playeriB <- rbind(yobsB, newinitialY(postsamples2[[i]], J))

  mds1B[[i]] <- apply(playeriB, 1, mahalanobis, center = ypredB, cov = sigmaypredB)
}

mdps1B <- double(length = length(which(singularsB == 0)))
for(i in 1:length(which(singularsB == 0))){
  pctl <- ecdf(mds1B[[i]][-1])

  if(mds1B[[i]][1] <= min(mds1B[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1B[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1B[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1B[[i]][1] > max(mds1B[[i]][-1])){
      minpct <- pctl(max(mds1B[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1B[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1B[[i]][-1])[max(which(sort(mds1B[[i]][-1]) < mds1B[[i]][1]))])
      maxpct <- pctl(mds1B[[i]][1])
      if(minpct == 1){
        L <- length(mds1B[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1B[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1B[i] <- runif(1, minpct, maxpct)
}

mds1C <- vector("list", length(which(singularsC == 0)))
for(i in 1:length(which(singularsC == 0))){
  yobsC <- ly_true2C[i,]
  ypredC <- colMeans(newinitialY(sampdists2C[[i]], J))
  sigmaypredC <- cov(newinitialY(sampdists2C[[i]], J))

  playeriC <- rbind(yobsC, newinitialY(sampdists2C[[i]], J))

  mds1C[[i]] <- apply(playeriC, 1, mahalanobis, center = ypredC, cov = sigmaypredC)
}

mdps1C <- double(length = length(which(singularsC == 0)))
for(i in 1:length(which(singularsC == 0))){
  pctl <- ecdf(mds1C[[i]][-1])

  if(mds1C[[i]][1] <= min(mds1C[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1C[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1C[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1C[[i]][1] > max(mds1C[[i]][-1])){
      minpct <- pctl(max(mds1C[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1C[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1C[[i]][-1])[max(which(sort(mds1C[[i]][-1]) < mds1C[[i]][1]))])
      maxpct <- pctl(mds1C[[i]][1])
      if(minpct == 1){
        L <- length(mds1C[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1C[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1C[i] <- runif(1, minpct, maxpct)
}

# Convert to Standard Normal and Make Plots #
# MLOGIT Model #
mdrs1A <- qnorm(mdps1A)

png("Figure13a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1A, col.lines = "red", line = "none", main = "MLOGIT Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure13b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1A, main = "MLOGIT Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure13.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1A, col.lines = "red", line = "none", main = "MLOGIT Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1A, main = "MLOGIT Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# MLN Model #
mdrs1B <- qnorm(mdps1B)

png("Figure15a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1B, col.lines = "red", line = "none", main = "MLN Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure15b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1B, main = "MLN Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure15.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1B, col.lines = "red", line = "none", main = "MLN Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1B, main = "MLN Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# DM Model #
mdrs1C <- qnorm(mdps1C)

png("Figure14a.png", units = "in", width = 5, height = 5, res = 300)
car::qqPlot(mdrs1C, col.lines = "red", line = "none", main = "DM Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure14b.png", units = "in", width = 5, height = 5, res = 300)
plot(mdrs1C, main = "DM Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure14.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1C, col.lines = "red", line = "none", main = "DM Baseball Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1C, main = "DM Baseball Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# Kolmogorov-Smirnov Tests #
ks.test(mdrs1A, "pnorm", mean = 0, sd = 1)
ks.test(mdrs1B, "pnorm", mean = 0, sd = 1)
ks.test(mdrs1C, "pnorm", mean = 0, sd = 1)

####Bird Count Data#########################
##  Bird Count Data  ##
bbs_63014_saxapahaw <- read_excel("bbs_63014_saxapahaw.xlsx")

# Decide which birds to use #
set.seed(33)
choosebirds <- order(colSums(bbs_63014_saxapahaw[,2:83]), decreasing = T)[1:J]
woodpeckers <- c(27,44,59,63,65)
sparrows <- c(35,37,46,72,23)
fewbirds <- bbs_63014_saxapahaw[,c(sparrows, 84:89)]
J <- 5

head(fewbirds)
colMeans(fewbirds)
birdnames <- names(fewbirds)
names(fewbirds)[1:J] <- c(paste0("B", 1:J))

X <- as.matrix(fewbirds[,c(6:11)])
Xs <- scale(X)
new_fewbirds <- cbind(fewbirds, Xs)
names(new_fewbirds)[12:17] <- c("X1", "X2", "X3", "X4", "X5", "X6")

####Note 3##################################
##  Note  ##
# The model selection is omitted here, but was done using the MGLMtune function with group penalty
# For example, among many other levels of interaction, we considered:
# MGLMtune(cbind(B1, B2, B3, B4, B5) ~ 1 + (X1 + X2 + X3 + X5 + X6)^2, data = new_fewbirds, dist = "MN", penalty = "group", ngridpt = 30, display = TRUE)
# which performs regularized estimation to determine which variables (with a group penalty) minimize BIC

# Fit final candidate models #
# Intercept
birdmlr <- MGLMreg(cbind(B1, B2, B3, B4, B5) ~ 1, data = new_fewbirds, dist = "MN")

# Includes MaxTempF, RainDays, MaxSusWindmph, and interaction between MinTempF and RainDays
birdmlr_full <- MGLMreg(cbind(B1, B2, B3, B4, B5) ~ 1 + X1 + X5 + X6 + X2:X5, data = new_fewbirds, dist = "MN")

# Construct the sampling distribution of the fitted model for each individual #
sampdists <- list(length = 22)
for(i in 1:22){
  sampdists[[i]] <- t(rmultinom(n = 10000, size = rowSums(new_fewbirds[,1:J])[i], prob = birdmlr@fitted[i,]))
}

sampdists_full <- list(length = 22)   #make the list to populate
for(i in 1:22){
  sampdists_full[[i]] <- t(rmultinom(n = 10000, size = rowSums(new_fewbirds[,1:J])[i], prob = birdmlr_full@fitted[i,]))
}

# Observed Log-Odds #
ly_true <- newinitialY(new_fewbirds[,1:J], base = J)
head(ly_true)

# Check for sampling distributions without variability #
singulars <- unlist(lapply(sampdists, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singulars) > 0){
  ly_true2 <- ly_true[-which(singulars == 1),]
  sampdists2 <- sampdists[-which(singulars == 1)]
} else{
  ly_true2 <- ly_true
  sampdists2 <- sampdists
}

singulars_full <- unlist(lapply(sampdists_full, function(x) ifelse(any(colMeans(x) == 0), 1, 0)))

if(sum(singulars_full) > 0){
  ly_true2_full <- ly_true[-which(singulars_full == 1),]
  sampdists2_full <- sampdists_full[-which(singulars_full == 1)]
} else{
  ly_true2_full <- ly_true
  sampdists2_full <- sampdists_full
}

# Calculate Mahalanobis distances #
mds1 <- vector("list", length(which(singulars == 0)))
for(i in 1:length(which(singulars == 0))){
  yobs <- ly_true2[i,]
  ypred <- colMeans(newinitialY(sampdists2[[i]], J))
  sigmaypred <- cov(newinitialY(sampdists2[[i]], J))

  playeri <- rbind(yobs, newinitialY(sampdists2[[i]], J))

  mds1[[i]] <- apply(playeri, 1, mahalanobis, center = ypred, cov = sigmaypred)
}

mdps1 <- double(length = length(which(singulars == 0)))
for(i in 1:length(which(singulars == 0))){
  pctl <- ecdf(mds1[[i]][-1])

  if(mds1[[i]][1] <= min(mds1[[i]][-1])){
    minpct <- 0
    maxpct <- pctl(min(mds1[[i]][-1]))
    if(maxpct == 0){
      L <- length(mds1[[i]][-1])
      maxpct <- (L-1)/L
    }
  } else{
    if(mds1[[i]][1] > max(mds1[[i]][-1])){
      minpct <- pctl(max(mds1[[i]][-1]))
      maxpct <- 1
      if(minpct == 1){
        L <- length(mds1[[i]][-1])
        minpct <- (L-1)/L
      }
    } else{
      minpct <- pctl(sort(mds1[[i]][-1])[max(which(sort(mds1[[i]][-1]) < mds1[[i]][1]))])
      maxpct <- pctl(mds1[[i]][1])
      if(minpct == 1){
        L <- length(mds1[[i]][-1])
        minpct <- (L-1)/L
      }
      if(maxpct == 0){
        L <- length(mds1[[i]][-1])
        maxpct <- (L-1)/L
      }
    }
  }

  mdps1[i] <- runif(1, minpct, maxpct)
}

mds1_full <- vector("list", length(which(singulars_full == 0)))
for(i in 1:length(which(singulars_full == 0))){
  yobs_full <- ly_true2_full[i,]
  ypred_full <- colMeans(newinitialY(sampdists2_full[[i]], J))
  sigmaypred_full <- cov(newinitialY(sampdists2_full[[i]], J))

  playeri_full <- rbind(yobs_full, newinitialY(sampdists2_full[[i]], J))

  mds1_full[[i]] <- apply(playeri_full, 1, mahalanobis, center = ypred_full, cov = sigmaypred_full)
}

mdps1_full <- double(length = length(which(singulars_full == 0)))
for(i in 1:length(which(singulars_full == 0))){
  pctl_full <- ecdf(mds1_full[[i]][-1])

  if(mds1_full[[i]][1] <= min(mds1_full[[i]][-1])){
    minpct_full <- 0
    maxpct_full <- pctl_full(min(mds1_full[[i]][-1]))
    if(maxpct_full == 0){
      L <- length(mds1_full[[i]][-1])
      maxpct_full <- (L-1)/L
    }
  } else{
    if(mds1_full[[i]][1] > max(mds1_full[[i]][-1])){
      minpct_full <- pctl_full(max(mds1_full[[i]][-1]))
      maxpct_full <- 1
      if(minpct_full == 1){
        L <- length(mds1_full[[i]][-1])
        minpct_full <- (L-1)/L
      }
    } else{
      minpct_full <- pctl_full(sort(mds1_full[[i]][-1])[max(which(sort(mds1_full[[i]][-1]) < mds1_full[[i]][1]))])
      maxpct_full <- pctl_full(mds1_full[[i]][1])
      if(minpct_full == 1){
        L <- length(mds1_full[[i]][-1])
        minpct_full <- (L-1)/L
      }
      if(maxpct_full == 0){
        L <- length(mds1_full[[i]][-1])
        maxpct_full <- (L-1)/L
      }
    }
  }

  mdps1_full[i] <- runif(1, minpct_full, maxpct_full)
}

# Convert to Standard Normal and Make Plots #
# Without Covariates #
mdrs1 <- qnorm(mdps1)

png("Figure16a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1, col.lines = "red", line = "none", main = "Null MLOGIT Bird Count Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure16b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1, main = "Null MLOGIT Bird Count Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure16.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1, col.lines = "red", line = "none", main = "Null MLOGIT Bird Count Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1, main = "Null MLOGIT Bird Count Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# With Covariates #
mdrs1_full <- qnorm(mdps1_full)

png("Figure17a.png", units = "in", width = 5, height = 5, res = 600)
car::qqPlot(mdrs1_full, col.lines = "red", line = "none", main = "MLOGIT w/ X Bird Count Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
dev.off()
png("Figure17b.png", units = "in", width = 5, height = 5, res = 600)
plot(mdrs1_full, main = "MLOGIT w/ X Bird Count Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

png("Figure17.png", units = "in", width = 12, height = 6, res = 600)
par(mfrow=c(1,2), cex=1.5)
car::qqPlot(mdrs1_full, col.lines = "red", line = "none", main = "MLOGIT w/ X Bird Count Model", ylab = "MD2 Residuals", xlab = "Theoretical Quantiles", envelope = c(style = "lines"))
abline(0, 1, col = "red", lwd = 2)
plot(mdrs1_full, main = "MLOGIT w/ X Bird Count Model", ylab = "MD2 Residuals")
abline(h = 0, col = "red")
dev.off()

# Kolmogorov-Smirnov Tests #
ks.test(mdrs1, "pnorm", mean = 0, sd = 1)
ks.test(mdrs1_full, "pnorm", mean = 0, sd = 1)
