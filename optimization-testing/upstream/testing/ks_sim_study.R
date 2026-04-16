# Kolmogorov-Smirnov Tests
require(iterators)
require(doParallel)

#First, do a simulation study with 100 observations but different simulated data sets
start.simstudy <- Sys.time()
sims <- 100

registerDoParallel(cores = 4)

ks.log <- foreach(s=1:sims, .combine = rbind) %dopar% {
  #Load Packages (must be done within foreach)
  require(mvnfast) #for multivariate normal sampling
  require(MCMCpack) #for inverse wishart sampling
  require(EnvStats) #for uniform qqPlot
  require(car) #for confidence envelope in qqPlot
  require(MGLM) #for multinomial logistic regression with counts
  
  #Observed Log-Odds Function
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
  
  #MLN Functions
  #Functions and algorithms
  {
    #Inverse of matrix using Cholesky Decomposition (faster than solve())
    cholinv <- function(x){
      chol2inv(chol(x))
    }
    
    #For the Static Beta Proposal
    pstartoy <- function(pstarvec){
      y <- numeric(length=length(pstarvec))
      for(i in 1:length(pstarvec)){
        y[i] <- log(pstarvec[i]/(1-pstarvec[i]))
      }
      return(y)
    }
    
    ytopstar <- function(yvec){
      pstar <- numeric(length=length(yvec))
      for(i in 1:length(yvec)){
        pstar[i] <- exp(yvec[i])/(1 + exp(yvec[i]))
      }
      return(pstar)
    }
    
    betapropdist <- function(WMu_vec, Sigma){
      k <- ncol(Sigma)
      w_vec <- WMu_vec[1:(k+1)]
      mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
      result <- numeric(length=k)
      for(i in 1:k){
        alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
        if(alpha < 0){
          alpha <- alpha + 1
        }
        beta <- alpha*exp(-mu_vec[i])
        alpha_star <- w_vec[i] + alpha
        beta_star <- w_vec[k+1] + beta
        result[i] <- rbeta(1, alpha_star, beta_star)
      }
      return(result)
    }
    
    betaloglike <- function(WMuPstar_vec, Sigma){
      k <- ncol(Sigma)
      w_vec <- WMuPstar_vec[1:(k+1)]
      mu_vec <- WMuPstar_vec[(k+2):(2*k+1)]
      pstar_vec <- WMuPstar_vec[(2*k+2):length(WMuPstar_vec)]
      loglike <- 0
      logjac <- 0
      for(i in 1:k){
        alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
        if(alpha < 0){
          alpha <- alpha + 1
        }
        beta <- alpha*exp(-mu_vec[i])
        alpha_star <- w_vec[i] + alpha
        beta_star <- w_vec[k+1] + beta
        loglike <- loglike + dbeta(pstar_vec[i], alpha_star, beta_star, log = T)
        logjac <- logjac + log(pstar_vec[i]) + log(1 - pstar_vec[i])
      }
      return(loglike + logjac)
    }
    
    #For the Static Normal Approx to Beta Proposal
    normbetapropdist <- function(WMu_vec, Sigma){
      k <- ncol(Sigma)
      w_vec <- WMu_vec[1:(k+1)]
      mu_vec <- WMu_vec[(k+2):length(WMu_vec)]
      result <- numeric(length=k)
      for(i in 1:k){
        alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
        beta <- alpha*exp(-mu_vec[i])
        alpha_star <- w_vec[i] + alpha
        beta_star <- w_vec[k+1] + beta
        muprop <- digamma(alpha_star) - digamma(beta_star)
        sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
        result[i] <- rnorm(1, muprop, sigprop)
      }
      return(result)
    }
    
    normbetaloglike <- function(WMuY_vec, Sigma){
      k <- ncol(Sigma)
      w_vec <- WMuY_vec[1:(k+1)]
      mu_vec <- WMuY_vec[(k+2):(2*k+1)]
      y_vec <- WMuY_vec[(2*k+2):length(WMuY_vec)]
      result <- 0
      for(i in 1:k){
        alpha <- ((1 + exp(mu_vec[i])) / Sigma[i,i]) - (exp(mu_vec[i]) / (1 + exp(mu_vec[i])))
        beta <- alpha*exp(-mu_vec[i])
        alpha_star <- w_vec[i] + alpha
        beta_star <- w_vec[k+1] + beta
        muprop <- digamma(alpha_star) - digamma(beta_star)
        sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
        result <- result + dnorm(y_vec[i], muprop, sigprop, log = T)
      }
      return(result)
    }
    
    #Quickly convert Y to Pi
    ytopi <- function(vec){
      c(exp(vec)/(1 + sum(exp(vec))), 1/(1 + sum(exp(vec))))
    }
    
    mixgibbs <- function(W, X, Z, randeff = 1, base, iter = 1000, mhiter = 1, mhscale = 1, proposal = c("norm", "beta", "normbeta"), block = TRUE, r.seed = 42, trY = NULL, trBeta = NULL, trPsi = NULL, trSigma = NULL, trPhi = NULL, startru = FALSE, fixY = FALSE, fixBeta = FALSE, fixPsi = FALSE, fixSigma = FALSE, fixPhi = FALSE, whichYstart = c("mean", "inity"), mixed = TRUE, test = FALSE){
      
      #Set seed for replicability
      set.seed(r.seed)
      
      #Mixed Model
      if(mixed == TRUE){
        
        #Make sure base is the last category
        W <- cbind(W[,-base], W[,base])
        base <- ncol(W)
        
        #Some dimensions
        p <- ncol(X) #number of covariates, including intercept
        d <- ncol(W) - 1 #number of multivariate categories
        N <- nrow(X) #number of total observations
        q <- randeff #number of random player effects (default is simply an intercept)
        m <- ncol(Z)/randeff #number of players
        
        #Algorithm uses a wide X
        X <- t(X) #should be p by N
        
        #Get initial estimates of the log-odds
        if(fixY == FALSE){
          Y <- t(newinitialY(W, base)) #should be d by N
        } else{
          Y <- t(trY)
        }
        
        #Beta set-up
        if(fixBeta == FALSE){
          prior_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
          betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
          betatrack[1,] <- as.vector(prior_beta)
          if(startru == TRUE){
            prior_beta <- trBeta
            betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
            betatrack[1,] <- as.vector(prior_beta)
          }
        } else{
          prior_beta <- trBeta
          betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
          betatrack[1,] <- as.vector(prior_beta)
        }
        
        #Psi set-up
        if(fixPsi == FALSE){
          #prior_psi <- apply(rbind(colSums(Z), c(0, colSums(Z)[-m])), 2, function(x) colMeans(t(Y - tcrossprod(prior_beta, t(X)))[(x[2]+1):(x[1]+x[2]),]) )
          prior_psi <- matrix(0, nrow = d, ncol = q*m)
          psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
          psitrack[1,] <- as.vector(prior_psi)
          if(startru == TRUE){
            prior_psi <- t(trPsi) #should be d by q*m
            psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
            psitrack[1,] <- as.vector(prior_psi)
          }
        } else{
          prior_psi <- t(trPsi)
          psitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_psi))
          psitrack[1,] <- as.vector(prior_psi)
        }
        
        #Sigma set-up
        #Hyperparameters for prior distribution of Sigma
        v_1 <- d
        Lambda_1 <- diag(d) + matrix(1, nrow = d)%*%matrix(1, ncol = d)
        
        if(fixSigma == FALSE){
          prior_sigma <- diag(d)
          sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
          sigmatrack[1,] <- as.vector(prior_sigma)
          if(startru == TRUE){
            prior_sigma <- trSigma
            sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
            sigmatrack[1,] <- as.vector(prior_sigma)
          }
        } else{
          prior_sigma <- trSigma
          sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
          sigmatrack[1,] <- as.vector(prior_sigma)
        }
        
        #Phi set-up
        #Hyperparameters for prior distribution of Phi
        v_2 <- d*q
        Lambda_2 <- diag(d*q) + matrix(1, nrow = d*q)%*%matrix(1, ncol = d*q)
        
        if(fixPhi == FALSE){
          prior_phi <- diag(d*q)
          phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
          phitrack[1,] <- as.vector(prior_phi)
          if(startru == TRUE){
            prior_phi <- trPhi
            phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
            phitrack[1,] <- as.vector(prior_phi)
          }
        } else{
          prior_phi <- trPhi
          phitrack <- matrix(0, nrow = iter + 1, ncol = length(prior_phi))
          phitrack[1,] <- as.vector(prior_phi)
        }
        
        #Get initial log-odds means, and then set starting values of Y to that (except doing that leads to horrible MH acceptance ratio when moving to more categories)
        prior_mu <- tcrossprod(prior_beta, t(X)) + prior_psi%*%as.matrix(t(Z))
        
        #Starting values for Y
        if(fixY == FALSE){
          if(whichYstart == "mean"){
            Y <- prior_mu
          }
        } else{
          Y <- t(trY)
        }
        
        ytrack <- matrix(0, nrow = iter + 1, ncol = d*N)
        ytrack[1,] <- as.vector(Y)
        
        #Certain values will stay fixed throughout the Gibbs chain
        v_1n <- v_1 + N - d*p #posterior df of Sigma
        v_2n <- v_2 + m #posterior df of Phi
        S_xx <- tcrossprod(X) #X sums of squares matrix
        
        #Finally, also want to keep track of MH acceptance probabilities for proposing Y
        MHacceptancelist <- matrix(0, nrow = iter, ncol = mhiter)
        
        #Full posterior log-likelihood and error term tracking can be recovered from the resulting output matrices, and to save computation time are ommitted from the main function
        
        #The MH-within-Gibbs Algorithm
        for(i in 1:iter){
          
          #        if(i %in% seq(100, iter, round(iter/100))){
          #          print(i)
          #        }
          
          #Metropolis-Hastings step
          #Y|W, X, Z, Beta, Psi, Sigma, Phi
          
          if(fixY == FALSE){
            for(j in 1:mhiter){
              Y <- t(Y) #put it long ways for proposing
              Mu <- t(prior_mu) #put it long ways for proposing
              wmu <- cbind(W, Mu) #for the proposal function
              wmuy <- cbind(wmu, Y)
              
              if(proposal == "norm"){
                Y_new <- Y + rmvn(N, integer(length = d), mhscale*prior_sigma) #propose centered at old Y
              }
              
              if(proposal == "beta"){
                Pstar <- t(apply(X = Y, MARGIN = 1, FUN = ytopstar)) #old Pstar
                Pstar_new <- t(apply(X = wmu, MARGIN = 1, FUN = betapropdist, Sigma = mhscale*prior_sigma))
                Y_new <- t(apply(X = Pstar_new, MARGIN = 1, FUN = pstartoy)) #newY based on proposed probabilities
                wmuy_new <- cbind(wmu, Y_new)
                wmuPstar <- cbind(wmu, Pstar)
                wmuPstar_new <- cbind(wmu, Pstar_new)
              }
              
              if(proposal == "normbeta"){
                Y_new <- t(apply(X = wmu, MARGIN = 1, FUN = normbetapropdist, Sigma = mhscale*prior_sigma)) #Fixed Normal proposal based on logit(Beta) expected values and variances
                wmuy_new <- cbind(wmu, Y_new)
              }
              
              sexpY <- rowSums(exp(Y))
              sexpYn <- rowSums(exp(Y_new))
              
              newnormpart <- .5*tcrossprod(tcrossprod((Y_new - Mu) , t(cholinv(prior_sigma))), (Y_new - Mu))
              oldnormpart <- .5*tcrossprod(tcrossprod((Y - Mu) , t(cholinv(prior_sigma))), (Y - Mu))
              newloglike <- rowSums(as.matrix(W[,1:d]*Y_new[,1:d])) - rowSums(W*(log1p(sexpYn))) - newnormpart[cbind(1:nrow(newnormpart), 1:nrow(newnormpart))]
              oldloglike <- rowSums(as.matrix(W[,1:d]*Y[,1:d])) - rowSums(W*(log1p(sexpY))) - oldnormpart[cbind(1:nrow(oldnormpart), 1:nrow(oldnormpart))]
              
              rm(newnormpart, oldnormpart) #remove from memory to help efficiency
              
              if(proposal == "norm"){
                #Random walk is symmetric, Metropolis, not MH
                logptolike <- 0
                logotplike <- 0
                
                rm(wmu, wmuy)
              }
              
              if(proposal == "beta"){
                logptolike <- apply(X = wmuPstar, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to old p_star under proposal distribution, including Jacobian
                logotplike <- apply(X = wmuPstar_new, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to proposed p_star under proposal, including Jacobian
                
                rm(wmu, wmuy, wmuy_new, wmuPstar, wmuPstar_new)
              }
              
              if(proposal == "normbeta"){
                logptolike <- apply(X = wmuy, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
                logotplike <- apply(X = wmuy_new, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
                
                rm(wmu, wmuy, wmuy_new)
              }
              
              ratio <- newloglike - oldloglike + logptolike - logotplike
              temp <- runif(N)
              difr <- as.numeric((ratio - log(temp)) >= 0)
              MHacceptancelist[i,j] <- mean(difr)
              Y <- difr*Y_new + (1 - difr)*Y #accept or reject samples to get draw
              Y <- t(Y) #put it back for the rest of the Gibbs
            }
            
            ytrack[i+1,] <- as.vector(Y)
          } else{
            #To Troubleshoot: keep Y fixed
            ytrack[i+1,] <- as.vector(Y)
          }
          
          #Using sampled Y, create necessary sums of squares matrices and OLS mean estimates for:
          #individual update of Beta
          if(fixBeta == FALSE){
            S_yx <- tcrossprod((Y - prior_psi%*%as.matrix(t(Z))), X)
            beta_hat <- tcrossprod(S_yx, t(cholinv(S_xx)))
            
            #block update of Beta
            margS_yx <- tcrossprod(Y, X)
            margbeta_hat <- tcrossprod(margS_yx, t(cholinv(S_xx)))
          }
          
          #Block update of Beta and Psi is more efficient than individual updates
          
          if(block == TRUE){
            
            #Marginal Update of Beta, not conditional on Psi
            #Beta|Y, X, Z, Sigma, Phi
            if(fixBeta == FALSE){
              sig_beta_inv <- matrix(0, ncol = p*d, nrow = p*d)
              for(k1 in 1:m){
                Zi <- as.matrix(t(Z[which(t(Z)[(k1*q - (q - 1)),]==1),(k1*q - (q - 1)):(k1*q)]))
                n_i <- ncol(Zi)
                Xi <- X[,which(t(Z)[(k1*q - (q - 1)),]==1)]
                if(nrow(X) == 1){
                  Xi <- matrix(Xi, nrow = 1)
                }
                sig_beta_inv <- sig_beta_inv + kronecker(tcrossprod(Xi), cholinv(n_i*prior_phi + prior_sigma))
              }
              
              sig_beta <- cholinv(sig_beta_inv)
              mu_beta <- margbeta_hat
              
              post_beta_vec <- rmvn(1, mu_beta, sig_beta)
              prior_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
              betatrack[i+1,] <- as.vector(prior_beta)
            } else{
              #To Troubleshoot: keep Beta fixed
              post_beta <- prior_beta
              betatrack[i+1,] <- as.vector(prior_beta)
            }
            
            #Psi| Beta, Y, X, Z, Sigma, Phi
            if(fixPsi == FALSE){
              post_psi <- matrix(0, nrow = d*q, ncol = m)
              for(k in 1:m){
                Xi <- X[,which(t(Z)[(k*q - (q - 1)),]==1)]
                Yi <- Y[,which(t(Z)[(k*q - (q - 1)),]==1)]
                Zi <- as.matrix(t(Z[which(t(Z)[(k*q - (q - 1)),]==1),(k*q - (q - 1)):(k*q)]))
                Ui <- cholinv(cholinv(prior_phi) + kronecker(cholinv(prior_sigma) , as.matrix(Zi%*%t(Zi))))
                
                sigZi <- kronecker(cholinv(prior_sigma), Zi)
                errvec <- matrix(as.vector(t(Yi - prior_beta%*%Xi)), byrow = T, ncol = 1)
                psi_hat <- Ui%*%sigZi%*%errvec
                
                post_psivec <- rmvn(1, as.vector(psi_hat), Ui)
                post_psi[,k] <- post_psivec
              }
              prior_psi_star <- post_psi #need it in alternate form for sampling Phi (d*q by m)
              post_psi <- matrix(as.vector(post_psi), nrow = d) #for the rest of the time (d by q*m)
              psitrack[i+1,] <- as.vector(post_psi)
              prior_psi <- post_psi
            } else{
              #To troubleshoot: keep Psi fixed
              prior_psi_star <- matrix(as.vector(prior_psi), ncol = m) #Not entirely sure this will work with multiple random effects
              post_psi <- prior_psi
              psitrack[i+1,] <- as.vector(post_psi)
            }
            
          }
          
          #Individual updates of Beta and Psi change the order a bit, to improve efficiency
          #Psi| Beta, Y, X, Z, Sigma, Phi
          
          if(block == FALSE){
            
            #Same as in block update
            if(fixPsi == FALSE){
              post_psi <- matrix(0, nrow = d*q, ncol = m)
              for(k in 1:m){
                Xi <- X[,which(t(Z)[(k*q - (q - 1)),]==1)]
                Yi <- Y[,which(t(Z)[(k*q - (q - 1)),]==1)]
                Zi <- as.matrix(t(Z[which(t(Z)[(k*q - (q - 1)),]==1),(k*q - (q - 1)):(k*q)]))
                Ui <- cholinv(cholinv(prior_phi) + kronecker(cholinv(prior_sigma) , as.matrix(Zi%*%t(Zi))))
                
                
                sigZi <- kronecker(cholinv(prior_sigma), Zi)
                errvec <- matrix(as.vector(t(Yi - prior_beta%*%Xi)), byrow = T, ncol = 1)
                psi_hat <- Ui%*%sigZi%*%errvec
                
                post_psivec <- rmvn(1, as.vector(psi_hat), Ui)
                post_psi[,k] <- post_psivec
              }
              prior_psi_star <- post_psi
              post_psi <- matrix(as.vector(post_psi), nrow = d)
              psitrack[i+1,] <- as.vector(post_psi)
              prior_psi <- post_psi
            } else{
              #To troubleshoot: keep Psi fixed
              prior_psi_star <- matrix(as.vector(prior_psi), ncol = m) #Not entirely sure this will work with multiple random effects
              post_psi <- prior_psi
              psitrack[i+1,] <- as.vector(post_psi)
            }
            
          }
          
          #Phi|Psi
          
          if(fixPhi == FALSE){
            Lambda_2n <- Lambda_2 + tcrossprod(prior_psi_star)
            post_phi <- riwish(v = v_2n, S = Lambda_2n)
            phitrack[i+1,] <- as.vector(post_phi)
            prior_phi <- post_phi
          } else{
            #To troubleshoot: keep Phi fixed
            post_phi <- prior_phi
            phitrack[i+1,] <- as.vector(prior_phi)
          }
          
          if(fixBeta == FALSE){
            #Y given X sums of squares matrix for Sigma constructed using beta_hat and updated Psi
            if(block == FALSE){
              S_ygx <- tcrossprod((Y - tcrossprod(beta_hat, t(X)) - prior_psi%*%as.matrix(t(Z))))
            } else{
              S_ygx <- tcrossprod((Y - tcrossprod(prior_beta, t(X)) - prior_psi%*%as.matrix(t(Z))))
            }
          } else{
            S_ygx <- tcrossprod((Y - tcrossprod(prior_beta, t(X)) - prior_psi%*%as.matrix(t(Z))))
          }
          
          
          #Sigma| Y, X, Beta, Psi
          
          if(fixSigma == FALSE){
            Lambda_1n <- Lambda_1 + S_ygx
            post_sigma <- riwish(v = v_1n, S = Lambda_1n)
            sigmatrack[i+1,] <- as.vector(post_sigma)
            prior_sigma <- post_sigma
          } else{
            #To Troubleshoot: keep Sigma fixed
            post_sigma <- prior_sigma
            sigmatrack[i+1,] <- as.vector(post_sigma)
          }
          
          #Beta| Y, X, Z, Psi, Sigma
          
          if(block == FALSE){
            
            if(fixBeta == FALSE){
              post_beta_vec_cov <- kronecker(cholinv(S_xx), prior_sigma)
              post_beta_0 <- beta_hat
              post_beta_0_vec <- as.vector(post_beta_0)
              post_beta_vec <- rmvn(n = 1, mu = post_beta_0_vec, sigma = post_beta_vec_cov)
              post_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
              betatrack[i+1,] <- as.vector(post_beta)
              prior_beta <- post_beta
            } else{
              #To Troubleshoot: keep Beta fixed
              post_beta <- prior_beta
              betatrack[i+1,] <- as.vector(post_beta)
            }
            
          }
          
          if(fixBeta == FALSE){
            rm(S_yx, margS_yx, beta_hat, margbeta_hat) #clean up memory
          }
          
          #Update log-odds means, Mu
          prior_mu <- tcrossprod(prior_beta, t(X)) + prior_psi%*%as.matrix(t(Z))
          
        }
        
        paths <- list(Y = ytrack, Sigma = sigmatrack, Beta = betatrack, Psi = psitrack, Phi = phitrack, MHRatio = MHacceptancelist)
        return(paths)
        
      }
      
      #Fixed Model
      if(mixed == FALSE){
        
        #Some dimensions
        p <- ncol(X) #number of covariates, including intercept
        d <- ncol(W) - 1 #number of multivariate categories
        N <- nrow(X) #number of total observations
        
        #Algorithm uses a wide X
        X <- t(X) #should be p by N
        
        #Get initial estimates of the log-odds
        if(fixY == FALSE){
          Y <- t(newinitialY(W, base)) #should be d by N
        } else{
          Y <- t(trY)
        }
        
        #Beta set-up
        if(fixBeta == FALSE){
          prior_beta <- t(tcrossprod(cholinv(tcrossprod(X)), tcrossprod(Y, X)))
          betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
          betatrack[1,] <- as.vector(prior_beta)
          if(startru == TRUE){
            prior_beta <- trBeta
            betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
            betatrack[1,] <- as.vector(prior_beta)
          }
        } else{
          prior_beta <- trBeta
          betatrack <- matrix(0, nrow = iter + 1, ncol = d*p)
          betatrack[1,] <- as.vector(prior_beta)
        }
        
        #Sigma set-up
        #Hyperparameters for prior distribution of Sigma
        v_1 <- d + 1
        Lambda_1 <- diag(d)
        
        if(fixSigma == FALSE){
          prior_sigma <- diag(d)
          sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
          sigmatrack[1,] <- as.vector(prior_sigma)
          if(startru == TRUE){
            prior_sigma <- trSigma
            sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
            sigmatrack[1,] <- as.vector(prior_sigma)
          }
        } else{
          prior_sigma <- trSigma
          sigmatrack <- matrix(0, nrow = iter + 1, ncol = length(prior_sigma))
          sigmatrack[1,] <- as.vector(prior_sigma)
        }
        
        #Get initial log-odds means, and then set starting values of Y to that (except doing that leads to horrible MH acceptance ratio when moving to more categories)
        prior_mu <- tcrossprod(prior_beta, t(X))
        
        #Starting values for Y
        if(fixY == FALSE){
          if(whichYstart == "mean"){
            Y <- prior_mu
          }
        } else{
          Y <- t(trY)
        }
        
        ytrack <- matrix(0, nrow = iter + 1, ncol = d*N)
        ytrack[1,] <- as.vector(Y)
        
        #Certain values will stay fixed throughout the Gibbs chain
        v_1n <- v_1 + N #posterior df of Sigma
        S_xx <- tcrossprod(X) #X sums of squares matrix
        
        #Finally, also want to keep track of MH acceptance probabilities for proposing Y
        MHacceptancelist <- matrix(0, nrow = iter, ncol = mhiter)
        
        #Full posterior log-likelihood and error term tracking can be recovered from the resulting output matrices, and to save computation time are ommitted from the main function
        
        #The MH-within-Gibbs Algorithm
        for(i in 1:iter){
          
          #Metropolis-Hastings step
          #Y|W, X, Z, Beta, Psi, Sigma, Phi
          
          if(fixY == FALSE){
            for(j in 1:mhiter){
              Y <- t(Y) #put it long ways for proposing
              Mu <- t(prior_mu) #put it long ways for proposing
              wmu <- cbind(W, Mu) #for the proposal function
              wmuy <- cbind(wmu, Y)
              
              if(proposal == "norm"){
                Y_new <- Y + rmvn(N, integer(length = d), mhscale*prior_sigma) #propose centered at old Y
              }
              
              if(proposal == "beta"){
                Pstar <- t(apply(X = Y, MARGIN = 1, FUN = ytopstar)) #old Pstar
                Pstar_new <- t(apply(X = wmu, MARGIN = 1, FUN = betapropdist, Sigma = mhscale*prior_sigma))
                Y_new <- t(apply(X = Pstar_new, MARGIN = 1, FUN = pstartoy)) #newY based on proposed probabilities
                wmuy_new <- cbind(wmu, Y_new)
                wmuPstar <- cbind(wmu, Pstar)
                wmuPstar_new <- cbind(wmu, Pstar_new)
              }
              
              if(proposal == "normbeta"){
                Y_new <- t(apply(X = wmu, MARGIN = 1, FUN = normbetapropdist, Sigma = mhscale*prior_sigma)) #Fixed Normal proposal based on logit(Beta) expected values and variances
                Y_new <- ifelse(Y_new > 500, 500, Y_new)
                wmuy_new <- cbind(wmu, Y_new)
              }
              
              sexpY <- rowSums(exp(Y))
              sexpYn <- rowSums(exp(Y_new))
              
              newnormpart <- .5*tcrossprod(tcrossprod((Y_new - Mu) , t(cholinv(prior_sigma))), (Y_new - Mu))
              oldnormpart <- .5*tcrossprod(tcrossprod((Y - Mu) , t(cholinv(prior_sigma))), (Y - Mu))
              newloglike <- rowSums(as.matrix(W[,1:d]*Y_new[,1:d])) - rowSums(W*(log1p(sexpYn))) - newnormpart[cbind(1:nrow(newnormpart), 1:nrow(newnormpart))]
              oldloglike <- rowSums(as.matrix(W[,1:d]*Y[,1:d])) - rowSums(W*(log1p(sexpY))) - oldnormpart[cbind(1:nrow(oldnormpart), 1:nrow(oldnormpart))]
              
              rm(newnormpart, oldnormpart) #remove from memory to help efficiency
              
              if(proposal == "norm"){
                #Random walk is symmetric, Metropolis, not MH
                logptolike <- 0
                logotplike <- 0
                
                rm(wmu, wmuy)
              }
              
              if(proposal == "beta"){
                logptolike <- apply(X = wmuPstar, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to old p_star under proposal distribution, including Jacobian
                logotplike <- apply(X = wmuPstar_new, MARGIN = 1, FUN = betaloglike, Sigma = mhscale*prior_sigma) #going to proposed p_star under proposal, including Jacobian
                
                rm(wmu, wmuy, wmuy_new, wmuPstar, wmuPstar_new)
              }
              
              if(proposal == "normbeta"){
                logptolike <- apply(X = wmuy, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
                logotplike <- apply(X = wmuy_new, MARGIN = 1, FUN = normbetaloglike, Sigma = mhscale*prior_sigma)
                
                rm(wmu, wmuy, wmuy_new)
              }
              
              ratio <- newloglike - oldloglike + logptolike - logotplike
              temp <- runif(N)
              difr <- as.numeric((ratio - log(temp)) >= 0)
              MHacceptancelist[i,j] <- mean(difr)
              Y <- difr*Y_new + (1 - difr)*Y #accept or reject samples to get draw
              Y <- t(Y) #put it back for the rest of the Gibbs
            }
            
            ytrack[i+1,] <- as.vector(Y)
          } else{
            #To Troubleshoot: keep Y fixed
            ytrack[i+1,] <- as.vector(Y)
          }
          
          #Using sampled Y, create necessary sums of squares matrices and OLS mean estimates for:
          #individual update of Beta
          if(fixBeta == FALSE){
            S_yx <- tcrossprod(Y, X)
            beta_hat <- tcrossprod(S_yx, t(cholinv(S_xx)))
          }
          
          #If Fixed, there should be no block update
          
          if(fixBeta == FALSE){
            #Y given X sums of squares matrix for Sigma constructed using beta_hat and updated Psi
            S_ygx <- tcrossprod(Y - tcrossprod(beta_hat, t(X)))
          } else{
            S_ygx <- tcrossprod(Y - tcrossprod(prior_beta, t(X)))
          }
          
          #Sigma| Y, X, Beta, Psi
          
          if(fixSigma == FALSE){
            Lambda_1n <- cholinv(Lambda_1) + S_ygx
            post_sigma <- riwish(v = v_1n, S = Lambda_1n)
            sigmatrack[i+1,] <- as.vector(post_sigma)
            prior_sigma <- post_sigma
          } else{
            #To Troubleshoot: keep Sigma fixed
            post_sigma <- prior_sigma
            sigmatrack[i+1,] <- as.vector(post_sigma)
          }
          
          #Beta| Y, X, Z, Psi, Sigma
          
          if(fixBeta == FALSE){
            post_beta_vec_cov <- kronecker(cholinv(S_xx), prior_sigma)
            post_beta_0 <- beta_hat
            post_beta_0_vec <- as.vector(post_beta_0)
            post_beta_vec <- rmvn(n = 1, mu = post_beta_0_vec, sigma = post_beta_vec_cov)
            post_beta <- matrix(post_beta_vec, nrow = d) #should be d by p
            betatrack[i+1,] <- as.vector(post_beta)
            prior_beta <- post_beta
          } else{
            #To Troubleshoot: keep Beta fixed
            post_beta <- prior_beta
            betatrack[i+1,] <- as.vector(post_beta)
          }
          
          if(fixBeta == FALSE){
            rm(S_yx, beta_hat) #clean up memory
          }
          
          #Update log-odds means, Mu
          prior_mu <- tcrossprod(prior_beta, t(X))
          
          #This if statement can be moved around for troubleshooting at different points
          if(test == TRUE){
            #if(i == 114){
            #browser()         
            print(i)
            #}
          }
          
        }
        
        paths <- list(Y = ytrack, Sigma = sigmatrack, Beta = betatrack, MHRatio = MHacceptancelist)
        return(paths)
        
      }
      
      
    }
    
    posteriorsampling.givenparams <- function(Exposure, X, dimW, sigma_mean, beta_mean, samples = 10000){
      
      temp_mu <- X%*%beta_mean
      
      np_y_pred <- rmvn(samples, temp_mu, sigma_mean)
      np_pi_pred <- cbind(exp(np_y_pred)/(1 + rowSums(exp(np_y_pred))), 1/(1 + rowSums(exp(np_y_pred))))
      np_w_pred <- t(apply(np_pi_pred, 1, rmultinom, n = 1, size = Exposure))
      
      #result <- list(predy = np_y_pred, predpi = np_pi_pred, predw = np_w_pred)
      return(np_w_pred)
    }
    
  }
  
  set.seed(s)
  #Create data sets
  {
    #M-Logit with quadratic term
    m <- 1000 #observations
    X <- sort(rnorm(m))   #covariate
    Xmat <- cbind(1, X, X^2)   #design matrix
    n <- 200   #exposure
    J <- 3   #number of  categories
    
#    B <- cbind(c(-1, .5, .25), c(-.75, -1, -.5))  #parameters
#    XB <- Xmat%*%B
#    Pi <- cbind(exp(XB)/(1 + rowSums(exp(XB))), 1/(1 + rowSums(exp(XB))))   #true probabilities
#    Y <- t(apply(Pi, 1, rmultinom, n = 1, size = n))   #simulate data
    
    #Create data frame
#    ml <- as.data.frame(cbind(Y, Xmat))
    
    #MLN (overdispersed) with no quadratic
    Xmat2 <- cbind(1, X)
    B2 <- cbind(c(-1, .5), c(-.75, -1))
#    XB2 <- Xmat2%*%B2 + rmvn(n = m, mu = numeric(length = J-1), sigma = diag(J-1))  #true log-odds, including additional variability
#    Pi2 <- cbind(exp(XB2)/(1 + rowSums(exp(XB2))), 1/(1 + rowSums(exp(XB2))))   #true probabilities
#    Y2 <- t(apply(Pi2, 1, rmultinom, n = 1, size = n))   #simulate data
    
    #Create data frame
#    ml2 <- as.data.frame(cbind(Y2, Xmat2))
    
    #DM (overdispersed/negative correlations) with no quadratic
    B3 <- rbind(c(-.09, .12, .61), c(.25, -.15, .05))
    XB3 <- Xmat2%*%B3
    alpha <- exp(XB3) #true dirichlet parameter
    Y3.2 <- t(apply(alpha, 1, rdirmn, n = 1, size = n)) #simulate data
    
    #Create data frame
    ml3 <- as.data.frame(cbind(Y3.2, Xmat2))
    
    #MLN (overdispersed) with no quadratic (and positive correlations)
    epsilon2 <- rmvn(n = m, mu = numeric(length = J-1), sigma = matrix(c(1,-.9,-.9,4), byrow=T, ncol=2))
    XB4 <- Xmat2%*%B2 + epsilon2  #true log-odds, including additional variability
    Pi4 <- cbind(exp(XB4)/(1 + rowSums(exp(XB4))), 1/(1 + rowSums(exp(XB4))))   #true probabilities
    #cor(Pi4)
    Y4 <- t(apply(Pi4, 1, rmultinom, n = 1, size = n))   #simulate data
    
    #Create data frame
    ml4 <- as.data.frame(cbind(Y4, Xmat2))
    
    #Observed Log-Odds
#    ly_true <- newinitialY(Y, base = J)
#    ly2_true <- newinitialY(Y2, base = J)
    ly3_true <- newinitialY(Y3.2, base = J)
    ly4_true <- newinitialY(Y4, base = J)
  }
  
  #Model Fitting
#  correct_mlogit <- MGLMreg(cbind(V1, V2, V3) ~ X + V6, data = ml, dist = "MN")
#  incorrect_mlogit <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml, dist = "MN")
#  correct_mln <- mixgibbs(W = Y2, X = Xmat2, Z = diag(1), base = J, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE, r.seed = s)
#  incorrect_mln <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml2, dist = "MN")
  correct_dm <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml3, dist = "DM")
  incorrect_dm <- mixgibbs(W = Y3.2, X = Xmat2, Z = diag(1), base = J,  mhscale = 1, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)
#  incorrect_dm <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml3, dist = "MN") #FIT DM WITH MN
  correct_mln_pos <- mixgibbs(W = Y4, X = Xmat2, Z = diag(1), base = J, mhscale = 1, iter = 1000, proposal = "normbeta", whichYstart = "inity", mixed = FALSE)
  incorrect_mln_dm <- MGLMreg(cbind(V1, V2, V3) ~ X, data = ml4, dist = "DM")
  
  #Construct the sampling distribution of the fitted models for each individual
  {
    #Correct MLOGIT
#    sampdists_cmlogit <- list(length = m)   #make the list to populate
#    for(i in 1:m){
#      sampdists_cmlogit[[i]] <- t(rmultinom(n = 1000, size = n, prob = correct_mlogit@fitted[i,]))   #draw a bunch of times from fitted model; I verified that the fitted values from the fit are the same as the probabilities from inverse log-ratio of XB_hat
#    }
    
    #Incorrect MLOGIT
#    sampdists_imlogit <- list(length = m)
#    for(i in 1:m){
#      sampdists_imlogit[[i]] <- t(rmultinom(n = 1000, size = n, prob = incorrect_mlogit@fitted[i,]))
#    }
    
    #Correct MLN
#    burnthin <- seq(401, 1001, 2)
    
#    posteriormean.y <- matrix(colMeans(correct_mln$Y[burnthin,]), ncol = J-1, byrow = T)
#    posteriormean.beta <- matrix(colMeans(correct_mln$Beta[burnthin,]), ncol = J-1, byrow = T)
#    posteriormean.sigma <- matrix(colMeans(correct_mln$Sigma[burnthin,]), ncol = J-1, byrow = T)
    
#    Y_samples <- correct_mln$Y[burnthin,]
#    Beta_samples <- correct_mln$Beta[burnthin,]
#    Sigma_samples <- correct_mln$Sigma[burnthin,]
    
#    sampdists_cmln <- vector(mode = "list", length = m)
#    for(i in 1:m){
#      sampdists_cmln[[i]] <- posteriorsampling.givenparams(Exposure = sum(Y2[i,]), X = t(as.matrix(Xmat2[i,])), dimW = J, sigma_mean = posteriormean.sigma, beta_mean = posteriormean.beta)
#    }
    
    #Incorrect MLN
#    sampdists_imln <- list(length = m)
#    for(i in 1:m){
#      sampdists_imln[[i]] <- t(rmultinom(n = 1000, size = n, prob = incorrect_mln@fitted[i,]))
#    }
    
    #Correct DM
    pred_alpha <- exp(Xmat2%*%correct_dm@coefficients)
    sampdists_cdm <- list(length = m)
    for(i in 1:m){
      sampdists_cdm[[i]] <- rdirmn(n = 10000, size = n, alpha = pred_alpha[i,])
    }
    
    #Incorrect DM (MLN)
    burnthin <- seq(401, 1001, 2)

    posteriormean.y2 <- matrix(colMeans(incorrect_dm$Y[burnthin,]), ncol = J-1, byrow = T)
    posteriormean.beta2 <- matrix(colMeans(incorrect_dm$Beta[burnthin,]), ncol = J-1, byrow = T)
    posteriormean.sigma2 <- matrix(colMeans(incorrect_dm$Sigma[burnthin,]), ncol = J-1, byrow = T)

    Y_samples2 <- incorrect_dm$Y[burnthin,]
    Beta_samples2 <- incorrect_dm$Beta[burnthin,]
    Sigma_samples2 <- incorrect_dm$Sigma[burnthin,]

    sampdists_idm <- vector(mode = "list", length = m)
    for(i in 1:m){
      sampdists_idm[[i]] <- posteriorsampling.givenparams(Exposure = sum(Y3.2[i,]), X = t(as.matrix(Xmat2[i,])), dimW = J, sigma_mean = posteriormean.sigma2, beta_mean = posteriormean.beta2)
    }
    
    #Incorrect DM (MN)
    # sampdists_idm <- list(length = m)
    # for(i in 1:m){
    #   sampdists_idm[[i]] <- t(rmultinom(n = 1000, size = n, prob = incorrect_dm@fitted[i,]))
    # }
    
    #Correct MLN+
    burnthin <- seq(401, 1001, 2)

    posteriormean.y3 <- matrix(colMeans(correct_mln_pos$Y[burnthin,]), ncol = J-1, byrow = T)
    posteriormean.beta3 <- matrix(colMeans(correct_mln_pos$Beta[burnthin,]), ncol = J-1, byrow = T)
    posteriormean.sigma3 <- matrix(colMeans(correct_mln_pos$Sigma[burnthin,]), ncol = J-1, byrow = T)

    Y_samples3 <- correct_mln_pos$Y[burnthin,]
    Beta_samples3 <- correct_mln_pos$Beta[burnthin,]
    Sigma_samples3 <- correct_mln_pos$Sigma[burnthin,]

    sampdists_cmlnp <- vector(mode = "list", length = m)
    for(i in 1:m){
      sampdists_cmlnp[[i]] <- posteriorsampling.givenparams(Exposure = sum(Y4[i,]), X = t(as.matrix(Xmat2[i,])), dimW = J, sigma_mean = posteriormean.sigma3, beta_mean = posteriormean.beta3)
    }
    
    #Incorrect MLN+ (DM)
    pred_alpha_mlndm <- exp(Xmat2%*%incorrect_mln_dm@coefficients)
    sampdists_imlndm <- list(length = m)
    for(i in 1:m){
      sampdists_imlndm[[i]] <- rdirmn(n = 10000, size = n, alpha = pred_alpha_mlndm[i,])
    }
    
    
  }
  
  #Determining which ones are estimable
  {
    #Correct MLOGIT
#    singulars_A <- unlist(lapply(sampdists_cmlogit, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
#    if(sum(singulars_A) > 0){
#      ly_true_A <- ly_true[-which(singulars_A == 1),]
#      sampdists_A <- sampdists_cmlogit[-which(singulars_A == 1)]
#    } else{
#      ly_true_A <- ly_true
#      sampdists_A <- sampdists_cmlogit
#    }
    
    #Incorrect MLOGIT (MN NO QUAD)
#    singulars_B <- unlist(lapply(sampdists_imlogit, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
#    if(sum(singulars_B) > 0){
#      ly_true_B <- ly_true[-which(singulars_B == 1),]
#      sampdists_B <- sampdists_imlogit[-which(singulars_B == 1)]
#    } else{
#      ly_true_B <- ly_true
#      sampdists_B <- sampdists_imlogit
#    }
    
    #Correct MLN
#    singulars_C <- unlist(lapply(sampdists_cmln, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
#    if(sum(singulars_C) > 0){
#      ly_true_C <- ly2_true[-which(singulars_C == 1),]
#      sampdists_C <- sampdists_cmln[-which(singulars_C == 1)]
#    } else{
#      ly_true_C <- ly2_true
#      sampdists_C <- sampdists_cmln
#    }
    
    #Incorrect MLN (MN)
    # singulars_D <- unlist(lapply(sampdists_imln, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    # 
    # #If there aren't any, no need to change things
    # if(sum(singulars_D) > 0){
    #   ly_true_D <- ly2_true[-which(singulars_D == 1),]
    #   sampdists_D <- sampdists_imln[-which(singulars_D == 1)]
    # } else{
    #   ly_true_D <- ly2_true
    #   sampdists_D <- sampdists_imln
    # }
    
    #Correct DM
    singulars_E <- unlist(lapply(sampdists_cdm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
    if(sum(singulars_E) > 0){
      ly_true_E <- ly3_true[-which(singulars_E == 1),]
      sampdists_E <- sampdists_cdm[-which(singulars_E == 1)]
    } else{
      ly_true_E <- ly3_true
      sampdists_E <- sampdists_cdm
    }
    
    #Incorrect DM (Either MN OR MLN)
    singulars_F <- unlist(lapply(sampdists_idm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
    if(sum(singulars_F) > 0){
      ly_true_F <- ly3_true[-which(singulars_F == 1),]
      sampdists_F <- sampdists_idm[-which(singulars_F == 1)]
    } else{
      ly_true_F <- ly3_true
      sampdists_F <- sampdists_idm
    }
    
    #Correct MLN+
    singulars_G <- unlist(lapply(sampdists_cmlnp, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify

    #If there aren't any, no need to change things
    if(sum(singulars_G) > 0){
      ly_true_G <- ly4_true[-which(singulars_G == 1),]
      sampdists_G <- sampdists_cmlnp[-which(singulars_G == 1)]
    } else{
      ly_true_G <- ly4_true
      sampdists_G <- sampdists_cmlnp
    }
    
    #Incorrect MLN+ (DM)
    singulars_H <- unlist(lapply(sampdists_imlndm, function(x) ifelse(any(diag(cov(x)) == 0), 1, 0)))   #identify
    
    #If there aren't any, no need to change things
    if(sum(singulars_H) > 0){
      ly_true_H <- ly4_true[-which(singulars_H == 1),]
      sampdists_H <- sampdists_imlndm[-which(singulars_H == 1)]
    } else{
      ly_true_H <- ly4_true
      sampdists_H <- sampdists_imlndm
    }
    
  }
  
  #Calculate the squared Mahalanobis distances
  {
    #Correct MLOGIT
#    mds_A <- vector("list", length(which(singulars_A == 0)))
#    for(i in 1:length(which(singulars_A == 0))){
#      yobs_A <- ly_true_A[i,]
#      ypred_A <- colMeans(newinitialY(sampdists_A[[i]], J))
#      sigmaypred_A <- cov(newinitialY(sampdists_A[[i]], J))
      
#      playeri_A <- rbind(yobs_A, newinitialY(sampdists_A[[i]], J))
      
#      mds_A[[i]] <- apply(playeri_A, 1, mahalanobis, center = ypred_A, cov = sigmaypred_A)
#    }
    
    #Incorrect MLOGIT
#    mds_B <- vector("list", length(which(singulars_B == 0)))
#    for(i in 1:length(which(singulars_B == 0))){
#      yobs_B <- ly_true_B[i,]
#      ypred_B <- colMeans(newinitialY(sampdists_B[[i]], J))
#      sigmaypred_B <- cov(newinitialY(sampdists_B[[i]], J))
      
#      playeri_B <- rbind(yobs_B, newinitialY(sampdists_B[[i]], J))
      
#      mds_B[[i]] <- apply(playeri_B, 1, mahalanobis, center = ypred_B, cov = sigmaypred_B)
#    }
    
    #Correct MLN
#    mds_C <- vector("list", length(which(singulars_C == 0)))
#    for(i in 1:length(which(singulars_C == 0))){
#      yobs_C <- ly_true_C[i,]
#      ypred_C <- colMeans(newinitialY(sampdists_C[[i]], J))
#      sigmaypred_C <- cov(newinitialY(sampdists_C[[i]], J))
      
#      playeri_C <- rbind(yobs_C, newinitialY(sampdists_C[[i]], J))
      
#      mds_C[[i]] <- apply(playeri_C, 1, mahalanobis, center = ypred_C, cov = sigmaypred_C)
#    }
    
    #Incorrect MLN
#    mds_D <- vector("list", length(which(singulars_D == 0)))
#    for(i in 1:length(which(singulars_D == 0))){
#      yobs_D <- ly_true_D[i,]
#      ypred_D <- colMeans(newinitialY(sampdists_D[[i]], J))
#      sigmaypred_D <- cov(newinitialY(sampdists_D[[i]], J))
      
#      playeri_D <- rbind(yobs_D, newinitialY(sampdists_D[[i]], J))
      
#      mds_D[[i]] <- apply(playeri_D, 1, mahalanobis, center = ypred_D, cov = sigmaypred_D)
#    }
    
    #Correct DM
    mds_E <- vector("list", length(which(singulars_E == 0)))
    for(i in 1:length(which(singulars_E == 0))){
      yobs_E <- ly_true_E[i,]
      ypred_E <- colMeans(newinitialY(sampdists_E[[i]], J))
      sigmaypred_E <- cov(newinitialY(sampdists_E[[i]], J))
      
      playeri_E <- rbind(yobs_E, newinitialY(sampdists_E[[i]], J))
      
      mds_E[[i]] <- apply(playeri_E, 1, mahalanobis, center = ypred_E, cov = sigmaypred_E)
    }
    
    #Incorrect DM
    mds_F <- vector("list", length(which(singulars_F == 0)))
    for(i in 1:length(which(singulars_F == 0))){
      yobs_F <- ly_true_F[i,]
      ypred_F <- colMeans(newinitialY(sampdists_F[[i]], J))
      sigmaypred_F <- cov(newinitialY(sampdists_F[[i]], J))
      
      playeri_F <- rbind(yobs_F, newinitialY(sampdists_F[[i]], J))
      
      mds_F[[i]] <- apply(playeri_F, 1, mahalanobis, center = ypred_F, cov = sigmaypred_F)
    }
    
    #Correct MLN+
    mds_G <- vector("list", length(which(singulars_G == 0)))
    for(i in 1:length(which(singulars_G == 0))){
      yobs_G <- ly_true_G[i,]
      ypred_G <- colMeans(newinitialY(sampdists_G[[i]], J))
      sigmaypred_G <- cov(newinitialY(sampdists_G[[i]], J))

      playeri_G <- rbind(yobs_G, newinitialY(sampdists_G[[i]], J))

      mds_G[[i]] <- apply(playeri_G, 1, mahalanobis, center = ypred_G, cov = sigmaypred_G)
    }
    
    #Incorrect MLN+
    mds_H <- vector("list", length(which(singulars_H == 0)))
    for(i in 1:length(which(singulars_H == 0))){
      yobs_H <- ly_true_H[i,]
      ypred_H <- colMeans(newinitialY(sampdists_H[[i]], J))
      sigmaypred_H <- cov(newinitialY(sampdists_H[[i]], J))
      
      playeri_H <- rbind(yobs_H, newinitialY(sampdists_H[[i]], J))
      
      mds_H[[i]] <- apply(playeri_H, 1, mahalanobis, center = ypred_H, cov = sigmaypred_H)
    }
    
    
  }
  
  #Jittering within the empirical distribution of the Squared Mahalanobis Distances to get the single "uniform" residual
  {
    #Correct MLOGIT
#    mdps_A <- double(length = length(which(singulars_A == 0)))
#    for(i in 1:length(which(singulars_A == 0))){
#      pctl <- ecdf(mds_A[[i]][-1])
      
#      if(mds_A[[i]][1] <= min(mds_A[[i]][-1])){
#        minpct <- 0
#        maxpct <- pctl(min(mds_A[[i]][-1]))
#        if(maxpct == 0){
#          L <- length(mds_A[[i]][-1])
#          maxpct <- (L-1)/L
#        }
#      } else{
#        if(mds_A[[i]][1] > max(mds_A[[i]][-1])){
#          minpct <- pctl(max(mds_A[[i]][-1]))
#          maxpct <- 1
#          if(minpct == 1){
#            L <- length(mds_A[[i]][-1])
#            minpct <- (L-1)/L
#          }
#        } else{
#          minpct <- pctl(sort(mds_A[[i]][-1])[max(which(sort(mds_A[[i]][-1]) < mds_A[[i]][1]))])
#          maxpct <- pctl(mds_A[[i]][1])
#          if(minpct == 1){
#            L <- length(mds_A[[i]][-1])
#            minpct <- (L-1)/L
#          }
#          if(maxpct == 0){
#            L <- length(mds_A[[i]][-1])
#            maxpct <- (L-1)/L
#          }
#        }
#      }
#      
#      mdps_A[i] <- runif(1, minpct, maxpct)
#    }
    
    #Incorrect MLOGIT
    # mdps_B <- double(length = length(which(singulars_B == 0)))
    # for(i in 1:length(which(singulars_B == 0))){
    #   pctl <- ecdf(mds_B[[i]][-1])
    #   
    #   if(mds_B[[i]][1] <= min(mds_B[[i]][-1])){
    #     minpct <- 0
    #     maxpct <- pctl(min(mds_B[[i]][-1]))
    #     if(maxpct == 0){
    #       L <- length(mds_B[[i]][-1])
    #       maxpct <- (L-1)/L
    #     }
    #   } else{
    #     if(mds_B[[i]][1] > max(mds_B[[i]][-1])){
    #       minpct <- pctl(max(mds_B[[i]][-1]))
    #       maxpct <- 1
    #       if(minpct == 1){
    #         L <- length(mds_B[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #     } else{
    #       minpct <- pctl(sort(mds_B[[i]][-1])[max(which(sort(mds_B[[i]][-1]) < mds_B[[i]][1]))])
    #       maxpct <- pctl(mds_B[[i]][1])
    #       if(minpct == 1){
    #         L <- length(mds_B[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #       if(maxpct == 0){
    #         L <- length(mds_B[[i]][-1])
    #         maxpct <- (L-1)/L
    #       }
    #     }
    #   }
    #   
    #   mdps_B[i] <- runif(1, minpct, maxpct)
    # }
    
    #Correct MLN
    # mdps_C <- double(length = length(which(singulars_C == 0)))
    # for(i in 1:length(which(singulars_C == 0))){
    #   pctl <- ecdf(mds_C[[i]][-1])
    #   
    #   if(mds_C[[i]][1] <= min(mds_C[[i]][-1])){
    #     minpct <- 0
    #     maxpct <- pctl(min(mds_C[[i]][-1]))
    #     if(maxpct == 0){
    #       L <- length(mds_C[[i]][-1])
    #       maxpct <- (L-1)/L
    #     }
    #   } else{
    #     if(mds_C[[i]][1] > max(mds_C[[i]][-1])){
    #       minpct <- pctl(max(mds_C[[i]][-1]))
    #       maxpct <- 1
    #       if(minpct == 1){
    #         L <- length(mds_C[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #     } else{
    #       minpct <- pctl(sort(mds_C[[i]][-1])[max(which(sort(mds_C[[i]][-1]) < mds_C[[i]][1]))])
    #       maxpct <- pctl(mds_C[[i]][1])
    #       if(minpct == 1){
    #         L <- length(mds_C[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #       if(maxpct == 0){
    #         L <- length(mds_C[[i]][-1])
    #         maxpct <- (L-1)/L
    #       }
    #     }
    #   }
    #   
    #   mdps_C[i] <- runif(1, minpct, maxpct)
    # }
    
    #Incorrect MLN
    # mdps_D <- double(length = length(which(singulars_D == 0)))
    # for(i in 1:length(which(singulars_D == 0))){
    #   pctl <- ecdf(mds_D[[i]][-1])
    #   
    #   if(mds_D[[i]][1] <= min(mds_D[[i]][-1])){
    #     minpct <- 0
    #     maxpct <- pctl(min(mds_D[[i]][-1]))
    #     if(maxpct == 0){
    #       L <- length(mds_D[[i]][-1])
    #       maxpct <- (L-1)/L
    #     }
    #   } else{
    #     if(mds_D[[i]][1] > max(mds_D[[i]][-1])){
    #       minpct <- pctl(max(mds_D[[i]][-1]))
    #       maxpct <- 1
    #       if(minpct == 1){
    #         L <- length(mds_D[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #     } else{
    #       minpct <- pctl(sort(mds_D[[i]][-1])[max(which(sort(mds_D[[i]][-1]) < mds_D[[i]][1]))])
    #       maxpct <- pctl(mds_D[[i]][1])
    #       if(minpct == 1){
    #         L <- length(mds_D[[i]][-1])
    #         minpct <- (L-1)/L
    #       }
    #       if(maxpct == 0){
    #         L <- length(mds_D[[i]][-1])
    #         maxpct <- (L-1)/L
    #       }
    #     }
    #   }
    #   
    #   mdps_D[i] <- runif(1, minpct, maxpct)
    # }
    
    #Correct DM
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
    
    #Incorrect DM
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
    
    #Correct MLN+
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
    
    #Incorrect MLN+
    mdps_H <- double(length = length(which(singulars_H == 0)))
    for(i in 1:length(which(singulars_H == 0))){
      pctl <- ecdf(mds_H[[i]][-1])
      
      if(mds_H[[i]][1] <= min(mds_H[[i]][-1])){
        minpct <- 0
        maxpct <- pctl(min(mds_H[[i]][-1]))
        if(maxpct == 0){
          L <- length(mds_H[[i]][-1])
          maxpct <- (L-1)/L
        }
      } else{
        if(mds_H[[i]][1] > max(mds_H[[i]][-1])){
          minpct <- pctl(max(mds_H[[i]][-1]))
          maxpct <- 1
          if(minpct == 1){
            L <- length(mds_H[[i]][-1])
            minpct <- (L-1)/L
          }
        } else{
          minpct <- pctl(sort(mds_H[[i]][-1])[max(which(sort(mds_H[[i]][-1]) < mds_H[[i]][1]))])
          maxpct <- pctl(mds_H[[i]][1])
          if(minpct == 1){
            L <- length(mds_H[[i]][-1])
            minpct <- (L-1)/L
          }
          if(maxpct == 0){
            L <- length(mds_H[[i]][-1])
            maxpct <- (L-1)/L
          }
        }
      }
      
      mdps_H[i] <- runif(1, minpct, maxpct)
    }
    
  }
  
  #Convert to Standard Normal
  {
    #Correct MLOGIT
#    mdrs_A <- qnorm(mdps_A)
    
    #Incorrect MLOGIT
#    mdrs_B <- qnorm(mdps_B)
    
    #Correct MLN
#    mdrs_C <- qnorm(mdps_C)
    
    #Incorrect MLN
#    mdrs_D <- qnorm(mdps_D)
    
    #Correct DM
    mdrs_E <- qnorm(mdps_E)
    
    #Incorrect DM
    mdrs_F <- qnorm(mdps_F)
    
    #Correct MLN+
    mdrs_G <- qnorm(mdps_G)
    
    #Incorrect MLN+
    mdrs_H <- qnorm(mdps_H)
  }
  
  options(warn = - 1)
  #Kolmogorov-Smirnov Test (to compare it to standard normal)
  #And Deviance for the MLOGIT models, but not the MLN models (since can't calculate deviance of correctly fit model)
  {
    #MLOGIT models
#    kA <- ks.test(mdrs_A, "pnorm", mean = 0, sd = 1)
#    kB <- ks.test(mdrs_B, "pnorm", mean = 0, sd = 1)
    
    #MLN models
#    kC <- ks.test(mdrs_C, "pnorm", mean = 0, sd = 1)
#    kD <- ks.test(mdrs_D, "pnorm", mean = 0, sd = 1)
    
    #DM models
    kE <- ks.test(mdrs_E, "pnorm", mean = 0, sd = 1)
    kF <- ks.test(mdrs_F, "pnorm", mean = 0, sd = 1)
    
    #MLN+ models
    kG <- ks.test(mdrs_G, "pnorm", mean = 0, sd = 1)
    kH <- ks.test(mdrs_H, "pnorm", mean = 0, sd = 1)
  }
  options(warn = 0)
  
  
#  c(kA$p.value, kB$p.value, kC$p.value, kD$p.value, kE$p.value, kF$p.value, kG$p.value, kH$p.value, kA$statistic, kB$statistic, kC$statistic, kD$statistic, kE$statistic, kF$statistic, kG$statistic, kH$statistic)
  c(kE$p.value, kF$p.value, kG$p.value, kH$p.value, kE$statistic, kF$statistic, kG$statistic, kH$statistic)
}

end.simstudy <- Sys.time()

end.simstudy - start.simstudy

head(ks.log)
#colMeans(ks.log > .05)
colMeans(ks.log[,1:4] > .05)
mean(ks.log[,5] > ks.log[,6])
mean(ks.log[,7] > ks.log[,8])
