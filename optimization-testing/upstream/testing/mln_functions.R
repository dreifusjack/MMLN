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
            difr <- ifelse(abs(sexpY) == Inf | abs(sexpYn) == Inf, 0, difr)
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
            difr <- ifelse(abs(sexpY) == Inf | abs(sexpYn) == Inf, 0, difr)
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
        
        #This if statement can be moved around for troubleshooting at different points
        if(test == TRUE){
          if(i == 14){
            browser()         
          }
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
