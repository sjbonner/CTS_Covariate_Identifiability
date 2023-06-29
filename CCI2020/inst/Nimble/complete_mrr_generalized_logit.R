model{

    ##### Likelihood #####

    for(i in 1:nrelease){ ## i -- releases
        for(j in 1:delta[i]){
            ## Impute covariate values
            mu.z.tmp[i,j] <- Z[i,j] + mu.z[release[i]+j-1]
            Z[i,j+1] ~ dnorm(mu.z.tmp[i,j],tau.z)
        
            ## Survival probabilities
            logit(phi_tmp[i,j]) <- beta.phi[1] + beta.phi[2] * Z[i,j]
            phi[i,j] <- lower + (upper - lower) *  phi_tmp[i,j]
       }
    }

  ## Likelihood contributions
  ## 1) Recaptures
  for(i in 1:nrecapture){
    Prob[i] <- drecap(delta[i], phi[i,1:delta[i]], p)
  }
  
  ## 2) Recoveries
  for(i in 1:nrecovery){
    Prob[nrecapture + i] <- drecov(delta[nrecapture + i], phi[nrecapture + i,1:delta[nrecapture + i]], p, lambda)
  }
  
  ## 3) Evasions
  for(i in 1:(nrelease - nrecapture - nrecovery)){
    Prob[nrecapture + nrecovery + i] <- devade(delta[nrecapture + nrecovery + i], phi[nrecapture + nrecovery + i,1:delta[nrecapture + nrecovery + i]], p, lambda)
  }
  
  ## Likelihood
  for(i in 1:nrelease){
    dummy[i] ~ dbern(Prob[i])
  }
   ##### Priors #####

    ## Covariate model
    for(t in 1:(nocc-1)){
        mu.z[t] ~ dnorm(mu.z.hyper[1], mu.z.hyper[2])
    }

    sigma.z ~ T(dt(0, sigma.z.hyper[1], sigma.z.hyper[2]), 0,)
    tau.z <- 1/(sigma.z * sigma.z)

  ## Survival probability
  for(k in 1:2){
    beta.phi[k] ~ dt(beta.phi.hyper[k,1], beta.phi.hyper[k,2], beta.phi.hyper[k,3])
  }

  ## Scale parameter
  lower ~ dunif(0,1)
  upper ~ dunif(lower,1)

  ## Recovery probability
  lambda ~ dunif(lambda.hyper[1],lambda.hyper[2])
  
  ## Capture probability
  p ~ dunif(p.hyper[1],p.hyper[2])
}

    

        
                  
