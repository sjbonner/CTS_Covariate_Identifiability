model{

    ##### Likelihood #####

    for(i in 1:nrelease){ ## i -- releases
        for(j in 1:drecap[i]){
            ## Impute covariate values
            mu.z.tmp[i,j] <- Z[i,j] + mu.z[release[i]+j-1]
            Z[i,j+1] ~ dnorm(mu.z.tmp[i,j],tau.z)
        
            ## Survival probabilities
            logit(phi_tmp[i,j]) <- beta.phi[1] + beta.phi[2] * Z[i,j]
            phi[i,j] <- lower + (upper - lower) *  phi_tmp[i,j]

            ## Capture probabilities
            q[i,j+1] <- 1-p
        }

        q[i,1] <- 1
    }

    ## Likelihood contributions
    ## 1) Recaptured individuals
    for(i in 1:nrecapture){
        Prob[i] <- prod(phi[i,1:drecap[i]]) * # Survival 
            prod(q[i,1:drecap[i]]) *          # Missed captures
                p                             # Recapture

        dummy[i] ~ dbern(Prob[i])
    }

    ## 2) Not recaptured individuals
    for(i in 1:(nrelease-nrecapture)){

      Prob.tmp[i,1] <- phi[nrecapture+i,1] * p
        
        for(j in 2:drecap[nrecapture+i]){
            Prob.tmp[i,j] <-  Prob.tmp[i,j-1] / p *
                q[nrecapture+i,j]*phi[nrecapture+i,j] * p
        }
        
        Prob[nrecapture+i] <- 1-sum(Prob.tmp[i,1:drecap[nrecapture+i]])

        dummy[nrecapture+i] ~ dbern(Prob[nrecapture+i])
    }

    ##### Priors #####

    ## Covariate model
    for(t in 1:(nocc-1)){
        mu.z[t] ~ dnorm(0,.0001)
    }
    tau.z ~ dgamma(.001,.001)
    sigma.z <- 1/sqrt(tau.z)

  ## Survival probability
  for(k in 1:2){
    beta.phi[k] ~ dnorm(beta.phi.hyper[k,1], beta.phi.hyper[k,2])
  }

  ## Scale parameter
  lower <- 0
  upper ~ dunif(lower,1)

  ## Capture probability
  p ~ dunif(p.hyper[1],p.hyper[2])
}

    

        
                  
