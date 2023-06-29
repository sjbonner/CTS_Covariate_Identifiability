model{

    ##### Likelihood #####

    for(i in 1:nrelease){ ## i -- releases
      ## Survival probability
      logit(phi_tmp[i]) <- inprod(beta.phi[1:nZphi], Zphi[i,1:nZphi])
      
      phi[i] <- lower + (upper - lower) *  phi_tmp[i]

      ## Capture probability
      logit(p[i]) <- inprod(beta.p[1:nZp], Zp[i, 1:nZp])
      
      ## Binomial probabilities
      Probs[i,2] <- phi[i] * p[i] # Recapture
      Probs[i,1] <- 1 - Probs[i,2] # Not recapture

      ## Event
      outcome[i] ~ dcat(Probs[i,1:2])
    }

##### Priors #####
  ## Survival probability
  for(k in 1:nZphi){
    beta.phi[k] ~ dnorm(0,.01)
  }

  ## Scale parameter
  lower <- 0
  upper <- 1
  
  ## Capture probability
  for(k in 1:nZp){
    beta.p[k] ~ dnorm(0,.01)
  }
}

    

        
                  
