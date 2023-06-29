model{

    ##### Likelihood #####

    for(i in 1:nrelease){ ## i -- releases
      ## Survival probability
      logit(phi_tmp[i]) <- beta.phi[1] + beta.phi[2] * z[i]
      phi[i] <- lower + (upper - lower) *  phi_tmp[i]
      
      ## Binomial probabilities
      Probs[i,2] <- phi[i] * p # Recapture
      Probs[i,1] <- 1 - Probs[i,2] # Not recapture

      ## Event
      outcome[i] ~ dcat(Probs[i,1:2])
    }

    ##### Priors #####
  ## Survival probability
  for(k in 1:2){
    beta.phi[k] ~ dt(beta.phi.hyper[k,1], beta.phi.hyper[k,2], beta.phi.hyper[k,3])
  }

  ## Scale parameter
  lower ~ dunif(0,1)
  upper ~ dunif(lower,1)
  
  ## Capture probability
  p ~ dunif(p.hyper[1],p.hyper[2])
}

    

        
                  
