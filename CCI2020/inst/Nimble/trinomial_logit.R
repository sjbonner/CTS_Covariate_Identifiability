model{

  ##### Likelihood #####
  
  for(i in 1:nrelease){ ## i -- releases
    ## Survival probability
    logit(phi_tmp[i]) <- beta.phi[1] + beta.phi[2] * z[i]
    phi[i] <- lower + (upper - lower) *  phi_tmp[i]
    
    ## Trinomial probabilities
    Probs[i,2] <- phi[i] * p # Recapture
    Probs[i,3] <- (1 - phi[i]) * lambda # Recovery
    Probs[i,1] <- 1 - Probs[i,2] - Probs[i,3] # Neither
    
    ## Event
    outcome[i] ~ dcat(Probs[i,1:3])
  }
  
  ##### Priors #####
  ## Survival probability
  for(k in 1:2){
    beta.phi[k] ~ dt(beta.phi.hyper[k,1], beta.phi.hyper[k,2], beta.phi.hyper[k,3])
  }

  ## Scale parameter
  lower <- 0
  upper <- 1
  
  ## Capture probability
  p ~ dunif(p.hyper[1],p.hyper[2])

  ## Recovery probability
  lambda ~ dunif(lambda.hyper[1],lambda.hyper[2])
}

    

        
                  
