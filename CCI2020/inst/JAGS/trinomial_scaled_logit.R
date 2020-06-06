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
    beta.phi[1] ~ dnorm(0,.01)
    beta.phi[2] ~ dnorm(0,.01)

  ## Scale parameter
  lower <- 0
  upper ~ dunif(0,1)
  
  ## Capture probability
  p ~ dunif(0,1)

  ## Recovery probability
  lambda ~ dunif(0,1)
}

    

        
                  
