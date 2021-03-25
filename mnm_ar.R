library(R2jags)
# Function to simulate data for the MNM-AR Model
# R = number of sites
# T = number of sampling occasions
# S = number of species
# k = number of years

# prob = S-dimensional vector containing probability of detection - can be either "small" or "large"
# prob = "small" generates probability values in the range 0.1-0.4. Otherwise p values are in the range 0.5-0.9

# abundance = size of latent species abundance - can be either "small" or "large"
# abundance = "small" generates lambda values in the approx. range 0-50. Otherwise lambda values are approx. 0-700



simulateData_ar<-function(S,R,T,K, prob, abundance){
  Y<-array(NA, dim=c(R,T,S,K), dimnames=list(NULL, NULL, paste("species", 1:S), NULL))
  N<-lambda<-array(NA, dim=c(R,S,K), dimnames=list(NULL, paste("species", 1:S), NULL))    
  gamma<-vector(length=S)
  
  
  # Probability of detection
  set.seed(9875)
  if(prob=="small") p <- runif(S, 0.1, 0.4)
  else p <- runif(S,0.5,0.9)
  
  # Normally dist. autocorrelation coefficient
  muPhi<-0
  sigmaPhi<-0.2
  set.seed(9875)
  phi<-rnorm(S, muPhi, sigmaPhi)
  
  
  # MVN random effect
  Omega<-diag(1, nrow=S, ncol=S) # Scale matrix for wishart distribution
  set.seed(9875)
  Sigma<-clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation<-cov2cor(Sigma)
  alpha<-rep(ifelse(abundance=="small", 2, 4),S)

  # Standard Deviations
  sigma<-vector(length=S)
  for(i in 1:S){ 
    sigma[i]<-sqrt(Sigma[i,i])
  }
  
  
  set.seed(NULL) 
  a <- mvtnorm::rmvnorm(R, mean=alpha, sigma=Sigma)
  
  
  for(i in 1:R){
    for(s in 1:S){
      lambda[i,s,1]<-exp(a[i,s])        
      N[i,s,1]<-rpois(1, lambda[i,s,1])
      
      for(k in 2:K){
        lambda[i,s,k]<-exp(a[i,s]+phi[s]*log(N[i,s,k-1]+1))        
        N[i,s,k]<-rpois(1, lambda[i,s,k])
      }
    }
  }
  
  for(i in 1:R){
    for(t in 1:T){
      for(s in 1:S){
        for(k in 1:K){
          Y[i,t,s,k]<-rbinom(1, size=N[i,s,k], prob=p[s])
        }
      }
    }
  }
  
  ylist<-list("Y"=Y, 
              "N"=N, 
              "R"=R, 
              "T"=T, 
              "S"=S,
              "K"=K,
              "p"=p, 
              "Sigma"=Sigma, 
              "lambda"=lambda, 
              "a"=a, 
              "alpha"=alpha, 
              "correlation"=correlation,
              "phi"=phi, 
              "muPhi"=muPhi,
              "varPhi"=sigmaPhi)
  return(ylist)
}


modelRun_ar<-function(R,T,S,K,Y){
  cat(  '
model {
  for(s in 1:S){
      logit(probability[s])<-gamma[s]
      
      
    for (i in 1:R) { 
      log(lambda[i,s,1])<-a[i,s]
      N[i,s,1]~dpois(lambda[i,s,1])
      
      for(k in 2:K){
        log(lambda[i,s,k]) <- a[i,s] + phi[s]*log(N[i,s,k-1]+1)
         N[i,s,k] ~ dpois(lambda[i,s,k])
      }
    }
  }

        # Loop over time points
        for(i in 1:R){
          for(s in 1:S){
            for(t in 1:T){   
              for(k in 1:K){
                Y[i,t,s,k] ~ dbin(probability[s], N[i,s,k])
              }
            }
        }
     }
   
  # Random effects a with mean vector mu, and variance-covariance matrix cov
  for(i in 1:R){
    a[i,1:S] ~ dmnorm(alpha[], precision[,])
  }
  

  # Wishart prior on precision with df=S+1 and diagonal matrix Omega
  df<-S+1
  precision[1:S,1:S] ~ dwish(Omega[,], df)
  covariance[1:S,1:S]<-inverse(precision[,])


 # Correlations and Standard deviations
  for (s in 1:S){    
    sigma[s] <- sqrt(covariance[s,s])
    for (s1 in 1:S){
      cor[s,s1]<-covariance[s,s1]/sqrt(covariance[s,s]*covariance[s1, s1])
    }
  }


  # Priors
  for(s in 1:S){  
    alpha[s] ~ dnorm(0.001,0.001)
    gamma[s] ~ dnorm(0, 0.001)
    phi[s] ~ dnorm(muPhi,tauPhi)
  }
    
  muPhi~dnorm(0,0.001)
  varPhi~dunif(0,100)
  tauPhi<-pow(varPhi, -2)
}
', file={model_code <- tempfile()})
  
  # Initial Values for N - need to be larger than Y
  count.inits<-apply(Y, c(1,3), max)+1
  
  initfunction <- function(chain) {
    return(switch(chain,
                  "1" = list("N"=count.inits, 
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 1),
                  "2" = list("N"=count.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 2),
                  "3" = list("N"=count.inits, 
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 3),
                  "4" = list("N"=count.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 4)))
  }
  # Choose the parameters to watch - 
  model_parameters=c("alpha", "cor", "covariance", "N", "sigma","probability", "phi", "muPhi", "varPhi")
  
  model_data = list("R" = R, "Y"=Y, "T"=T, "S"=S, "K"=K,"Omega"=diag(1, S),"count.inits"=count.inits)
  model_run=NULL
  # Run the model
  list2env(model_data , envir=globalenv())
  # Run the model
  model_run = jags.parallel(data = names(model_data),
                            inits=initfunction,
                            parameters.to.save = model_parameters,
                            model.file = model_code,
                            n.chains = 4,
                            n.iter = 50000,
                            n.burnin = 10000,
                            n.thin = 10)
  
  modList<-list("summary"=as.data.frame(model_run$BUGSoutput$summary),
                "estimatedCorrelation"=model_run$BUGSoutput$mean$cor,
                "estimatedN"=model_run$BUGSoutput$mean$N,
                "estimatedProbability"=model_run$BUGSoutput$mean$probability, 
                "estimatedCovariance"=model_run$BUGSoutput$mean$covariance,
                "estimatedAlpha"=model_run$BUGSoutput$mean$alpha,
                "estimatedphi"=model_run$BUGSoutput$mean$phi,
                "estimatedMuPhi"=model_run$BUGSoutput$mean$muPhi,
                "estimatedVarPhi"=model_run$BUGSoutput$mean$varPhi,
                "estimatedVariance"=model_run$BUGSoutput$mean$sigma)
  
  return(modList)
}

mnm_ar<-function(combinations, ncombinations, ndatasets){
  data1<-vector("list", length=ncombinations)
  
  for(i in 1:ncombinations){
    set.seed(i)
    data1[[i]]<-vector("list", length=ndatasets)
    R<-combinations[i,'R']
    T<-combinations[i,'T']
    S<-combinations[i,'S']
    K<-combinations[i,'K']
    prob<-combinations[i, 'probability']
    abundance<-combinations[i,'abundance']
    
    
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               simulateData_ar(R=R,
                                               T=T,
                                               S=S,
                                               K=K,
                                               prob=prob,
                                               abundance=abundance))
    }
    
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               modelRun_ar(R=R,
                                        T=T,
                                        S=S,
                                        K=K,
                                        Y=data1[[i]][[dataset]][["Y"]]))
      print(paste("dataset (", i, ",", dataset, ") complete", sep=""))
    }
    
  }
  return(data1)
}

# Choose number of datasets, sites, sampling occasions, species and years
ndatasets<-5
R<-100
T<-10
S<-10
K<-5
prob<-"large"
abundance<-"large"
combinations<-expand.grid("R"=R,"T"=T,"S"=S,"K"=K,"probability"=prob, "abundance"=abundance)
ncombinations<-dim(combinations)[1]
x<-mnm_ar(ndatasets=ndatasets, ncombinations=ncombinations, combinations=combinations)
