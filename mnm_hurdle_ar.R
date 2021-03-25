library(R2jags)


simulateData_hurdle_ar<-function(S,R,T,K, prob, abundance, zeros){
  Y<-array(NA, dim=c(R,T,S,K), dimnames=list(NULL, NULL, paste("species", 1:S), NULL))
  N<-lambda<-array(NA, dim=c(R,S,K), dimnames=list(NULL, paste("species", 1:S), NULL))    

  
  # Zero-counts
  if(zeros=="small") theta <- 0.2
  else theta <- 0.7
  zeroInflation <- array(rbinom(R*S*K, size=1, prob=1-theta), dim=c(R,S,K))
  
  
  # Probability of detection
  if(prob=="small") p <- runif(S, 0.1, 0.4)
  else p <- runif(S,0.5,0.9)
  
  
  # Uniform dist. autocorrelation coefficient
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
      N[i,s,1]<-ifelse(zeroInflation[i,s,1]==0,
                       0,
                       extraDistr::rtpois(n=1, lambda=lambda[i,s,1],a=0))
      
      for(k in 2:K){
        lambda[i,s,k]<-exp(a[i,s]+phi[s]*log(N[i,s,k-1]+1))        
        N[i,s,k]<-ifelse(zeroInflation[i,s,k]==0,
                         0,
                         extraDistr::rtpois(n=1, lambda=lambda[i,s,k],a=0))
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
              "sigma"=sigma,
              "theta"=theta,
              "zeros"=zeros,
              "lambda"=lambda, 
              "a"=a, 
              "alpha"=alpha, 
              "correlation"=correlation,
              "phi"=phi, 
              "muPhi"=muPhi,
              "varPhi"=sigmaPhi)
  return(ylist)
}


modelRun_hurdle_ar<-function(R,T,S,K,Y){
  cat(  '
model {
  for(s in 1:S){
      logit(probability[s])<-gamma[s]
      
      
    for (i in 1:R) { 
      x[i,s,1] ~ dbern(1-theta)
      log(lambda[i,s,1])<-a[i,s]
      count[i,s,1] ~ dpois(lambda[i,s,1]+1E-10)T(1,)  
      N[i,s,1]<-ifelse(x[i,s,1]==0, 0, count[i,s,1])

      for(k in 2:K){
        x[i,s,k] ~ dbern(1-theta)
        log(lambda[i,s,k]) <- a[i,s] + phi[s]*log(N[i,s,k-1]+1)
        count[i,s,k] ~ dpois(lambda[i,s,k]+1E-1)T(1,)  
        N[i,s,k]<-ifelse(x[i,s,k]==0, 0, count[i,s,k])
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
    phi[s] ~ dnorm(muPhi,tauPhi)
    gamma[s]~dnorm(0, 0.001)

  }
  
  theta~dbeta(1,1)

  muPhi~dnorm(0,0.001)
  varPhi~dunif(0,100)
  tauPhi<-pow(varPhi, -2)
}
', file={model_code <- tempfile()})
  
  # Initial Values for N - need to be larger than the largest Y value for that species at that site
  count.inits<-apply(Y, c(1,3,4), max)+1
  x.inits<-apply(Y, c(1,3,4),function(z) ifelse(any(z)>0, 1, 0))
  
  initfunction <- function(chain) {
    return(switch(chain,
                  "1" = list("count"=count.inits, 
                             "x"=x.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 1),
                  "2" = list("count"=count.inits,
                             "alpha"=rnorm(S),
                             "x"=x.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 2),
                  "3" = list("count"=count.inits, 
                             "alpha"=rnorm(S), 
                             "x"=x.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 3),
                  "4" = list("count"=count.inits,
                             "alpha"=rnorm(S),
                             "x"=x.inits,
                             .RNG.name = "base::Super-Duper",
                             .RNG.seed = 4)))
  }
  # Choose the parameters to watch - 
  model_parameters=c("alpha", "cor", "covariance", "N", "sigma","probability", "phi", "muPhi", "varPhi", "theta", "x")
  
  model_data = list("R" = R, "Y"=Y, "T"=T, "S"=S, "K"=K,"Omega"=diag(1, S),"count.inits"=count.inits, "x.inits"=x.inits)
  model_run=NULL
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
                "estimatedVariance"=model_run$BUGSoutput$mean$sigma,
                "estimatedTheta"=model_run$BUGSoutput$mean$theta)
  
  return(modList)
}


# Wrapper function that calls on function to simulate data and function to fit model in JAGS
mnm_hurdle_ar<-function(combinations, ncombinations, ndatasets){
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
    zeros<-combinations[i,'zeros']
    
    
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               simulateData_hurdle_ar(R=R,
                                               T=T,
                                               S=S,
                                               K=K,
                                               prob=prob,
                                               abundance=abundance,
                                               zeros=zeros))
    }
    
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               modelRun_hurdle_ar(R=R,
                                        T=T,
                                        S=S,
                                        K=K,
                                        Y=data1[[i]][[dataset]][["Y"]]))
      print(paste("dataset (", i, ",", dataset, ") complete", sep=""))
    }
    
  }
  return(data1)
}

# Choose number of datasets, sites, sampling occasions, species and years,along with small or large probability of detection, abundance and probability of obtaining zero
ndatasets<-10
R<-100
T<-10
S<-10
K<-5
prob<-"small"
abundance<-"large"
zeros<-"small"
combinations<-expand.grid("R"=R,"T"=T,"S"=S,"K"=K,"probability"=prob, "abundance"=abundance, "zeros"=zeros)
ncombinations<-dim(combinations)[1]
x<-mnm_hurdle_ar(ndatasets=ndatasets, ncombinations=ncombinations, combinations=combinations)

