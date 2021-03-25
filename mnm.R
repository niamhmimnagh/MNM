# Function to simulate multi-species correlated abundance data
# R = number of sites
# T = number of sampling occasions
# S = number of species

# prob = S-dimensional vector containing probability of detection - can be either "small" or "large"
# prob = "small" generates probability values in the range 0.1-0.4. Otherwise p values are in the range 0.5-0.9

# abundance = size of latent species abundance - can be either "small" or "large"
# abundance = "small" generates lambda values in the approx. range 0-50. Otherwise lambda values are approx. 0-700

library(R2jags)

simulateData<-function(S,R,T, prob, abundance){
  N<-lambda<-matrix(ncol=S, nrow=R)
  Y <-  array(NA, dim = c(R,T,S), dimnames = list(NULL, NULL, paste("species", 1:S)))
  
  
  # Probability of detection
  set.seed(9875)
  if(prob=="small") p <- runif(S, 0.1, 0.4)
  else p <- runif(S,0.5,0.9)
  
  # Mean and Covariance matrix for MVN variable a
  set.seed(9875)
  Sigma<-clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation<-cov2cor(Sigma)
  alpha<-rep(ifelse(abundance=="small", 2, 4),S)  
  
  set.seed(NULL)
  a<-mvtnorm::rmvnorm(R, mean=alpha, sigma=Sigma)
  lambda<-exp(a)
  
  ## observed and true abundances
  for(i in 1:R){
    for(s in 1:S){
      N[i,s]<-rpois(1, lambda[i,s])
      for(t in 1:T){
        Y[i,t,s]<-rbinom(1, N[i,s], p[s])
      }
    }
  }
  
  ylist<-list("Y"=Y, 
              "N"=N, 
              "R"=R, 
              "T"=T, 
              "S"=S, 
              "p"=p, 
              "Sigma"=Sigma, 
              "lambda"=lambda, 
              "a"=a, 
              "alpha"=alpha, 
              "correlation"=correlation)
  return(ylist)
}


# Function to fit MNM model to data simulated in above function
modelRun<-function(R,T,S,Y, correlation, N, p){
  cat('
model {
  # Likelihood
  for(s in 1:S){
    logit(probability[s])<-gamma[s]
 
    for (i in 1:R) {
      N[i,s] ~ dpois(lambda[i,s])
      log(lambda[i,s]) <- a[i,s]

      for(t in 1:T){
      Y[i,t,s] ~ dbin(probability[s],N[i,s])
      }
    }
  }
   
  # Random effects a with mean vector mu, and variance-covariance matrix cov
  for(i in 1:R){
    a[i,1:S] ~ dmnorm(alpha, precision[,])
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
    gamma[s] ~ dnorm(0, 0.001)    
    alpha[s] ~ dnorm(0, 0.001)
  }
}
', file={model_code <- tempfile()})
  
  # Initial Values
  # Initial values for N must be greater than Y
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
  model_parameters =  c("alpha", "cor", "sigma", "N", "probability", "analyticCorrelation")
  model_data = list("R" = R, "Y"=Y, "T"=T, "S"=S, "Omega"=diag(S), "count.inits"=count.inits)
  model_run=NULL
  model_runtime<-system.time({model_run =jags.parallel(data = names(model_data),
                                                       inits=initfunction,
                                                       parameters.to.save = model_parameters,
                                                       model.file = model_code,
                                                       n.chains = 4,
                                                       n.iter = 50000,
                                                       n.burnin = 10000,
                                                       n.thin = 10,
                                                       n.cluster=4)})
  
  
  
  modList<-list("summary"=as.data.frame(model_run$BUGSoutput$summary),
                "estimatedCorrelation"=model_run$BUGSoutput$mean$cor,
                "estimatedN"=model_run$BUGSoutput$mean$N,
                "analyticCorrelations"=model_run$BUGSoutput$mean$analyticCorrelation,
                "estimatedProbability"=model_run$BUGSoutput$mean$probability, 
                "estimatedVariance"=model_run$BUGSoutput$mean$sigma, 
                "estimatedAlpha"=model_run$BUGSoutput$mean$alpha)
  return(modList)
}


# Wrapper function which simulates data and fits model for MNM model
# combinations = dataframe containing 1 row per combination of parameters and 1 column per parameter
# ncombinations = number of rows in combinations dataframe
# ndatasets = number of datasets to simulate


mnm<-function(combinations, ncombinations, ndatasets){
  data1<-vector("list", length=ncombinations)
  
  for(i in 1:ncombinations){
    data1[[i]]<-vector("list", length=ndatasets)
    R<-combinations[i,'R']
    T<-combinations[i,'T']
    S<-combinations[i,'S']
    prob<-combinations[i, 'probability']
    abundance<-combinations[i,'abundance']
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               simulateData(R=R,
                                            T=T,
                                            S=S,
                                            prob=prob,
                                            abundance=abundance))
    }
    
    
    
    for(dataset in 1:ndatasets){
      data1[[i]][[dataset]]<-c(data1[[i]][[dataset]],
                               modelRun(R=R,
                                        T=T,
                                        S=S,
                                        Y=data1[[i]][[dataset]][["Y"]],
                                        N=data1[[i]][[dataset]][["N"]]))
      print(paste("dataset (", i, ",", dataset, ") complete", sep=""))
    }
    
  }
  return(data1)
}


# Choose number of datasets, sites, sampling occasions, species and run the model
ndatasets<-20
R<-100
T<-10
S<-10
prob<-"large"
abundance<-"small"
combinations<-expand.grid("R"=R,"T"=T,"S"=S,"probability"=prob, "abundance"=abundance)
ncombinations<-dim(combinations)[1]
x<-mnm(ndatasets=ndatasets, ncombinations=ncombinations, combinations=combinations)
