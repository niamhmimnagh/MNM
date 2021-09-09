# This file provides R code for the analytic correlations derived in Appendix A
# Author: Niamh Mimnagh

## MNM Model
analyticCorrelationMNM<-function(covariance_matrix, mu){
  analyticCovariance<-matrix(nrow=nrow(covariance_matrix), ncol=ncol(covariance_matrix))
  S<-nrow(covariance_matrix)

for(s in 1:(S-1)){
  for(s.prime in (s+1):S){
    analyticCovariance[s,s.prime]<-exp(mu[s]+0.5*(covariance_matrix[s,s]))*exp(mu[s.prime]+0.5*(covariance_matrix[s.prime, s.prime]))*(exp(covariance_matrix[s,s.prime])-1)
    analyticCovariance[s.prime, s]<-analyticCovariance[s,s.prime]
    }
}

for(s in 1:S){
  analyticCovariance[s,s]<-exp(mu[s]+0.5*(covariance_matrix[s,s]))*exp(mu[s]+0.5*(covariance_matrix[s, s]))*(exp(covariance_matrix[s,s])-1)
}

analyticVariance<-vector()
for(s in 1:S){
  analyticVariance[s]<-analyticCovariance[s,s]
}

analyticCorrelation<-matrix(nrow=S, ncol=S)
for(s in 1:(S-1)){
  for(s.prime in (s+1):S){
    analyticCorrelation[s,s.prime]<-analyticCovariance[s,s.prime]/(sqrt(analyticVariance[s])*sqrt(analyticVariance[s.prime]))
    analyticCorrelation[s.prime, s]<-analyticCorrelation[s,s.prime]
  }
}

for(s in 1:S){
  analyticCorrelation[s,s]<-analyticCovariance[s,s]/(sqrt(analyticVariance[s])*sqrt(analyticVariance[s]))
}
return(analyticCorrelation)
}


## AR Model
analyticCorrelationMNMAR<-function(covariance_matrix_a, mu_a, mu_phi, var_phi, N){
  mu_lambda<-array(dim=c(R,S,K))
  for(i in 1:R){
    for(s in 1:S){
      mu_lambda[i,s,1]<-mu_a[s]
      for(k in 2:K){
        mu_lambda[i,s,k]<-mu_a[s]+log(N[i,s,k-1]+1)*mu_phi[s]
      }
    }
  }
  
  Sigma_lambda<-array(dim=c(S,S,R,K))
  for(i in 1:R){
    Sigma_lambda[,,i,]<-covariance_matrix_a
  }
  
  estVar<-diag(var_phi,S)
  
  for(i in 1:R){
    for(k in 2:K){
      Sigma_lambda[,,i,k]<-Sigma_lambda[,,i,k]+log(N[i,,k-1]+1)^2*estVar[,]
    }
  }
  
  analyticCovariance<-analyticCorrelation<-array(dim=c(S,S,R,K))
  for(i in 1:R){
    for(s in 1:(S-1)){
      for(s.prime in (s+1):S){
        for(k in 1:K){
          analyticCovariance[s,s.prime,i,k]<-exp(mu_lambda[i,s,k]+0.5*(Sigma_lambda[s,s,i,k]))*exp(mu_lambda[i,s.prime,k]+0.5*(Sigma_lambda[s.prime, s.prime,i,k]))*(exp(Sigma_lambda[s,s.prime,i,k])-1)
          analyticCovariance[s.prime, s,i,k]<-analyticCovariance[s,s.prime,i,k]
        }
      }
    }
  }
  for(i in 1:R){
    for(s in 1:S){
      for(k in 1:K){
        analyticCovariance[s,s,i,k]<-exp(mu_lambda[i,s,k]+0.5*(Sigma_lambda[s,s,i,k]))*exp(mu_lambda[i,s,k]+0.5*(Sigma_lambda[s, s,i,k]))*(exp(Sigma_lambda[s,s,i,k])-1)
      }
    }
  }
  
  analyticVariance<-array(dim=c(S,R,K))
  for(i in 1:R){
    for(s in 1:S){
      for(k in 1:K){
        analyticVariance[s,i,k]<-analyticCovariance[s,s,i,k]
      }
    }
  }
  for(i in 1:R){
    for(s in 1:(S-1)){
      for(s.prime in (s+1):S){
        for(k in 1:K){
          analyticCorrelation[s,s.prime,i,k]<-analyticCovariance[s,s.prime,i,k]/(sqrt(analyticVariance[s,i,k])*sqrt(analyticVariance[s.prime,i,k]))
          analyticCorrelation[s.prime, s,i,k]<-analyticCorrelation[s,s.prime,i,k]
        }
      }
    }
  }
  
  for(i in 1:R){
    for(s in 1:S){
      for(k in 1:K){
        analyticCorrelation[s,s,i,k]<-analyticCovariance[s,s,i,k]/(sqrt(analyticVariance[s,i,k])*sqrt(analyticVariance[s,i,k]))
      }
    }
  }
  return(analyticCorrelation)
}

## Hurdle Model
analyticCorrelationHurdle<-function(theta, mu_a, covariance_matrix_a){
  analyticCov<-analyticCorrelation<-Cov_lambda<-matrix(ncol=S, nrow=S)
  analyticVariance<-vector(length=S)
  
  mu_l<-exp(mu_a+0.5*diag(covariance_matrix_a))
  sd_l<-exp(2*mu_a+diag(covariance_matrix_a))*(exp(diag(covariance_matrix_a))-1)
  
  for(s in 1:(S-1)){
    for(s.prime in (s+1):S){
      Cov_lambda[s,s.prime]=exp(mu_a[s]+mu_a[s.prime]+0.5*(covariance_matrix_a[s,s]+covariance_matrix_a[s.prime,s.prime]))*(exp(covariance_matrix_a[s,s.prime])-1)
      Cov_lambda[s.prime,s]<-Cov_lambda[s,s.prime]
      }
    }
    
  for(s in 1:S){
    Cov_lambda[s,s]<-sd_l[s]
    }

  for(s in 1:(S-1)){
    for(s.prime in (s+1):S){
       analyticCov[s,s.prime]<-(((1-theta)*mu_l[s])/(1-exp(-mu_l[s])))*(((1-theta)*mu_l[s.prime])/(1-exp(-mu_l[s.prime])))+
        0.5*sd_l[s]*(((exp(-mu_l[s])*mu_l[s.prime])*(1-theta)^2*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2))/((1-exp(-mu_l[s.prime]))*(1-exp(-mu_l[s]))^3))+
        0.5*sd_l[s.prime]*(((exp(-mu_l[s.prime])*mu_l[s])*(1-theta)^2*(exp(-mu_l[s.prime])*mu_l[s.prime]+mu_l[s.prime]+2*exp(-mu_l[s.prime])-2))/((1-exp(-mu_l[s]))*(1-exp(-mu_l[s.prime]))^3))+
        ((Cov_lambda[s,s.prime])*(1-exp(-mu_l[s])-exp(-mu_l[s])*mu_l[s])*(1-theta)^2*(1-exp(-mu_l[s.prime])-exp(-mu_l[s.prime])*mu_l[s.prime])/((1-exp(-mu_l[s]))^2*(1-exp(-mu_l[s.prime]))^2))-
        ((((1-theta)*mu_l[s]/(1-exp(-mu_l[s])))+(sd_l[s]/2)*(exp(-mu_l[s])*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2)*(1-theta))/((1-exp(-mu_l[s]))^3))*
        ((1-theta)*mu_l[s.prime]/(1-exp(-mu_l[s.prime])))+(sd_l[s.prime]/2)*(exp(-mu_l[s.prime])*(exp(-mu_l[s.prime])*mu_l[s.prime]+mu_l[s.prime]+2*exp(-mu_l[s.prime])-2)*(1-theta))/((1-exp(-mu_l[s.prime]))^3))
       analyticCov[s.prime,s]<-analyticCov[s,s.prime]
      }
    }

  for(s in 1:S){
      analyticCov[s,s]<-(((1-theta)*mu_l[s])/(1-exp(-mu_l[s])))*(((1-theta)*mu_l[s])/(1-exp(-mu_l[s])))+
        0.5*sd_l[s]*(((exp(-mu_l[s])*mu_l[s])*(1-theta)^2*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2))/((1-exp(-mu_l[s]))*(1-exp(-mu_l[s]))^3))+
        0.5*sd_l[s]*(((exp(-mu_l[s])*mu_l[s])*(1-theta)^2*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2))/((1-exp(-mu_l[s]))*(1-exp(-mu_l[s]))^3))+
        ((Cov_lambda[s,s])*(1-exp(-mu_l[s])-exp(-mu_l[s])*mu_l[s])*(1-theta)^2*(1-exp(-mu_l[s])-exp(-mu_l[s])*mu_l[s])/((1-exp(-mu_l[s]))^2*(1-exp(-mu_l[s]))^2))-
        ((((1-theta)*mu_l[s]/(1-exp(-mu_l[s])))+(sd_l[s]/2)*(exp(-mu_l[s])*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2)*(1-theta))/((1-exp(-mu_l[s]))^3))*
        ((1-theta)*mu_l[s]/(1-exp(-mu_l[s])))+(sd_l[s]/2)*(exp(-mu_l[s])*(exp(-mu_l[s])*mu_l[s]+mu_l[s]+2*exp(-mu_l[s])-2)*(1-theta))/((1-exp(-mu_l[s]))^3))
      analyticVariance[s]<-analyticCov[s,s]
      }

    for(s in 1:(S-1)){
      for(s.prime in (s+1):S){
          analyticCorrelation[s,s.prime]<-analyticCov[s,s.prime]/(sqrt(analyticVariance[s])*sqrt(analyticVariance[s.prime]))
          analyticCorrelation[s.prime, s]<-analyticCorrelation[s,s.prime]
        }
      }

    for(s in 1:S){
        analyticCorrelation[s,s]<-analyticCov[s,s]/(sqrt(analyticVariance[s])*sqrt(analyticVariance[s]))
      }
  return(analyticCorrelation)
}

## Hurdle-AR Model
analyticCorrelationARHurdle<-function(theta,mu_a, mu_phi, var_phi, covariance_matrix_a,N){
  mu_AR<-array(dim=c(R,S,K))
  for(i in 1:R){
    for(s in 1:S){
      mu_AR[i,s,1]<-mu_a[s]
      for(k in 2:K){
        mu_AR[i,s,k]<-mu_a[s]+log(N[i,s,k-1]+1)*mu_phi[s] 
      }
    }
  }
  
  Sigma_AR<-array(dim=c(S,S,R,K))
  for(i in 1:R){
    Sigma_AR[,,i,]<-covariance_matrix_a
  }
  
  Sigma_phi<-diag(var_phi,S)
  
  for(i in 1:R){
    for(k in 2:K){
      Sigma_AR[,,i,k]<-Sigma_AR[,,i,k]+log(N[i,,k-1]+1)^2*Sigma_phi[,]
    }
  }
  
  Cov_lambda<-analyticCovariance<-analyticCorrelation<-array(dim=c(S,S,R,K))
  mu_lambda<-sd_lambda<-analyticVariance<-array(dim=c(R,S,K))
  
  for(i in 1:R){
    for(s in 1:S){
      for(k in 1:K){
        mu_lambda[i,s,k]<-exp(mu_AR[i,s,k]+0.5*Sigma_AR[s,s,i,k])
        sd_lambda[i,s,k]<-exp(2*mu_AR[i,s,k] + Sigma_AR[s,s,i,k])*(exp(Sigma_AR[s,s,i,k])-1)
      }
    }
  }
  
  for(i in 1:R){
    for(k in 1:K){
      for(s in 1:(S-1)){
        for(s.prime in (s+1):S){
          Cov_lambda[s,s.prime,i,k]=exp(mu_AR[i,s,k]+mu_AR[i,s.prime,k]+0.5*(Sigma_AR[s,s,i,k]+Sigma_AR[s.prime,s.prime,i,k]))*(exp(Sigma_AR[s,s.prime,i,k])-1)
          Cov_lambda[s.prime,s,i,k]<-Cov_lambda[s,s.prime,i,k]
        }
      }
      
      for(s in 1:S){
        Cov_lambda[s,s,i,k]<-sd_lambda[i,s,k]
      }
    }
  }
  
  for(i in 1:R){
    for(k in 1:K){
      for(s in 1:(S-1)){
        for(s.prime in (s+1):S){
          analyticCovariance[s,s.prime,i,k]<-(((1-theta)*mu_lambda[i,s,k])/(1-exp(-mu_lambda[i,s,k])))*(((1-theta)*mu_lambda[i,s.prime,k])/(1-exp(-mu_lambda[i,s.prime,k])))+
            0.5*sd_lambda[i,s,k]*(((exp(-mu_lambda[i,s,k])*mu_lambda[i,s.prime,k])*(1-theta)^2*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2))/((1-exp(-mu_lambda[i,s.prime,k]))*(1-exp(-mu_lambda[i,s,k]))^3))+
            0.5*sd_lambda[i,s.prime,k]*(((exp(-mu_lambda[i,s.prime,k])*mu_lambda[i,s,k])*(1-theta)^2*(exp(-mu_lambda[i,s.prime,k])*mu_lambda[i,s.prime,k]+mu_lambda[i,s.prime,k]+2*exp(-mu_lambda[i,s.prime,k])-2))/((1-exp(-mu_lambda[i,s,k]))*(1-exp(-mu_lambda[i,s.prime,k]))^3))+
            ((Cov_lambda[s,s.prime,i,k])*(1-exp(-mu_lambda[i,s,k])-exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k])*(1-theta)^2*(1-exp(-mu_lambda[i,s.prime,k])-exp(-mu_lambda[i,s.prime,k])*mu_lambda[i,s.prime,k])/((1-exp(-mu_lambda[i,s,k]))^2*(1-exp(-mu_lambda[i,s.prime,k]))^2))-
            ((((1-theta)*mu_lambda[i,s,k]/(1-exp(-mu_lambda[i,s,k])))+(sd_lambda[i,s,k]/2)*(exp(-mu_lambda[i,s,k])*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2)*(1-theta))/((1-exp(-mu_lambda[i,s,k]))^3))*
            ((1-theta)*mu_lambda[i,s.prime,k]/(1-exp(-mu_lambda[i,s.prime,k])))+(sd_lambda[i,s.prime,k]/2)*(exp(-mu_lambda[i,s.prime,k])*(exp(-mu_lambda[i,s.prime,k])*mu_lambda[i,s.prime,k]+mu_lambda[i,s.prime,k]+2*exp(-mu_lambda[i,s.prime,k])-2)*(1-theta))/((1-exp(-mu_lambda[i,s.prime,k]))^3))
          analyticCovariance[s.prime,s,i,k]<-analyticCovariance[s,s.prime,i,k]
        }
      }
      
      for(s in 1:S){
        analyticCovariance[s,s,i,k]<-(((1-theta)*mu_lambda[i,s,k])/(1-exp(-mu_lambda[i,s,k])))*(((1-theta)*mu_lambda[i,s,k])/(1-exp(-mu_lambda[i,s,k])))+
          0.5*sd_lambda[i,s,k]*(((exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k])*(1-theta)^2*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2))/((1-exp(-mu_lambda[i,s,k]))*(1-exp(-mu_lambda[i,s,k]))^3))+
          0.5*sd_lambda[i,s,k]*(((exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k])*(1-theta)^2*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2))/((1-exp(-mu_lambda[i,s,k]))*(1-exp(-mu_lambda[i,s,k]))^3))+
          ((Cov_lambda[s,s,i,k])*(1-exp(-mu_lambda[i,s,k])-exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k])*(1-theta)^2*(1-exp(-mu_lambda[i,s,k])-exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k])/((1-exp(-mu_lambda[i,s,k]))^2*(1-exp(-mu_lambda[i,s,k]))^2))-
          ((((1-theta)*mu_lambda[i,s,k]/(1-exp(-mu_lambda[i,s,k])))+(sd_lambda[i,s,k]/2)*(exp(-mu_lambda[i,s,k])*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2)*(1-theta))/((1-exp(-mu_lambda[i,s,k]))^3))*
          ((1-theta)*mu_lambda[i,s,k]/(1-exp(-mu_lambda[i,s,k])))+(sd_lambda[i,s,k]/2)*(exp(-mu_lambda[i,s,k])*(exp(-mu_lambda[i,s,k])*mu_lambda[i,s,k]+mu_lambda[i,s,k]+2*exp(-mu_lambda[i,s,k])-2)*(1-theta))/((1-exp(-mu_lambda[i,s,k]))^3))
        analyticVariance[i,s,k]<-analyticCovariance[s,s,i,k]
      }  
    }
  }
  
  
  for(i in 1:R){
    for(k in 1:K){
      for(s in 1:(S-1)){
        for(s.prime in (s+1):S){
          analyticCorrelation[s,s.prime,i,k]<-analyticCovariance[s,s.prime,i,k]/(sqrt(analyticVariance[i,s,k])*sqrt(analyticVariance[i,s.prime,k]))
          analyticCorrelation[s.prime, s,i,k]<-analyticCorrelation[s,s.prime,i,k]
        }
      }
      for(s in 1:S){
        analyticCorrelation[s,s,i,k]<-analyticCovariance[s,s,i,k]/(sqrt(analyticVariance[i,s,k])*sqrt(analyticVariance[i,s,k]))
      }  
    }
  }
  return(analyticCorrelation)
}
