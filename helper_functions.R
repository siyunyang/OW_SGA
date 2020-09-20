################################################
## 2) Helper functions ##########################
#######A few helper functions:###################
################################################

calc_ate <- function(A, Y, weights){
  #Calculate average treatment effect # 
  # Inputs:
  # A: a vector of treatment indicator
  # Y: a vecor of the outcomes
  # weights: a vector of weights
  if (anyNA(weights)) {return (NA)} else
    ate <- sum((A*weights*Y))/sum((weights*A)) - sum(((1-A)*weights*Y))/sum((weights*(1-A)))
  return(ate)
}

abs_stand_diff <- function(x_j, z, w){
  # Calculate absolute standardized difference # 
  # Inputs:
  # x_j: a vecor of the covariate
  # z: a vector of treatment indicator
  # w: a vector of weights
  if (anyNA(w)) {return (NA)} else
    delta = abs(sum(x_j*z*w)/sum(z*w) - sum(x_j*(1-z)*w)/sum((1-z)*w))
    tstat <- delta/sqrt((var(x_j[which(z==1)])+var(x_j[which(z==0)]))/2)
  return (tstat)
}


ow <- function(Z, ps){
  return( Z*(1-ps) + (1-Z)*ps)
}

ipw <- function(Z, ps){
  return( Z/ps + (1-Z)/(1-ps))
}

vifw <- function(S, z, w){
  vif <- NULL
  #overall vifw
  vifw_overall<- (sum(w^2*z)/sum(w*z)^2 +sum(w^2*(1-z))/(sum(w*(1-z)))^2)/(1/sum(z)+1/(length(z)-sum(z)))
  for(r in 1:ncol(S))
  { 
    for(g in 1:2)
    {
      find_g <- which(S[, r] == (g-1))
      n1 <- sum(z[find_g]==1)
      n0 <- sum(z[find_g]==0)
      vifw <- (sum(w[find_g]^2*z[find_g])/sum(w[find_g]*z[find_g])^2 +sum(w[find_g]^2*(1-z[find_g]))/(sum(w[find_g]*(1-z[find_g])))^2)/(1/n1+1/n0)
      vif <-c(vif, vifw)
      # names_col <- c(names_col, paste("group", paste(r,g, sep="_"), sep="-"))
      
    }
  }
    return(c(vifw_overall,vif))
}


###########################################################
#######  4)check_balance  ##########################
###########################################################
check_balance <- function(S, Z, X, weights){
  # This function evaulate balance overall and within subgroups
  # Inputs:
  # S - a n by n_S matrix of subgroup indicators
  # Z - a vector of treatment indicators
  # X - a n by n_X matrix ofcovariates
  # weights - a vector of weights
  if (anyNA(weights)) {
    overall_asd <- rep(NA, ncol(X))
    groups_asd <- matrix(NA, ncol(X), ncol(S)*2)
    
    return (list(groups_asd=groups_asd, overall_asd=overall_asd))} else
      #Calculate overall balance across covariates
      overall_asd <- apply(X, 2, abs_stand_diff, Z, weights)
    
    #Calculate balance per subgroup across covariates
    groups_asd <- c()
    names_col <- c()
    
    for(r in 1:ncol(S)){ 
      for(g in 1:length(unique(S[,r]))){
        #ASD
        #X <- cbind(X,S[,-r])
        find_g <- which(S[,r]==(g-1))
        g_asd <- apply(X[find_g, ], 2, abs_stand_diff, Z[find_g], weights[find_g])
        groups_asd <- cbind(groups_asd, g_asd)
        names_col <- c(names_col, paste("group", paste(r,g, sep="_"), sep="-"))
      }
    }
    
    colnames(groups_asd) <- names_col
    return(list(groups_asd=groups_asd, overall_asd=overall_asd))
}


###########################################################
#######  5)Estimate ATE and SATE ##########################
####### Calculate estimation error ####################################
###########################################################
cal_esterror <- function(S, Z, weights, Y, ATE, SATE)
  # This function estimates overall and within subgroups causal effect 
  # and calculate the esterror for ate and sate
  
  # Inputs:
  # S - a n by n_S matrix of subgroup indicators
  # Z - a vector of treatment indicators
  # weights - a vector of weights
  # Y - a vector of outcomes
  # ATE - true ate
  # SATE - a 2 by ncol(S) matrix of true SATE
{
  #Calculate ATE overall
  ate_overall <- sum((Z*weights*Y))/sum((weights*Z)) - sum(((1-Z)*weights*Y))/sum((weights*(1-Z)))
  #ate_unweighted <- sum((Z*Y))/sum((Z)) - sum(((1-Z)*Y))/sum(((1-Z)))
  #ate_all <- c(ate_overall, ate_unweighted)
  ate_all <- ate_overall
  
  #estimated sate 
  sate <- matrix(0, nrow=2, ncol=ncol(S))
  
  for(r in 1:ncol(S))
  { 
    for(g in 1:2)
    {
      find_g <- which(S[, r] == (g-1))
      sate[g, r] <- sum((Z*weights*Y)[find_g])/sum((weights*Z)[find_g]) - sum(((1-Z)*weights*Y)[find_g])/sum((weights*(1-Z))[find_g])
      
    }
  }
  esterror_ate <- ate_all - ATE
  esterror_sate <- c(sate - SATE)
  esterror <-  matrix(c(esterror_ate, esterror_sate),5,1)
  #rownames(esterror) <- c("Overall","V1=0", "V1=1", "V2=0", "V2=1")
  return(esterror)
}


