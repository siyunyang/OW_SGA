
load.lib<-c("ranger", "glmnet", "gbm", "twang", "BART", "mboost", "crossmatch", "ggplot2",
            "rpart", "mvtnorm", "xgboost", "caret", "parallel", "doParallel")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE, repos = "http://cran.us.r-project.org")
sapply(load.lib, require, character=TRUE)
# source functions
dir <- "/dscrhome/sy89/SGA/new_simulation/"


source(paste0(dir,"datagen.R"))
source(paste0(dir,"helper_functions.R"))
source(paste0(dir,"PS_model.R"))
load(paste0(dir,"true_effect.Rdata"))
############################################################
####  Wrap function to generate data and #################
####    perform analysis ###################################
############################################################

run_sim_analysis <- function(j, n, n_X, R, S_test, factor1, factor2, perc_conf, taus){
  
  ## Arguments:
  # j - number of simulation scenario
  # n - number of individuals
  # n_X = number of covariates
  # R = # of true subgroup variables when generate the data
  # S_test = # of subgroups to be tested in the model fit process
  # factor1 = controls strength of alpha_x for main effects in the PS model
  # factor2 = controls strenght of interaction terms in the PS model
  # perc_conf - percent of n_X that are confounding
  # taus = effect sizes of interaction terms in the outcome model, the length should match R
  #####################################################################
  
  ## 1) Generate data
  data1 = sim_data(n, n_X, R, taus, factor1, factor2, perc_conf)
  #subgroup variable
  S = data1$S
  colnames(S) <- paste("S",1:ncol(S),sep="")
  Z = as.matrix(data1$Z) #treatment
  dat <- data.frame(PS = data1$pZ, Treatment= as.factor(data1$Z))
 
  Y = data1$Y
  main_X = cbind(data1$just_X, S)
  true_X <- data1$true_X
  just_X <- data1$just_X
  colnames(just_X) <- paste("X",1:ncol(just_X),sep="")
  data <- data.frame(Treatment= as.factor(data1$Z), just_X, S, Y )
  ## 2) Fit PS models
  #Create model equation for all possible interactions between X and R_test
  if (R == S_test) {
    full_X <- data1$prop_X[, -1]
  }  else if (R < S_test ){
    #random select R_nulls from the categorical variables
    #S_test = min(n_X/2-R, S_test)
    R_null_idx = sample(1:(n_X/2-R) , min((n_X/2-R), S_test - R))
    new_X = just_X[ , - R_null_idx]
    new_S = cbind(S, just_X[, R_null_idx])
    model_eqn <- as.formula(paste("~new_X", paste(paste(paste("new_S[,", 1:ncol(new_S), sep=""), "]*new_X", sep=""), collapse=" + "), sep="+"))
    full_X <- model.matrix(model_eqn) #this includes intercept
    #Fix the colnames
    colnames(full_X) <- c("Intercept", paste("X", (1:ncol(new_X)), sep="_"), paste("S", 1:ncol(new_S), sep="_"),
                        c(sapply(paste("S", 1:ncol(new_S), sep="_"), function(x) paste(x,paste("X", (1:ncol(new_X)) , sep="_"), sep="*"))))
    full_X <- full_X[, -1]
    }  else if  (R > S_test ) {
    new_S = S[, 1: S_test]
    new_X = cbind(just_X, S[, (S_test+1) :ncol(S)])
    
    model_eqn <- as.formula(paste("~new_X", paste(paste(paste("new_S[,", 1:S_test, sep=""), "]*new_X", sep=""), collapse=" + "), sep="+"))
    full_X <- model.matrix(model_eqn) #this includes intercept
    #Fix the colnames
    colnames(full_X) <- c("Intercept", paste("X", (1:(n_X- S_test)), sep="_"), paste("S", 1:S_test, sep="_"),
                          c(sapply(paste("S", 1:S_test, sep="_"), function(x) paste(x,paste("X", (1:(n_X-S_test)) , sep="_"), sep="*"))))
    full_X <- full_X[, -1]}
  
  
    model_returns <- fit_models(Z= Z, fullX= full_X, mainX= main_X, trueX= true_X, S = S, method=c("m1","m2","m3","m4","m5","m7","m8","m9"))
  ############ Now compare across models ############
  ############ balance check ########################
  # number of models
  nmod <- ncol(model_returns$all_weights_ow)
  # number of covariates
  num_allX <- ncol(just_X)
  
# initiate storage matrix
  GASD_OW <- array(NA, dim=c(n_X-R, 2*R, nmod))
  GASD_IPW <- array(NA, c(n_X-R,2*R, nmod))
  OASD_OW <- matrix(NA, n_X-R, nmod)
  OASD_IPW <- matrix(NA, n_X-R, nmod)
    
  for (k in 1:nmod){
  b_ow <- check_balance(S , Z , just_X , model_returns$all_weights_ow[,k])
  b_ipw <- check_balance (S , Z , just_X , model_returns$all_weights_ipw[,k])
  
  GASD_OW[, , k] <- b_ow$groups_asd
  GASD_IPW[, ,k] <- b_ipw$groups_asd
  
  OASD_OW[,k] <- b_ow$overall_asd
  OASD_IPW[,k] <- b_ipw$overall_asd
  }
 
  ##3) Calcualte bias
  # supply true causal effects from the "true_effect.Rdata"

  ATE <- ATEALL[j]
  SATE <- matrix(c(SATE0ALL[j],SATE0ALL[j],SATE1ALL[j],SATE1ALL[j] ), ncol= 2, nrow= 2, byrow = T)
  ATO <- ATOALL[j]
  SATO <- matrix(c(SATO0ALL[j],SATO0ALL[j],SATO1ALL[j],SATO1ALL[j] ), ncol= 2, nrow= 2, byrow = T)
  
  # initiate storage matrix for estimation error
  ERROR_OW <- matrix(NA, 2*R+1, nmod)
  ERROR_IPW <- matrix(NA, 2*R+1, nmod)
  
  for (k in 1:nmod){
    #estimation error
    ERROR_OW[,k] <- cal_esterror(S , Z , model_returns$all_weights_ow[,k], Y, ATO, SATO)
    ERROR_IPW[,k] <- cal_esterror(S , Z , model_returns$all_weights_ipw[,k], Y, ATE, SATE)
  }
  rownames(ERROR_OW) <- c("Overall","V1=0", "V1=1", "V2=0", "V2=1")
  rownames(ERROR_IPW) <- c("Overall","V1=0", "V1=1", "V2=0", "V2=1")
  
  return(list(GASD_OW=GASD_OW, GASD_IPW=GASD_IPW, OASD_OW= OASD_OW, OASD_IPW=OASD_IPW, ERROR_OW=ERROR_OW, ERROR_IPW=ERROR_IPW, conv_stat=conv_stat))
}


run <- function(j, taus, nsim){
  res <- list()
  res= foreach(i=1: nsim,  .export=c('run_sim_analysis', 'sim_data', 'fit_models','check_balance', 'cal_esterror','all_params', 'ow', 'ipw', 'calc_ate', 'abs_stand_diff' ), .packages=c("ranger", "glmnet", "randomForest","gbm", "twang", "BART","rpart", "mvtnorm", "xgboost"))%dopar% 
    {
      tmp= run_sim_analysis(j, all_params[j,"ns"], all_params[j,"n_Xs"],
                            all_params[j,"Rs"], all_params[j,"S_tests"],
                            all_params[j,"factor1s"], all_params[j,"factor2s"],
                            all_params[j,"perc_confs"], taus, ATEALL, ATOALL,SATE0ALL,SATE1ALL,SATO0ALL,SATO1ALL)
      return(tmp)
      
    }
  # average across nsim data sets
  BIAS_OW <- apply(simplify2array(lapply(res, function(x) x$ERROR_OW)),1:2,mean)
  BIAS_IPW <- apply(simplify2array(lapply(res, function(x) x$ERROR_IPW)),1:2,mean)
  
  MSE_OW <- apply(simplify2array(lapply(res, function(x) x$ERROR_OW^2)),1:2,mean)
  MSE_IPW <- apply(simplify2array(lapply(res, function(x) x$ERROR_IPW^2)),1:2,mean)
  
  OASD_OW <- apply(simplify2array(lapply(res, function(x) x$OASD_OW)),1:2,mean)
  OASD_IPW <- apply(simplify2array(lapply(res, function(x) x$OASD_IPW)),1:2,mean)
  
  GASD_OW <- apply(simplify2array(lapply(res, function(x) x$GASD_OW)),1:3,mean)
  GASD_IPW <- apply(simplify2array(lapply(res, function(x) x$GASD_IPW)),1:3,mean)
 
  return(list(GASD_OW=GASD_OW, GASD_IPW=GASD_IPW, OASD_OW= OASD_OW, OASD_IPW=OASD_IPW, BIAS_OW=BIAS_OW, BIAS_IPW=BIAS_IPW, MSE_OW=MSE_OW, MSE_IPW=MSE_IPW))
}


############################################################
#### 7) Start Simulation ###################################
############################################################

all_params <- expand.grid(ns = c(3000),
                          n_Xs = c(20, 50),
                          Rs = c(2),
                          S_tests = c(2, 4, 10),
                          factor1s = c(1, 1.25, 1.5),
                          factor2s = c(0.25,0.5, 0.75),
                          perc_confs = c(.25, 0.75),
                          taus = c(0, 0.5))
all_results <- list()
nsim <- 100

cl <- makeCluster(10)
registerDoParallel(cl)


for(j in 1:nrow(all_params)){
  if (j%%50 ==0) {print(j); save(all_results, file="simulation_part.Rdata")}

  taus = rep( all_params[j,"taus"], all_params[j,"Rs"] )
  all_results[[j]] <-run(j,taus,nsim)

}

stopCluster(cl)
save(all_results, file="simulation_comp.Rdata")
