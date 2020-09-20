# update 2020/3/16
################################################################
####### #3) Fit propensity score function ######################
################################################################


fit_models <- function(Z, fullX, mainX, trueX, S, method=c("m1","m2","m3","m4","m5","m6","m7","m8")){
  ## This functino fits propensity score models based on input data, 
  ## and returns estimated PS, IPW, and OW
  
  ## Inputs: 
  #### Z - treatment,
  #### fullX - all interactions with X and S, 
  #### mainX - all main effects of covariates and subgroup variables, 
  #### trueX - X with nonzero coefficients,
  #### S - subgroup variables
  #### method - m1:  true model, m2: logistic Main effects, m3: LASSO 
  ####          m4: post-LASSO, m5: RF_main, m6: RF_All, m7: GBM, m8: BART, m9: logistic full interaction
  ## Returns: propensuity scores, resulting overlap indices, and weights
  
  PS <- NULL
  ### Model 1 : Logistic main X  ####
  if("m1" %in% method){
  model_eqn <- as.formula("Z ~ mainX ")
  pscore_mod <- glm(model_eqn, family="binomial")
  pscore <- predict(pscore_mod, type="response")
  PS <- cbind(PS, pscore)
  }
  
  ###  Model 2 :  True effects ####
  if("m2" %in% method){
  pscore_mod2 <- glm(Z~trueX, family="binomial")
  pscore <- predict(pscore_mod2, type="response")
  PS <- cbind(PS, pscore)
  }
 
  ####  Model3: Lasso - with higher order terms###
  ### use penalty.factor to keep the main effect
  if("m3" %in% method){
  pscore_mod3 <- cv.glmnet(y=Z, x=fullX, penalty.factor=c(rep(0, ncol(mainX)), rep(1, ncol(fullX)-ncol(mainX))), family="binomial", maxit=50000)
  pscore <- predict(pscore_mod3, newx=fullX, type="response", s='lambda.min')
  PS <- cbind(PS, pscore)
  }

 
  ######  Model 4:  Relaxed Lasso - Find which variables are selected by lasso ###
  if("m4" %in% method){
  nonzero_coef <- rownames(coef(pscore_mod3, s='lambda.min'))[which(coef(pscore_mod3, s='lambda.min')!=0)][-1]
  if(length(nonzero_coef)==0){
    print("no coef")
    pscore_mod4 <- glm(Z~1, family="binomial")
  }else{
    pscore_mod4 <- glm(Z~., data=data.frame(fullX[,nonzero_coef]), family="binomial")
    
  }
  pscore <- pscore_mod4$fitted.values
  PS <- cbind(PS, pscore)
  }
 
  #M5: RF_main 
  if("m5" %in% method){
  comb_dat <- data.frame(Z=Z, mainX )
  pscore_mod5 <- ranger(factor(Z) ~ ., data = comb_dat, num.trees = 1000, replace = T, verbose = FALSE, probability = TRUE,
                        num.threads =10, importance = "none", seed = 1234)
  # outbag prediction
  pscore <- pscore_mod5$predictions[,2]
  PS <- cbind(PS, pscore)
  
}
  ### Model6: Random Forest with all interactions
  if("m6" %in% method){
    
  comb_dat2 <- data.frame(Z=Z, fullX)
  pscore_mod6 <- ranger(factor(Z) ~ ., data = comb_dat2, num.trees = 1000, replace = T, verbose = FALSE, probability = TRUE,
                        num.threads = 10, importance = "none", seed = 1234)
  
  # outbag prediction
  pscore <- pscore_mod6$predictions[,2]
  PS <- cbind(PS, pscore)
  
  }
  
  ### Model7: gbm
  if("m7" %in% method){
    comb_dat <- data.frame(Z=Z, mainX )
    pscore_mod7 <- ps(Z ~ ., data = comb_dat, n.trees=5000, interaction.depth=2, shrinkage=0.01,
                      perm.test.iters=0, stop.method=c("es.mean"), estimand = "ATE",  verbose=FALSE)
    pscore <- pscore_mod7$ps
    PS <- cbind(PS, pscore)
    
  }
  
 
  
  ### Model8: bart
  if("m8" %in% method){
  pscore_mod8 <- pbart( x.train=mainX, y.train=Z, printevery=1000L)
  tmp <- pnorm(pscore_mod8$yhat.train)
  ## average over the posterior samples
  pscore <- apply(tmp, 2, mean)
  PS <- cbind(PS, pscore)
  
  }
  ### Model 9 : Logistic full interaction  ####
  if("m9" %in% method){
    model_eqn9 <- as.formula("Z ~ fullX ")
    #pscore_mod9 <- glm(model_eqn9, family="binomial")
    pscore_mod9 <-try( glm(model_eqn9, family="binomial"))
    conv <- inherits(pscore_mod9,"warnings")
    pscore9 <- predict(pscore_mod9, type="response")
    PS <- cbind(PS, pscore9)
  }

  #Estimate weights
  OW <- NULL
  IPW <- NULL
  for (k in 1:length(method)){
  ows <- ow(Z, PS[,k])
  ipws <- ipw(Z, PS[,k])
  OW <- cbind(OW, ows)
  IPW <- cbind(IPW, ipws)
  }

  return(list(all_pscores= PS, all_weights_ow= OW, all_weights_ipw= IPW, conv =conv))
}
