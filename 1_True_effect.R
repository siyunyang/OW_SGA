
library(gridExtra)
library(ggplot2)
library(mvtnorm)

#set working directory
#dir<-


source(paste0(dir,"datagen.R"))
source(paste0(dir,"helper_functions.R"))

theme_set(theme_bw(base_size = 16))
# simulation scenarios
all_params <- expand.grid(ns = c(1000000),
                          n_Xs = c(20, 50),
                          Rs = c(2),
                          S_tests = c(2, 4, 10),
                          factor1s = c(1, 1.25, 1.5),
                          factor2s = c(0.25,0.5, 0.75),
                          perc_confs = c(.25, 0.75),
                          taus = c(0, 0.5))
set.seed(123)
ATOALL <- ATEALL <- SATO0ALL <- SATO1ALL <- SATE0ALL <- SATE1ALL<- NULL

for(j in 1:nrow(all_params)){
taus = rep( all_params[j,"taus"], all_params[j,"Rs"] )
data1 <- sim_data(all_params[j,"ns"], all_params[j,"n_Xs"],
                  all_params[j,"Rs"], taus,
                  all_params[j,"factor1s"], all_params[j,"factor2s"],
                  all_params[j,"perc_confs"] )
dat <- data.frame(PS = data1$pZ, S= data1$S, Treatment= as.factor(data1$Z), Z= data1$Z, Y=data1$Y)

# estimate ATO
true_ow <- ow(dat$Z, dat$PS)
true_ipw <- ipw(dat$Z, dat$PS)
dat$true_ow <- true_ow
dat$true_ipw <- true_ipw

ATO <- calc_ate(dat$Z, dat$Y, dat$true_ow)
ATE <- calc_ate(dat$Z, dat$Y, dat$true_ipw)
ATOALL <- c(ATOALL,ATO) # -0.67 when taus = c(0.5,0.5); -1 when taus = c(0,0)
ATEALL <- c(ATEALL, ATE) # -0.75 when taus = c(0.5,0.5); -1 when taus = c(0,0)

# SATO
sdat0 <- subset(dat, dat$S.1 == 0 )
sdat1 <- subset(dat, dat$S.1 == 1 )

SATO0 <- calc_ate(sdat0$Z, sdat0$Y, sdat0$true_ow)
SATO1 <- calc_ate(sdat1$Z, sdat1$Y, sdat1$true_ow)
SATO0ALL <- c(SATO0ALL,SATO0 )
SATO1ALL <- c(SATO1ALL,SATO1) 

SATE0 <- calc_ate(sdat0$Z, sdat0$Y, sdat0$true_ipw)
SATE1 <- calc_ate(sdat1$Z, sdat1$Y, sdat1$true_ipw)
SATE0ALL <- c(SATE0ALL,SATE0 )
SATE1ALL <- c(SATE1ALL,SATE1) 
}
save(ATOALL ,ATEALL , SATO0ALL , SATO1ALL, SATE0ALL, SATE1ALL, file=paste("true_effect",".Rdata",sep=""))
