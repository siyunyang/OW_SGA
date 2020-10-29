

# source functions
# set working directory
dir <- 

source(paste0(dir,"helper_functions.R"))
source(paste0(dir,"/Connect_S_func.R"))

#########################################
######## Illustration 1 #################

data <- read.csv("ConnectS_data.csv")
# Fit a logistic main effect PS model 
pscore_mod <- glm(Treatment~., family="binomial", data =data)
pscore <- predict(pscore_mod, type="response")

#Estimate IPW weights
ipw <- ipw(data$Treatment, pscore)

w <- ipw
cov <- data[,2:19]
# subgroups of interest
subg <- select(data, S1,S2,X1,X3,X7,X8)
ncov <- dim(cov)[2]
path <- "/desktop/Connect-S.png"
Z <- data$Treatment

### example call
### plot will be saved to the path 
connectS(subg , Z , cov , w, path, width=10, height=6)



#########################################
######## Illustration 2 #################


ufdata <- read.csv(paste(dir,"selected.csv",sep=""))
asmd <- data.frame(ASMD = abs(ufdata$w_std_i1),
                   Subgroups =  ufdata$subg_new_c,
                   VarName =  ufdata$vbles
) 
nv <- length(unique(ufdata$subg_new_c))
ncov <- length(unique(ufdata$vbles))
vif <- round(ufdata$vi_ipw1[1:nv],1)
ss <-  ufdata$ss[1:nv]

# specify orders of the subgroup and covariates
asmd$VarName <- factor(asmd$VarName, levels =unique(asmd$VarName))
asmd$Subgroups <- factor(asmd$Subgroups, levels =unique(asmd$Subgroups))
path <- "desktop/Connect-S2.png"

connectS(subg , Z , cov , w=NULL, path, width=10, height=6,ss, vif, asmd)

### all
ufdata <- read.csv("all.csv")
asmd <- data.frame(ASMD = abs(ufdata$w_std_i1),
                   Subgroups =  ufdata$subg_c,
                   VarName =  ufdata$vbles
) 
nv <- length(unique(ufdata$subg_c))
ncov <- length(unique(ufdata$vbles))
vif <- round(ufdata$vi_ipw1[1:nv],2)
ss <-  ufdata$ss[1:nv]

# specify orders of the subgroup and covariates
asmd$VarName <- factor(asmd$VarName, levels =unique(asmd$VarName))

asmd$Subgroups <- factor(asmd$Subgroups, levels =unique(asmd$Subgroups))
path <- "/desktop/Connect-S3.png"

connectS(subg , Z , cov , w=NULL, path, width=13, height=16,ss, vif, asmd)
