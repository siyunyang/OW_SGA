
library("ggplot2")
library(plyr)
library(dplyr)
library(mvtnorm)


# source functions
# set working directory
# dir <- 

source(paste0(dir,"helper_functions.R"))

#################################################
### R function for Connect-S Plot ###############
 
connectS <- function(subg , Z , cov , w=NULL, path, width, height, ss=NULL, vif=NULL, asmd=NULL){
# subg: an array or dataframe which contains the subgrouping variables by column 
#       categorical variables need to be converted to factors 
# z : a numeric 0,1 vector that specifies treatment assignment ( 0 for control, 1 for treated)
# cov: a numeric array or dataframe which contains the covariates.
#       categorical variables need to be converted to dummy variables 
# w: a numeric vector that specifies weights for each unit
# path : specifies the path to save the plot
# width : width of the graph
# height: height of the graph
  library(ggplot2)
  library(plyr)
  library(dplyr)
#  if weight is not provided, need to provide summary statistics
  
if  (is.null(w)){
  mydata <- asmd
  nsubg <- ss
  VIFW <- vif
  nv <- length(unique(asmd$Subgroups))
  ncov <- length(unique(asmd$VarName))
} else{ # if weight is provided

bal_return <- check_balance(subg , Z , cov , w)

GASD <- bal_return$groups_asd
nlevels <- bal_return$nlevels
# subgroup sample size
 nsubg <- bal_return$nsubg
 # number of subgroups
 nv <- length(nsubg)
 
# calculate VIF  
# weights approximation

VIFW <- round( vifw(subg, Z, w)[-1],2)
  

  
 mydata <- data.frame(ASMD = c(GASD),
                           Subgroups = rep(colnames(GASD), each=nrow(GASD)),
                           VarName = rep(rownames(GASD), nv)
                           ) 

}


toplot<- mutate(mydata, lcolor = ifelse(ASMD <=0.1, "0",
                                            ifelse((ASMD <=0.15), "1",
                                                   ifelse((ASMD <=0.25), "2","3"))))
  
group.colors <- c("0" = "White", "1" = "grey85", "2" ="grey50","3"="black","NA"=NA)

g1 <- ggplot(toplot, aes(x =VarName , y =Subgroups , fill=factor(lcolor))) +
  geom_dotplot(binaxis = "y",binwidth = 0.75, stackdir = "center")+labs(x="Confounder", y="Subgroup"  )+
  scale_fill_manual(values=group.colors,name="ASMD",  labels=c("<0.1", "0.1-0.15", "0.15-0.25", ">0.25","NA"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text( angle=-45,vjust = -0.5, hjust=0.3),axis.text = element_text( size = 16 ),
        axis.title = element_text( size = 20, face = "bold" ),legend.position="bottom",legend.title = element_text(size=18),
        legend.text=element_text(size=18), legend.key.size = unit(2,"line"))+
  annotate("text", x = rep(ncov+1.25,nv), y = seq(1, nv,1), label = format(nsubg,zero.print = T),size = 6)+
  annotate("text", rep(ncov+3,nv), y = seq(1, nv,1), label = format(VIFW,zero.print = T),size = 6)+
  annotate("text", x = c(ncov+1,ncov+2), y = rep(nv+1,2), label = c("Size", "VIF"),size = 6)+
 # geom_hline(yintercept=cumsum(nlevels)+0.5)+
  scale_y_discrete(expand=expand_scale(mult = c(0, 0.15)))  + scale_x_discrete(expand=expand_scale(mult = c(0, 0.17)))

splot <- g1+ theme(panel.border = element_rect(linetype = "solid", fill = NA))
ggsave(path,splot, width = width, height = height,dpi = 300 )

return(splot)
}

