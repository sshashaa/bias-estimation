###############################################
###############################################
##########    Data Generating   ###############
############# Mar 13 2023 #####################
############# Kimia Vahdat ####################
###############################################
###############################################

packages <- c("caret", "dplyr", "GA", "doParallel","extraDistr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

############# Libraries #####################
library(caret)
library(GA)
library(dplyr)
library(doParallel)
library(extraDistr)
library(data.table)
library(MASS)
library(purrr)
##### End loading packages##################

dataGen5= function(p=5, n=100, Corr=TRUE){
  mus= rdunif(p, 2*p, 1)
  covMat= matrix(0, nrow = p, ncol = p)
  covTmp = matrix(runif(p*p, -5, 5), nrow = p, ncol = p)
  covMat = covTmp %*% t(covTmp)
  vars = diag(covMat)
  if(!Corr){
    covMat=matrix(0, nrow = p, ncol = p)
    diag(covMat)=vars
  }
  X= as.data.frame(mvrnorm(n, mus, covMat))
  colnames(X) = c(paste0(rep("X",p),c(1:p)))
  X[,c(1:5)]=X[,c(1:5)]+ replicate(5, rnorm(n, 1, 1))
  y = X[,1]+X[,2]*3-X[,3]*5+X[,4]*X[,1]-X[,5] +rnorm(n,0,2) # true formula
  df = as.data.frame(cbind(X,y))
  colnames(df) = c(names(X),"y")
  return(df)
}

dataGen10= function(p=10, n=100, Corr=TRUE){
  mus= rdunif(p, 2*p, 1)
  covMat= matrix(0, nrow = p, ncol = p)
  covTmp = matrix(runif(p*p, -5, 5), nrow = p, ncol = p)
  covMat = covTmp %*% t(covTmp)
  vars = diag(covMat)
  if(!Corr){
    covMat=matrix(0, nrow = p, ncol = p)
    diag(covMat)=vars
  }
  X= as.data.frame(mvrnorm(n, mus, covMat))
  colnames(X) = c(paste0(rep("X",p),c(1:p)))
  X[,c(1:10)]=X[,c(1:10)]+ replicate(10, rnorm(n, 1, 1))
  y = X[,1]+X[,2]*3-X[,3]*5+X[,4]*X[,1]-X[,5]^2+X[,6]*2+X[,7]*X[,8]/2-X[,9]*X[,10]+X[,10]*5+rnorm(n,0,2) # true formula
  df = as.data.frame(cbind(X,y))
  colnames(df) = c(names(X),"y")
  return(df)
}

setwd("/home/kvahdat/Data/ch4/")
### Generating 3 datasets ###
set.seed(1)
df1= dataGen5(p=15, n=300, Corr = FALSE)
# summary(lm(y~., df1)) # misses X5 
write.csv(df1, "./SmallDF.csv", row.names = FALSE)

set.seed(1)
df2= dataGen10(p=30, n=300, Corr = FALSE)
# summary(lm(y~., df2)) # picks x16 extra
write.csv(df2, "./LargeDF.csv", row.names = FALSE)

set.seed(10)
df3= dataGen5(p=15, n=300, Corr = TRUE)
# summary(lm(y~., df3)) 
write.csv(df3, "./SmallDF_Corr.csv", row.names = FALSE)

set.seed(10)
df4= dataGen10(p=30, n=300, Corr = TRUE)
# summary(lm(y~., df4))
write.csv(df4, "./LargelDF_Corr.csv", row.names = FALSE)