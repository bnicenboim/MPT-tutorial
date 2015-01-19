library(MASS)
library(Hmisc)

gen_pair_clustering_MPT <- function(            

                   rho_subj=NULL
                   ,
                   rho_item=NULL
                   ,
                   sdev_subj=NULL
                   ,
                   sdev_item=NULL
                   ,
                   param_MPT=NULL
                   ,
                   N_subj=40
                   ,
                   N_item=20
       
                   ) {


N_par <- 3
       
if(is.null(param_MPT)) param_MPT <- c(.6,.8,.9)
param_MPT <- qnorm(param_MPT)
if(is.null(sdev_subj)) sdev_subj  <- c(0.3,0.3,0.3)
if(is.null(sdev_item)) sdev_item  <- c(0.03,0.03,0.03)
ms <- matrix(rep(.3,N_par*N_par),N_par,N_par)
diag(ms) <-1
mi <- matrix(rep(.03,N_par*N_par),N_par,N_par)
diag(mi) <-1
if(is.null(rho_subj)) rho_subj <-   ms
if(is.null(rho_item)) rho_item <-  mi


#stdevs is the vector that contains the standard deviations of the var
b_subj <- sdev_subj %*% t(sdev_subj)  
# b is an 5*5 matrixstdevs_item whose generic term is stdev[i]*stdev[j]
Sigma_par_subj <- b_subj * rho_subj  #variance covariance matrix for subj

#stdevs is the vector that contains the standard deviations of the var
b_item <- sdev_item %*% t(sdev_item)  
Sigma_par_item <- b_item * rho_item 


#displacement parameters by subj and by item
raneff_subj <- mvrnorm(n = N_subj, rep(0,N_par), Sigma_par_subj)
raneff_item <- mvrnorm(n = N_item, rep(0,N_par), Sigma_par_item)


raneff_subj <- data.frame(u=raneff_subj)

raneff_item <- data.frame(w=raneff_item)


data <- data.frame(subj =(rep(seq(1:N_subj),each=N_item)),
                    item =(rep(seq(1:N_item),times=N_subj)),
                    raneff_subj[rep(1:nrow(raneff_subj),each=N_item),],
                    raneff_item[rep(1:nrow(raneff_item),times=N_subj),]
)



#Parameters for each observation in the probability space:

data$c <- with(data,pnorm(param_MPT[1] + u.1 + w.1))
data$r <- with(data,pnorm(param_MPT[2] + u.2 + w.2))
data$u <- with(data,pnorm(param_MPT[3] + u.3 + w.3))

PC2wc <- with(data,c * r)
PC2wn <- with(data,(1-c) * u^2)
PC1w  <- with(data,2*u *(1-c) *(1-u))
PC0w  <- with(data,c * (1-r) +(1-c)*(1-u)^2)

mProb <- matrix(c(PC2wc,PC2wn,PC1w,PC0w),ncol=4)

data$cat <- rMultinom(mProb,1)

data$catname <- c("C2wc","C2wn","C1w","C0w")[data$cat]

dataout <- data[c("subj","item","cat", "catname")]


  return(dataout)
}





gen_pair_clustering_MPT_bad <- function(            

                   rho_subj=NULL
                   ,
                   rho_item=NULL
                   ,
                   sdev_subj=NULL
                   ,
                   sdev_item=NULL
                   ,
                   param_MPT=NULL
                   ,
                   N_subj=40
                   ,
                   N_item=20
       
                   ) {


N_par <- 3
       
if(is.null(param_MPT)) param_MPT <- c(.6,.8,.9)
param_MPT <- qnorm(param_MPT)

if(is.null(sdev_subj)) sdev_subj  <- c(0.3,0.3,0.3)
if(is.null(sdev_item)) sdev_item  <- c(0.03,0.03,0.03)
ms <- matrix(rep(.3,N_par*N_par),N_par,N_par)
diag(ms) <-1
mi <- matrix(rep(.03,N_par*N_par),N_par,N_par)
diag(mi) <-1
if(is.null(rho_subj)) rho_subj <-   ms
if(is.null(rho_item)) rho_item <-  mi


#stdevs is the vector that contains the standard deviations of the var
b_subj <- sdev_subj %*% t(sdev_subj)  
# b is an 5*5 matrixstdevs_item whose generic term is stdev[i]*stdev[j]
Sigma_par_subj <- b_subj * rho_subj  #variance covariance matrix for subj

#stdevs is the vector that contains the standard deviations of the var
b_item <- sdev_item %*% t(sdev_item)  
Sigma_par_item <- b_item * rho_item 


#displacement parameters by subj and by item
raneff_subj <- mvrnorm(n = N_subj, rep(0,N_par), Sigma_par_subj)
raneff_item <- mvrnorm(n = N_item, rep(0,N_par), Sigma_par_item)


raneff_subj <- data.frame(u=raneff_subj)

raneff_item <- data.frame(w=raneff_item)


data <- data.frame(subj =(rep(seq(1:N_subj),each=N_item)),
                    item =(rep(seq(1:N_item),times=N_subj)),
                    raneff_subj[rep(1:nrow(raneff_subj),each=N_item),],
                    raneff_item[rep(1:nrow(raneff_item),times=N_subj),]
)



#Parameters for each observation in the probability space:

data$c <- with(data,pnorm(param_MPT[1] + u.1 + w.1))
data$r <- with(data,pnorm(param_MPT[2] + u.2 + w.2))
data$u <- with(data,pnorm(param_MPT[3] + u.3 + w.3))

PC2wc <- with(data,c )
PC2wn <- with(data,(1-c) * u*.1)
PC1w  <- with(data,(1-c) *(u*(1-.1)+.1* (1-u)))
PC0w  <- with(data,(1-c)*(1-u)*(1-.1))

mProb <- matrix(c(PC2wc,PC2wn,PC1w,PC0w),ncol=4)

data$cat <- rMultinom(mProb,1)

data$catname <- c("C2wc","C2wn","C1w","C0w")[data$cat]

dataout <- data[c("subj","item","cat", "catname")]


  return(dataout)
}






gen_pair_clustering_MPT_WMC <- function(            

                   rho_subj=NULL
                   ,
                   rho_item=NULL
                   ,
                   sdev_subj=NULL
                   ,
                   sdev_item=NULL
                   ,
                   param_MPT=NULL
                   ,
                   N_subj=40
                   ,
                   N_item=20
       
                   ) {


N_par <- 3
       
if(is.null(param_MPT)) param_MPT <- c(.6,.8,.9)
param_MPT <- qnorm(param_MPT)
if(is.null(sdev_subj)) sdev_subj  <- c(0.3,0.3,0.3)
if(is.null(sdev_item)) sdev_item  <- c(0.03,0.03,0.03)
ms <- matrix(rep(.3,N_par*N_par),N_par,N_par)
diag(ms) <-1
mi <- matrix(rep(.03,N_par*N_par),N_par,N_par)
diag(mi) <-1
if(is.null(rho_subj)) rho_subj <-   ms
if(is.null(rho_item)) rho_item <-  mi


#stdevs is the vector that contains the standard deviations of the var
b_subj <- sdev_subj %*% t(sdev_subj)  
# b is an 5*5 matrixstdevs_item whose generic term is stdev[i]*stdev[j]
Sigma_par_subj <- b_subj * rho_subj  #variance covariance matrix for subj

#stdevs is the vector that contains the standard deviations of the var
b_item <- sdev_item %*% t(sdev_item)  
Sigma_par_item <- b_item * rho_item 


#displacement parameters by subj and by item
raneff_subj <- mvrnorm(n = N_subj, rep(0,N_par), Sigma_par_subj)
raneff_item <- mvrnorm(n = N_item, rep(0,N_par), Sigma_par_item)


raneff_subj <- data.frame(u=raneff_subj)

raneff_item <- data.frame(w=raneff_item)


data <- data.frame(subj =(rep(seq(1:N_subj),each=N_item)),
                    item =(rep(seq(1:N_item),times=N_subj)),
                    raneff_subj[rep(1:nrow(raneff_subj),each=N_item),],
                    raneff_item[rep(1:nrow(raneff_item),times=N_subj),],
                    pcu = rep(round(pnorm(rnorm(N_subj,0,1)),2),each=N_item)

)



#Parameters for each observation in the probability space:

data$c <- with(data,pnorm(param_MPT[1]  + u.1 + w.1))
data$r <- with(data,pnorm(param_MPT[2] + u.2 + w.2))
data$u <- with(data,pnorm(param_MPT[3]+scale(pcu) + u.3 + w.3))

PC2wc <- with(data,c * r)
PC2wn <- with(data,(1-c) * u^2)
PC1w  <- with(data,2*u *(1-c) *(1-u))
PC0w  <- with(data,c * (1-r) +(1-c)*(1-u)^2)

mProb <- matrix(c(PC2wc,PC2wn,PC1w,PC0w),ncol=4)

data$cat <- rMultinom(mProb,1)

data$catname <- c("C2wc","C2wn","C1w","C0w")[data$cat]

dataout <- data[c("subj","item","cat", "catname","pcu")]


  return(dataout)
}

