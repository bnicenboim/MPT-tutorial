rm(list=ls())
set.seed(42);

source("gen_MPT_Matzke.r")
library(rstan)
library(reshape2)
library(ggplot2)




model_0 <- "
// Multinomial Processing Tree
data { 
  int<lower=1> n;
  int<lower=0,upper=n> k[4];
}
parameters {
  real<lower=0,upper=1> c;
  real<lower=0,upper=1> r;
  real<lower=0,upper=1> u;
} 
transformed parameters {
  simplex[4] theta;
  
  // MPT Category Probabilities for Word Pairs
  theta[1] <- c * r;  //C2wc
  theta[2] <- (1 - c) * u ^ 2; //C2wn
  theta[3] <- (1 - c) * 2 * u * (1 - u); //C1w
  theta[4] <- c * (1 - r) + (1 - c) * (1 - u) ^ 2; //C0w
}
model {
  // Priors
  c ~ beta(1, 1);  // can be removed
  r ~ beta(1, 1);  // can be removed
  u ~ beta(1, 1);  // can be removed

  // Data
  k ~ multinomial(theta);
}
generated quantities{
	int predq[4];
  predq <- multinomial_rng(theta, n);
}
"
stan_nhmpt <- stan_model(model_code=model_0)



#############################################################33
#############################################################33
### Example 1. Data from a low variance population


d1 <- gen_pair_clustering_MPT(param_MPT=c(.6,.5,.9),
                        N_subj=40,
                        N_item=20,
                        sdev_subj  <- c(0.03,0.03,0.03),
						sdev_item  <- c(0.03,0.03,0.03)
                        )


#The data has to be arranged into a list.
head(d1)
d1$catname<-factor(d1$catname,levels =c("C2wc","C2wn","C1w","C0w"))
d1_f <- dcast(d1,subj+item~catname,length)
# First we sort have each category in one column.
head(d1_f)

# Then we collapse the data ignoring subjects and items:
k <- colSums(d1_f[c("C2wc","C2wn","C1w","C0w")])
n <- sum(k)
list_d1 <- list(k=k, n=n) # To be passed on to Stan
list_d1

#We fit the data
samples_MPT_1 <- sampling(stan_nhmpt,    
                data=list_d1, 
                iter=500,
                chains=4
                )


samples_MPT_1

#Compare predq with list_d1$k
list_d1


### Plots
c1 <- rstan::extract(samples_MPT_1)$c
r1 <- rstan::extract(samples_MPT_1)$r
u1 <- rstan::extract(samples_MPT_1)$u
dposterior <- data.frame(parameter=rep(c("c","r","u"),each=length(c1) ))
dposterior$posterior <- c(c1,r1,u1)

ggplot(dposterior,aes(posterior,group=parameter,fill=parameter)) +geom_density(alpha=.3)+ geom_vline(aes(xintercept= c(by(dposterior$posterior,dposterior$parameter, median)) ), color="red", linetype="dashed", size=.3)






#############################################################33
#############################################################33
### Example 2. Data from a low variance population, but where our MPT model does not fit the data

#We generate bad data first:
d2 <- gen_pair_clustering_MPT_bad(
                        N_subj=40,
                        N_item=20,
                        sdev_subj  <- c(0.03,0.03,0.03),
						sdev_item  <- c(0.03,0.03,0.03)
                        )


#We arrange the data as before:
head(d2)
d2$catname<-factor(d2$catname,levels =c("C2wc","C2wn","C1w","C0w"))
d2_f <- dcast(d2,subj+item~catname,length)
head(d2_f)

# We ignore subjects and items
k <- colSums(d2_f[c("C2wc","C2wn","C1w","C0w")])
n <- sum(k)
list_d2 <- list(k=k, n=n) # To be passed on to Stan

samples_MPT_2 <- sampling(stan_nhmpt,    
                data=list_d2, 
                iter=500,
                chains=4
                )

samples_MPT_2

#Compare predq of samples_MPT_2 with list_d2$k
list_d2


#############################################################33
#############################################################33
### Example 3. Data from a high variance population

d3 <- gen_pair_clustering_MPT(param_MPT=c(.6,.5,.9),
                        N_subj=40,
                        N_item=20,
                        sdev_subj  <- c(1.2,1.2,1.2),
						sdev_item  <- c(0.03,0.03,0.03),
						rho_subj <- matrix(c( 1,.6,.6,
		  									 .6,1,.6,
		  									 .6,.6,1),nrow=3,byrow=T)

                        )


#The has to be arranged into a list.
head(d3)
library(reshape2)
d3$catname<-factor(d3$catname,levels =c("C2wc","C2wn","C1w","C0w"))
d3_f <- dcast(d3,subj+item~catname,length)
# First we sort have each category in one column.
head(d3_f)

# Then we collapse the data ignoring subjects and items:
k <- colSums(d3_f[c("C2wc","C2wn","C1w","C0w")])
n <- sum(k)
list_d3 <- list(k=k, n=n) # To be passed on to Stan
list_d3

#We fit the data
samples_MPT_3 <- sampling(stan_nhmpt,    
                data=list_d3, 
                iter=500,
                chains=4
                )


samples_MPT_3

#Compare predq with list_d3$k
list_d3


### Plots
c1 <- rstan::extract(samples_MPT_3)$c
r1 <- rstan::extract(samples_MPT_3)$r
u1 <- rstan::extract(samples_MPT_3)$u
dposterior <- data.frame(parameter=rep(c("c","r","u"),each=length(c1) ))
dposterior$posterior <- c(c1,r1,u1)

(p3<-ggplot(dposterior,aes(posterior,group=parameter,fill=parameter)) +geom_density(alpha=.3)+ geom_vline(aes(xintercept= c(by(dposterior$posterior,dposterior$parameter, median)) ), color="red", linetype="dashed", size=.3)
)

