#rm(list=ls()) 
model_h <- "
// Multinomial Processing Tree with Latent Traits
data { 
	int<lower=1> nsubjs; 
  int<lower=1> nitems; 
  int<lower=1> n; 
	int<lower=1, upper=n> subj[n];
  int<lower=0,upper=20> k[nsubjs,4];

}
transformed data {
  
  int<lower=1> npar; 
  npar <- 3;
	
}
parameters {
  vector[npar] mus_hat;
  vector<lower=0>[npar] sigma_s;
  cholesky_factor_corr[npar] L_s;
  matrix[npar,nsubjs] z_s;

} 
transformed parameters {
  simplex[4] theta[nsubjs];
  vector<lower=0,upper=1>[nsubjs] c;
  vector<lower=0,upper=1>[nsubjs] r;
  vector<lower=0,upper=1>[nsubjs] u;
  matrix[nsubjs,npar] s;
  


  s <- (diag_pre_multiply(sigma_s,L_s) * z_s)'; // subj random effects'

		
	for (i in 1:nsubjs) {
		
		// Probitize Parameters c, r, and u 
		c[i] <- Phi(mus_hat[1] +  s[subj[i],1]);
		r[i] <- Phi(mus_hat[2] +  s[subj[i],2]);
		u[i] <- Phi(mus_hat[3] +  s[subj[i],3]);
		
		// MPT Category Probabilities for Word Pairs
		theta[i,1] <- c[i] * r[i];
		theta[i,2] <- (1 - c[i]) * (u[i]) ^ 2;
		theta[i,3] <- (1 - c[i]) * 2 * u[i] * (1 - u[i]);
		theta[i,4] <- c[i] * (1 - r[i]) + (1 - c[i]) * (1 - u[i]) ^ 2;
	}

}
model {
  // Priors
  mus_hat ~ normal(0,1);
  sigma_s ~ normal(0,2);
  L_s ~ lkj_corr_cholesky(4.0); 
  to_vector(z_s) ~ normal(0,1);


  for (i in 1:nsubjs)
    k[i] ~ multinomial(theta[i]);
}
generated quantities {
  real<lower=0,upper=1> muc;
  real<lower=0,upper=1> mur;
  real<lower=0,upper=1> muu;
#  vector<lower=0>[4] predq[n];
  matrix[npar,npar] Cor_s;

  // Post-Processing Means, Standard Deviations, Correlations
	muc <- Phi(mus_hat[1]);
	mur <- Phi(mus_hat[2]);
	muu <- Phi(mus_hat[3]);
  Cor_s <- L_s * L_s';

#This doesn't work, no idea why 
    # for (i in 1:n)
    #    predq[i] <- multinomial_rng(theta[i],nitems);

}"


stan_hmpt <- stan_model(model_code=model_h)



########################################3
#########################################
##  Example 4. Model that takes into account "subjects"
d3_bysubj <- dcast(d3,subj~catname,length)
head(d3_bysubj)

list_d4 <- list(nsubjs=length(unique(d3_bysubj$subj)), 
                nitems=sum(d3_bysubj[1,2:5]),  #sum of responses for one line is the number of items
                n =nrow(d3_bysubj),
                subj = d3_bysubj$subj,
                k= as.matrix(d3_bysubj[2:5])
)

samples_MPT_4 <- sampling(stan_hmpt,    
                data=list_d4, 
                iter=500,
                chains=4
                )

#The output of samples_MPT_4 is too long now

print(samples_MPT_4,pars=c("mus_hat", 
"muc","mur","muu","Cor_s","sigma_s"),digits_summary = 2)

#The fit is not perfect but much better.


### Plots
c1 <- rstan::extract(samples_MPT_4)$muc
r1 <- rstan::extract(samples_MPT_4)$mur
u1 <- rstan::extract(samples_MPT_4)$muu
dposterior2 <- data.frame(parameter=rep(c("c","r","u"),each=length(c1) ))
dposterior2$posterior <- c(c1,r1,u1)

dposterior$type <-"aggregated"
dposterior2$type <-"hierarchical"

dposterior <- rbind(dposterior,dposterior2)

(p4<-ggplot(dposterior,aes(posterior,group=parameter,fill=parameter)) +geom_density(alpha=.3)+ geom_vline(aes(xintercept= with(dposterior[dposterior$type=="hierarchical",],c(by(posterior,parameter, median))) ), color="red", linetype="dashed", size=.3)+ facet_grid(type~.)
)









