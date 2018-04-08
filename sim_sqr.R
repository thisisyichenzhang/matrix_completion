source("./R_solver_change.R")

#gen_0<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 0) # SNR INF
#gen_0.01<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 0.01) # SNR 300
#gen_0.1<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 0.1) # SNR 30
#gen_0.5<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 0.5) # SNR 6
#gen_1<-gen_data(n = 50,rank = 2,percent_missing = 0.1,sigma = 2) # SNR 3.0
#gen_10<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 10) # SNR 0.3
#gen_100<-gen_data(n = 50,rank = 2,percent_missing = 0.01,sigma = 100) # SNR 0.03

# find the optimal lambda
#results$lambdas[which.min(results$errors)]

# data frame for save
# if existing use the following command instead 
# load("./Results/sqr_results.Rdata")
sim_results<-data.frame(n=numeric(),
                        rank=integer(),
                        percent_missing=numeric(),
                        sigma=numeric(),
                        SNR = numeric(),
                        lambda=numeric(),
                        error=numeric())
# seed for reproducibility
set.seed(666)
# number of simulations per set 
k=1000
# results for each set of training values
results<-numeric(k)
# training values
n_vec<-c(50,100)
rank_vec<-c(2,3,4,5)
percent_missing_vec<-c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
sigma_vec<-c((0:9)*0.01,0.1,0.5,1,10)
lambda_vec<-c(seq(0,1,0.1))

for(n in n_vec){
  for(rank in rank_vec){
    for(percent_missing in percent_missing_vec){
      for(sigma in sigma_vec){
        sim_dat<-gen_data(n = n,rank = rank,percent_missing = percent_missing, sigma = sigma) 
        for(lambda in lambda_vec){
          for(i in 1:k){
            results[i]<-cv_matrixcmplt_ggd_sqr(theta = sim_dat$mu,
                                            lambda_sequence = lambda,
                                            data = sim_dat$Y,
                                            mask = sim_dat$mask,
                                            true = sim_dat$mu,
                                            thresh = 1e-6,
                                            nfold = 10)$error
          }
          sim_results<-rbind(sim_results,
                             data.frame(n=n,
                                        rank=rank,
                                        percent_missing=percent_missing,
                                        sigma=sigma,
                                        SNR=sim_dat$SNR,
                                        lambda=lambda,
                                        error=mean(results)))
          save(sim_results,file = "./Results/sqr_results.Rdata")
        }
      }
    }
  }
}


#library(data.table)
#sim_results<-as.data.table(sim_results)
#sim_results[,.(optimal.lambda=lambda[which.min(error)]),by=SNR]
