library(ggplot2)
library(data.table)
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_lasso_n50_results.Rdata")
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_lasso_n75_results.Rdata")
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_lasso_n100_results.Rdata")
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_sqrtlasso_n50_results.Rdata")
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_sqrtlasso_n75_results.Rdata")
load("/Users/Isaac/Dropbox/server run/matrix-complete/Results/oracle_sqrtlasso_n100_results.Rdata")

oracle_lasso<-rbind(oracle_lasso_n50,oracle_lasso_n75,oracle_lasso_n100)

oracle_lasso<-as.data.table(oracle_lasso)
data<-oracle_lasso
data<-data[,.(optimal.lambda=lambda[which.min(error)]),by=.(n,rank,percent_missing,sigma)]

oracle_srqtlasso<-rbind(oracle_srqtlasso_n50,oracle_srqtlasso_n75,oracle_srqtlasso_n100)
data5<-oracle_srqtlasso
data_5<-data5[,.(optimal.lambda=lambda[which.min(error)]),by=.(n,rank,percent_missing,sigma)]

oracle_srqtlasso_n50<-as.data.table(oracle_srqtlasso_n50)
data_2<-oracle_srqtlasso_n50
data_2<-data_2[,.(optimal.lambda=lambda[which.min(error)]),by=.(n,rank,percent_missing,sigma)]

oracle_srqtlasso_n75<-as.data.table(oracle_srqtlasso_n75)
data_3<-oracle_srqtlasso_n75
data_3<-data_3[,.(optimal.lambda=lambda[which.min(error)]),by=.(n,rank,percent_missing,sigma)]

oracle_srqtlasso_n100<-as.data.table(oracle_srqtlasso_n100)
data_4<-oracle_srqtlasso_n100
data_4<-data_4[,.(optimal.lambda=lambda[which.min(error)]),by=.(n,rank,percent_missing,sigma)]

#data_2
ggplot(data = data_2[sigma != 0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() +
  geom_point(data = data[n==50 & sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_line(data = data[n==50& sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  ggtitle("Lasso & Square-root Lasso n=50 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)

ggplot(data = data_2[sigma !=0, ],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() + ylim(0,1)+
  ggtitle("Square-root Lasso MC \nn=50 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)


#data_3
ggplot(data = data_3[sigma != 0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() +
  geom_point(data = data[n==75 & sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_line(data = data[n==75& sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  ggtitle("Lasso & Square-root Lasso n=75 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)

ggplot(data = data_3[sigma !=0, ],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() + ylim(0,1)+
  ggtitle("Square-root Lasso MC \nn=75 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)


#data_4
ggplot(data = data_4[sigma != 0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() +
  geom_point(data = data[n==100 & sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_line(data = data[n==100& sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  ggtitle("Lasso & Square-root Lasso n=100 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)

ggplot(data = data_4[sigma !=0, ],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line() + ylim(0,1)+
  ggtitle("Square-root Lasso MC \nn=100 \nFacet by missing proportions \nColor by rank") +facet_wrap(~percent_missing)



#data


ggplot(data = data[percent_missing==0.01,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ggtitle("Lasso\nfacet by rank\ncolor by n\npercent of missing=1%") +facet_wrap(~rank)

ggplot(data = data[percent_missing==0.05,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ggtitle("Lasso\nfacet by rank\ncolor by n\npercent_missing=5%") +facet_wrap(~rank)

ggplot(data = data[percent_missing==0.10,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ggtitle("Lasso\nfacet by rank\ncolor by n\npercent_missing=10%") +facet_wrap(~rank)

ggplot(data = data[n==50 & sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line()+
  ggtitle("Lasso\nfacet by missing\ncolor by rank\nn=50") +facet_wrap(~percent_missing)

ggplot(data = data[n==50 &sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(percent_missing)))+
  geom_point() +
  geom_line()+
  ggtitle("Lasso\nfacet by rank\ncolor by missing\nn=50") +facet_wrap(~rank)

#data_5
ggplot(data = data_5[percent_missing==0.01 &sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ylim(0,1)+
  ggtitle("Square-root Lasso\nfacet by rank\ncolor by n\npercent of missing=1%") +facet_wrap(~rank)

ggplot(data = data_5[percent_missing==0.05 &sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ylim(0,1)+
  ggtitle("Square-root Lasso\nfacet by rank\ncolor by n\npercent_missing=5%") +facet_wrap(~rank)

ggplot(data = data_5[percent_missing==0.1 &sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(n)))+
  geom_point() +
  geom_line()+
  ylim(0,1)+
  ggtitle("Square-root Lasso\nfacet by rank\ncolor by n\npercent_missing=10%") +facet_wrap(~rank)

ggplot(data = data_5[n==50 & sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(rank)))+
  geom_point() +
  geom_line()+
  ylim(0,1)+
  ggtitle("Square-root Lasso\nfacet by missing\ncolor by rank\nn=50") +facet_wrap(~percent_missing)

ggplot(data = data_5[n==50 &sigma !=0,],aes(x = sigma, y = optimal.lambda,color=as.factor(percent_missing)))+
  geom_point() +
  geom_line()+
  ylim(0.1,0.3)+
  ggtitle("Square-root Lasso\nfacet by rank\ncolor by missing\nn=50") +facet_wrap(~rank)
