require(onemap)
require(Rcpp)
require(qtl)
sourceCpp("cpp/twopt_est_out.cpp")
sourceCpp("cpp/twopt_est_f2.cpp")
sourceCpp("cpp/twopt_est_bc.cpp")
############################
geno<-NULL
nmar<-1000
nind<-250
for(i in 1:nmar)
{
  j<-sample(1:3,1)
  if(j==1)
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.25,.5,.25)))
  else if(j==2)
    geno<-cbind(geno,sample(c(1,5), size=nind, replace=TRUE, prob = c(.25,.75)))
  else if(j==3)
    geno<-cbind(geno,sample(c(3,4), size=nind, replace=TRUE, prob = c(.25,.75)))
}
z<-system.time(y<-est_rf_f2(x=geno, n = nind))
speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000, col="red", lwd=2, xlab = "number of markers", ylab = "time in seconds")
############################
geno<-NULL
segrega.type<-sample(c(1:7), size=1000, replace=TRUE)
for(i in segrega.type)
{
  if(i==1)
    geno<-cbind(geno,sample(1:4, size=250, replace=TRUE))
  else if(i==2)
    geno<-cbind(geno,sample(1:3, size=250, replace=TRUE, prob = c(.5,.25,.25)))
  else if(i==3)
    geno<-cbind(geno,sample(1:3, size=250, replace=TRUE, prob = c(.5,.25,.25)))
  else if(i==4)
    geno<-cbind(geno,sample(1:3, size=250, replace=TRUE, prob = c(.25,.5,.25)))
  else if(i==5)
    geno<-cbind(geno,sample(1:2, size=250, replace=TRUE, prob = c(.75,.25)))
  else if(i==6)
    geno<-cbind(geno,sample(1:2, size=250, replace=TRUE, prob = c(.5,.5)))
  else if(i==7)
    geno<-cbind(geno,sample(1:2, size=250, replace=TRUE, prob = c(.5,.5)))
}
z<-system.time(y<-est_rf_out(x=geno, segreg_type = segrega.type, n = 250))
speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000, col="blue", lwd=2, add=TRUE)
############################
geno<-NULL
nmar<-1000
nind<-250
for(i in 1:nmar)
{
  geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE))
}
z<-system.time(y<-est_rf_bc(x=geno, n_ind = nind))

speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000,add = TRUE, col="green", lwd=2)
############################
legend("topleft", legend = c("F2", "Outcross", "Backcross"),
       lty = 1:2, xjust = 1, yjust = 1,
        col=c(2,4,3))

