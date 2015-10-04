require(onemap)
require(Rcpp)
data("fake.bc.onemap")
system.time(x<-rf.2pts(fake.bc.onemap,verbose = FALSE))
fake.bc.onemap$geno.mmk
sourceCpp("cpp/twopt_est_bc.cpp")
system.time(y<-est_rf_bc(x=fake.bc.onemap$geno, n_ind = fake.bc.onemap$n.ind))

geno<-NULL
nmar<-1000
nind<-250
for(i in 1:nmar)
{
    geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE))
}
z<-system.time(y<-est_rf_bc(x=geno, n_ind = nind))

speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000,add = TRUE, col=2)
speed(20000)

