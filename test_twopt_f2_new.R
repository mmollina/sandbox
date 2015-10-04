require(onemap)
require(Rcpp)
data("fake.f2.onemap")
system.time(x<-rf.2pts(fake.f2.onemap,verbose = FALSE))
sourceCpp("cpp/twopt_est_f2.cpp")
system.time(y<-est_rf_f2(x=fake.f2.onemap$geno, n = fake.f2.onemap$n.ind))

geno<-NULL
nmar<-1000
nind<-250
for(i in 1:nmar)
{
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.25,.5,.25)))
}
z<-system.time(y<-est_rf_f2(x=geno, n = nind))

speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000,add = TRUE, col=3)
speed(20000)

