require(onemap)
require(Rcpp)
require(qtl)

dat1 <- read.cross("mm", file=file, mapfile="fake.f2.onemap.map")
data(fake.f2.onemap)
fake.f2.onemap$geno.mmk[[1]][[1]][is.na(fake.f2.onemap$geno.mmk[[1]][[1]])]<-0
dat1<-est.rf(dat1, tol =  0.00001)
r.qtl<-dat1$rf
colnames(r.qtl)
diag(r.qtl)<-NA

sourceCpp("cpp/twopt_est_f2.cpp")
sourceCpp("cpp/twopt_est_f2_backup.cpp")

z1<-z2<-numeric(100)
for(i in 1:100)
{
 z1[i]<-system.time(y1<-est_rf_f2(x=as.numeric(fake.f2.onemap$geno.mmk[[1]][,colnames(r.qtl)]), n = fake.f2.onemap$n.ind))[3]
 z2[i]<-system.time(y2<-est_rf_f2_backup(x=as.numeric(fake.f2.onemap$geno.mmk[[1]][,colnames(r.qtl)]), n = fake.f2.onemap$n.ind))[3]
}

summary(z1)
summary(z2)

round(r.qtl-y1,5)
round(r.qtl-y2,5)








geno<-NULL
nmar<-1000
nind<-250
for(i in 1:nmar)
{
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.25,.5,.25)))
}
z<-system.time(y<-est_rf_f2_new(x=geno, n = nind))

speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000,add = TRUE, col=4, lwd=2)
speed(20000)

