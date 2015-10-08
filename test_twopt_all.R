require(onemap)
require(Rcpp)
require(qtl)
source("general_simulation.R")
sourceCpp("cpp/twopt_est_out.cpp")
sourceCpp("cpp/twopt_est_f2.cpp")
sourceCpp("cpp/twopt_est_bc.cpp")


z.bc<-system.time(y.bc<-est_rf_bc(x=dat.bc, n = n.ind))
mat.bc<-as.matrix(as.dist(y.bc, upper = TRUE))
speed.bc<-function(w) choose(w,2) * z.bc[3]/choose(n.mar,2)

##F2
dat.f2<-sim.ch.f2(n.ind, n.mrk, ch.len, dom43 = 20, dom51 =20, missing = 10)
dat.temp<-dat.f2
dat.temp[dat.temp==0]<-NA
check.type<-function(x)
{
    if(sum(is.na(match(na.omit(unique(x)), c(1,2,3))))==0) return(1)
    else if (sum(is.na(match(na.omit(unique(x)), c(4,3))))==0) return(2)
    else if (sum(is.na(match(na.omit(unique(x)), c(5,1))))==0) return(3)
    else return(4)
}
type<-apply(dat.temp, 2, check.type)
z.f2<-system.time(y.f2<-est_rf_f2(x=dat.f2, type = type, n = n.ind))
mat.f2<-as.matrix(as.dist(y.f2, upper = TRUE))
speed.f2<-function(w) choose(w,2) * z.f2[3]/choose(n.mar,2)

##Backcross
dat.out<-sim.ch.out(n.ind, n.mrk, ch.len, mis, prob=c(1,1,1,1,1,1,1))
z.out<-system.time(y.out<-est_rf_out(x=dat.out[[1]], segreg_type =  dat.out[[2]], n = n.ind))
mat.out<-as.matrix(as.dist(y.out[[1]], upper = TRUE))
speed.out<-function(w) choose(w,2) * z.bc[3]/choose(n.mar,2)


image(mat.bc, axes=FALSE, col=rainbow(n=500, start=min(mat.bc,na.rm=TRUE)*1.3, end=max(mat.bc,na.rm=TRUE)*1.3))
image(mat.f2, axes=FALSE, col=rainbow(n=500, start=min(mat.f2,na.rm=TRUE)*1.3, end=max(mat.f2,na.rm=TRUE)*1.3))
image(mat.out[1:100,1:100], axes=FALSE, col=rainbow(n=500, start=min(mat.out,na.rm=TRUE)*1.3, end=max(mat.out,na.rm=TRUE)*1.3))
curve(speed, 10, 10000, col="red", lwd=2, xlab = "number of markers", ylab = "time in seconds")





############################
nmar<-1000
nind<-250
geno<-NULL
segrega.type<-sample(c(1:7), size=nmar, replace=TRUE)
for(i in segrega.type)
{
  if(i==1)
    geno<-cbind(geno,sample(1:4, size=nind, replace=TRUE))
  else if(i==2)
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.5,.25,.25)))
  else if(i==3)
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.5,.25,.25)))
  else if(i==4)
    geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.25,.5,.25)))
  else if(i==5)
    geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE, prob = c(.75,.25)))
  else if(i==6)
    geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE, prob = c(.5,.5)))
  else if(i==7)
    geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE, prob = c(.5,.5)))
}
z<-system.time(y<-est_rf_out(x=geno, segreg_type = segrega.type, n = nind))
speed<-function(w) choose(w,2) * z[3]/choose(nmar,2)
curve(speed, 10, 10000, col="red", lwd=2, xlab = "number of markers", ylab = "time in seconds")
############################
type<-NULL
geno<-NULL
for(i in 1:nmar)
{
    j<-sample(1:3,1)
    type<-c(type,j)
    if(j==1)
      geno<-cbind(geno,sample(1:3, size=nind, replace=TRUE, prob = c(.25,.5,.25)))
    else if(j==2)
      geno<-cbind(geno,sample(c(3,4), size=nind, replace=TRUE, prob = c(.25,.75)))
    else if(j==3)
      geno<-cbind(geno,sample(c(1,5), size=nind, replace=TRUE, prob = c(.25,.75)))
}
z<-system.time(y<-est_rf_f2(x=geno, type=type, n_ind = nind))
speed<-function(w) choose(w,2) * z[3]/choose(nmar,2)
curve(speed, 10, 10000, col="orange", lwd=2, add=TRUE)

############################
geno<-NULL

for(i in 1:nmar)
{
  geno<-cbind(geno,sample(1:2, size=nind, replace=TRUE))
}
z<-system.time(y<-est_rf_bc(x=geno, n_ind = nind))
speed<-function(w) choose(w,2) * z[3]/choose(nmar,2)
curve(speed, 10, 10000,add = TRUE, col="blue", lwd=2)
############################
legend("topleft", legend = c("Outcross","F2",  "Backcross"),
       lty = 1, xjust = 1, yjust = 1,
        col=c("red", "orange", "blue"))

