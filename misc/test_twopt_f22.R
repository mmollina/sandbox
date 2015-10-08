require(onemap)
require(Rcpp)
require(qtl)

file<-paste(system.file("example",package="onemap"),"fake.f2.onemap.raw", sep="/")
dat1 <- read.cross("mm", file=file, mapfile="fake.f2.onemap.map")
data(fake.f2.onemap)
fake.f2.onemap$geno.mmk[[1]][[1]][is.na(fake.f2.onemap$geno.mmk[[1]][[1]])]<-0
dat1<-est.rf(dat1, tol =  10e-10)
r.qtl<-dat1$rf
colnames(r.qtl)
diag(r.qtl)<-NA

geno<-fake.f2.onemap$geno.mmk[[1]][,colnames(r.qtl)]
geno[is.na(geno)]<-0

check.type<-function(x){
 a<-paste(sort(unique(x[!x==0])), collapse = "")
 if(a=="123") return(1)
 else if(a=="34") return(2)
 else if (a=="15") return(3)
 else return(4)
}

sourceCpp("cpp/twopt_est_f2.cpp")
type<-apply(geno, 2, check.type)
type
y1<-est_rf_f2(x=as.numeric(geno), type = type,  n = fake.f2.onemap$n.ind)

image(round(r.qtl-y1,2))

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

