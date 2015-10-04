require(onemap)
require(Rcpp)
data(example.out)
x<-rf.2pts(example.out,verbose = FALSE)
s<-example.out$segr.type
names(s)<-colnames(example.out$geno)

##----------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------
##M1  M2  M3  M4  M5  M6  M7  M8  M9 M10 M11 M12 M13 M14 M15 M16 M17 M18 M19 M20 M21 M22 M23 M24 M25 M26 M27 M28 M29 M30 
##4   7   6   1   7   4   7   4   6   7   7   1   5   1   1   7   3   1   2   1   7   6   5   4   3   1   6   1   6   4 
##----------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------

sourceCpp("cpp/twopt_est_out.cpp")
y<-est_rf_out_new(x=example.out$geno, segreg_type = example.out$segr.type.num, n = example.out$n.ind)
a<-print(x, mrk1="M12", mrk2 = "M14")
(b<-t(sapply(y, function(x) c(x[14,12], x[12,14]))))
round(a-b,5)

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
z<-system.time(y<-est_rf_out_new(x=geno, segreg_type = segrega.type, n = 250))

speed<-function(nmar) choose(nmar,2) * z[3]/choose(1000,2)
curve(speed, 10, 10000)
speed(20000)

