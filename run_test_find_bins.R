source("utils.R")
sourceCpp("find_bins.cpp")
ch.len<-300
n.mrk<-20000
r<-mf.k(ch.len/n.mrk)
n.ind<-500
dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
dat<-dat+1
dat[sample(1:length(dat), length(dat)*.05)]<-0
dat<-dat[,sample(1:n.mrk)]

system.time(res.3<-getbins(dat))
fun(arg1, arg2, arg3)
res <- mclapply(vec, fun, arg1=x, arg2=y, arg3=z, mc.cores=n.cluster)





dat<-dat[,sample(1:n.mrk)]
system.time(res.4<-getbins(dat))
