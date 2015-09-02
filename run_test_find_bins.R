source("utils.R")
sourceCpp("find_bins.cpp")

##Firts test 500 individuals, 50000 markers
##ch.len<-300
##n.mrk<-50000
##r<-mf.k(ch.len/n.mrk)
##n.ind<-2000
##dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
##dat<-dat+1
##dat[sample(1:length(dat), length(dat)*.05)]<-0
##colnames(dat)<-paste("M", 1:n.mrk, sep="")

##456.6 seconds
##system.time(res.1<-getbins(dat))

## 1397.5 seconds
##dat<-dat[,sample(1:n.mrk)]
##system.time(res.2<-getbins(dat))


##Firts test 500 individuals, 50000 markers
ch.len<-3000
n.mrk<-250000
r<-mf.k(ch.len/n.mrk)
n.ind<-2000
dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
dat<-dat+1
dat[sample(1:length(dat), length(dat)*.05)]<-0
colnames(dat)<-paste("M", 1:n.mrk, sep="")
save.image("big_example_250000_mrk_2000_ind.RData")

system.time(res.3<-getbins(dat))

dat<-dat[,sample(1:n.mrk)]
system.time(res.4<-getbins(dat))
