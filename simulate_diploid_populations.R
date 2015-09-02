#ch.len<-3000
#n.mrk<-250000
#r<-mf.k(ch.len/n.mrk)
#n.ind<-2000
#dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
#dat<-dat+1
#dat[sample(1:length(dat), length(dat)*.05)]<-NA
#colnames(dat)<-paste("M", 1:n.mrk, sep="")
#dat<-dat[,sample(1:n.mrk)]
#save.image("big_example_250000_mrk_2000_ind.RData")


#dat.back<-dat
#pdf("no_bins.pdf")
#par(bg="gray")
#image(dat.back, col=c(2,3,4))
#dev.off()

#system.time(bins.slow<-find.bins(dat, n.cpus=4, exact=FALSE))
#pdf("bins.slow.pdf")
#par(bg="gray")
#image(dat.back[,sort(names(bins.slow))], col=c(2,3,4))
#dev.off()

#dat[is.na(dat)]<-0
#system.time(bins.fast<-find.bins(dat[1:200,1:31250], exact=FALSE))
#pdf("bins.fast.pdf")
#par(bg="gray")
#image(dat.back[,sort(names(bins.fast))], col=c(2,3,4))
#dev.off()

#save.image("bins.RData")


#n.cpu<-8
#n.ind<-nrow(dat)
#n.mrk<-ncol(dat)
#n.sam<-200
#v.mrk<-1:n.mrk
#v.mrk<-split(v.mrk, ceiling(seq_along(v.mrk)/ceiling(n.mrk/n.cpu)))
#ind<-sample(1:n.ind, n.sam)


#find.bins.par<-function(x, dat, ind, exact)
#    {
#        bins.temp<-find.bins(dat[ind,x], exact=FALSE)
#        bins.temp
#    }

    
#fun(arg1, arg2, arg3)

#res <- mclapply(vec, fun, arg1=x, arg2=y, arg3=z, mc.cores=n.cluster)


source("utils.R")




bla<-check.bin.R(dat)
ble<-find.bins.slow.no.cpp.code(dat, n.cpus=4, exact=TRUE)




ch.len<-300
n.mrk<-10000
r<-mf.k(ch.len/n.mrk)
n.ind<-1000
dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
dat<-dat+1
colnames(dat)<-paste("M", 1:n.mrk, sep="")
dat<-dat[,sample(1:n.mrk)]

system.time(bli<-find.bins(dat, TRUE))
v.bli<-numeric(n.mrk)
for(i in 1: length(bli))
    v.bli[match(names(bli[[i]]),colnames(dat))]<-i


all(blo==v.bli)
