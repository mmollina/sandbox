require(snow)
##Kosambi mapping function
mf.k<-function (d) 
    0.5 * tanh(d/50)

##Inverse of Kosambi mapping function
imf.k<-function (r)
{
    r[r >= 0.5] <- 0.5 - 1e-14
    50 * atanh(2 * r)
}

##Simulates one diploid gamete
sim.gam<-function(n.mrk, r, init=0)
    {
        g<-rep(init, n.mrk)
        for(i in 2:n.mrk)
            g[i]<-sample(c(g[i-1],0^g[i-1]), size=1, prob=c(1-r, r))
        g
    }

##Simulates a F2 population for one chomossome 
sim.ch.f2<-function(n.ind, n.mrk, ch.len)
    {
        dat<-matrix(0,n.ind, n.mrk)
        pb <- txtProgressBar(min = 1, max = n.ind, style = 3)
        for(i in 1:n.ind)
            {
                setTxtProgressBar(pb, i)
                init<-sample(c(0,1), size=1)
                dat[i,]<-sim.gam(n.mrk,r,init)+sim.gam(n.mrk,r,init)
            }
        close(pb)
        dat
    }

##Find clusters
find.bins<-function(dat, n.cpus, exact=FALSE)
    {
        n.mrk<-nrow(dat)
        cl <- makeCluster(n.cpus)
        on.exit(stopCluster(cl))
        mis<-parCapply(cl, dat, function(x) sum(is.na(x)))
        bins<-vector("list", 1)
        bt.mrk<-character()
        count<-1
        while(length(ncol(dat)) > 0)
            {
                cat("\n bin number: ", count, " --- remaining markers:", ncol(dat))
                a<-dat[,1]
                clusterExport(cl,"a")
                mrk.a<-colnames(dat)[1]
                if(class(dat[,-1])=="numeric")
                    {
                        M<-matrix(dat[,-1], nrow = n.mrk)
                        if(exact)
                            aa <- apply(M, 2, function(x) identical(x,a))
                        else
                            aa <- apply(M, 2, function(x) all(x==a, na.rm=TRUE))
                    }
                else
                    {
                        if(exact)
                            aa <- parCapply(cl, dat[,-1], function(x) identical(x,a))
                        else
                            aa <- parCapply(cl, dat[,-1], function(x) all(x==a, na.rm=TRUE))
                    }
                bins[[count]]<-mis[c(mrk.a,names(which(aa)))]
                bt.mrk[count]<-names(which.min(bins[[count]]))
                b<-match(names(bins[[count]]), colnames(dat))
                count=count+1
                dat<-dat[,-b]
            }
        names(bins)<-bt.mrk
        bins
    }



ch.len<-200
n.mrk<-10000
r<-mf.k(ch.len/n.mrk)
n.ind<-100
dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
dat[sample(1:length(dat), length(dat)*.05)]<-NA
colnames(dat)<-paste("M", 1:n.mrk, sep="")
dat<-dat[,sample(1:n.mrk)]

dat.back<-dat
par(bg="gray")
image(dat.back, col=c(2,3,4))

system.time(
    bins<-find.bins(dat, n.cpus=4, exact=FALSE)
    )

x11()
par(bg="gray")
image(dat.back[,names(bins)])
 
