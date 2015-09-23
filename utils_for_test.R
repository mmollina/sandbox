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
sim.ch.f2<-function(n.ind, n.mrk, ch.len, r)
    {
        dat<-matrix(0,n.ind, n.mrk)
        pb <- txtProgressBar(min = 1, max = n.ind, style = 3)
        for(i in 1:n.ind)
            {
                setTxtProgressBar(pb, i)
                dat[i,]<-sim.gam(n.mrk,r,sample(c(0,1), size=1))+sim.gam(n.mrk,r,sample(c(0,1), size=1))
            }
        close(pb)
        dat
}






##Find clusters
find.bins.slow.no.cpp.code<-function(w, n.cpus, exact=FALSE)
    {
        w.temp<-w
        n.mrk<-ncol(w.temp)
        cl <- makeCluster(n.cpus)
        on.exit(stopCluster(cl))
        mis<-parCapply(cl, w.temp, function(x) sum(is.na(x)))
        bins<-vector("list", 1)
        bt.mrk<-character()
        count<-1
        while(length(ncol(w.temp)) > 0 && ncol(w.temp) > 0)
            {
#                cat("\n bin number: ", count, " --- remaining markers:", ncol(w.temp))
                a<-w.temp[,1]
                #clusterExport(cl,"a")
                mrk.a<-colnames(w.temp)[1]
                if(class(w.temp[,-1])=="numeric")
                    {
                        M<-matrix(w.temp[,-1], nrow = n.mrk)
                        if(exact)
                            aa <- apply(M, 2, function(x,a) identical(x,a), a=a)
                        else
                            aa <- apply(M, 2, function(x,a) all(x==a, na.rm=TRUE), a=a)
                    }
                else
                    {
                        if(exact)
                            aa <- parCapply(cl, w.temp[,-1], function(x,a) identical(x,a), a=a)
                        else
                            aa <- parCapply(cl, w.temp[,-1], function(x,a) all(x==a, na.rm=TRUE), a=a)
                    }
                bins[[count]]<-mis[c(mrk.a,names(which(aa)))]
                bt.mrk[count]<-names(which.min(bins[[count]]))
                b<-match(names(bins[[count]]), colnames(w.temp))
                count=count+1
                #print(-b)
                w.temp<-w.temp[,-b]
            }
        names(bins)<-bt.mrk
        bins
    }


find.bins.par<-function(x, dat, ind, exact)
    {
        bins.temp<-find.bins(dat[ind,x], exact=FALSE)
        bins.temp
    }


check.occourance<-function(v,x)
    {
        for(i in 1:length(v))
            if(v[i]==x)
                return(i)
        return(-1)
    }
check.bin.R<-function(dat)
    {
        b.vec<-1
        b<-rep(1,ncol(dat))
        for(i in 1:ncol(dat))
            {
                for(j in 1:length(b.vec))
                    {
                        flag<-0
                         l<-check.occourance(b, b.vec[j])
                        for(k in 1:nrow(dat))
                            {

                                if(dat[k,i]!=dat[k, l])
                                    {
                                        flag=1
                                        break
                                    }
                            }
                        if(flag==0)
                            {
                                b[i]<-b.vec[j]
                                break
                            }
                    }
                if(flag==1)
                    {
                        b[i]<-max(b.vec)+1
                        b.vec<-c(b.vec, b[i])
                    }
            }
        return(b)
    }
