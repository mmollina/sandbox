require(snow)

##Rcpp example
library(Rcpp)
library(inline)

code <- '
Rcpp::NumericMatrix geno(genoR);
int exact = Rcpp::as<int>(exactR);
std::vector<int> v(1);
int n_mar = geno.ncol();
int n_ind = geno.nrow();
if(exact==1)
  {
    for(int i = 1; i < n_mar; i++)
      {
	int flag=0;
	for(int j = 0; j < n_ind; j++)
	  {
	    if(geno(j,0)!=geno(j,i))
	      {
		flag=1;
		break;
	      }
	  }
	if(flag==0)
	  v.push_back(i);
      }
  }
 else
  {
    for(int i = 1; i < n_mar; i++)
      {
	int flag=0;
	for(int j = 0; j < n_ind; j++)
	  {
	    if(geno(j,0)!=geno(j,i) && geno(j,0)!=0 && geno(j,i)!=0)
	      {
		flag=1;
		break;
	      }
	  }
	if(flag==0)
	  v.push_back(i);
      }
  }
return(wrap(v));
'

comp.vec <- cxxfunction(signature(genoR="numeric", exactR="numeric"),
                         plugin="Rcpp",
                         body=code)

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
find.bins<-function(w, n.cpus, exact=FALSE)
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
                cat("\n bin number: ", count " --- remaining markers:", ncol(w.temp))
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



##Find clusters
find.bins.fast<-function(w, n.cpus, exact=FALSE)
    {
        w.temp<-w
        n.mrk<-ncol(w.temp)
        #cl <- makeCluster(n.cpus)
        ##on.exit(stopCluster(cl))
        mis<-apply(w.temp, 2, function(x) sum(x==0))
        #stopCluster(cl)
        bins<-vector("list", 1)
        bt.mrk<-character()
        count<-1
        while(length(ncol(w.temp)) > 0 && ncol(w.temp) > 0)
            {
                #cat("\n bin number: ", count, " --- remaining markers:", ncol(w.temp))
                cat(".")
                if(exact)
                    aa <-comp.vec(w.temp, 1)
                else
                    aa <-comp.vec(w.temp, 0)
                bins[[count]]<-mis[aa+1]
                bt.mrk[count]<-names(which.min(bins[[count]]))
                count=count+1
                #print(-(aa+1))
                w.temp<-w.temp[,-(aa+1)]
                mis<-mis[-(aa+1)]
            }
        names(bins)<-bt.mrk
        bins
    }




ch.len<-300
n.mrk<-10000
r<-mf.k(ch.len/n.mrk)
n.ind<-100
dat<-sim.ch.f2(n.ind, n.mrk, ch.len)
dat<-dat+1
dat[sample(1:length(dat), length(dat)*.05)]<-NA
colnames(dat)<-paste("M", 1:n.mrk, sep="")
dat<-dat[,sample(1:n.mrk)]

dat.back<-dat
pdf("no_bins.pdf")
par(bg="gray")
image(dat.back, col=c(2,3,4))
dev.off()

system.time(bins.slow<-find.bins(dat, n.cpus=4, exact=FALSE))
pdf("bins.slow.pdf")
par(bg="gray")
image(dat.back[,sort(names(bins.slow))], col=c(2,3,4))
dev.off()


dat[is.na(dat)]<-0
system.time(bins.fast<-find.bins.fast(dat, n.cpus=4, exact=FALSE))
pdf("bins.fast.pdf")
par(bg="gray")
image(dat.back[,sort(names(bins.fast))], col=c(2,3,4))
dev.off()

save.image("bins.RData")



