source("utils_for_test.R")
require(onemap)

n.chom<-5
ch.len<-250
n.mrk<-1000
r<-mf.k(ch.len/n.mrk)
n.ind<-250

geno<-NULL
for(i in 1:n.chom){
  dat<-sim.ch.f2(n.ind, n.mrk, ch.len, r)
  dat<-dat+1
  dat[sample(1:length(dat), length(dat)*.05)]<-0
  colnames(dat)<-paste("ch_", i, "M", 1:n.mrk, sep="")
  geno<-cbind(geno, dat)
}
#geno<-geno[,sample(1:ncol(geno))]
data("fake.f2.onemap")
nm<-names(fake.f2.onemap)
fake.big.data.f2<-structure(vector(mode="list", length(nm)), class=class(fake.f2.onemap))
names(fake.big.data.f2)<-nm
fake.big.data.f2$geno<-geno
temp<-geno
temp[temp==0]<-NA
fake.big.data.f2$geno.mmk<-list(geno=temp, type="f2")
fake.big.data.f2$n.ind<-nrow(geno)
fake.big.data.f2$n.mar<-ncol(geno)
fake.big.data.f2$segr.type<-rep("B3.7", ncol(geno))
fake.big.data.f2$segr.type.num<-rep(4, ncol(geno))
fake.big.data.f2$phase<-rep(1, ncol(geno))
fake.big.data.f2$input<-"none.txt"
fake.big.data.f2$n.phe<-0
fake.big.data.f2$pheno<-NA
fake.big.data.f2
plot(fake.big.data.f2)

system.time(x<-find.bins(fake.big.data.f2, exact=FALSE))

x[[1]][1:3]
new.data<-create.data.bins(fake.big.data.f2, x)
plot(new.data)

twopts.f2 <- rf.2pts(new.data)

input.seq<-make.seq(twopts.f2, "all")
markers <- length(input.seq$seq.num)
## create reconmbination fraction matrix
r <- matrix(NA,markers,markers)
for(i in 1:(markers-1)) {
  for(j in (i+1):markers) {
    big <- pmax.int(input.seq$seq.num[i],input.seq$seq.num[j])
    small <- pmin.int(input.seq$seq.num[i],input.seq$seq.num[j])
    temp <- get(input.seq$twopt)$analysis[acum(big-2)+small,,]
    ## check if any assignment meets the criteria
    relevant <- which(temp[,2] > (max(temp[,2])-0.005)) # maximum LOD scores
    phases <- relevant[which((temp[relevant,1] <= max.rf & temp[relevant,2] >= LOD))]
    if(length(phases) == 0) r[i,j] <- r[j,i] <- 0.5
    else r[i,j] <- r[j,i] <- temp[phases[1],1]
  }
}
image(r)



load(url("https://github.com/mmollina/onemap/blob/master/fake.big.data.f2.RData?raw=true"))
fake.big.data.f2
(bins<-find.bins(fake.big.data.f2, exact=FALSE))
(new.data<-create.data.bins(fake.big.data.f2, bins))
new.data$geno[new.data$geno==0]<-NA


