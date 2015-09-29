require(onemap)
require(Rcpp)
sourceCpp("two_pt_test.cpp")
data(fake.f2.onemap)
twopt<-rf.2pts(fake.f2.onemap)
lg<-group(make.seq(twopt, "all"))

##"pre-allocate" an empty list of length lg$n.groups (3, in this case)
maps.list<-vector("list", lg$n.groups)
for(i in 1:lg$n.groups){
  ##create linkage group i
  LG.cur <- make.seq(lg,i)
  ##ordering
  map.cur<-record(LG.cur)
  ##assign the map of the i-th group to the maps.list
  maps.list[[i]]<-map.cur
}


##write maps.list to "fake.f2.onemap.map" file
write.map(maps.list, "fake.f2.onemap.map")
require(qtl)
sourceCpp("two_pt_test.cpp")
file<-paste(system.file("example",package="onemap"),"fake.f2.onemap.raw", sep="/")

for(i1 in 1:10)
{
dat1 <- read.cross("mm", file=file, mapfile="fake.f2.onemap.map")
data(fake.f2.onemap)
fake.f2.onemap$geno.mmk[[1]][[1]][is.na(fake.f2.onemap$geno.mmk[[1]][[1]])]<-0
dat1<-est.rf(dat1)
which(fake.f2.onemap$segr.type=="B3.7")
(r.new<-round(est_rf(as.numeric(fake.f2.onemap$geno.mmk[[1]][, c(1,30)])),5))
r.qtl<-(round(dat1$rf,5))
r.qtl["M30","M1"]
cat("\n------------------------------------\n")
cat(r.qtl["M30","M1"]==r.new[2])
temp<-dat1$geno$`1`$data[,"M30"]
id<-which(temp!=3 & !is.na(temp))
dat1$geno$`1`$data[sample(id,20),"M30"]<-4
temp<-dat1$geno$`1`$data[,"M30"]
id<-which(temp!=1 & temp!=4 & !is.na(temp))
dat1$geno$`1`$data[sample(id,20),"M30"]<-5
fake.f2.onemap$geno.mmk[[1]][, 30]<-dat1$geno$`1`$data[,"M30"]
fake.f2.onemap$geno.mmk[[1]][is.na(fake.f2.onemap$geno.mmk[[1]][, 30]), 30]<-0
dat1<-est.rf(dat1)
r.qtl<-(round(dat1$rf,5))
(r.new<-round(est_rf(as.numeric(fake.f2.onemap$geno.mmk[[1]][, c(1,30)])),5))
cat(r.qtl["M30","M1"]==r.new[2])
temp<-dat1$geno$`1`$data[,"M1"]
id<-which(temp!=3 & !is.na(temp))
dat1$geno$`1`$data[sample(id,20),"M1"]<-4
temp<-dat1$geno$`1`$data[,"M1"]
id<-which(temp!=1 & temp!=4 & !is.na(temp))
dat1$geno$`1`$data[sample(id,20),"M1"]<-5
fake.f2.onemap$geno.mmk[[1]][, 1]<-dat1$geno$`1`$data[,"M1"]
fake.f2.onemap$geno.mmk[[1]][is.na(fake.f2.onemap$geno.mmk[[1]][, 1]), 1]<-0
dat1<-est.rf(dat1)
r.qtl<-(round(dat1$rf,5))
(r.new<-round(est_rf(as.numeric(fake.f2.onemap$geno.mmk[[1]][, c(1,30)])),5))
cat(r.qtl["M30","M1"]==r.new[2], "\n------------------------------------\n")
}


require(onemap)
require(Rcpp)
load(url("https://github.com/mmollina/onemap/blob/master/fake.big.data.f2.RData?raw=true"))
(bins<-find.bins(fake.big.data.f2, exact=FALSE))
(new.data<-create.data.bins(fake.big.data.f2, bins))

dim(geno<-fake.big.data.f2$geno.mmk[[1]])
geno.t<-geno[,substring(colnames(geno),4,4)=="1"]
colnames(geno.t)
#geno.s<-geno.t[,sample(1:ncol(geno.t), 5000)]
geno.s<-geno.t
#geno.s<-geno.s[,order(as.numeric(substring(colnames(geno.s), 6)))]
dim(geno.s)
sourceCpp("two_pt_test_updated.cpp")
system.time(r.new2<-est_rf(as.numeric(geno.s), 250))






