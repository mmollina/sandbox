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
system.time(r.new2<-est_rf_f2(as.numeric(geno.s), 250))

require(onemap)

sourceCpp("twopt_est_out.cpp")
data("example.out")
id<-which(example.out$segr.type.num==1 | example.out$segr.type.num==4 )
system.time(twopt <- rf.2pts(example.out))
system.time(ble<-est_rf_out(x=example.out$geno, segreg_type = example.out$segr.type.num, n = example.out$n.ind))
print(twopt, mrk1="M4", mrk2 = "M6")
t(sapply(ble, function(x) c(x[6,4], x[4,6])))
for(i in 1:4){ dimnames(ble[[i]])<-list(example.out$segr.type, example.out$segr.type)}
ble[[1]]



require(Rcpp)
require(onemap)
data("example.out")
x<-rf.2pts(example.out)
print(x, "M4", "M19")
which(example.out$segr.type.num==1)

geno<-NULL
segrega.type<-sample(c(1,4), size=100, replace=TRUE)
for(i in segrega.type)
    {
        if(i==1)
            geno<-cbind(geno,sample(1:4, size=100, replace=TRUE))
        else
            geno<-cbind(geno,sample(1:3, size=100, replace=TRUE))
    }

segrega.type[c(4,19)]<-c(1,2)
geno[,c(4,19)]<-example.out$geno[,c(4,19)]

bli[[1]][19,4]

sourceCpp("twopt_est_out.cpp")
t1<-system.time(ble<-est_rf_out(x=geno, segreg_type = segrega.type, 250))
sourceCpp("twopt_est_out_update.cpp")
t2<-system.time(bli<-est_rf_out_new(x=geno, segreg_type = segrega.type, 250))
t1;t2
identical(lapply(ble, round, 3),lapply(bli, round,3))

all(round(ble[[1]][1:10,1:10],3) == round(bli[[1]][1:10,1:10],3))

ble[[1]][1:10,1:10] - bli[[1]][1:10,1:10]


blo<-NULL
for(i in 1:499)
{
  for( j in (i+1):500)
  {
    cat(i, "--", j, "\n")
    if(segrega.type[i]==1 &&  segrega.type[j]==4)
    {
      blo<-c(blo, all(round(bli[[1]][j,i]+bli[[4]][j,i],5)==1L,
          round(bli[[2]][j,i]+bli[[3]][j,i],5)==1L,
          round(bli[[2]][i,j],5)==round(bli[[3]][i,j],5),
          round(bli[[2]][i,j],5)==round(bli[[3]][i,j],5)))
    }
  }
}

all(blo)


require(onemap)
require(Rcpp)
sourceCpp("twopt_est_out_update.cpp")
data(example.out)
system.time(for(i in 1:100) ble<-est_rf_out_new(x=example.out$geno, segreg_type = example.out$segr.type.num, n = example.out$n.ind))
system.time(for(i in 1:100) x<-rf.2pts(example.out,verbose = FALSE))

print(x, mrk1="M4", mrk2 = "M19")
t(sapply(ble, function(x) c(x[19,4], x[4,19])))
print(x, mrk1="M4", mrk2 = "M17")
t(sapply(ble, function(x) c(x[17,4], x[4,17])))
print(x, mrk1="M13", mrk2 = "M14")
t(sapply(ble, function(x) c(x[14,13], x[13,14])))
print(x, mrk1="M3", mrk2 = "M4")
t(sapply(ble, function(x) c(x[4,3], x[3,4])))
print(x, mrk1="M4", mrk2 = "M5")
t(sapply(ble, function(x) c(x[5,4], x[4,5])))
temp<-example.out$geno[,26]
temp[temp==1 | temp==2]<-1
temp[temp==3]<-2
temp[temp==4]<-3
example.out$geno[,26]<-temp
example.out$segr.type[26]<-"B2.6"
example.out$segr.type.num[26]<-2
print(x, mrk1="M19", mrk2 = "M26")
t(sapply(ble, function(x) c(x[26,19], x[19,26])))
print(x, mrk1="M19", mrk2 = "M25")
t(sapply(ble, function(x) c(x[25,19], x[19,25])))
print(x, mrk1="M1", mrk2 = "M19")
t(sapply(ble, function(x) c(x[19,1], x[1,19])))
print(x, mrk1="M13", mrk2 = "M19")
t(sapply(ble, function(x) c(x[19,13], x[13,19])))






