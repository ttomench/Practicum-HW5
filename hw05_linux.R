rm(list=ls())
#setwd("Z:/Practicum-HW5")
setwd("~/Documents/Practicum")
#setwd("~/Documents/Practicum/Practicum_HW5")
load("predictors.Rdata")

install.packages("sn","fields","mvtnorm")
library(sn)
library(fields)
library(mvtnorm)

cols=1:16  	#Columns (i.e. levels) of the temperature profiles to be used in the 
years=as.numeric(substr(dates,1,4))
months=as.numeric(substr(dates,5,6))
all.months= as.numeric(substr(dates[date.ind],5,6))
all.years = as.numeric(substr(dates[date.ind],1,4))

prior.probs=array(0,dim=c(length(stations),length(unique(months)),4))

for(i in 1:length(stations)){
  
  if(i%%100==0){print(i)}	
  
  for(j in 1:length(unique(months))){
    mon=sort(unique(months))[j]
    #Finding the right stations	
    station.i=which(station.ind==i)
    #Finding the right months
    #month.labels=which(months==mon)
    #month.rows=which(date.ind%in%month.labels)
    month.rows=which(all.months==mon)
    #Getting the right stations AND months
    rows.needed=intersect(station.i,month.rows)
    
    rain.nn=length(which(ptype[rows.needed]=="RA"))
    snow.nn=length(which(ptype[rows.needed]=="SN"))
    pellet.nn=length(which(ptype[rows.needed]=="IP"))
    ice.nn=length(which(ptype[rows.needed]=="FZRA"))
    
    prior.probs[i,j,1:4]=c(rain.nn,snow.nn,pellet.nn,ice.nn)/length(rows.needed)		
  }
}


##### Brier Score Calculator #####
BS.fun <- function(prob.hats) {
  classes=c("RA","SN","IP","FZRA")
  BS=0
  for(i in 1:4){
    
    matches=which(prob.hats[,5]==classes[i])
    o.ik=rep(0,length(prob.hats[,5]))
    o.ik[matches]=rep(1,length(matches))
    
    p.ik=prob.hats[,i]
    
    BS=BS+sum((p.ik-o.ik)^2,na.rm=T)
    
    
  }
  return(BS/length(prob.hats[,5]))
}


param = c(1,10)
##### Optimizing a & b #####
ab.BSS <- function(param){
  print("Inner Loop")
  a = param[1]
  b = param[2]
  
  if (a < 0 || b < 0 ){return(0)}
  
  prob.hats=data.frame(matrix(0,nrow=sum(train.nn),ncol=5))
  prob.hats.ref = data.frame(matrix(0,nrow=sum(train.nn),ncol=5))

  
  
  train.rows = train.rows[[i]]
  test.rows = test.rows[[i]]
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows]=="RA")
  snow.rows=which(ptype[train.rows]=="SN")
  pellet.rows=which(ptype[train.rows]=="IP")
  ice.rows=which(ptype[train.rows]=="FZRA")

  cov.reg.rain=a*cov.train[[1]][[i]]+b*I
  cov.reg.snow=a*cov.train[[2]][[i]]+b*I
  cov.reg.pellet=a*cov.train[[3]][[i]]+b*I
  cov.reg.freeze=a*cov.train[[4]][[i]]+b*I
  
  print("ab loop:")
  print(i)
  
  for(j in 1:train.nn[i]){
    if(j%%1000==0){print(j)}
    ind=ind+1
    
    station.j=station.ind[train.rows[[i]][j]]
    mon.j=months[date.ind[train.rows[[i]][j]]]
    mon.col=which(sort(unique(months))==mon.j)
    
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[1]][,i], cov.reg.rain)
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[2]][,i], cov.reg.snow)
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[3]][,i], cov.reg.pellet)
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[4]][,i], cov.reg.freeze)
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    
    
    prob.hats[ind,1:4]=collection/sum(collection)
    prob.hats[ind,5]=ptype[test.rows[[i]][j]]
    
    prob.hats.ref[ind,1:4]=pi.smk
    prob.hats.ref[ind,5]=ptype[test.rows[[i]][j]]
    }
  
  
  BS[i]= BS.fun(prob.hats)
  BS.ref[i]= BS.fun(prob.hats.ref)
  
  BSS[i] = 1-(BS[i]/BS.ref[i])
  
  return(-BSS[i])
}


########################################################
##Find the total number of testing profiles and all of their indices
########################################################
test.nn=array()
ALL.testing.rows=NULL
for(i in 1:12){
  test.years=2000+i
  test.labels.start=min(which(all.months >= 9 & all.years==test.years ))
  test.labels.end=max(which( all.months <= 5 & all.years==test.years+1 ))
  test.rows = test.labels.start:test.labels.end
  test.nn[i]=test.labels.end - test.labels.start+1
  ALL.testing.rows=c(ALL.testing.rows,test.rows)}
sum(test.nn)



train.nn=array()
test.nn=array()

mean.train=list()
mean.train[[1]]=matrix(0,nrow=16,ncol=12)
mean.train[[2]]=matrix(0,nrow=16,ncol=12)
mean.train[[3]]=matrix(0,nrow=16,ncol=12)
mean.train[[4]]=matrix(0,nrow=16,ncol=12)
names(mean.train)=c("rain","snow","pellets","ice")

cov.train=list()
cov.train[[1]]=list()
cov.train[[2]]=list()
cov.train[[3]]=list()
cov.train[[4]]=list()
names(cov.train)=c("rain","snow","pellets","ice")

train.rows = list()
test.rows = list()

ind=0


for(i in 1:12){
  
  ##UPDATE THESE
  
  train.years=1996:2001+i-1
  test.years=2000+i
  
  print(i)
  
  train.labels.start = min(which( all.months >= 9 & all.years >=train.years[1]))
  train.labels.end = max(which( all.months <= 5 & all.years<=train.years[6]))
  test.labels.start=min(which(all.months >= 9 & all.years==test.years ))
  test.labels.end=max(which( all.months <= 5 & all.years==test.years+1 ))
  
  train.nn[i]=train.labels.end-train.labels.start+1
  test.nn[i]=test.labels.end - test.labels.start+1
  
  train.rows[[i]] = train.labels.start:train.labels.end
  test.rows[[i]] = test.labels.start:test.labels.end
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows[[i]]]=="RA")
  snow.rows=which(ptype[train.rows[[i]]]=="SN")
  pellet.rows=which(ptype[train.rows[[i]]]=="IP")
  ice.rows=which(ptype[train.rows[[i]]]=="FZRA")
  
  mean.train[[1]][,i]=apply(Twb.prof[train.rows[[i]][rain.rows],cols],2,mean)
  mean.train[[2]][,i]=apply(Twb.prof[train.rows[[i]][snow.rows],cols],2,mean)
  mean.train[[3]][,i]=apply(Twb.prof[train.rows[[i]][pellet.rows],cols],2,mean)
  mean.train[[4]][,i]=apply(Twb.prof[train.rows[[i]][ice.rows],cols],2,mean)
  
  cov.train[[1]][[i]]=cov(Twb.prof[train.rows[[i]][rain.rows],cols])
  cov.train[[2]][[i]]=cov(Twb.prof[train.rows[[i]][snow.rows],cols])
  cov.train[[3]][[i]]=cov(Twb.prof[train.rows[[i]][pellet.rows],cols])
  cov.train[[4]][[i]]=cov(Twb.prof[train.rows[[i]][ice.rows],cols])
  
}




prob.hats=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")

BS = array()
BS.ref = array()
BSS = array()

cov.reg = list()
cov.reg[[1]]=list() 
cov.reg[[2]]=list()
cov.reg[[3]]=list()
cov.reg[[4]]=list()
I = diag(1,length(cols))

abstart = matrix(0,13,2)
abstart[1,] = c(1,10)
print("Main Loop:")

for(i in 1:12){
  
  abstart[(i+1),] <- optim(abstart[i,],ab.BSS)$par
  a = abstart[i+1,1]
  b = abstart[i+1,2]
  
  cov.reg[[1]][[i]] = a*cov.train[[1]][[i]] + b*I 
  cov.reg[[2]][[i]] = a*cov.train[[2]][[i]] + b*I
  cov.reg[[3]][[i]] = a*cov.train[[3]][[i]] + b*I
  cov.reg[[4]][[i]] = a*cov.train[[4]][[i]] + b*I
  print("ML: ")
  print(i)
  
  for(j in 1:test.nn[i]){
    if(j%%1000==0){print(j)}
    ind=ind+1
    
    
    
    
    station.j=station.ind[test.rows[[i]][j]]
    #mon.j=months[date.ind[test.rows[[i]][j]]]
    mon.j=all.months[test.rows[[i]][j]]
    mon.col=which(sort(unique(months))==mon.j)
    
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[test.rows[[i]][j],cols], mean.train[[1]][,i], cov.reg[[1]][[i]])
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[test.rows[[i]][j],cols], mean.train[[2]][,i], cov.reg[[2]][[i]])
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[test.rows[[i]][j],cols], mean.train[[3]][,i], cov.reg[[3]][[i]])
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[test.rows[[i]][j],cols], mean.train[[4]][,i], cov.reg[[4]][[i]])
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    
    
    prob.hats[ind,1:4]=collection/sum(collection)
    prob.hats[ind,5]=ptype[test.rows[[i]][j]]
    
    #prob.hats.clim[ind,1:4]=pi.smk
    #prob.hats.clim[ind,5]=ptype[test.rows[[i]][j]]
  }
  
  write.table(prob.hats,file="regularizatoin_classification.txt")
  
}


