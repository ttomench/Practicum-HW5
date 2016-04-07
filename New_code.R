
## Statistics Practicum: Spring 2016
## Hmwk #5 Partial Solutions

rm(list=ls())
#setwd("Z:/Practicum-HW5")
setwd("~/Documents/Practicum")
load("predictors.Rdata")
library(sn)
library(fields)
library(mvtnorm)


# dates:        all dates for which data are available
# stations:     names of all stations for which data are available
# lat, lon:     longitude and latitude coordinate for each of these stations
# elev:         station elevation (m)
# ptype:        all cases where one of the four precipitation types of interest was reported
# Twb.prof:     the corresponding vertical profiles of wetbulb temperature (0m to 3000m above the surface in steps of 100m)
# station.ind:  for each case, the index of the station where it was reported
# date.ind:     for each case, the index of the date on which it was reported

#############################
##Answer From Last Assignment
#############################
cols=1:16  	#Columns (i.e. levels) of the temperature profiles to be used in the 
years=as.numeric(substr(dates,1,4))
months=as.numeric(substr(dates,5,6))
all.months= as.numeric(substr(dates[date.ind],5,6))
all.years = as.numeric(substr(dates[date.ind],1,4))

# write.table(prior.probs,file="prior_probs.txt")
prior.probs<-read.table("prior_probs.txt")
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

########################################################
##Baseline (Assignment #4) Classification Model
########################################################
prob.hats=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")

prob.hats.clim = data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats.clim)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")


prob.hats.train=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats.train)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")

prob.hats.clim.train = data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats.clim.train)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")

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

######################################################
#Computing probabilities of observations belonging to
#each of the 4 groups
######################################################


for(i in 1:12){
  print(i)
  for(j in 1:train.nn[i]){
    if(j%%1000==0){print(j)}
    ind=ind+1
    
    station.j=station.ind[train.rows[[i]][j]]
    mon.j=months[date.ind[train.rows[[i]][j]]]
    mon.col=which(sort(unique(months))==mon.j)
    
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[train.nn[[i]][j],cols], mean.train[[1]][,i], cov.train[[1]][[i]])
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[train.nn[[i]][j],cols], mean.train[[2]][,i], cov.train[[2]][[i]])
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[train.nn[[i]][j],cols], mean.train[[3]][,i], cov.train[[3]][[i]])
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[train.nn[[i]][j],cols], mean.train[[4]][,i], cov.train[[4]][[i]])
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    
    
    prob.hats[ind,1:4]=collection/sum(collection)
    prob.hats[ind,5]=ptype[test.rows[[i]][j]]
    
    prob.hats.clim[ind,1:4]=pi.smk
    prob.hats.clim[ind,5]=ptype[test.rows[[i]][j]]
  }
}

write.table(prob.hats,file="baseline_classification_train.txt")
write.table(prob.hats.clim, file="climatology_classification_train.txt")

BS=BS.fun(prob.hats);BS
BS.ref=BS.fun(prob.hats.clim);BS.ref
BSS=1-(BS/BS.ref);BSS