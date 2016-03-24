
## Statistics Practicum: Spring 2016
## Hmwk #5 Partial Solutions

rm(list=ls())
setwd("Z:/Practicum-HW5")
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
############################
##Computing Prior Probabilities
############################
prior.probs=array(0,dim=c(length(stations),length(unique(months)),4))

for(i in 1:length(stations)){
  
  print(i)	
  
  for(j in 1:length(unique(months))){
    mon=sort(unique(months))[j]
    #Finding the right stations	
    station.i=which(station.ind==i)
    #Finding the right months
    month.labels=which(months==mon)
    month.rows=which(date.ind%in%month.labels)
    #Getting the right stations AND months
    rows.needed=intersect(station.i,month.rows)
    
    rain.nn=length(which(ptype[rows.needed]=="RA"))
    snow.nn=length(which(ptype[rows.needed]=="SN"))
    pellet.nn=length(which(ptype[rows.needed]=="IP"))
    ice.nn=length(which(ptype[rows.needed]=="FZRA"))
    
    prior.probs[i,j,1:4]=c(rain.nn,snow.nn,pellet.nn,ice.nn)/length(rows.needed)		
  }
}

########################################################
##Find the total number of testing profiles and all of their indices
########################################################
test.nn=array()
ALL.testing.rows=NULL
for(i in 1:12){
  test.years=2000+i
  test.labels=which(years==test.years)
  test.rows=which(date.ind%in%test.labels)
  test.nn[i]=length(test.rows)
  ALL.testing.rows=c(ALL.testing.rows,test.rows)}


########################################################
##Baseline (Assignment #4) Classification Model
########################################################
prob.hats=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")

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
  
  train.rows = train.labels.start:train.labels.end
  test.rows = test.labels.start:test.labels.end
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows]=="RA")
  snow.rows=which(ptype[train.rows]=="SN")
  pellet.rows=which(ptype[train.rows]=="IP")
  ice.rows=which(ptype[train.rows]=="FZRA")
  
  mean.train[[1]][,i]=apply(Twb.prof[train.rows[rain.rows],cols],2,mean)
  mean.train[[2]][,i]=apply(Twb.prof[train.rows[snow.rows],cols],2,mean)
  mean.train[[3]][,i]=apply(Twb.prof[train.rows[pellet.rows],cols],2,mean)
  mean.train[[4]][,i]=apply(Twb.prof[train.rows[ice.rows],cols],2,mean)
  
  cov.train[[1]][[i]]=cov(Twb.prof[train.rows[rain.rows],cols])
  cov.train[[2]][[i]]=cov(Twb.prof[train.rows[snow.rows],cols])
  cov.train[[3]][[i]]=cov(Twb.prof[train.rows[pellet.rows],cols])
  cov.train[[4]][[i]]=cov(Twb.prof[train.rows[ice.rows],cols])
  
  


}

######################################################
#Computing probabilities of observations belonging to
#each of the 4 groups
######################################################
for(i in 1:12){
  train.years=1996:2001+i-1
  test.years=2000+i

  train.labels.start = min(which( all.months >= 9 & all.years >=train.years[1]))
  train.labels.end = max(which( all.months <= 5 & all.years<=train.years[6]))
  test.labels.start=min(which(all.months >= 9 & all.years==test.years ))
  test.labels.end=max(which( all.months <= 5 & all.years==test.years+1 ))

  train.rows = train.labels.start:train.labels.end
  test.rows = test.labels.start:test.labels.end
  for(j in 1:test.nn[i]){
    #if(j%%1000==0){print(j)}
    ind=ind+1

    station.j=station.ind[test.rows[j]]
    mon.j=months[date.ind[test.rows[j]]]
    mon.col=which(sort(unique(months))==mon.j)

    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[1]][,i], cov.train[[1]][[i]])
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[2]][,i], cov.train[[2]][[i]])
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[3]][,i], cov.train[[3]][[i]])
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[4]][,i], cov.train[[4]][[i]])
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    prob.hats[ind,1:4]=collection/sum(collection)
    prob.hats[ind,5]=ptype[test.rows[j]]
  }
}




write.table(prob.hats,file="baseline_classification.txt")
prob.hats = read.table(file="baseline_classification.txt", header=T)
prob.hats.base = read.table(file="base.txt", header=T)
prob.hats.clim = read.table(file="climatology.txt", header=T)
prob.hats.noaa = read.table(file="michael.txt", header=T)

########################################################
##Adjusting for different columns  
########################################################
prob.hats.noaa = prob.hats.noaa[,c(2,1,3,4)]
prob.hats.noaa = cbind(prob.hats.noaa,prob.hats.base[,5])


########################################################
##Forecast Evaluation---START HERE
########################################################

prob.hats.baseline = prob.hats.base
classes=c("RA","SN","IP","FZRA")
BS=0
BS.ref=0
for(i in 1:4){
  
  matches=which(prob.hats.baseline[,5]==classes[i])
  o.ik=rep(0,length(prob.hats.baseline[,5]))
  o.ik[matches]=rep(1,length(matches))
  
  p.ik=prob.hats.baseline[,i]
  p.clim.ik=prob.hats.clim[,i]
  
  BS=BS+sum((p.ik-o.ik)^2,na.rm=T)
  BS.ref = BS.ref+sum((p.clim.ik-o.ik)^2,na.rm=T)
  
}

BS
BS/length(prob.hats.baseline[,5])
#0.2217187

BS.ref
BS.ref/length(prob.hats.clim[,5])
#0.2185553



#Deterministic classification assessment
observed=prob.hats.base[,5]
observed[which(prob.hats.base[,5]=="RA")]=rep(1,length(which(prob.hats.base[,5]=="RA")))
observed[which(prob.hats.base[,5]=="SN")]=rep(2,length(which(prob.hats.base[,5]=="SN")))
observed[which(prob.hats.base[,5]=="IP")]=rep(3,length(which(prob.hats.base[,5]=="IP")))
observed[which(prob.hats.base[,5]=="FZRA")]=rep(4,length(which(prob.hats.base[,5]=="FZRA")))
observed=as.numeric(observed)
summary(observed)
table(observed)
table(prob.hats.base[,5])

prob.class=as.matrix(prob.hats.base[,1:4])
hard.class=as.integer(apply(prob.class,1,which.max))

library(s20x)

CX=crosstabs(hard.class~observed)
round(CX$whole.props,4)
sum(diag(CX$whole.props)) #0.8638785, => 86.4% of the profiles are correctly classified



# observed
# class      1      2      		3      	4
# 1 0.6297 0.0614 0.0012 0.0026
# 2 0.0360 0.2280 0.0009 0.0040
# 3 0.0023 0.0021 0.0010 0.0006
# 4 0.0153 0.0093 0.0005 0.0053



###################
##Problem #2
###################
years=as.numeric(substr(dates,1,4))
months=as.numeric(substr(dates,5,6))

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

for(i in 1:12){
  train.years=1996:2000+i-1
  test.years=2000+i
  
  #print(i)
  #print(train.years)
  #print(test.years)
  
  train.labels=which(years>=train.years[1]& years<=train.years[5])
  test.labels=which(years==test.years)
  
  train.rows=which(date.ind%in%train.labels)
  test.rows=which(date.ind%in%test.labels)
  
  train.nn[i]=length(train.rows)
  test.nn[i]=length(test.rows)
  
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows]=="RA")
  snow.rows=which(ptype[train.rows]=="SN")
  pellet.rows=which(ptype[train.rows]=="IP")
  ice.rows=which(ptype[train.rows]=="FZRA")
  
  mean.train[[1]][,i]=apply(Twb.prof[train.rows[rain.rows],1:16],2,mean)
  mean.train[[2]][,i]=apply(Twb.prof[train.rows[snow.rows],1:16],2,mean)
  mean.train[[3]][,i]=apply(Twb.prof[train.rows[pellet.rows],1:16],2,mean)
  mean.train[[4]][,i]=apply(Twb.prof[train.rows[ice.rows],1:16],2,mean)
  
  cov.train[[1]][[i]]=cov(Twb.prof[train.rows[rain.rows],1:16])
  cov.train[[2]][[i]]=cov(Twb.prof[train.rows[snow.rows],1:16])
  cov.train[[3]][[i]]=cov(Twb.prof[train.rows[pellet.rows],1:16])
  cov.train[[4]][[i]]=cov(Twb.prof[train.rows[ice.rows],1:16])
  
  
  #######################################################
  ##Plotting the covariances for each precip type
  #######################################################
  #pdf(file=paste("Figures/covariance_training_", i, ".pdf",sep=""), width=18, height=5)
  par(mfrow=c(1,4))
  image.plot(1:16,1:16,t(cov.train[[1]][[i]][16:1,]), xaxt="n",yaxt="n",xlab="",ylab="")
  title(paste("Rain Cov, Training ", i, sep=""), cex.main=2.5)
  
  image.plot(1:16,1:16,t(cov.train[[2]][[i]][16:1,]), xaxt="n",yaxt="n",xlab="",ylab="")
  title(paste("Snow Cov, Training ", i, sep=""), cex.main=2.5)
  
  image.plot(1:16,1:16,t(cov.train[[3]][[i]][16:1,]),xaxt="n",yaxt="n",xlab="",ylab="")
  title(paste("Pellets Cov, Training ", i, sep=""), cex.main=2.5)
  
  image.plot(1:16,1:16,t(cov.train[[4]][[i]][16:1,]),xaxt="n",yaxt="n",xlab="",ylab="")
  title(paste("Freeing Rain Cov, Training ", i, sep=""), cex.main=2.5)
  
  #dev.off()
  
}

#pdf(file=paste("Figures/mean_prof_rain", ".pdf",sep=""), width=6, height=6)
xrange=c(260,285)
#Plots of Training Rain Mean Profiles
i=1
plot(mean.train[[1]][1:16,i],0:15,xlab="Temperature",ylab="Height Above Ground",type="n",xlim=xrange)
for(i in 1:12){
  lines(mean.train[[1]][1:16,i],0:15,col=rgb(seq(.5,1,len=12)[i],0,0,1))}
title("Mean Rain Training Profiles",cex.main=2)
abline(v=273.15,col="gray",lwd=3,lty=2)
#dev.off()


#pdf(file=paste("Figures/mean_prof_snow", ".pdf",sep=""), width=6, height=6)
#Plots of Training Snow Mean Profiles
i=1
plot(mean.train[[2]][1:16,i],0:15,xlab="Temperature",ylab="Height Above Ground",type="n",xlim=xrange)
for(i in 1:12){
  lines(mean.train[[2]][1:16,i],0:15,col=rgb(0,0,seq(.5,1,len=12)[i],1))}
title("Mean Snow Training Profiles",cex.main=2)
abline(v=273.15,col="gray",lwd=3,lty=2)
#dev.off()


#pdf(file=paste("Figures/mean_prof_pellets", ".pdf",sep=""), width=6, height=6)
#Plots of Training Ice Pellets Mean Profiles
i=1
plot(mean.train[[3]][1:16,i],0:15,xlab="Temperature",ylab="Height Above Ground",type="n",xlim=xrange)
for(i in 1:12){
  lines(mean.train[[3]][1:16,i],0:15,col=rgb(0,seq(.5,1,len=12)[i],0,1))}
title("Mean Ice Pellets Training Profiles",cex.main=2)
abline(v=273.15,col="gray",lwd=3,lty=2)
#dev.off()

#pdf(file=paste("Figures/mean_prof_freeze", ".pdf",sep=""), width=6, height=6)
#Plots of Training Freezing Rain Mean Profiles
i=1
plot(mean.train[[4]][1:16,i],0:15,xlab="Temperature",ylab="Height Above Ground",type="n",xlim=xrange)
for(i in 1:12){
  lines(mean.train[[4]][1:16,i],0:15,col=gray(seq(.5,1,len=12)[i]))}
title("Mean Freezing Rain Training Profiles",cex.main=2)
abline(v=273.15,col="gray",lwd=3,lty=2)
#dev.off()




###################
##Problem #3
###################


#i=6 is the 2001-2005 period
i=6

show.means=cbind(mean.train[[1]][,i], mean.train[[2]][,i], mean.train[[3]][,i], mean.train[[4]][,i])
colnames(show.means)=c("rain","snow","pellets","ice")
round(show.means,2)

cov.train[[1]][[i]]
cov.train[[2]][[i]]
cov.train[[3]][[i]]
cov.train[[4]][[i]]

image.plot(1:16,1:16,t(cov.train[[1]][[i]][16:1,]), main="Rain Covariance",xaxt="n",yaxt="n",xlab="",ylab="")

image.plot(1:16,1:16,t(cov.train[[2]][[i]][16:1,]), main="Snow Covariance",xaxt="n",yaxt="n",xlab="",ylab="")

image.plot(1:16,1:16,t(cov.train[[3]][[i]][16:1,]), main="Pellets Covariance",xaxt="n",yaxt="n",xlab="",ylab="")

image.plot(1:16,1:16,t(cov.train[[4]][[i]][16:1,]), main="Freezing Rain Covariance",xaxt="n",yaxt="n",xlab="",ylab="")

##Simple Example
#XX=matrix(c(1:16),nrow=4,byrow=T)
#image.plot(1:4,1:4,t(XX[4:1,]))








###################
##Problem #4
###################


month.names=c("Jan","Feb","Mar","Apr","May","Sep","Oct","Nov","Dec")
for(j in 1:length(unique(months))){
  mon=month.names[j]
  
  formap=matrix(0,nrow=length(stations),ncol=4)
  for(i in 1:length(stations)){
    formap[i,]=prior.probs[[i]][,j]	
  }
  
  pdf(file=paste("Figures/prior_probs_v2_", j, ".pdf",sep=""), width=18, height=5)
  
  par(mfrow=c(1,4))
  
  type="Rain"
  quilt.plot(lon-360,lat, formap[,1],zlim=c(0,1))
  map("world",add=T)
  map("state",add=T)
  title(paste("Proportion",type,mon,sep=" "),cex.main=2)
  
  type="Snow"
  quilt.plot(lon-360,lat, formap[,2],zlim=c(0,1))
  map("world",add=T)
  map("state",add=T)
  title(paste("Proportion",type,mon,sep=" "),cex.main=2)
  
  type="Ice Pellets"
  quilt.plot(lon-360,lat, formap[,3],zlim=c(0,.1))
  map("world",add=T)
  map("state",add=T)
  title(paste("Proportion",type,mon,sep=" "),cex.main=2)
  
  type="Freezing Rain"
  quilt.plot(lon-360,lat, formap[,4],zlim=c(0,.1))
  map("world",add=T)
  map("state",add=T)
  title(paste("Proportion",type,mon,sep=" "),cex.main=2)
  dev.off()
  
}	















