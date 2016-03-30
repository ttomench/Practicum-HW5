##############################
####### Functions File #######
##############################


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
  

##### Optimizing a & b #####
ab.BSS <- function(param) {
  print("Inner Loop")
  a = param[1]
  b = param[2]
  
  if (a < 0 || b < 0 ){return(0)}
  
  
  
  #train.rows = train.rows[[i]]
  #test.rows = test.rows[[i]]
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows[[i]]]=="RA")
  snow.rows=which(ptype[train.rows[[i]]]=="SN")
  pellet.rows=which(ptype[train.rows[[i]]]=="IP")
  ice.rows=which(ptype[train.rows[[i]]]=="FZRA")
  
  cov.reg.rain=a*cov.train[[1]][[i]]+b*I
  cov.reg.snow=a*cov.train[[2]][[i]]+b*I
  cov.reg.pellet=a*cov.train[[3]][[i]]+b*I
  cov.reg.freeze=a*cov.train[[4]][[i]]+b*I
  
  print("ab loop:");print(i)
  
  for(j in 1:train.nn[i]){
    
    if(j%%1000==0){print(paste("j: ",j))}
    ind=ind+1
    
    station.j=station.ind[train.rows[[i]][j]]
    #print("station.j: ");print(station.j)
    mon.j=months[date.ind[train.rows[[i]][j]]]
    mon.col=which(sort(unique(months))==mon.j)
    
    pi.smk=prior.probs[station.j,mon.col,]
    #print(pi.smk)
    pi.den.rain=(pi.smk[1]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[1]][,i], cov.reg.rain))
    pi.den.snow=(pi.smk[2]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[2]][,i], cov.reg.snow))
    pi.den.pellet=(pi.smk[3]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[3]][,i], cov.reg.pellet))
    pi.den.freeze=(pi.smk[4]*dmvnorm(Twb.prof[train.rows[[i]][j],cols], mean.train[[4]][,i], cov.reg.freeze))
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    
    prob.hats_1[ind,1:4]=collection/sum(collection)
    #prob.hats_1[ind,1:4] = as.numeric(prob.hats_1[ind,1:4])
    prob.hats_1[ind,5]=ptype[train.rows[[i]][j]]
    #print("ind: "); print(ind)
    
    prob.hats.ref[ind,1:4]=pi.smk
    prob.hats.ref[ind,5]=ptype[train.rows[[i]][j]]
    #if(j==tail(train.rows[[i]],1)){print("prob.hats_1: "); print(tail(prob.hats_1))}
    
  }
  
  #prob.hats.temp = prob.hats_1
  BS.temp[i]=BS.fun(prob.hats_1)
  
  BS.ref[i]= BS.fun(prob.hats.ref)
  BSS[i] = 1-(BS.temp[i]/BS.ref[i])
  return(-BSS[i])
}



library(verification)
reliability <- function(forecast.probs){
  
  classes=c("RA","SN","IP","FZRA")
  names = c("Rain","Snow", "Pellets","Freeze")
  color = c(rainbow(12)[5],rainbow(12)[9],rainbow(12)[3],rainbow(12)[1])

  o.ik = matrix(0,length(forecast.probs[,5]),4)
  for(i in 1:4) {
    matches=which(forecast.probs[,5]==classes[i])
    o.ik[matches,i]=rep(1,length(matches))
  }
  plot(forecast.probs[,1],o.ik[,1],col=color[1], pch="b", xlab = "Forcast Probabilities", ylab = "Observed frequencies")
  lines(forecast.probs[,2],o.ik[,2], col= color[2], pch="b")
  lines(forecast.probs[,3],o.ik[,3], col= color[3], pch="b")
  lines(forecast.probs[,4],o.ik[,4], col= color[4], pch="b")
  freqs= c(1,50000)
  
  par(fig=c(0,0.2,0.76,1), new=TRUE)
  hist(forecast.probs[,1], breaks=20, col=color[1], main = names[1], ylim=freqs)
  par(fig=c(0,0.2,0.5,.74), new=TRUE)
  hist(forecast.probs[,2], breaks=20, col=color[2], main = names[2], ylim=freqs)
  par(fig=c(0.7,0.9,0.26,0.5), new=TRUE)
  hist(forecast.probs[,3], breaks=20, col=color[3], main = names[3], ylim=freqs)
  par(fig=c(0.7,0.9,0.75,1), new=TRUE)
  hist(forecast.probs[,4], breaks=20, col=color[4], main = names[4], ylim=freqs)
  
  abline(0,1)
}








