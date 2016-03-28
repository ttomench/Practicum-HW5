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
ab.BSS <- function(param,i,j){
  a = param[1]
  b = param[2]
  
  if (a < 0 || b < 0 ){return(0)}
  I = diag(1,length(cols))
  BSS = 0
  prob.hats=data.frame(matrix(0,nrow=sum(train.nn),ncol=5))
  prob.hats.ref = data.frame(matrix(0,nrow=sum(train.nn),ncol=5))
  
  #for(i in 1:12){
    
    
    train.rows = train.rows[[i]]
    test.rows = test.rows[[i]]
    #######################################################
    ##Computing means and covariances for each precip type
    #######################################################
    rain.rows=which(ptype[train.rows]=="RA")
    snow.rows=which(ptype[train.rows]=="SN")
    pellet.rows=which(ptype[train.rows]=="IP")
    ice.rows=which(ptype[train.rows]=="FZRA")
    
    # mean.train[[1]][,i]=apply(Twb.prof[train.rows[rain.rows],cols],2,mean)
    # mean.train[[2]][,i]=apply(Twb.prof[train.rows[snow.rows],cols],2,mean)
    # mean.train[[3]][,i]=apply(Twb.prof[train.rows[pellet.rows],cols],2,mean)
    # mean.train[[4]][,i]=apply(Twb.prof[train.rows[ice.rows],cols],2,mean)
    
    cov.reg.rain=a*cov.train[[1]][[i]]+b*I
    cov.reg.snow=a*cov.train[[2]][[i]]+b*I
    cov.reg.pellet=a*cov.train[[3]][[i]]+b*I
    cov.reg.freeze=a*cov.train[[4]][[i]]+b*I
    
    print(i)
    #for(j in 1:test.nn[i]){
      if(j%%1000==0){print(j)}
      ind=ind+1
      
      station.j=station.ind[test.rows[[i]][j]]
      mon.j=months[date.ind[test.rows[[i]][j]]]
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
    #}
    
  
  #}
  
  BS= BS.fun(prob.hats)
  BS.ref= BS.fun(prob.hats.ref)
  
  BSS = 1-(BS/BS.ref)
  
  return(-BSS)
}


reliability <- function(forcast.probs){
  classes=c("RA","SN","IP","FZRA")
  names = c("Rain","Snow", "Pellets","Freeze")
  color = c(rainbow(12)[5],rainbow(12)[9],rainbow(12)[3],rainbow(12)[1])

  o.ik = matrix(0,length(forcast.probs[,5]),4)
  for(i in 1:4) {
    matches=which(forcast.probs[,5]==classes[i])
    o.ik[matches,i]=rep(1,length(matches))
  }
  plot(forcast.probs[,1],o.ik[,1],col=color[1], pch="b", xlab = "Forcast Probabilities", ylab = "Observed frequencies")
  lines(forcast.probs[,2],o.ik[,2], col= color[2], pch="b")
  lines(forcast.probs[,3],o.ik[,3], col= color[3], pch="b")
  lines(forcast.probs[,4],o.ik[,4], col= color[4], pch="b")
  freqs= c(1,50000)
  
  par(fig=c(0,0.2,0.76,1), new=TRUE)
  hist(forcast.probs[,1], breaks=20, col=color[1], main = names[1], ylim=freqs)
  par(fig=c(0,0.2,0.5,.74), new=TRUE)
  hist(forcast.probs[,2], breaks=20, col=color[2], main = names[2], ylim=freqs)
  par(fig=c(0.7,0.9,0.26,0.5), new=TRUE)
  hist(forcast.probs[,3], breaks=20, col=color[3], main = names[3], ylim=freqs)
  par(fig=c(0.7,0.9,0.75,1), new=TRUE)
  hist(forcast.probs[,4], breaks=20, col=color[4], main = names[4], ylim=freqs)
  
  abline(0,1)
}








