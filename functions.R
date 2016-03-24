##############################
####### Functions File #######
##############################


##### Brier Score Calculator #####
BS <- function(prob.hats) {
  classes=c("RA","SN","IP","FZRA")
  BS=0
  for(i in 1:4){
    
    matches=which(prob.hats[,5]==classes[i])
    o.ik=rep(0,length(prob.hats[,5]))
    o.ik[matches]=rep(1,length(matches))
    
    p.ik=prob.hats[,i]
    
    BS=BS+sum((p.ik-o.ik)^2,na.rm=T)
    
  }
  
}
  
  
##### Optimizing a & b #####
ab.BSS <- function(param){
  a = param[1]
  b = param[2]
  
  if (a < 0 || b < 0 ){return(0)}
  
  
  
}