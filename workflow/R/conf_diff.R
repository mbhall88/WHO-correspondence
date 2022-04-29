.getNP<-function(v1, z=1.96){
  v =sort(v1)
  p=v[2];lower=v[1]; upper=v[3]
  err1= abs(upper-p)
  N1= z^2/err1^2 * p*(1-p)
  err1= abs(lower-p)
  N2 = z^2/err1^2 * p*(1-p)
  N = (N1+N2)/2
 res =  as.list(c(v[2],N))
 names(res) = c("p","N")
 res
}



##v1 is c(lower,p,upper)
#v2 is c(lower,p,upper)
.getDiff<-function(v1,v2,z=1.96){
   s1 = .getNP(v1,z=z)
   s2 = .getNP(v2,z=z)
   diff = z* sqrt(s1$p*(1-s1$p)/s1$N  + s2$p*(1-s2$p)/s2$N )
   mid = abs(s2$p-s1$p)
v3=   c(mid-diff,mid, mid+diff)
  v3
  
}

.calcP<-function(v3, z=1.96,x=0){
  diff = (v3[3] - v3[1])/2
  SE = (diff)/z
  #mid/SE
  mid = v3[2]
  zscore = (mid-x)/SE
  2* pnorm(-abs(zscore))
}

##example 
v1 = c(84.9, 80.5,88.4)/100
v2=c(77.3 ,74.9,79.5)/100
v3=.getDiff(v1,v2)
pv=.calcP(v3)
sprintf("%5.3g",pv)
#.getN(p,l,u)





