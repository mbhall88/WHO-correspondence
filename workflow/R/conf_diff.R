.getNP<-function(v1, z=1.96){
  v =sort(v1)
  #p=(v[1]+v[3])/2; # 
  p = v[2]
  lower=v[1]; upper=v[3]
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
#second minus first
.getDiff<-function(v1,v2,z=1.96){
   s1 = .getNP(v1,z=z)
   s2 = .getNP(v2,z=z)
   diff = z* sqrt(s1$p*(1-s1$p)/s1$N  + s2$p*(1-s2$p)/s2$N )
   mid = s2$p-s1$p
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

setwd("../pvalues")

tab=read.csv("who-results.csv")
types = c("sensitivity","specificity")
walker_res=lapply(types,function(x){
 t =  apply(tab[,grep(x,names(tab))],c(1,2),function(y)y/100)
 dimnames(t)[[1]] = tolower(tab[,1])
 dimnames(t)[[2]] = paste(1:3,"Walker")
 t
})

sn_who=read.csv("sn_who.csv")
sp_who=read.csv("sp_who.csv")
who_res = lapply(list(sensitivity=sn_who, specificity=sp_who), function(t){
  t1 = t[,-(1:2)]
  dimnames(t1)[[1]] = tolower(t[,1])
  dimnames(t1)[[2]] = paste(1:3,"WHO")
  t1
})
names(who_res)=types
names(walker_res) = types

.merge<-function(tab1,tab2){
  mi = match(dimnames(tab1)[[1]],dimnames(tab2)[[1]])
  tab_all = cbind(tab1,tab2[mi,])
  tab_all
}

tabs = list()
for(i in 1:length(walker_res)){
  tabs[[i]] = .merge(walker_res[[i]], who_res[[i]])
}
names(tabs) = names(walker_res)

.calcAllP<-function(tab,comparison){
  cols1 = grep(comparison[[1]], names(tab))
  cols2 = grep(comparison[[2]], names(tab))
  pvs=apply(tab,1, function(v){
    v3=.getDiff(v[cols1],v[cols2])
    pv=.calcP(v3)
    c(pv, 100*v3[2])
  })
  pvs
}

.calcAllP1<-function(tabs,comparison){
  outfile=paste(c("pval",comparison,"csv"),collapse=".")
  res = lapply(tabs, .calcAllP, comparison)
  df = data.frame(cbind(t(res[[1]]), t(res[[2]])))
names(df) = unlist(lapply(names(tabs), function(x) paste(x,c("pv","delta (%)"))))
 df[,1] = gsub(" ","",sprintf("%3.2g",df[,1]))
 df[,3] = gsub(" ","",sprintf("%3.2g",df[,3]))
 df[,2] =round(df[,2]*10)/10
 df[,4] =round(df[,4]*10)/10
 drugs = dimnames(df)[[1]]
 write.table(cbind(drugs,df),quote=F,row.names=F,sep="\t",file=outfile)
 invisible(df)
}

all_res=.calcAllP1(tabs,comparison=c("Walker","WHO"))


allres1 = read.csv("results.csv")
.conv<-function(str){
  as.numeric(strsplit(gsub("-"," ",gsub("[%())]","",str))," ")[[1]])/100
  
}
df2 =data.frame( t(data.frame(apply(res1,1, function(x){
  c(.conv(x[5]),.conv(x[6]))
}))))
names(df2) = c(paste("sensitivity",1:3),paste("specificity",1:3))

df3=cbind(res1[,1:2],df2)
dimnames(df3)[[1]]=1:nrow(df3)
#df3
tab_all = lapply(types, function(x) {
  t =  df3[,c(1,2,grep(x,names(df3)))]
  t$Panel = factor(t$Panel)
  t1 = lapply(levels(t$Panel), function(y){
    t[t$Panel==y,]
  })
  names(t1) = levels(t$Panel)
  t2 =t1[[1]][-(1:2)]
  dimnames(t2)[[1]] = t1[[1]][,1]
  nmes = c( paste(1:3,names(t1)[[1]]))
  for(j in 2:length(t1)){
    t2 = cbind(t2,t1[[j]][-(1:2)])
    nmes = c(nmes,paste(1:3,names(t1)[[j]]))
  }
  names(t2) = nmes
  t2
})
names(tab_all) = (types)

all_res1=.calcAllP1(tab_all,comparison=c("WHO","Combined"))
all_res1=.calcAllP1(tab_all,comparison=c("Mykrobe","Combined"))
all_res1=.calcAllP1(tab_all,comparison=c("Mykrobe","WHO"))
