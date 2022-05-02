##FUNCTIONS

.getNP<-function(v1, z=1.96){
 
  #p=(v[1]+v[3])/2; # 
  p = v[2]
#  if(length(v==0)) return (NULL)
#  if(v[2]==0) v[2] = (v[1]+v[3])/2  #to avoid NA
  lower=v[1]; upper=v[3]
  
  if(TRUE){
   err= (upper-lower)/2
  #if(p<0.0000001) p = (upper+lower)/2
    N= (z^2/err^2) * p*(1-p)
   # if(N==0){
    #  p = (lower +upper)/2
    #  N= (z^2/err^2) * p*(1-p)
    #}  
  }else{
    err1= abs(upper-p)
    N1= z^2/err1^2 * p*(1-p)
    err1= abs(lower-p)
    N2 = z^2/err1^2 * p*(1-p)
    N = (N1+N2)/2
  }
 res =  as.list(c(v[2],N))
 names(res) = c("p","N")
 res
}



##v1 is c(lower,p,upper)
#v2 is c(lower,p,upper)
#second minus first
#adjusted wald, see https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/diffprop.htm
.getDiff<-function(v1_,v2_,z=1.96){
  v1 =sort(v1_,na.last=T)
  v2 = sort(v2_,na.last=T)
  if(length(which(is.na(v1)))>0 || length(which(is.na(v2)))) return (c(NA,NA,NA))
  if(!is.na(v1[2]) && (v1[2]==0 || v1[2] ==1) )v1[2] = (v1[1]+v1[3])/2.0
  if(!is.na(v2[2]) && (v2[2]==0 || v2[2] ==1) )v2[2] = (v2[1]+v2[3])/2.0
  
   s1 = .getNP(v1,z=z)
   s2 = .getNP(v2,z=z)
   diff = z* sqrt(s1$p*(1-s1$p)/(s1$N+2)  + s2$p*(1-s2$p)/(s2$N+2) )
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


.merge<-function(tab1,tab2){
  mi = match(dimnames(tab1)[[1]],dimnames(tab2)[[1]])
  tab_all = cbind(tab1,tab2[mi,])
  tab_all
}

.conv<-function(str){
  as.numeric(strsplit(gsub("-"," ",gsub("[%())]","",str))," ")[[1]])/100
}




.calcAllP<-function(tab,comparison,CI_in=0.95, CI_out=0.95){
  
  cols1 = grep(comparison[[1]], names(tab))
  cols2 = grep(comparison[[2]], names(tab))
  z_in = qnorm((1-CI_in)/2, lower.tail=F)
  z_out= qnorm((1-CI_out)/2, lower.tail=F)
  pvs=apply(tab,1, function(v){
    v3=.getDiff(v[cols1],v[cols2], z = z_in)
    pv=.calcP(v3,z=z_out)
    c(pv, 100*v3[2], 100*v3[1],100*v3[3])
  })
  pvs
}
##CI_in is assumed confidence levels of tables
##CI_out is desired confidence interval for result
.calcAllP1<-function(tabs,comparison, CI_in = 0.95, CI_out=0.95){
  outfile=paste(c("pval",comparison,"csv"),collapse=".")
  res = lapply(tabs, .calcAllP, comparison, CI_in, CI_out)
  df = data.frame(cbind(t(res[[1]]), t(res[[2]])))
  names(df) = unlist(lapply(names(tabs), function(x) paste(x,c("pv","delta (%)","delta lower(%)", "delta upper(%)"))))
  df
}

.conv1<-function(v){
  v[1] =  round(v[1]*10)/10
  v[2] = round(v[2]*10)/10
  v[3] = round(v[3]*10)/10
  paste0(v[1]," (",v[2],",",v[3],")")
}
.formatTable<-function(df){
  
  
  df[,1] = gsub(" ","",sprintf("%3.2g",df[,1]))
  df[,5] = gsub(" ","",sprintf("%3.2g",df[,5]))
  df[,2] =apply(df[,2:4],1,.conv1)
  df[,6] =apply(df[,6:8],1,.conv1)
  Drug = dimnames(df)[[1]]
  apply(cbind(Drug, df[,c(1,2,5,6)]),c(1,2), function(x) if(length(grep("NA",x))>0) "-" else x)
  #write.table(cbind(drugs,df),quote=F,row.names=F,sep="\t",file=outfile)
 
}


##example 
#v1 = c(84.9, 80.5,88.4)/100
#v2=c(77.3 ,74.9,79.5)/100
#v3=.getDiff(v1,v2)
#pv=.calcP(v3)
#sprintf("%5.3g",pv)
#.getN(p,l,u)

setwd("../pvalues")

##READ AND REFORMAT
.readWalker<-function(file="who-results.csv",  types = c("sensitivity","specificity")){
  tab=read.csv("who-results.csv")

  walker_res=lapply(types,function(x){
   t =  apply(tab[,grep(x,names(tab))],c(1,2),function(y)y/100)
   dimnames(t)[[1]] = tolower(tab[,1])
   dimnames(t)[[2]] = paste(1:3,"Walker")
   t
  })
  names(walker_res) = types
  walker_res
}


.readWHO<-function(f1 = "sn_who.csv",f2 ="sp_who.csv", types = c("sensitivity","specificity") ){
    sn_who=read.csv(f1)
    sp_who=read.csv(f2)
    who_res = lapply(list(sensitivity=sn_who, specificity=sp_who), function(t){
      t1 = t[,-(1:2)]
      dimnames(t1)[[1]] = tolower(t[,1])
      dimnames(t1)[[2]] = paste(1:3,"WHO")
      t1
    })
    names(who_res)=types
    who_res
}



.mergeWalkerWHO<-function(walker_res, who_res){
  tabs = list()
  for(i in 1:length(walker_res)){
    tabs[[i]] = .merge(walker_res[[i]], who_res[[i]])
  }
  names(tabs) = names(walker_res)
  tabs
}



#FOR COMPARING COMBINATIONS
.readResults1<-function(file="results.csv", types=c("sensitivity","specificity")){
  allres1 = read.csv(file)
  df2 =data.frame( t(data.frame(apply(allres1,1, function(x){
    c(.conv(x[5]),.conv(x[6]))
  }))))
  names(df2) = c(paste("sensitivity",1:3),paste("specificity",1:3))
  df3=cbind(allres1[,1:2],df2)
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
  tab_all
}



.writeTable<-function(tab,file, sep){
 for(j in 1:nrow(tab)){
   writeLines( paste(tab[j,],collapse=" & "),con=file,sep=sep)
   
 }
}

.writeLatexTable<-function(formatted,outfile="latex.txt"){
 
  out=file(outfile,"wt")
  sep=" eol \n"
  nmes1 = unlist(lapply(names(formatted[[1]]), function(x) paste0("\\textbf{",x,"}")))
  writeLines( paste(nmes1,collapse=" & "),con=out,sep=sep)
  line1="\\multicolumn{_n_}{c}{\\textbf{_entry_}}"
for(i in 1:length(formatted)){
 line2=gsub("_entry_", names(formatted)[[i]],gsub("_n_", ncol(formatted[[i]]), line1))
  writeLines(line2, con=out,sep=sep)
  .writeTable(formatted[[i]], file=out,sep=sep)
}
  close(out)  
}


walker_res = .readWalker("who-results.csv",types=c("sensitivity","specificity"))
who_res = .readWHO("sn_who.csv", "sp_who.csv")
tabs = .mergeWalkerWHO(walker_res,who_res)
tab_all=.readResults1("results.csv",types=c("sensitivity","specificity"))



all_res=list(.calcAllP1(tabs,comparison=c("Walker","WHO")))
names(all_res)="Walker et al vs WHO"
comparisons = list(c("WHO","Combined"),c("Mykrobe","Combined"),c("Mykrobe","WHO"))
all_res1= lapply(comparisons, function(x) .calcAllP1(tab_all,comparison=x))
names(all_res1) = unlist(lapply(comparisons, function(x) paste(rev(x),collapse=" vs ")))
all_res2 = c(all_res, all_res1)

formatted=lapply(all_res2, .formatTable)
mi =  match(tolower(dimnames(formatted[[2]])[[1]]),dimnames(formatted[[1]])[[1]])
formatted[[1]] = formatted[[1]][mi,]
formatted[[1]][,1] = formatted[[2]][,1]
dimnames(formatted[[1]])[[1]] = dimnames(formatted[[2]])[[1]]


.writeLatexTable(formatted,"pvals3.tex")
closeAllConnections()
