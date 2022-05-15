rnaNNLS = function(cnt,pb,normsd=TRUE){
  if(normsd){
    w = apply(cnt,1,sd)
    cnt = sweep(cnt,1,w,'/')
    pb = sweep(pb,1,w,'/')
  }

  r = apply(cnt,2,function(x){
    coef(nnls(pb,x))
  })
  r = t(r)
  colnames(r) = colnames(pb)
  r
}


rnaNNLS.by.ann = function(v,ann,normsd=TRUE,min.gcounts=500){
  x = as.matrix(v@assays$Spatial@counts)
  tot = apply(x,1,sum)
  x = x[tot>=min.gcounts,]
  pb = calcMeanCols(x[,as.logical(v$is.tissue)],ann[as.logical(v$is.tissue)])
  list(nnls=rnaNNLS(x,pb,normsd=normsd),pb=pb)
}


dist2ann = function(v,ann,weights=NULL,fun=min,norm2spot.dist=TRUE){
  dist = as.matrix(dist(v@images$slice1@coordinates[,c('imagerow','imagecol')]))
  #dist = as.matrix(dist(v@images$slice1@coordinates[,c('row','col')],m='manhattan'))
  if(norm2spot.dist)
    dist = dist/min(dist[dist>0])
  as = unique(ann[as.logical(v$is.tissue)])
  r = sapply(as,function(a){
    d =dist[,ann==a,drop=FALSE]
    if(!is.null(weights))
      w = weights[ann==a,a]
    sapply(1:nrow(d),function(i){
      if(!is.null(weights))
        fun(d[i,],w)
      else
        fun(d[i,])
      })
    })
  rownames(r) = rownames(dist)
  r
}


mean.by.bin = function(x,y,bins=NULL,...){
  if(is.null(bins))
    bins = seq(...)
  i = findInterval(x,bins,all.inside = TRUE)
  x   = sapply(split(x,i),mean)
  ysd = sapply(split(y,i),sd)
  y   = sapply(split(y,i),mean)
  cnt = table(i)
  i = as.character(1:(length(bins)-1))
  r = data.frame(x=x[i],y=y[i],ysd=ysd[i],cnt=as.numeric(cnt[i]),row.names = i)
  r[!is.na(r$y),]
}


getDiffFun = function(x,w){
  function(a)sum(a*w*exp(-(a*x)^2))/sqrt(pi)
}

estimateDiffusion = function(x,w,obs,lower=0,upper=0.5,init=0,method='L-BFGS-B'){
  f = getDiffFun(x,w)
  optim(init,function(a)(f(a)-obs)^2,lower = lower,upper = upper,method = 'L-BFGS-B')$par
}



plotWithMeanTrend = function(x,y,col='black',new=TRUE,from=min(x),to=max(x),by=1,lwd=3,min.obs=5,plot.points=TRUE,plot=TRUE,...){
  f = !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
  x = x[f]
  y = y[f]

  m = mean.by.bin(x,y,from=from,to=to,by=by)
  m$ysd = m$ysd/sqrt(m$cnt)
  m = m[!is.infinite(m$x) & m$cnt>=min.obs,]
  if(!plot) return(m)
  if(new)
    plot(x,y,col=col,type=ifelse(plot.points,'p','n'),...)
  else if(plot.points)
    points(x,y,col=col,...)

  plotArea(m$x,cbind(m$y,m$ysd),t='b',lwd=lwd,col=col,new = FALSE)
  invisible(m)
}

getWeightedMeanCols = function(nls,anscol){
  z = nls$nnls
  z = sweep(z,2,apply(nls$pb[,colnames(z)],2,sum),'*')
  z = sweep(z,1,apply(z,1,sum),'/')
  z[is.na(z)] = 0
  colrgb = col2rgb(anscol[colnames(z)])
  z = z %*% t(colrgb)
  apply(z,1,function(x)rgb(x[1],x[2],x[3],maxColorValue = 255))
}

exogenousFraq = function(n){
  rna = sweep(n$nnls,2,apply(n$pb[,colnames(n$nnls)],2,sum),'*')
  total = apply(rna,1,sum)
  endo = rep(NA,length(total))
  for(a in colnames(rna))
    endo[n$annotation==a] = rna[n$annotation==a,a]
  1 - (endo / total)
}

empty2tissueRatio = function(v,d,max.dist=Inf){
  mean(v$nCount_Spatial[v$annotation=='not_tissue' & apply(d,1,min)<=max.dist])/mean(v$nCount_Spatial[v$annotation!='not_tissue'])
}
