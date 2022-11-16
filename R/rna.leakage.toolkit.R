removeExoRNA_ = function(m,ann,nnls){
  m = as.matrix(m)
  for(a in colnames(nnls$nnls)){
    pb = rowMeans(m[,ann ==a])
    f = ann != a
    m[,f] = m[,f] - matrix(pb,ncol=1) %*% matrix(nnls$nnls[f,a],nrow=1)
  }
  m = round(m)
  m[m<0] = 0
  m = as(m,'dgCMatrix')
  m
}


removeExoRNA = function(m,ann,is.tissue,n=1,normsd = T,min.gcounts=500,min.spots=10){
  # remove small ann types
  t = table(ann)
  small.anns = names(t)[t<min.spots]
  is.tissue = as.logical(is.tissue) & !(ann %in% small.anns)
  ann[ann %in% small.anns] = 'not_tissue'
  if(sum(is.tissue) < min.spots){
    warning("Too few tissue spots! Stop procesing.")
    return(NULL)
  }

  r = list(list(vis=m))
  cat('0; total=',sum(m),'; non-zero=',length(m@x),'\n',sep='')
  m = as.matrix(m)
  for(i in 1:n){
    nnls = rnaNNLS.by.ann(m,ann,is.tissue,normsd = normsd,min.gcounts=min.gcounts)
    m = removeExoRNA_(m,ann,nnls)
    cat(i,'; total=',sum(m),'; non-zero=',length(m@x),'\n',sep='')
    nnls$vis = m
    r[[i+1]] = nnls
  }
  r
}


rnaNNLS = function(cnt,pb,normsd=TRUE){
  require(nnls)
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


rnaNNLS.by.ann = function(v,ann,is.tissue,normsd=TRUE,min.gcounts=500){
  if(is(v,'Seurat')){
    x = as.matrix(v@assays$Spatial@counts)
    is.tissue = v$is.tissue
  }else
    x = (v)

  is.tissue = as.logical(is.tissue)

  tot = rowSums(x)
  x = x[tot>=min.gcounts,]
  pb = calcColSums(x[,is.tissue],ann[is.tissue],mean = TRUE)
  list(nnls=rnaNNLS(x,pb,normsd=normsd),pb=pb)
}


dist2ann = function(xy,ann,weights=NULL,fun=min,norm2spot.dist=TRUE){
  dist = as.matrix(dist(xy[,c('imagerow','imagecol')]))
  if(norm2spot.dist)
    dist = dist/min(dist[dist>0])
  as = unique(ann)
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
  colnames(r) = as
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
