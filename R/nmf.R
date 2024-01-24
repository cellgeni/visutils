#' Plots NMF clustering, silhouette and consensuse clustering
#'
#' @param coefs coef matrix from NMF
#' @param cons conensus matrix
#' @param clcols colors for clusters
#' @param max.cex dot size for dotPlot
#' @param colfun color function to be used in dotPlot
#' @param title title for whole plot
#' @param ylab.cex,xlab.cex sizes of axis labels
#'
#' @return list of clusters, colors and silhouette
#' @export
plotNMFCons = function(coefs,cons,clcols=NULL,max.cex = 4.5/14*8,
                       colfun=function(x)num2col(x,c('blue','gray','orange','violet','black'),minx = 0,maxx = 1),
                       title='',ylab.cex=1,xlab.cex=1){
  cons = cons[colnames(coefs),colnames(coefs)]
  cls = apply(coefs,2,which.max)
  if(is.null(clcols))
    clcols = RColorBrewer::brewer.pal(nrow(coefs),'Set3')

  o = names(cls)[order(cls)]
  layout(matrix(1:3,ncol=3),widths = c(0.9,0.3,1.8))
  parl = par(no.readonly=TRUE)
  par(mar=c(3,9,2.5,0),bty='n',oma=c(0,0,1,0),cex=1,tcl=-0.2,mgp=c(1.2,0.3,0),las=1)

  dotPlot(t(coefs[,o]),rowColours = cbind(clcols[cls[o]]),colColours = clcols,max.cex = max.cex,ylab.cex = ylab.cex,xlab.cex = xlab.cex,
          scaleWM=F,colfun=colfun,plot.legend=FALSE)

  slh = cluster::silhouette(cls,dmatrix=1-cons)
  rownames(slh) = names(cls)


  par(mar=c(3,0.2,2.5,0.3))
  b=barplot(slh[rev(o),3],horiz = T,width = 1,space = 0, yaxs = "i",ylim=c(-1.5,length(o)+0.5),names.arg = '',border=NA,col=clcols[cls[rev(o)]],xlab='silhouette',
            xlim=c(-1,1),xaxt='n')
  par(cex=0.6)
  axis(1)
  par(cex=1,xaxt='n',yaxt='n')
  dotPlot(cons[o,o],max.cex = 1.3,rowColours = cbind(clcols[cls[o]]),colColours = clcols[cls[o]],main='Consensus clustering matrix',plot.legend=FALSE)

  mtext(title,3,outer = TRUE,line = 0)
  par(parl)
  layout(matrix(1,ncol=1,nrow=1))
  invisible(list(clusters=cls,cols=clcols,slh=slh))
}

#' Get function to normalize NMF results
#'
#' same NMF deconvolution can be expressed in infinite number of ways (different by scaling of factors), one need to define either normalization or regularization to unambiguously choose one.
#'
#' @param type either 'no' for no normalization or 'max' for max-normalization of basis
#'
#' @return function that takes output of nmf function as input and calculates normalization coeficients
#' @export
getNMFNormFs = function(type){
  if(type == 'no')
    res = function(n){rep(1,ncol(basis(n)))}
  if(type == 'max')
    res = function(n){apply(basis(n),2,max)}
  res
}

#' Get clusterings from list of nmfs
#'
#' @param nmfs list of nmf outputs
#' @param normF normalization function (output of getNMFNormFs)
#'
#' @return list of clusterings
#' @export
getNMFcl = function(nmfs,normF){
  lapply(nmfs,function(n){
    c = coef(n)
    f = normF(n)
    c = sweep(c,1,f,'*')
    apply(c,2,which.max)
  })
}


#' Calculate consensus matrix for list of clusterings
#'
#' @param cs list of clusterings (output of getNMFcl)
#'
#' @return matrix n_celltypes*n_celltypes, with number of times given celltype pair was assigned to same cluster
#' @export
getConsensus = function(cs){
  c = cs[[1]]
  r = matrix(0,ncol=length(c),nrow=length(c),dimnames=list(names(c),names(c)))
  for(c in cs){
    for(i in 1:length(c))
      for(j in 1:length(c))
        r[i,j] = r[i,j] + (c[i]==c[j])
  }
  r
}



#' Assesses clusterings against distance matrix
#'
#' @param cls list of clusterings
#' @param d distance matrix (for example 1 - coincedence frequency)
#'
#' @return data.frame with silhouette summaries
getSilhouetteStat = function(cls,d){
  do.call(rbind,lapply(cls,function(cl){
    x = cluster::silhouette(cl,d)[,3]
    data.frame(min=min(x),mean=mean(x),median=median(x),n0 = sum(x<0))
  }))
}


nmf2cl = function(n,normF){
  c = coef(n)
  b = basis(n)
  f = normF(n)
  c = sweep(c,1,f,'*')
  b = sweep(b,2,f,'/')
  cl = apply(c,2,which.max)
  c = sweep(c,2,apply(c,2,sum),'/')
  list(nmf=n,coefn = c,basisn=b,cl=cl)
}


#' Finds best nmf run
#'
#' by 1) its agreement (mean Silhouette) with consensus matrix calculated based on all runs 2) deviance
#'
#' the one that agree best with coincedence matrix
#'
#' @param nmfs list of nmf results
#' @param normF normalization function (output of getNMFNormFs)
#'
#' @return list containing:
#' best nmf
#' normalized basis and coedicients (these are in addition per celltype sum-normalized)
#' celltype clustering
#' consensus matrix
#' Silhouette statistics
#' intex of best run
#'
#' @export
nmfGetBest = function(nmfs,normF){
  cls = getNMFcl(nmfs,normF)
  cns = getConsensus(cls)
  slh = getSilhouetteStat(cls,1-cns/length(nmfs))
  slh$dev = sapply(nmfs,deviance)
  i = order(-slh$mean,slh$dev)[1]
  n = nmfs[[i]]

  # assemble output object
  c = coef(n)
  b = basis(n)
  f = normF(n)
  c = sweep(c,1,f,'*')
  b = sweep(b,2,f,'/')
  cl = apply(c,2,which.max)

  c = sweep(c,2,apply(c,2,sum),'/')

  r = list(nmf=n,coefn = c,basisn=b,cl=cl)

  r$cons=cns
  r$slh.stat=slh
  r$besti=i
  r
}
