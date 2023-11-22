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

  dotPlot(t(coefs[,o]),rowColours = cbind(clcols[cls[o]]),colColours = clcols,max.cex = max.cex,ylab.cex = ylab.cex,xlab.cex = xlab.cex,scaleWM=F,colfun=colfun,plot.legend=FALSE)

  slh = cluster::silhouette(cls,dmatrix=1-cons)
  rownames(slh) = names(cls)


  par(mar=c(3,0.2,2.5,0.3))
  b=barplot(slh[rev(o),3],horiz = T,width = 1,space = 0, yaxs = "i",ylim=c(-1.5,length(o)+0.5),names.arg = '',border=NA,col=clcols[cls[rev(o)]],xlab='silhouette',xlim=c(-1,1),xaxt='n')
  par(cex=0.6)
  axis(1)
  par(cex=1,xaxt='n',yaxt='n')
  dotPlot(cons[o,o],max.cex = 1.3,rowColours = cbind(clcols[cls[o]]),colColours = clcols[cls[o]],main='Consensus clustering matrix',plot.legend=FALSE)

  mtext(title,3,outer = TRUE,line = 0)
  par(parl)
  layout(matrix(1,ncol=1,nrow=1))
  invisible(list(clusters=cls,cols=clcols,slh=slh))
}
