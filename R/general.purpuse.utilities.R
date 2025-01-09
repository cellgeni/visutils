#' String concatenation operator
#'
#' @export
`%p%` = paste0

#' Compare two lists of sets
#'
#' return matrix of set similarities defined by fun.
#'
#' Default is size of intersection divided by minimal size (overlap or Szymkiewicz–Simpson coefficient).
#'
#'
#' @param l1 list of sets to be compared
#' @param l2 list of sets to be compared (compare l1 with itself if l2 is not set)
#' @param fun function to calculate similarity or character. 'o' will result in intersect/min; 'j' in intersect/union
#'
#' @return matrix if length(l1),length(l2) size
#' @export
#'
#' @examples
#' compareSets(list(letters[1:6],1:5),list(letters[1:8],4:8,letters[4:10]),fun = 'j')
compareSets = function(l1,l2=l1,fun='o'){
  if(is.character(fun) && fun=='o')
    fun = function(x,y){length(intersect(x,y))/min(length(x),length(y))}
  if(is.character(fun) && fun=='j')
    fun = function(x,y){length(intersect(x,y))/length(union(x,y))}
  m = matrix(NA,nrow=length(l1),ncol=length(l2),dimnames = list(names(l1),names(l2)))
  for(i in 1:length(l1)){
    for(j in 1:length(l2)){
      m[i,j] = fun(l1[[i]],l2[[j]])
    }
  }
  m
}

#' log10(1+x)
#'
#' Same as log1p but with base equal to 10
#'
#' @param x numeric vector
#'
#' @return log10(1+x)
#' @export
#'
log10p1 = function(x){
  log1p(x)/log(10)
}

#' Cuts part of delimited text
#' analogous to bash cut
#'
#' @param x character vector to cut from
#' @param del delimiter
#' @param inx indexes of items to return
#' @param fixed logical, whether del is fixed (to be passed to strsplit)
#' @param collapse logical, whether to collapse selected items with same delimeter.
#' @param simplify logical, whether output should be simplified (to be passed to sapply).
#'
#' @return character vector or
#' @export
#'
#' @examples
#' splitSub(c('a,b','d,c'),',',2)
splitSub = function(x,del,inx,fixed=TRUE,collapse=TRUE,simplify=TRUE){
  r = sapply(strsplit(x,del,fixed=fixed),'[',inx,simplify = simplify & (!collapse | length(inx)==1))
  if(collapse & length(inx)>1)
    r = sapply(r,paste,collapse=del)
  r
}


#' Calculate weighted mean colour
#'
#' Weights are normalized per row by total
#'
#' @param cols vector of colors (any notation)
#' @param weights matrix with number of columns equal to the length of \code{cols} (in same order).
#'
#' @return vector of colours (length equal to nrow of weights)
#' @export
#'
weightedColourMeans = function(cols,weights){
  weights = sweep(weights,1,apply(weights,1,sum),'/')
  na = apply(is.na(weights),1,sum)>0
  weights[is.na(weights)] = 0
  colrgb = col2rgb(cols)
  z = weights %*% t(colrgb)
  r = apply(z,1,function(x)rgb(x[1],x[2],x[3],maxColorValue = 255))
  r[na] = NA
  r
}

#' 2D density plot
#'
#' Replacement of scatterplot to be used when number of points is to high and/or they highly overlap
#'
#' @param x,y point coordinates
#' @param xbins,ybins number of bins or bin borders. Everything outside of bins are removed
#' @param cols colours to form gradient
#' @param zfun transformation of frequency (try \code{\link{log1p}}, \code{\link{sqrt}}, etc). Identity by default
#' @param leg.title legend title
#' @param num.leg.tic desired number of tics in the legend
#' @param legend logical, should legend be plotted
#' @param trimZq frequencies outside of quantile(trimZq,1-trimZq) are set to the corresponding (closest) quantiles. Default is 0 (not trimming). See \code{\link{trimQ}}.
#' @param xlim,ylim,zlim,xlab,ylab graphical parameters
#' @param new create new plot (default) or add to existing
#' @param ... other arguments for plot function
#'
#' @return list bith bins and frequencies (invisible)
#' @export hist2D
hist2D = function(x,y,xbins=100,ybins=100,cols=c('white','gray','blue','orange','red'),zfun=identity,
                  leg.title='',num.leg.tic=NULL,legend=TRUE,trimZq=0,xlim=NULL,ylim=NULL,zlim=NULL,
                  new=TRUE,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),...){
  xlab;ylab;
  f = !is.na(x) & !is.infinite(x) & !is.na(y) & !is.infinite(y)
  if(!is.null(xlim)) f = f & x >= xlim[1] & x <= xlim[2]
  if(!is.null(ylim)) f = f & y >= ylim[1] & y <= ylim[2]
  x = x[f]
  y = y[f]
  if(length(xbins)==1) xbins = seq(min(x),max(x),length.out = xbins+1)
  if(length(ybins)==1) ybins = seq(min(y),max(y),length.out = ybins+1)

  f = x >= xbins[1] & x <= xbins[length(xbins)] & y >= ybins[1] & y <= ybins[length(ybins)]
  x = x[f]
  y = y[f]

  xb = findInterval(x,xbins,rightmost.closed = TRUE)
  yb = findInterval(y,ybins,rightmost.closed = TRUE)


  z = as.matrix(table(factor(yb,levels=1:(length(ybins)-1)),factor(xb,levels=1:(length(xbins)-1))))
  z = trimQ(z,trimZq)
  if(is.null(xlim)) xlim = range(x)
  if(is.null(ylim)) ylim = range(y)
  if(is.null(zlim)) zlim = range(z)
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  z2col=function(x)num2col(zfun(x),cols,minx=zfun(zlim[1]),maxx=zfun(zlim[2]))
  rect(rep(xbins[-length(xbins)],each =length(ybins)-1),
       rep(ybins[-length(ybins)],times=length(xbins)-1),
       rep(xbins[-1]            ,each =length(ybins)-1),
       rep(ybins[-1]            ,times=length(xbins)-1),
       col=z2col(z),border = NA)
  if(legend){
    z2col=function(x)num2col(x,cols)
    plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0,'npc','nfc'),grconvertY(1,'npc','nfc'),fullzlim = zlim,
                     zlim = range(z),zfun = zfun,z2col=z2col,leg=num.leg.tic,title=leg.title)
  }
  invisible(list(x=(xbins[-1]+xbins[-length(xbins)])/2,
                 y=(ybins[-1]+ybins[-length(ybins)])/2,
                 z=z,xbins=xbins,ybins=ybins))
}


#' Estimates density at specified 2D points
#'
#' similar to MASS::kde2d but calculates density for supplied points
#'
#' @param x,y data coordinates
#' @param kernel kernel to be used (dnorm by default)
#' @param approx logical, should density be estimated on subset of point. Approximate mode is much faster for large (length(x)>5000) datasets
#' @param random_seed random seed to chose subset in approximate mode
#' @param k subset size to use in approximate mode. 1000 is good starting point
#'
#' @return numeric vector with density estimates
#' @export
#'
#' @examples
#' x = rnorm(10000,mean = c(-300,300,300),sd=100)
#' y = rnorm(10000,x*c(-1,1),sd=100)
#' system.time(d1 <- pointKde2d(x,y,approx = F))
#' system.time(d2 <- pointKde2d(x,y,approx = T,k=100))
#' system.time(d3 <- pointKde2d(x,y,approx = T,k=1000))
#'
#' par(mfrow=c(2,2))
#' plot(x,y,col=num2col(d1),pch=19,cex=0.2)
#' plot(x,y,col=num2col(d2),pch=19,cex=0.2)
#' plot(x,y,col=num2col(d3),pch=19,cex=0.2)
#'
#' plot(d1,d3,pch=16,cex=0.2)
#' abline(a=0,b=1,col='red')
pointKde2d = function(x,y,kernel=dnorm,approx=length(x)>2000,random_seed=123,k=min(length(x),1000)){
  require(MASS)
  r = c()
  if(approx){
    set.seed(random_seed)
    neighs = sample(length(x),k)
    xn = x[neighs]
    yn = y[neighs]
  }else{
    xn = x
    yn = y
  }

  h = c(bandwidth.nrd.not0(xn), bandwidth.nrd.not0(yn))
  h = h/4
  xn = xn/h[1]
  yn = yn/h[2]
  x = x/h[1]
  y = y/h[2]

  for(i in 1:length(x)){
    xx = (x[i]-xn)
    yy = (y[i]-yn)
    r[i] = sum(kernel(xx)*kernel(yy))
  }
  r/(length(xn) * h[1L] * h[2L])
}

#' Identify optimal bandwidth via normal distribution
#'
#' ensure that result is not zero: removes the most frequent value in x if bandwidth.nrd return zero
#'
#' @param x a data vector
#'
#' @return A bandwidth on a scale suitable for the width argument of density.
#' @export
bandwidth.nrd.not0 = function (x) {
  repeat{
    if(length(x)<2) return(NA)
    r = bandwidth.nrd(x)
    if(r>0) return(r)
    t = sort(table(x),decreasing = T)
    if(t[1]/length(x) > 0.5){
      x = x[x != names(t)[1]]
    }
  }
}

#' Scatterplot wih point density shown by colour
#'
#' @param x,y point coordinates
#' @param pch,bty,log parameters of plot function
#' @param ... other parameters to be passed to plot
#'
#' @export
plotPointDensity = function(x,y,pch=16,bty='n',log='',...){
  x.=x
  y.=y
  pch=recycle(pch,length(x))
  if(grepl('x',log)){
    x. = log(x.)
  }
  if(grepl('y',log)){
    y. = log(y.)
  }
  c = pointKde2d(x.,y.)
  o = order(c)
  plot(x[o],y[o],col=num2col(c[o]),pch=pch[o],bty=bty,log=log,...)
}


#' Trim by quantiles
#'
#' Sets all values of x below qth quantile to the quantile. Sets all values of x above (1-q)th quantile to the quantile.
#'
#' @param x numeric vector
#' @param q quantile
#'
#' @return corrected x
#' @export
#'
trimQ = function(x,q){
  if(q==0) return(x)
  qq=quantile(x,sort(c(q,1-q)))
  x[x<=qq[1]] = qq[1]
  x[x>=qq[2]] = qq[2]
  x
}

#' Creates palette for character vector
#'
#' Sorts t numerically if possible, alphabetically otherwise.
#'
#' @param t input character vector
#' @param bpal RColorBrewer pallete to be used if number of colors allows
#' @param colfun pallete function to be used in number of colors 10 or more
#' @param palette logical, specifies whether pallete (color per unique value in t) or colour along t should be returned
#' @param random.seed seed to be set befor color generation (for stochastic colfuns)
#'
#' @return names colour vector
#' @examples
#' char2col(c('a','a','b'))
#' char2col(c(1,1,11,2))
#' char2col(c(T,F))
#' @export
char2col = function(t,bpal='Set1',colfun=randomcoloR::distinctColorPalette,palette=TRUE,random.seed=1234){
  set.seed(random.seed)
  require(randomcoloR)
  t = as.character(t)
  torig = t
  t = sort(unique(t))
  suppressWarnings({
    if(all(!is.na(as.numeric(t)))){
      t = as.character(sort(as.numeric(t)))
    }
  })
  if(length(t) <= RColorBrewer::brewer.pal.info[bpal,'maxcolors'])
    r=setNames(RColorBrewer::brewer.pal(max(3,length(t)),bpal),t)[1:length(t)]
  else if(length(t)<=2000){
    r=setNames(colfun(length(t)),t)
  }else{
    r=setNames(randomcoloR::randomColor(length(t)),t)
  }
  if(!palette)
    r = r[torig]
  r
}

#' Maps distance matrix into colour space
#'
#' use MDS (via cmdscale) to map objects into colour space
#'
#' @param d distance matrix to be used with MDS
#' @param use3D logical, specifies whether whole [0,1]^3 space should be used. About equally-bright 2D space is used otherwise (default)
#' @param orderBySim logical, whether objects should be ordered by similarity. Preserves original order otherwise (default)
#'
#' @return vector of colors
#' @export
#'
#' @examples
#' d=1-cor(t(mtcars))
#' c=getColoursByDistance(d,F)
#' plot(cmdscale(d),col=c,pch=19)
getColoursByDistance = function(d,use3D=FALSE,orderBySim=FALSE){
  if(use3D){
    mds = cmdscale(d,k=3)
    mds[,1] = scaleTo(mds[,1])
    mds[,2] = scaleTo(mds[,2])
    mds[,3] = scaleTo(mds[,3])
    col=rgb(mds[,1],mds[,2],mds[,3],maxColorValue = 1)
  }else{
    mds = cmdscale(d,k=2)
    mds[,1] = scaleTo(mds[,1])
    mds[,2] = scaleTo(mds[,2])
    col=rgb(mds[,1],mds[,2],scaleTo(-mds[,1]-mds[,2]),maxColorValue = 1)
  }
  col=setNames(col,rownames(mds))
  if(orderBySim){
    o = cmdscale(d,k=1)
    col = col[order(o[,1])]
  }
  col
}

#' Recycle vector v along i
#'
#' @param v vector to be recycled
#' @param i either number of values to be obtained or 1:number_of_values
#'
#' @return recycled vector
#' @export
recycle = function(v,i){
  if(length(i)==1)
    i = 1:i
  v[(i-1)%%length(v)+1]
}

#' Transforms numbers to colours
#'
#' Wrapper function for colorRamp. Transforms numbers to the gradient specified by col parameter
#'
#' @param d numbers to be transformed
#' @param col colours to form gradient (blue-red by default)
#' @param minx,maxx range of x values (min(x),max(x) by default)
#'
#' @return vector of colours (same length as d)
#' @export
num2col = function(d,col=c('blue','cyan','gray','orange','red'),minx=min(d,na.rm = TRUE),maxx=max(d,na.rm = TRUE)){
  if(sd(d,na.rm=TRUE)==0)
    return(rep(col[1],length(d)))
  d[d<minx] = minx
  d[d>maxx] = maxx
  bwr = colorRamp(col,alpha=TRUE)
  apply(bwr(scaleTo(d,0,1,minx=minx,maxx=maxx))/255,1,function(x){
    if(is.na(x[1])) return(NA)
    rgb(x[1],x[2],x[3],x[4])
  })
}

#' Tests whether value is colour
#'
#' @param x vector to be tested
#'
#' @return logical vector
#' @export
isColors <- function(x) {
  if(is.factor(x)) return(rep(FALSE,length(x)))
  x %in% colors() | (substr(x,1,1)=='#' & (nchar(x)==7 |  nchar(x)==9))
}

#' Plot heatmap with numbers
#'
#' Wraps image functions
#'
#' @param d matrix to be plotted
#' @param t text to be ploted on top of image (rounds d by default)
#' @param digits number of digits retained under rounding of d
#' @param text.col color for text
#' @param xaxlab,yaxlab axis labels, dimnames(d) by default
#' @param centerColors0 logical, specifies whether 0 should be considered as palette center (makes sense for correlation matrices)
#' @param las orientation of axis labels (see las in ?par)
#' @param text.xadj text adjastment by x
#' @param colAnns,rowAnns list (or matrix of columns) of column (row) annotation to be shown by color
#' @param colAnnsCols,rowAnnsCols - colors to be used for col (row) annotation. Defined by char2color if null.
#' @param legend.cex.at,legend.col.at values to be used in legend, set both to have two independent legends for size and colour
#' @param legend.cex.title,legend.col.title titles of legends
#' @param rowAnnWidth,colAnnWidth - size of colour annotations as fraction of plot area
#' @param annSpacer - spacer between heatmap and annotation as fraction of plot area
#' @param ... other options to be supplied to image
#'
#' @export
#' @examples
#' par(mar=c(10,6,1,1),bty='n')
#' imageWithText(mtcars[1:4,1:3],rowAnns=list(r1=c(1,1,2),r2=c('a','a','a')),rowAnnCols=list(r1=c('1'='red','2'='blue'),r2=c(a='green')),colAnns = list(c1=c(1,1,1,2),c2=c(2,1,2,1),c3=c('a','b','b','b')))
imageWithText = function(d,t=NULL,digits=2,text.col=NULL,xaxlab=rownames(d),yaxlab=colnames(d),centerColors0=FALSE,las=2,text.xadj=0.5,
                         colAnns=NULL,rowAnns=NULL,
                         colAnnCols=NULL,rowAnnCols=NULL,
                         col= hcl.colors(100, "YlOrRd", rev = TRUE),
                         rowAnnWidth=0.7,colAnnWidth=0.7,annSpacer=0.1,...){

  d = as.matrix(d)
  if(is.null(t))
    t = round(d,digits = digits)
  pars = list(...)
  if(is.null(pars$col)) pars$col= num2col(1:100)
  pars$col = col
  pars$x = 1:nrow(d)
  pars$y = 1:ncol(d)

  if(is.null(pars$xlab)) pars$xlab=''
  if(is.null(pars$ylab)) pars$ylab=''
  pars$z = d
  pars$xaxt='n'
  pars$yaxt='n'
  if(!is.null(col) & centerColors0){
    m = max(abs(d),na.rm=T)
    if(is.numeric(centerColors0))
      m = centerColors0
    pars$breaks = seq(-m,m,length.out = length(pars$col)+1)
  }
  do.call(image,pars)
  if(is.null(text.col))
    text.col = 'black'
  text(rep(pars$x,times=length(pars$y))+scaleTo(text.xadj,-0.45,0.45,minx=0,maxx=1),rep(pars$y,each=length(pars$x)),t,col=text.col,adj=c(text.xadj,0.5))

  mgp = par('mgp')
  if(!is.null(colAnns)){
    if(is.list(colAnns))
      n = length(colAnns)
    else
      n = ncol(colAnns)
    par(mgp=mgp+n*colAnnWidth+annSpacer)
    plotColorAnn(pars$x,colAnns,cex=colAnnWidth,horis=TRUE,col=colAnnCols,spacer = annSpacer)
  }

  if(!is.null(xaxlab)){
    axis(1,pars$x,xaxlab,las=las)
  }

  if(!is.null(rowAnns)){
    if(is.list(rowAnns))
      n = length(rowAnns)
    else
      n = ncol(rowAnns)
    par(mgp=mgp+n*rowAnnWidth+annSpacer)
    plotColorAnn(pars$y,rowAnns,cex=rowAnnWidth,horis=FALSE,col=rowAnnCols,spacer = annSpacer)
  }
  if(!is.null(yaxlab))
    axis(2,pars$y,yaxlab,las=las)

  par(mgp=mgp)
}

#' Returns with of line in user coordinates
#'
#' @param axis character, either 'x' or 'y'
#'
#' @return with of line in user coordinates
#' @export
#'
#' @examples
#' getLineWidth('x')
getLineWidth = function(axis){
  if(!(axis %in% c('x','y')))
    stop('axis shold be eitehr x or y')
  if(axis=='x')
    r = grconvertX(1:2,'lines','user')
  if(axis=='y')
    r = grconvertY(1:2,'lines','user')
  abs(r[1]-r[2])
}

#' Plots matrix as dotPlot
#'
#' Shows values in matix by point size and color. Size matrix (m) is not scaled by default, dot size is defined as m*max.cex
#'
#' @param m numeric matrix to be shown as dot size
#' @param mc numeric matrix to be shown as dot colour (uses m as default)
#' @param rfun function to calculate radius from matrix values. Use sqrt (default) to make area proportional to value
#' @param colfun function to transform values to colour gradient
#' @param grid logial, should grid be plotted
#' @param grid.lty lty of grid
#' @param grid.col line color of grid
#' @param max.cex max size of dots
#' @param xlab,ylab axis labels
#' @param ylab.cex magnification label for ylabs
#' @param colColours colour matrix to plot annotation for columns (matrix with nrow equal to ncol(m); each column of colColours is annotation)
#' @param rowColours colour matrix to plot annotation for row (matrix with nrow equal to nrow(m); each column of rowColours is annotation)
#' @param scaleWM logical, specifies wheter computed radiuses should be scaled into [0,1] interval (FALSE by default).
#' @param pch point character (19 by default)
#' @param plot.legend logical, whether legend should be plotted. Single legend will be plotted if m is identicall to mc.
#' @param legend.cex.at,legend.col.at values to be used in legend, set both to have two independent legends for size and colour
#' @param legend.cex.title,legend.col.title titles of legends
#' @param rowAnnWidth,colAnnWidth - size of colour annotations in user coordinates
#' @param ... other parameters to be passed to plot function
#'
#' @export
#' @examples
#' c = matrix(1:12,ncol=3)
#' par(mar=c(4,4,1,10),bty='n')
#' dotPlot(c/12,-c,max.cex = 3,colColours = cbind(col1=c('red','blue','red'),col2=c('green','green','magenta')),
#'         rowColours = cbind(rrr1=c('red','blue','blue','red'),rrr2=c('green','green','magenta','green')),
#'       legend.cex.title='size',legend.col.title='col',
#'       colAnnWidth = 0.5,
#'       rowAnnWidth = 0.5)
dotPlot = function(m,mc=m,rfun=sqrt,colfun=function(x)num2col(x,c('white','yellow','violet','black')),grid=TRUE,grid.lty=2,grid.col='gray',
                   max.cex=1,xlab='',ylab='',ylab.cex=1,xlab.cex=1,
                   colColours=NULL,rowColours=NULL,
                   rowAnnWidth=1,colAnnWidth=1,
                   scaleWM=FALSE,pch=19,plot.legend=TRUE,
                   legend.cex.at=NULL,legend.col.at=legend.cex.at,
                   legend.cex.title='',legend.col.title='',...){
  x = rep(1:ncol(m),each=nrow(m))
  y = rep(nrow(m):1,times=ncol(m))

  xlim=c(0,max(x)+1)
  ylim=c(0,max(y)+1)

  if(!is.null(colColours)){
    if(!is.array(colColours))
      colColours = matrix(colColours,ncol=1)
    ylim[1] = -ncol(colColours)*colAnnWidth
  }

  if(!is.null(rowColours)){
    if(!is.array(rowColours))
      rowColours = matrix(rowColours,ncol=1)
    xlim[1] = -ncol(rowColours)*rowAnnWidth
  }


  plot(1,t='n',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,yaxs='i',xaxs='i',...)
  if(grid){
    abline(v=1:ncol(m),lty=grid.lty,col=grid.col)
    abline(h=1:nrow(m),lty=grid.lty,col=grid.col)
  }
  wh = par('cin')
  # sym.size=max(c(grconvertX(wh[1],'in','user')-grconvertX(0,'in','user'),grconvertY(wh[2],'in','user')-grconvertY(0,'in','user')))
  # r = scaleTo(rfun(m))/sym.size*max.cex
  r = rfun(m)
  if(scaleWM)
    r = scaleTo(r)
  r = r*max.cex
  points(x,y,cex=r,col=colfun(mc),pch=pch)

  # add col annotation
  f = colAnnWidth
  if(!is.null(colColours)){
    for(i in 1:ncol(colColours))
      for(j in 1:nrow(colColours)){
        rect(j-0.5,0-i*f,j+0.5,f-i*f,border = NA,col=colColours[j,i])
      }

    if(!is.null(colnames(colColours))){
      text(nrow(colColours)+1,-(1:ncol(colColours))*f+f/2,colnames(colColours),adj=c(0,0.5),xpd=T)
    }
  }

  # add row annotation
  f = rowAnnWidth
  if(!is.null(rowColours)){
    for(i in 1:ncol(rowColours))
      for(j in 1:nrow(rowColours)){
        rect(0-i*f,nrow(rowColours)-j+1.5,f-i*f,nrow(rowColours)-j+0.5,border = NA,col=rowColours[j,i])
      }
    if(!is.null(colnames(rowColours))){
      text(-(1:ncol(rowColours))*f+f/2,0,colnames(rowColours),srt=90,adj=c(1,0.5),xpd=T)
    }
  }
  if(par('xaxt')=='s'){
    par.out = par(cex=xlab.cex)
    axis(1,1:ncol(m),colnames(m))
    do.call(par,par.out)
  }
  if(par('yaxt')=='s'){
    par.out = par(cex=ylab.cex)
    axis(2,nrow(m):1,rownames(m))
    do.call(par,par.out)
  }
  # legend
  legend.col.at
  if(plot.legend){
    x = grconvertX(1,'npc','user')
    y = grconvertY(1,'npc','user')

    if(is.null(legend.cex.at))
      legend.cex.at = round(seq(min(m),max(m),length.out = 5),digits = 4)
    if(is.null(legend.col.at))
      legend.col.at = round(seq(min(mc),max(mc),length.out = 5),digits = 4)


    legend.cex = rfun(legend.cex.at)
    if(scaleWM)
      legend.cex = scaleTo(legend.cex,minx = rfun(min(m)),maxx = rfun(max(m)))

    legend.cex = legend.cex * max.cex
    legend.col = colfun(c(min(mc),max(mc),legend.col.at))[-1:-2]

    cexx = legend.cex
    coll = 'black'
    if(all(m==mc))
      coll = legend.col
    labb = legend.cex.at
    l = legend(x,y,xpd=NA,pch=19,pt.cex=cexx,col=coll,legend = labb,title=legend.cex.title,bty=par('bty'))
    if(!all(m==mc)){
      cexx = 1
      coll = legend.col
      if(all(m==mc))
        coll = legend.col
      labb = legend.col.at
      legend(x,l$rect$top-l$rect$h,xpd=NA,pch=19,pt.cex=cexx,col=coll,legend = labb,title=legend.col.title,bty=par('bty'))
    }
  }
}

#' Adds color annotation to plot
#'
#' Plots annotation in margins, supposed to be used with image (or imageWithText)
#'
#' @param centers coordinates along dimension to be annotated (centers)
#' @param labs list of character vectors of labels or matrix of columns. Each item in the list results in one annotation. List names will be used to label annotation.
#' @param cex width (in lines) of each annotation (0.7 is default)
#' @param cols list  of named color vectors specifying colors for each label, should have same length as labs. Generated by char2col if NULL
#' @param horis logical, whether annotation should be horizontal or vertical
#' @param spacer spacer (in lines) between plot and annotation
#'
#' @export
#'
#' @examples
#' par(mgp=c(3,2,1.5))
#' image(1:5,1:2,matrix(1:10,ncol=2))
#' plotColorAnn(1:5,list(a1=c(1,1,2,3,3),a2=c('a','b','a,','b','a')),horis=T)
#' plotColorAnn(1:2,cbind(y1=c(1,2),y2=c('a','a')),horis=F)
plotColorAnn = function(centers,labs,cex=0.7,cols=NULL,horis = FALSE,spacer=0.1){
  xpd = par('xpd'=NA)
  if(is.array(labs))
    labs = as.data.frame(labs)
  labs = lapply(labs,as.character)
  line.with = getLineWidth(ifelse(horis,'y','x'))
  width.range = c(-cex*length(labs)-spacer,-spacer)*line.with
  if(horis)
    width.range = width.range + grconvertY(0,'npc','user')
  else
    width.range = width.range + grconvertX(0,'npc','user')

  if(is.null(cols)){
    cols = lapply(labs,char2col)
  }
  l = abs(centers[2]-centers[1])
  w = (width.range[2] - width.range[1])/length(labs)
  centers = centers - l/2

  for(i in seq_along(labs)){
    if(horis){
      rect(centers,width.range[1]+(i-1)*w,centers+l,width.range[1]+i*w,border=NA,col=cols[[i]][labs[[i]]])
      text(max(centers)+l,width.range[1]+(i-0.5)*w,names(labs)[i],adj=c(0,0.5),cex=cex)
    }else{
      rect(width.range[1]+(i-1)*w,centers,width.range[1]+i*w,centers+l,border=NA,col=cols[[i]][labs[[i]]])
      text(width.range[1]+(i-0.5)*w,min(centers),names(labs)[i],adj=c(1,0.5),srt=90,cex=cex)
    }
  }
  par(xpd)
}


#' Sums two colours with alpha channel
#'
#' @param a top colour (numeric vector of length four)
#' @param b bottom colour (numeric vector of length four)
#'
#' @return vector of length four with resultant colour
#' @export
overlayTwoColours = function(a,b){
  # a over b
  if(a[4] == 0 && b[4] == 0)
    return((a+b)/2)
  a0 = a[4] + b[4]*(1-a[4])
  r = c(sapply(1:3,function(i){
    (a[i]*a[4]+b[i]*b[4]*(1-a[4]))/a0
  }),a0)
  r
}

#' Calculate per-row mean colour according to opacity
#'
#' @param c text matrix with colours to be overlayed per row
#' @param reorderByOpacity logical, should colours be ordered by increased opacity prior to summing
#'
#' @return vector of colors
#' @export
overlayColours = function(c,reorderByOpacity=FALSE){
  apply(c,1,function(x){
    m = col2rgb(x,alpha = TRUE)/255
    if(reorderByOpacity)
      m = m[,order(m[4,])]
    r = m[,1]
    for(i in 2:ncol(m))
      r = overlayTwoColours(m[,i],r)
    rgb(r[1],r[2],r[3],r[4],maxColorValue = 1)
  })
}

#' Transforms colour to hex representation
#'
#' @param c character vector with colours
#' @param withAlpha
#'
#' @return color in hex form
#' @export
col2hex = function(c,withAlpha = TRUE){
  r = apply(col2rgb(c,alpha = TRUE),2,function(x){
    if(!withAlpha)
      x[4] = 255
    rgb(x[1],x[2],x[3],x[4],maxColorValue = 255)
  })
  if(!withAlpha)
    r = substr(r,1,7)
  r
}



#' Gradient legend
#'
#' Wrapper for plotColorLegend to plot legend for given range
#'
#' @param x0,x1,y0,y1 rectange coordinates of the legend in nfc coordinates
#' @param fullzlim numeric vector with two items. Full range to be considered.
#' @param zlim numeric vector with two items. Range to show, trimmed to fullzlim if wider.
#' @param zfun transformation for value (gradient will be drawn along transformed value). Identity, log1p, sqrt, ets.
#' @param z2col function to transform numbers to colors
#' @param N number of steps in the gradient
#' @param ntic desired number of tics
#' @param leg tics values, if NULL (default) estimated automatically
#' @param title legend title
#' @param title.adj legend title adj parameter to be passed to text function
#'
#' @export
plotColorLegend2 = function(x0=NULL,x1=NULL,y0=NULL,y1=NULL,zlim,fullzlim=zlim,zfun=identity,z2col=num2col,N=100,ntic=5,leg=NULL,title=NULL,title.adj=c(0,-0.5),horizontal=FALSE){
  if(is.null(x0)){
    if(horizontal){
      x0=grconvertX(0,'npc','nfc')
      x1=grconvertX(1,'npc','nfc')
      y0=grconvertY(0,'npc','nfc')
      y1=0
    }else{
      x0=grconvertX(1,'npc','nfc')
      x1=1
      y0=grconvertY(0.1,'npc','nfc')
      y1=grconvertY(0.9,'npc','nfc')
    }
  }
  if(zlim[1]<fullzlim[1])
    zlim[1]=fullzlim[1]
  if(zlim[2]>fullzlim[2])
    zlim[2]=fullzlim[2]
  if(zlim[1]>zlim[2])
    zlim[2] = zlim[1]
  # make tics
  if(is.null(leg)){
    ztic = seq(zlim[1],zlim[2],length.out = 1e5)
    ztict = zfun(ztic)
    ztticat = seq(zfun(zlim[1]),zfun(zlim[2]),length.out = ntic)
    leg = ztic[findInterval(ztticat,ztict)]
    leg[1] = zlim[1]
    leg[ntic] = zlim[2]

    # adjust tic step
    d = (10^floor(log10(leg[2]-leg[1])))/2
    leg = unique(d*round(leg / d))
    leg = leg[leg>=zlim[1]]
  }
  # make colors
  ztat = zfun(leg)
  ztcol = sort(unique(c(ztat,seq(zfun(zlim[1]),zfun(zlim[2]),length.out = N))))
  col = z2col(c(zfun(fullzlim),ztcol))[-(1:length(fullzlim))]
  at = match(ztat,ztcol)
  plotColorLegend(x0,x1,y0,y1,col,at=at,legend=leg,title=title,title.adj=title.adj,horizontal=horizontal)
}

#' Gradient legend
#'
#' @param x0,x1,y0,y1 rectange coordinates of the legend in nfc coordinates
#' @param col vector of colors (full gradient)
#' @param at indexes (of col vector) to place tics
#' @param legend tics labels
#' @param title legend title
#' @param title.adj legend title adj parameter to be passed to text function
#'
#' @export
#'
#' @examples
#' plot(1)
#' plotColorLegend(0.4,0.5,0.8,0.3,getPal(n=100),0:10*10,1:10)
plotColorLegend = function(x0,x1,y0,y1,col,at,legend,title=NULL,title.adj=c(0,-0.5),horizontal = FALSE){
  xpd = par(xpd=NA)
  if(horizontal){
    # to make marks at borders, boxes will be slightly larger
    coors = makeColorLenedCoordinates(x0,x1,length(col),at,grconvert=grconvertX)
    rect(coors$x[-length(coors$x)],grconvertY(y0,'nfc','user'),coors$x[-1],grconvertY(y0+(y1-y0)*0.25,'nfc','user'),col=col,border = NA)
  }else{
    coors = makeColorLenedCoordinates(y0,y1,length(col),at,grconvert=grconvertY)
    rect(grconvertX(x0,'nfc','user'),coors$x[-length(coors$x)],grconvertX(x0+(x1-x0)*0.25,'nfc','user'),coors$x[-1],col=col,border = NA)
    text(grconvertX(x0+(x1-x0)*0.3,'nfc','user'),coors$at,legend,adj=c(0,0.5))
    if(!is.null(title)){
      text(grconvertX(x0,'nfc','user'),coors$x[length(coors$x)],title,adj=title.adj)
    }
  }
  par(xpd=xpd)
}

makeColorLenedCoordinates = function(x0,x1,n,at,grconvert){
  w = (x1-x0)/(n - 1)/2
  x0 = x0 - w
  x1 = x1 + w
  x = seq(grconvert(x0,'nfc','user'),grconvert(x1,'nfc','user'),length.out = n+1)
  at = x[at]+(x[2]-x[1])/2
  list(x=x,at=at)
}

#' Merge multiple png files into pdf
#'
#' @param dir path to the folder with pngfiles (used if \code{fls} is NULL)
#' @param fls character vector of paths to png files
#' @param pdfout name of output pdf
#' @param ... other parameters for \code{\link{pdf}} function
#'
#' @export
mergePNG2PFD = function(dir=NULL,fls=NULL,pdfout,...){
  require(png)
  pdf(pdfout,...)

  if(is.null(fls))
    fls = sort(list.files(dir,full.names = T,pattern = "*.png"))

  for(f in fls) {
    print(f)
    pngRaster = readPNG(f)
    plot.new()
    par(xpd=NA)
    rasterImage(pngRaster, grconvertX(0,'ndc','user'), grconvertY(0,'ndc','user'),grconvertX(1,'ndc','user'),grconvertY(1,'ndc','user'))
  }
  dev.off()
}

#' Plots barplot with text around bar tops
#'
#' Barplot wrapper
#'
#' @param x heights of bars
#' @param t text to be plotted (\code{x} by default)
#' @param adj relative text position (see \code{\link{text}})
#' @param srt text rotation (see \code{\link{text}})
#' @param text.y text position (\code{x} by default)
#' @param ylim ylim
#' @param ... other barplot options
#'
#' @export
barplotWithText = function(x,t=x,adj=c(0.5,1.1),srt=0,text.y=x,ylim=c(0,max(x)),...){
  if(length(ylim)==1)
    ylim = c(0,max(x)*ylim)
  b=barplot(x,ylim=ylim,...)
  text(b,text.y,t,adj=adj,srt=srt)
}

#' Splits numetic vector into equally sized bins
#'
#' @param v vector to be splitted
#' @param n number of bins
#'
#' @return numeric vector win bin number (same length as \code{v})
#' @export
number2bin = function(v,n){
  o = order(v)
  j = 1
  for(i in 1:length(v)){
    if(j < i * n / length(v))
      j = j + 1
    v[o[i]] = j
  }
  v
}


#' Calculates row sums for matrix subsets
#'
#' @param d matrix
#' @param f character vector, factor to split matrix columns (length should be equal to nrow(d))
#' @param mean logical, calculate mean instead of sum
#'
#' @return matrix with number of rows equal to nrow(d) and number of columns equal to number of unique(f)
#' @export
calcColSums = function(d,f,mean=FALSE){
  # all levels are unique
  if(all(table(f)==1)){
    colnames(d) = f
    return(d)
  }
  # just one level
  if(length(unique(f))==1){
    r = Matrix::rowSums(d)
    r = matrix(r,ncol=1)
    rownames(r) = rownames(d)
    colnames(r) = f[1]
    return(r)
  }
  m = Matrix::sparse.model.matrix(~ f+0)
  r = d %*% m
  colnames(r) = sub('^f','',colnames(r))

  if(mean){
    t = as.numeric(table(f)[colnames(r)])
    r = sweep(r,2,t,'/')
  }
  if(!is(d,'sparseMatrix'))
    r = as.matrix(r)
  r
}

#' Add panel label
#'
#' Adds letter into top left corner of the plot
#'
#' @param l letter to be added
#' @param cex letter size
#' @param adj see adj parameter of \code{\link{text}}
#' @param ... other parameters for text function
#'
#' @export
plotPanelLetter = function(l,cex=1.2,adj=c(0,1.1),...){
  x=grconvertX(0,from='nfc',to='user')
  y=grconvertY(1,from='nfc',to='user')
  text(x=x,y=y,labels=l,adj=adj,font=2,cex=cex,xpd=NA)
}

#' Line with confidence interval
#'
#' CI is shown by area
#'
#' @param x x coordinates
#' @param p y coordinates. Matrix with two (mean and sd) or three (mean, lower CI bound, higher CI bound) columns
#' @param col line and area colour
#' @param sd.mult coefficient to multiply sd to get CI
#' @param new make ne plot (default is to add line on existed plot)
#' @param ylim,xlim see \code{\link{plot}}
#' @param area.transp alpha for CI area
#' @param type type of \code{\link{plot}}
#' @param area.den density of \code{\link{polygon}}
#' @param cilim numerical vector with two values, gives lower and apper values to truncate CI. NULL (to truncation) by default.
#' @param ...
#'
#' @export
plotArea = function(x,p,col,sd.mult=2,new=FALSE,ylim=NULL,xlim=range(x),area.transp=0.2,type='l',area.den=-1,cilim=NULL,...){
  #p should contain either mean and sd
  #or mean, lower and upper bounds
  o = order(x)
  x = x[o]
  p = p[o,,drop=F]
  na = apply(is.na(p),1,sum)==0
  x = x[na]
  p = p[na,,drop=F]
  if(ncol(p)==2)
    yp = c(p[,1]-p[,2]*sd.mult,rev(p[,1]+p[,2]*sd.mult))
  else
    yp = c(p[,2],rev(p[,3]))
  if(!is.null(cilim))
    yp = pmin(pmax(yp,cilim[1]),cilim[2])
  if(new){
    if(is.null(ylim))
      ylim = range(yp,na.rm=T)
    plot(1,t='n',xlim=xlim,ylim=ylim,...)
  }
  col.pol = col2rgb(col,alpha=TRUE)/255
  col.pol = rgb(col.pol[1],col.pol[2],col.pol[3],col.pol[4]*area.transp)
  polygon(c(x,rev(x)),yp,col=col.pol,border=NA,den=area.den)
  lines(x,p[,1],col=col,type=type,...)
}

#' Make Character Strings Unique
#'
#' Makes the elements of a character vector unique by appending sequence numbers to duplicates.
#' The difference from make.unique is that my.make.unique adds replicate number to all items, not to duplicates only.
#'
#' @param x a character vector
#' @param sep a character string used to separate a duplicate name from its sequence number.
#'
#' @return A character vector of same length as names with duplicates changed
#' @export
my.make.unique = function(x,sep='.'){
  ux = unique(x)
  for(u in ux){
    f = x == u
    x[f] = paste0(x[f],sep,1:sum(f))
  }
  x
}


#' Plots scatterplot with regression line
#'
#' Wrapper for pllot function
#'
#' @param x,y coordinates
#' @param cor.method pearson (default) or spearman
#' @param line.col colour of line
#' @param line.lwd line width
#' @param plot.ci logical, should CI be plotted
#' @param ci.transparency alpha of CI area
#' @param line.on.top logical, should line be on the top of points
#' @param new logical, create new plot (default) or add to existing
#' @param ... other arguments for plot function
#'
#' @export
plotLine = function(x,y=NULL,cor.method='pearson',line.col='red',line.lwd=1,plot.ci=FALSE,ci.transparency=0.3,line.on.top=TRUE,new=TRUE,...){
  if(is.null(y)){
    y = x[,2]
    x = x[,1]
  }
  if(new)
    plot(x,y,t='n',...)
  if(line.on.top)
    points(x,y,...)
  f = !is.na(x) & !is.na(y) & !is.infinite(x) &  !is.infinite(y)
  x = x[f]
  y = y[f]
  o = order(x)
  y = y[o]
  x = x[o]
  x = as.numeric(x)
  y = as.numeric(y)
  m = lm(y~x)
  ci=predict.lm(m,interval='confidence')
  if(plot.ci){
    c=col2rgb(line.col)
    polygon(c(x,rev(x)),c(ci[,2],rev(ci[,3])),border=NA,col=rgb(c[1],c[2],c[3],ci.transparency*255,maxColorValue = 255))
  }
  lines(x,ci[,1],col=line.col,lwd=line.lwd)
  if(!line.on.top)
    points(x,y,...)
  if(tolower(cor.method)==substr('spearman',1,nchar(cor.method))){
    x = rank(x)
    y = rank(y)
    cor.method = 'pearson'
  }
  c = cor.test(x,y,m=cor.method)
  ci = round(c$conf.int,2)
  leg=paste(toupper(substr(cor.method,1,1)),'CC=',round(c$estimate,2),' [',ci[1],',',ci[2],']\npv=',format(c$p.value,digits=2,scientific=TRUE),'\nN=',length(x),sep='')
  text(grconvertX(0.01,'npc','user'),grconvertY(0.99,'npc','user'),leg,col=line.col,adj=c(0,1),font=2)
}

#' Rename cluster by size
#'
#' Gives larger clusters smaller index
#'
#' @param x clustering, numeric vector with membership (as returned by cutree for instance)
#'
#' @return numeric cluster with new membership
#' @export
renameClustsBySize = function(x){
  t = table(x)
  n = names(t)
  o = order(t,decreasing=T)
  r = x
  for(i in 1:length(o))
    r[x == n[o[i]]] = i
  r
}

#' Scale numeric vector
#'
#' @param x numeric vector to be scaled
#' @param from,to target range
#' @param minx,maxx true min of "real x" (in case if given x is just subset of the real)
#' @param fraction fraction of target range to be used
#'
#' @return scaled value
#' @export
scaleTo = function(x,from=0,to=1,minx=min(x,na.rm=TRUE),maxx=max(x,na.rm=TRUE),fraction=1){
  x = (x-minx)/(maxx-minx)
  x*(to-from)*fraction + from + (to-from)*(1-fraction)/2
}

#' Capitalize first letter of input text
#'
#' @param t character vector
#'
#' @return t with first letter of each item capitalized
#' @export
first2Upper = function(t){
  paste0(toupper(substr(t,1,1)),substr(t,2,nchar(t)))
}
