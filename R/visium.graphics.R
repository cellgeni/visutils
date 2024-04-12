#' Plot Visium sample
#'
#' Main function to plot both histological and dimention reduction plots
#'
#' @param v Seurat object or data.frame with two columns (x and y coordinates) or list of polygons (list with x and y items) in case of type='tiles')
#' @param z z coordinate. Either numeric or categorical (factor/character)
#' @param cex relative size of symbols
#' @param type character, one of img, hex, rect, tiles, or xy, see Details.
#' @param border colour of symbol borders
#' @param z2col function to transform z to colous (if z is numeric) or one of viridis gradients (for example 'magma') or named vector of colours if z is character (names are levels of \code{z})
#' @param plot.legend logical
#' @param zlim numerical vector with two items, used to trim z
#' @param zfun function, transformation (log1p, sqrt, ect) of z (used only if z is numerical). Identity by default.
#' @param spot.filter logical vector, specifies which spots should be plotted
#' @param pch point type to be used with type='xy'
#' @param num.leg.tic numerical, desired number of tics in gradient legend
#' @param label.clusters label cluster on top of spots. Either logical, or character, if latter should provide names of each spot, it this case could be different from \code{z}
#' @param legend.args list of arguments to be passed to legend function (for categorical z), only title is used for numerical z.
#' @param randomize.points logical, should points be randomazed (makes sense if type is 'xy' and if spots overlap and ordered by z)
#' @param order.points.by.z logical, make sense to show values with high z on the top if type is 'xy'  (doublets for instance)
#' @param xaxt,yaxt type of axis (see \code{par}), no axis are drown by default
#' @param cluster.lab.adj adj parameterfor cluster labels (see  \code{text})
#' @param cluster.lab.cex size of cluster labels
#' @param cluster.lab.font font of cluster labels
#' @param cluster.lab2col named vector that maps cluster names to colors to be used as text colors. All labels are black if NULL (default). Lables not mentioned in here will be shown in black as well.
#' @param show.cluster.sizes logical, specifies whether number of cell per cluster should be shown in the legend
#' @param bg color to be used as background on xy/tiles plots
#' @param image.name name of image to plot, uses first if set to NULL. Use this parameter if v contains multiple visium objects.
#' @param pie.fracs if specified plots pies instead of simple cycles. Matrix with number of rows equal to the number of spots, and number of columns equal to pie pieces.
#' @param he.img.width integer, defines width (in pixels) the H&E figure should be resized to. Helps with speed of plotting and size of output files. Default is 400.
#' @param ... other arguments to be passed to graphical functions (see Details)
#'
#' @details Plots spots on top of H&E image, if type is 'img' (see \code{\link{plotVisiumImg}} for additional parameters),
#' plots hexagonal representation if type is 'hex' (see \code{\link{plotVisiumHex}}). Type 'xy' forces dimention redutcion plot (by \code{plot}),
#' coordinates either directly povided in \code{v} (as two-column data.frame) or extracted from umap slot of \code{v},
#' in the later case function attempts to use seurat_clusters if \code{z} is not specified.
#' Type 'tiles' plots each spot as polygon, for example from Voronoi diagram (see deldir::deldir). Polygons should be supplied as first argument
#'
#' @return data.frame with user spot coordinates
#' @export
plotVisium = function(v,z=NULL,cex=1,type='img',border=NA,z2col=num2col,plot.legend=TRUE,zlim=NULL,zfun = identity,spot.filter=NULL,pch=16,
                      num.leg.tic=NULL,label.clusters=FALSE,legend.args=list(),randomize.points=FALSE,order.points.by.z=FALSE,xaxt='n',yaxt='n',
                      cluster.lab.adj=c(0.5,0.5),cluster.lab.cex=1,cluster.lab.font=1,cluster.lab2col=NULL,show.cluster.sizes=FALSE,bg=NA,image.name=NULL,
                      pie.fracs=NULL,he.img.width=400,...){
  xy = NULL
  if('Seurat' %in% class(v)){
    if(type=='xy'){
      if(is.null(z))
        z = as.character(v$seurat_clusters)
      umap = names(v@reductions)
      umap = umap[grep('umap',umap,ignore.case = TRUE)]
      v = v@reductions[[umap]]@cell.embeddings
    }else if(length(v@images) > 0){
      image.name = getImageName(v,image.name,stopIfMoreThanOne=FALSE)
      xy = v@images[[image.name]]@coordinates
      if(length(v@images)>1){
        f = match(rownames(xy) , rownames(v@meta.data))
        N = ncol(v)
        if(length(label.clusters) > 1)
          label.clusters = label.clusters[f]
        z = recycle(z,N)[f]
        cex = recycle(cex,N)[f]
        border = recycle(border,N)[f]
        spot.filter = spot.filter[f]
        if(!is.null(pie.fracs)) pie.fracs = pie.fracs[f,]
      }
    }else{
      stop('img types can be used only with Seurat (Visium) object as input. For Dim plot set type="xy"')
    }
  }
  if(((is.data.frame(v) || is.array(v)) && type %in% c('rect','xy')) || (is.list(v) && type == 'tiles')){
    xy = v
  }
  if(is.null(xy)){
    stop('v should be either Seurat or data.frame or matrix list')
  }
  if(type %in% c('rect','xy') && !is.null(dim(xy)) && !all(colnames(xy) %in% c('x','y'))){
    colnames(xy) = c('x','y')
  }


  zorig = z
  if(all(is.na(z)))
    z = 'gray'
  # recycle
  if(!is.null(dim(xy))){
    N = nrow(xy)
  }else{
    N = length(xy) # for tiles
  }
  z = recycle(z,1:N)
  cex = recycle(cex,1:N)
  border= recycle(border,1:N)
  # subset spots
  if(!is.null(spot.filter)){
    if(length(label.clusters) == N)
      label.clusters = label.clusters[spot.filter]
    if(!is.null(dim(xy))){
      xy = xy[spot.filter,]
    }else{
      xy = xy[spot.filter] # for tiles
    }
    z = z[spot.filter]
    cex = cex[spot.filter]
    border = border[spot.filter]
    if(!is.null(pie.fracs)) pie.fracs = pie.fracs[spot.filter,]
    N = length(z)
  }
  # if z2col one of viridis gradientds
  if(length(z2col)==1 && is.character(z2col) && z2col %in% c('magma','inferno','plasma','cividis','rocket','mako','turbo','viridis')){
    require(viridis)
    ccc=viridis::viridis(100,option=z2col)
    z2col = function(x){
      num2col(x,ccc)
    }
  }
  # make color
  if(is.factor(z) | is.logical(z))
    z = as.character(z)
  if(is.numeric(z)){
    if(!is.null(zlim)){
      z[z>zlim[2]] = zlim[2]
      z[z<zlim[1]] = zlim[1]
    }else
      zlim=range(z,na.rm = TRUE)
    z = zfun(z)
    col = z2col(c(zfun(zlim),z))[-(1:length(zlim))]
  }else if(all(isColors(z))){
    col = z
  }else{
    if(!is.character(z2col)){
      z2col = char2col(z)
    }
    col = z2col[z]
  }
  # plot
  if(type=='img'){
    image = v@images[[image.name]]
    if('image' %in% slotNames(image))
      xy=plotVisiumImg(xy,image@image,image@scale.factors$lowres,image@spot.radius,cex=cex,col=col,border=border,xaxt=xaxt,yaxt=yaxt,pie.fracs=pie.fracs,he.img.width=he.img.width,...)
    else
      type = 'xy'
  }
  if(type=='hex'){
    xy=plotVisiumHex(xy,cex=cex,col=col,border=border,xaxt=xaxt,yaxt=yaxt,...)
  }
  if(type=='rect'){
    xy=plotVisiumRect(xy,cex=cex,col=col,border=border,xaxt=xaxt,yaxt=yaxt,...)
  }
  if(type=='xy'){
    if(randomize.points | order.points.by.z){
      if(randomize.points)
        o = sample(nrow(xy))
      else
        o = order(z)
      xy = xy[o,]
      col = col[o]
      cex = cex[o]
      z = z[o]
      if(length(label.clusters) == length(o))
        label.clusters = label.clusters[o]
    }
    plot(xy[,1:2],t='n',xaxt=xaxt,yaxt=yaxt,...)
    fillBackground(bg)
    points(xy[,1:2],cex=cex,col=col,pch=pch)
  }
  if(type == 'tiles'){
    plotTiles(col=col,tiles=xy,border=border,bg=bg,xaxt=xaxt,yaxt=yaxt,...)
  }
  #legend
  if(plot.legend){
    if(is.character(z2col)){
      # categorical
      if(is.null(legend.args$x)){
        legend.args$x = grconvertX(1,'npc','user')
        legend.args$y = grconvertY(1,'npc','user')
      }
      legend.args.def = list(xpd=NA,pch=19,col=z2col,legend=names(z2col),bty=par('bty'))
      for(n in names(legend.args.def))
        if(is.null(legend.args[[n]])) legend.args[[n]] = legend.args.def[[n]]
      if(show.cluster.sizes){
        clsize = table(col)[z2col]
        clsize[is.na(clsize)] = 0
        legend.args$legend = paste0(legend.args$legend,' (',clsize,')')
      }
      do.call(legend,legend.args)
    }else if(is.function(z2col) && !is.null(zorig) && is.numeric(zorig)){
      # numerical
      plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0.1,'npc','nfc'),grconvertY(0.9,'npc','nfc'),zlim,range(zorig,na.rm=TRUE),zfun,z2col,
                       leg=num.leg.tic,title=legend.args$title)
    }
  }
  if(is.character(label.clusters) || (label.clusters[1] & is.character(z))){
    if(is.character(label.clusters)){
      clusters = label.clusters
    }else
      clusters = z
    uclusters = unique(clusters)

    xl = sapply(split(xy[,1],clusters),mean,na.rm=T)[uclusters]
    yl = sapply(split(xy[,2],clusters),mean,na.rm=T)[uclusters]

    clcol = 'black'
    if(!is.null(cluster.lab2col)){
      clcol = cluster.lab2col[uclusters]
      clcol[is.na(clcol)] = 'black'
    }
    text(xl,yl,uclusters,adj = cluster.lab.adj,cex=cluster.lab.cex,font=cluster.lab.font,col=clcol)
  }
  invisible(list(xy=xy,z=z,col=col))
}

#' Plot Visium sample on top of H&E
#'
#' Normally this function shouldn't be called directly. Use plotVisium instead
#'
#' @param xy seu@images[[.]]@coordinates
#' @param img image (3D array)
#' @param scale.factor see @images[[.]]@scale.factors$lowres
#' @param spot.radius see @images[[.]]@spot.radius
#' @param cex spot size
#' @param col spot color
#' @param border color of spot border
#' @param spot.dist distance between spots (computed from xy by default)
#' @param img.alpha alpha level for H&E image
#' @param xlim,ylim,xlab,ylab parameters of \code{plot} function
#' @param symmetric.lims logical, specifies whether image crop should be square
#' @param pie.fracs if specified plots pies instead of simple cycles. Matrix with number of rows equal to the number of spots, and number of columns equal to pie pieces.
#' @param pie.cols colors to be used for pie pieces (ncol(pie.fracs) should be equal to length(pie.cols))
#' @param pie.min.frac all pieces with relative size less than \code{pie.min.frac} will be discared
#' @param he.img.width integer, defines width (in pixels) the H&E figure should be resized to. No resizing if NULL (default)
#' @param he.grayscale logical, specifies whether H&E image should be converted to grayscale
#' @param ... other parameters to be passed to \code{plot} function
#'
#' @return data.frame with user spot coordinates
#' @export
plotVisiumImg = function(xy,img,scale.factor,spot.radius,cex=1,col='red',border=NA,spot.dist=NULL,img.alpha=1,xlim=NULL,
                         ylim=NULL,symmetric.lims=TRUE,xlab='',ylab='',pie.fracs=NULL,pie.cols=NULL,pie.min.frac=0.05,
                         he.img.width=NULL,he.grayscale=FALSE,...){

  if(!is.null(he.img.width)){
    require(EBImage)
    coef = he.img.width/dim(img)[1]
    img = EBImage::resize(img,w=he.img.width)
    xy$imagerow = xy$imagerow*coef
    xy$imagecol = xy$imagecol*coef
  }
  if(he.grayscale){
    img = enhanceImage(img,wb = TRUE)
  }
  if(is.null(spot.dist)){
    #spot.dist = min(dist(xy[,c('imagerow','imagecol')]))*0.5
    spot.dist = spot.radius*max(dim(img))/scale.factor/2
  }
  xlim. = range(xy$imagecol*scale.factor)
  ylim. = range(nrow(img) - xy$imagerow*scale.factor)
  if(symmetric.lims){
    maxrange = max(xlim.[2]-xlim.[1],ylim.[2]-ylim.[1])
    xlim.[2] = xlim.[1] + maxrange
    ylim.[2] = ylim.[1] + maxrange
  }
  if(is.null(xlim)) xlim=xlim.
  if(is.null(ylim)) ylim=ylim.
  plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  rasterImage(1-(1-img)*img.alpha,1,1,ncol(img),nrow(img))

  f = cex>0
  if(any(cex>0) & is.null(pie.fracs))
    symbols(xy$imagecol[f]*scale.factor,
            nrow(img)-xy$imagerow[f]*scale.factor,
            circles=cex[f]*spot.dist*scale.factor,
            bg=col[f],inches = FALSE,fg=border[f],add = T)

  if(any(cex>0) & !is.null(pie.fracs)){
    pie.fracs = as.matrix(pie.fracs)
    pie.fracs = sweep(pie.fracs,1,apply(pie.fracs,1,sum),'/')
    pie.fracs[pie.fracs<pie.min.frac] = 0
    symbols.pie(x=xy$imagecol[f]*scale.factor,
                y=nrow(img)-xy$imagerow[f]*scale.factor,
                r=cex[f]*spot.dist*scale.factor,
                d=pie.fracs[f,],
                cols= pie.cols,border = border[f])
  }
  invisible(data.frame(x=xy$imagecol*scale.factor,y=nrow(img)-xy$imagerow*scale.factor))
}

#' Plot Visium as honeycombs
#'
#' Normally this function shouldn't be called directly. Use plotVisium instead
#'
#' @param xy seu@images[[.]]@coordinates
#' @param cex spot size
#' @param col spot color
#' @param border color of spot border (by default)
#' @param xlim,ylim,xlab,ylab parameters of \code{plot} function
#' @param hexstep size of spot hex, do not change unless you know what you are doing
#' @param transpose logical, whether to transpose the figure
#' @param ... other parameters to be passed to \code{plot} function
#'
#' @return data.frame with spot coordinates
#' @export
plotVisiumHex = function(xy,cex=1,col='red',border=col,xlab='Cols',ylab='Rows',xlim=c(min(xy$col)-1,max(xy$col)+1),
                         ylim=c(max(xy$row)+1,min(xy$row)-1),hexstep=0.25,transpose=FALSE,...){
  if(transpose){
    c = xy$row
    r = xy$col
    t = ylab
    ylab=xlab
    xlab=t
    xlim=c(max(c)+1,min(c)-1)
    ylim=c(min(r)-1,max(r)+1)
    rmask = c(-1,-1,0,1,1,0)
    cmask = c(-hexstep,hexstep,1-hexstep,hexstep,-hexstep,-1+hexstep)
  }else{
    c = xy$col
    r = xy$row
    cmask = c(-1,-1,0,1,1,0)
    rmask = c(-hexstep,hexstep,1-hexstep,hexstep,-hexstep,-1+hexstep)
  }
  cex = recycle(cex,1:length(r))
  col = recycle(col,1:length(r))

  plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)

  for(i in 1:length(r)){
    cexi = cex[i]
    if(cexi>0)
      polygon(c[i]+cmask*cexi,r[i]+rmask*cexi,col=col[(i-1) %% length(col)+1],border=border[(i-1) %% length(border)+1])
  }
  invisible(data.frame(x=c,y=r))
}

#' @export
plotTiles = function(col,tiles,border=NA,bg = NA,xlim=NULL,ylim=NULL,...){
  col = recycle(col,length(tiles))
  border = recycle(border,length(tiles))
  if(is.null(xlim))
    xlim = range(unlist(lapply(tiles,function(z)z$x)))
  if(is.null(ylim))
    ylim = range(unlist(lapply(tiles,function(z)z$y)))
  plot(1,t='n',xlim=xlim,ylim=ylim,...)
  fillBackground(bg)
  for(i in 1:length(tiles))
    polygon(tiles[[i]]$x,tiles[[i]]$y,col=col[i],border=border[i])
}



plotVisiumRect = function(xy,cex=1,col='red',border=NA,xlab='x',ylab='y',
                          xlim=c(min(xy$x)-0.5,max(xy$x)+0.5),
                          ylim=c(min(xy$y)+0.5,max(xy$y)-0.5),...){
  cex = recycle(cex,1:nrow(xy))
  col = recycle(col,1:nrow(xy))
  border = recycle(border,1:nrow(xy))

  plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)

  for(i in 1:nrow(xy)){
    cexi = cex[i]*0.5
    if(cexi>0)
      rect(xy$x[i]-cexi,xy$y[i]-cexi,xy$x[i]+cexi,xy$y[i]+cexi,col=col[i],border=border[i])
  }
  invisible(xy)
}

#' Plot multiple numerical values (gene expression or cell abundances) on H&E image
#'
#' Each spot is colored by weighted mean colour. Opacity of spot is proportional to maximal feature intensity.
#'
#' @param v Seurat object
#' @param z matrix with values (in columns) to be plotted
#' @param cols colors to be used for columns in z
#' @param zfun function to transform values in z (z^2 is default)
#' @param scale.per.colour logical, specifies whether each color should cover whole range (that is, should z be scaled per column)
#' @param min.opacity minimal spot opacity. Default is 0, that means that spots with low intensity of all features will be almost transparent. Set it higer if you want at least one feature to be visible in each spot.
#' @param title.adj legend title adj (to be passed to text function)
#' @param legend.ncol number of legend columns. Set to 0 to suppress legend plotting.
#' @param ... other parameters to be passed to plotVisium
#'
#' @return data.frame with user spot coordinates
#' @export
plotVisiumMultyColours = function(v,z,cols=NULL,zfun=function(x)x^2,scale.per.colour=TRUE,
                                  min.opacity=0,title.adj=c(0,-0.5),
                                  legend.ncol=1,...){
  z = as.matrix(z)
  zscaled = zfun(z)
  if(scale.per.colour){
    for(i in 1:ncol(zscaled)) zscaled[,i] = scaleTo(zscaled[,i])
  }else{
    zscaled = scaleTo(zscaled)
  }

  zscaled[is.na(zscaled)] = 0
  cols = col2hex(cols,withAlpha = FALSE)

  # opacity is proportional to untransformed z, maybe I'll need to add opacity transformation function, lets see
  opacity = z
  if(scale.per.colour){
    for(i in 1:ncol(opacity)) opacity[,i] = scaleTo(opacity[,i])
  }else{
    opacity = scaleTo(opacity)
  }
  # opacity is proportional to max opacity, maybe sum can be used instead.
  opacity = scaleTo(apply(opacity,1,max,na.rm=TRUE),from=min.opacity,to=255)
  opacity[is.na(opacity)] = 0


  col = weightedColourMeans(cols,zscaled)
  na = is.na(col)
  col[na] = '#000000' # just to make col2rgb below works, anyway all NA spots (ones with zero intensity in all chanels) will be set transparent
  col = rbind(col2rgb(col),opacity)
  col = apply(col,2,function(x)rgb(x[1],x[2],x[3],x[4],maxColorValue = 255))
  # to overwrite transparency for empty spots
  col[na] = '#00000000'

  res = plotVisium(v,col,...)


  # plot legends
  if(legend.ncol>0){
    legend.nrow = ceiling(length(cols) / legend.ncol)

    x0 = grconvertX(1,'npc','nfc')
    y0 = grconvertY(c(0,1),'npc','nfc')

    lw = grconvertY(1:2,'line','nfc')
    lw = max(lw)-min(lw)

    dx = (1-x0)/legend.ncol

    dy = (y0[2]-y0[1])/legend.nrow

    y0 = y0[2]

    for(i in 1:length(cols)){
      col1 = paste0(cols[i],'00')
      col2 = paste0(cols[i],'FF')
      r = (i-1) %% legend.nrow
      c = (i-1) %/% legend.nrow
      plotColorLegend2(x0+dx*c,
                       x0+dx*(c+1),
                       y0-dy*(r+1)+lw*0.8,
                       y0-dy*r-lw*0.8,
                       zlim = range(z[,i]),fullzlim = range(z[,i]),zfun=zfun,
                       z2col=function(x)num2col(x,c(col1,col2)),title=colnames(z)[i],title.adj = title.adj)
    }
  }
  invisible(res)
}

#' Plots scatter pie
#'
#' Plots set of pie charts in specified positions
#'
#' @param x,y positions of pies to plot
#' @param r radiuses of pies (recycled along x)
#' @param d matrix of pie piece sizes. Row sums are normalized to 1. Nrow(d) should be equal to lenght(x).
#' @param cols colours of pieces. Ncol(d) should be equal to length(cols)
#' @param border colour of pie borders, NA (default) for no borders.
#' @param alpha0 start angle
#'
#' @export symbols.pie
#'
#' @examples
#' plot(1,t='n',xlim=c(0,5),ylim=c(0,10))
#' d = matrix(runif(6*4),ncol=4)
#' symbols.pie(c(1,2,3,1,2,3),c(1,1,1,1,1,1),1/2,d,1:4)
symbols.pie = function(x,y,r,d,cols,border=NA,alpha0=0){
  r = recycle(r,1:length(x))
  d = sweep(d,1,apply(d,1,sum),'/')*2*pi
  for(i in 1:nrow(d)){
    if(any(is.na(d[i])) || any(is.nan(d[i])) || any(is.infinite(d[i]))) next
    a = alpha0
    for(j in 1:ncol(d)){
      if(d[i,j]>0){
        as = c(seq(a,a+d[i,j],by=pi/18),a+d[i,j])
        polygon(c(x[i],x[i]+r[i]*cos(as),x[i]),c(y[i],y[i]+r[i]*sin(as),y[i]),col = cols[j],border = NA)
        a = a+d[i,j]
      }
    }
    if(!is.na(border[i])){
      as = seq(0,2*pi,length.out = 36)
      polygon(x[i]+r[i]*cos(as),y[i]+r[i]*sin(as),col = NA,border = border[i])
    }
  }
}

fillBackground = function(bg,...){
  if(!is.na(bg)){
    rect(grconvertX(0,'nfc','user'),grconvertY(0,'nfc','user'),
         grconvertX(1,'nfc','user'),grconvertY(1,'nfc','user'),col = bg,border = NA,...)
  }
}
