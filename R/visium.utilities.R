#' Crop image in Seuart Visium object
#'
#' Retain square image that contains all spots (listed in the object). Spot coordinates adjasted accordingly.
#'
#' @param v Seurat object
#'
#' @return
#' @export
cropVisiumImage = function(v){
  c = v@images$slice1@coordinates
  scalefactors = v@images$slice1@scale.factors
  img = v@images$slice1@image

  c[,4:5] = c[,4:5]*scalefactors$lowres
  r = min(dist(c[,4:5]))*0.5 + 1

  rr = range(c$imagerow)
  rc = range(c$imagecol)

  cc = round(mean(rc))
  cr = round(mean(rr))

  ir = ceiling(max(rr[2]/2 - rr[1]/2,rc[2]/2 - rc[1]/2)) + ceiling(r)
  #irr = ceiling(rr[2]/2 - rr[1]/2 + r)
  #irc = ceiling(rc[2]/2 - rc[1]/2 + r)
  irr = irc = ir

  img = img[(cr-irr):(cr+irr),(cc-irc):(cc+irc),]
  c$imagerow = c$imagerow - (cr-irr-1)
  c$imagecol = c$imagecol - (cc-irc-1)
  c[,4:5] = c[,4:5]/scalefactors$lowres
  # c[,4] = round(c[,4])
  # c[,5] = round(c[,5])
  v@images$slice1@coordinates = c
  v@images$slice1@image = img
  v
}

#' Plots scatter pie
#'
#' Plots set of pie charts in specified positions
#'
#' @param x,y positions of pies to plot
#' @param r radiuses of pies (recycled along x)
#' @param d matrix of pie piece sizes. Row sums are normolized to 1. Nrow(d) should be equal to lenght(x).
#' @param cols colours of pieces. Ncol(d) should be equal to length(cols)
#' @param border colour of pie borders, NA (default) for no borders.
#' @param alpha0 start angle
#'
#' @return
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

#' Load spaceranger results
#'
#' Wrapper for \code{\link[Seurat]{Load10X_Spatial}}
#'
#' @param data.dir path to the folder with h5 files
#' @param filter.matrix logical, specifies whether only tissue spots should be loaded
#' @param ens_id logical, specifies whetehr Ensambl IDs should be used to identify genes (instead of gene names)
#' @param ... other parameters to \code{\link[Seurat]{Load10X_Spatial}}
#'
#' @return Seurat object
#' @export
myLoad10X_Spatial = function(data.dir,filter.matrix=TRUE,ens_id=TRUE,...){
  #TODO: rename)
  #TODO: if(dir.exists(paste0(path,f,'/outs')))
  d = Load10X_Spatial(data.dir,ifelse(filter.matrix,'filtered_feature_bc_matrix.h5','raw_feature_bc_matrix.h5'),filter.matrix=filter.matrix,...)
  if(ens_id){
    gids = read.table(paste0(data.dir,ifelse(filter.matrix,'/filtered_feature_bc_matrix/features.tsv.gz','/raw_feature_bc_matrix/features.tsv.gz')))
    d@assays$Spatial@counts@Dimnames[[1]] = gids$V1
    d@assays$Spatial@data@Dimnames[[1]] = gids$V1
    rownames(d@assays$Spatial@meta.features) = gids$V1
    d[['Spatial']][['name']] = gids$V2
    d[['Spatial']][['ensid']] = gids$V1
    #  d@assays$Spatial@scale.data@Dimnames[[1]] = gids$V1
  }
  if(!filter.matrix){
    f = read.csv(paste0(data.dir,'/spatial/tissue_positions_list.csv'),header = F,row.names = 1)
    d[['is.tissue']] = f[colnames(d),1]
  }
  d
}

#' Loads scanpy h5ad file with Visium data into Seurat object
#'
#' @param filename character, path to h5ad file
#'
#' @return Seurat object
#' @export
myLoadH5AD_Spatial = function (filename){
  require(Matrix)
  require(rhdf5)
  a = rhdf5::H5Fopen(filename)
  ll = sapply(a$obs,length)
  obs = as.data.frame(a$obs[ll==max(ll)],check.names=F)
  rownames(obs) = a$obs[["_index"]]

  var = as.data.frame(a$var[-1:-2],check.names=F)
  rownames(var) = a$var[["_index"]]

  m = a$X
  mtx = sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = c(nrow(var),nrow(obs)))
  rownames(mtx) = rownames(var)
  colnames(mtx) = rownames(obs)
  object <- CreateSeuratObject(counts = mtx, assay = 'Spatial')

  image <- aperm(a$uns[[1]][[1]]$images$hires,3:1)

  scale.factors <- a$uns[[1]][[1]]$scalefactors
  # looks like in stomics both images are hires actually (at least they were identical for example I tried)
  scale.factors$tissue_lowres_scalef = scale.factors$tissue_hires_scalef
  tissue.positions = cbind(obs[,c('in_tissue','array_row','array_col')], t(a$obsm$spatial)[,2:1])
  colnames(tissue.positions) = c("tissue", "row", "col", "imagerow", "imagecol")
  rownames(tissue.positions) = rownames(obs)

  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius/max(dim(x = image))
  image = new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                              fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef,
                                                                              scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, spot.radius = spot.radius)


  image <- image[Cells(x = object)]
  DefaultAssay(object = image) = 'Spatial'
  object[['slice1']] = image
  rhdf5::H5Fclose(a)
  object@meta.data = cbind(object@meta.data,obs)
  object@assays$Spatial@meta.features = var
  return(object)
}



#' Loads scanpy h5ad file with multiple Visium data into list of Seurat object
#'
#' @param filename character, path to h5ad file
#' @param library_id_field name of adata.obs that specifies visium library. 'library_id' by default
#'
#' @return list of Seurat object
#' @export
myLoadH5AD_Spatials = function (filename,library_id_field='library_id'){
  require(Matrix)
  require(rhdf5)
  a = rhdf5::H5Fopen(filename)
  ll = sapply(a$obs,length)
  obs = as.data.frame(a$obs[ll==max(ll)],check.names=F)
  rownames(obs) = a$obs[["_index"]]
  for(fn in names(a$obs[['__categories']])){
    obs[[fn]] = a$obs[['__categories']][[fn]][obs[[fn]]+1]
  }


  ll = sapply(a$var,length)
  var = as.data.frame(a$var[ll==max(ll)],check.names=F)

  for(fn in names(a$var[['__categories']])){
    var[[fn]] = a$var[['__categories']][[fn]][var[[fn]]+1]
  }
  if(!is.null(var$gene_ids))
    rownames(var) = var$gene_ids
  else if (!is.null(var$SYMBOL))
    rownames(var) = var$SYMBOL

  var = as.data.frame(a$var[-1:-2],check.names=F)
  rownames(var) = a$var[["_index"]]

  m = a$X
  mtx = sparseMatrix(i=m$indices+1, p=m$indptr,x = as.numeric(m$data),dims = c(nrow(var),nrow(obs)))
  rownames(mtx) = rownames(var)
  colnames(mtx) = rownames(obs)

  res = list()
  for(lid in unique(obs[[library_id_field]])){
    f = obs[[library_id_field]] == lid
    obs_ = obs[f,]
    if(!is.null(obs_$barcode)){
      rownames(obs_) = obs_$barcode
    }
    mtx_ = mtx[,which(f)]
    colnames(mtx_) = rownames(obs_)

    object <- CreateSeuratObject(counts = mtx_, assay = 'Spatial')

    image <- aperm(a$uns$spatial[[lid]]$images$hires,3:1)

    scale.factors <- a$uns$spatial[[lid]]$scalefactors
    # looks like in stomics both images are hires actually (at least they were identical for example I tried)
    scale.factors$tissue_lowres_scalef = scale.factors$tissue_hires_scalef
    tissue.positions = cbind(obs[f,c('in_tissue','array_row','array_col')], t(a$obsm$spatial[,f])[,2:1])
    colnames(tissue.positions) = c("tissue", "row", "col", "imagerow", "imagecol")



    rownames(tissue.positions) = rownames(obs_)

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    spot.radius <- unnormalized.radius/max(dim(x = image))
    image = new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                                                                fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef,
                                                                                scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, spot.radius = spot.radius)


    image <- image[Cells(x = object)]
    DefaultAssay(object = image) = 'Spatial'
    object[['slice1']] = image

    object@meta.data = cbind(object@meta.data,obs_)
    rownames(object@meta.data) = rownames(obs_)

    object@assays$Spatial@meta.features = var
    res[[lid]] = object
  }

  rhdf5::H5Fclose(a)
  return(res)

}

#' Adjast image brightnes and contrast
#'
#' @param p image (3d numeric array)
#' @param wb logical, specifies whether output image should be transformed to grayscale
#' @param pow power of transformation
#' @param qs quantiles to trim. Numerical vector with two items. Trims all values outside of specified quantile range.
#' @param trim01 logical, wspecifies whether pixels with zero (black) and maximal (white) intensity should be trimed ahead of quantile trimming.
#'
#' @return image (3d numeric array)
#' @export
enhanceImage = function(p,wb=FALSE,qs=NULL,trim01 = TRUE){
  pm = apply(p,1:2,max)
  if(!is.null(qs)){
    f = pm == 0 | pm == 1
    if(trim01)
      qs = quantile(pm[!f],probs = qs)
    else
      qs = quantile(pm,probs = qs)
    pm = (pm-qs[1])/(qs[2]-qs[1])
    pm[pm>1] = 1
    pm[pm<0] = 0
  }
  f = apply(p,1:2,max)/pm
  f[is.na(f)] = 1
  p = sweep(p,1:2, f,'/')
  if(wb){
    z=apply(p,1:2,mean)
    for(j in 1:3)
      p[,,j] = z
  }
  p
}


#' Adjast image brightnes and contrast
#'
#' @param p image (3d numeric array)
#' @param wb logical, specifies whether output image should be transformed to grayscale
#' @param pow power of transformation
#' @param qs quantiles to trim. Numerical vector with two items. Trims all values outside of specified quantile range.
#'
#' @details it is old versioin of enhanceImage
#'
#' @return image (3d numeric array)
#' @export
enhanceImage_ = function(p,wb=FALSE,pow=1,qs=NULL){
  p = (1-p)^pow
  pm = apply(p,1:2,max)
  if(!is.null(qs)){
    qs = quantile(pm,probs = qs)
    pm = (pm-qs[1])/(qs[2]-qs[1])
    pm[pm>1] = 1
    pm[pm<0] = 0
  }
  f = apply(p,1:2,max)/pm
  f[is.na(f)] = 1
  p = sweep(p,1:2, f,'/')
  p = 1-p
  if(wb){
    z=apply(p,1:2,mean)
    for(j in 1:3)
      p[,,j] = z
  }
  p
}

#' Plot multiple numerical values (gene expression or cell abundancies) on H&E image
#'
#' Values are plotted as a sum of color gradients. Be carefull, result depends on the order of features.
#'
#' @param v Seurat object
#' @param x matrix with values (in columns) to be plotted
#' @param cols colour to be used for x columns
#' @param log.pc numeric, pseudocount to be added before log-transformation. No transformation is applied if NA.
#' @param scale.per.colour logical, specifies whether each colour should cover whole range (that is, should x be sclaed per column)
#' @param reorderByOpacity logical, specifes whether colours should be ordered by increasing opacity prior to summing
#' @param title.adj legend title adj (to be passed to text function)
#' @param bg color to use as spot background. NULL (default) for transparent background.
#' @param legend.ncol number of legend columns. Set to 0 to suppress legend plotting.
#' @param ... other parameters to be passed to plotVisium
#'
#' @return
#' @export
#'
#' @examples
plotVisiumMultyColours = function(v,x,cols=NULL,log.pc=NA,scale.per.colour=TRUE,reorderByOpacity=FALSE,title.adj=c(0,-0.5),bg='#FFFFFFFF',legend.ncol=1,...){
  xs = x
  if(!is.na(log.pc))
    xs = log(x+log.pc)
  if(scale.per.colour){
    for(i in 1:ncol(xs)) xs[,i] = scaleTo(xs[,i])
  }else{
    xs = scaleTo(xs)
  }
  xs[is.na(xs)] = 0
  cols = col2hex(cols,withAlpha = FALSE)
  col = sapply(1:ncol(xs),function(i)num2col(xs[,i],paste0(cols[i],c('00','FF')),minx = 0,maxx = 1))
  col = overlayColours(col,reorderByOpacity = reorderByOpacity)
  if(!is.null(bg))
    col = overlayColours(cbind(bg,col),reorderByOpacity = FALSE)

  r = plotVisium(v,col,...)


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
                       zlim = range(x[,i]),fullzlim = range(x[,i]),zfun=ifelse(log,base::log,identity),
                       z2col=function(x)num2col(x,c(col1,col2)),title=colnames(x)[i],title.adj = title.adj)
    }
  }
  invisible(r)
}


#' Plot Visium sample
#'
#' Main function to plot both histological and dimention reduction plots
#'
#' @param v Seurat object or data.frame with two columns (x and y coordinates)
#' @param z z coordinate. Either numeric or categorical (factor/character)
#' @param cex relative size of symbols
#' @param type character, one of img, hex, or xy, see Details.
#' @param border colour of symbol borders
#' @param z2col function to transform z to colous (if z is numeric) or named vector of colours if z is character (names are levels of \code{z})
#' @param plot.legend logical
#' @param zlim numerical vector with two items, used to trim z
#' @param zfun function, transformation (log1p, sqrt, ect) of z (used only if z is numerical). Identity by default.
#' @param spot.filter logical vector, specifies which spots should be plotted
#' @param num.leg.tic numerical, desired number of tics in gradient legend
#' @param label.clusters label cluster on top of spots. Either logical, or character, if latter should provide names of each spot, it this case could be different from \code{z}
#' @param legend.args list of arguments to be passed to legend function (for categorical z), only title is used for numerical z.
#' @param randomize.points logical, should points be randomazed (makes sense if type is 'xy' and if spots overlap and ordered by z)
#' @param order.points.by.z logical, make sense to show values with high z on the top if type is 'xy'  (doublets for instance)
#' @param xaxt,yaxt type of axis (see \code{par}), no axis are drown by default
#' @param cluster.lab.adj adj parameterfor cluster labels (see  \code{text})
#' @param cluster.lab.cex size of cluster labels
#' @param ... other arguments to be passed to graphical functions (see Details)
#'
#' @details Plots spots on top of H&E image, if type is 'img' (see \code{\link{plotVisiumImg}} for additional parameters),
#' plots hexagonal representation if type is 'hex' (see \code{\link{plotVisiumHex}}). Type 'xy' forces dimention redutcion plot (by \code{plot}),
#' coordinates either directly rpovided in \code{v} (as two-column data.frame) or extracted from umap slot of \code{v},
#' in the later case function attempts to use seurat_clusters if \code{z} is not specified
#'
#' @return
#' @export
plotVisium = function(v,z=NA,cex=1,type='img',border=NA,z2col=num2col,plot.legend=TRUE,zlim=NULL,zfun = identity,spot.filter=NULL,
                      num.leg.tic=NULL,label.clusters=FALSE,legend.args=list(),randomize.points=FALSE,order.points.by.z=FALSE,xaxt='n',yaxt='n',cluster.lab.adj=c(0.5,0.5),cluster.lab.cex=1,...){
  if('Seurat' %in% class(v) & type=='xy'){
    if(is.na(z[1]))
      z = as.character(v$seurat_clusters)
    v = v@reductions$umap@cell.embeddings
  }
  if('Seurat' %in% class(v) && length(v@images) > 0)
    xy = v@images$slice1@coordinates
  else{
    if(type %in% c('img')) stop('img types can be used only with Seurat (Visium) object as input. For Dim plot set type="xy"')
    xy = v
    if(!all(colnames(xy) %in% c('x','y')) & type != 'hex')
      colnames(xy) = c('x','y')
  }
  if(all(is.na(z)))
    z = 'gray'
  # recycle
  z = recycle(z,1:nrow(xy))
  cex = recycle(cex,1:nrow(xy))
  border= recycle(border,1:nrow(xy))
  # subset spots
  if(!is.null(spot.filter)){
    xy = xy[spot.filter,]
    z = z[spot.filter]
    cex = cex[spot.filter]
    border = border[spot.filter]
    if(length(label.clusters) == length(spot.filter))
      label.clusters = label.clusters[spot.filter]
  }
  # make color
  if(is.factor(z))
    z = as.character(z)
  if(is.numeric(z)){
    if(!is.null(zlim)){
      z[z>zlim[2]] = zlim[2]
      z[z<zlim[1]] = zlim[1]
    }else
      zlim=range(z,na.rm = TRUE)
    zorig = z
    z = zfun(z)
    col = z2col(c(zfun(zlim),z))[-(1:length(zlim))]
  }else if(all(isColors(z))){
    plot.legend = FALSE
    col = z
  }else{
    if(!is.character(z2col)){
      z2col = char2col(z)
    }
    col = z2col[z]
  }
  # plot
  if(type=='img'){
    xy=plotVisiumImg(xy,v@images$slice1@image,v@images$slice1@scale.factors$lowres,cex=cex,col=col,border=border,xaxt=xaxt,yaxt=yaxt,...)
  }
  if(type=='hex'){
    xy=plotVisiumHex(xy,cex=cex,col=col,border=border,xaxt=xaxt,yaxt=yaxt,...)
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
    plot(xy[,1:2],cex=cex,col=col,xaxt=xaxt,yaxt=yaxt,...)
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
      do.call(legend,legend.args)
    }else{
      # numerical
      plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0.1,'npc','nfc'),grconvertY(0.9,'npc','nfc'),zlim,range(zorig,na.rm=TRUE),zfun,z2col,leg=num.leg.tic,title=legend.args$title)
    }
  }
  if(is.character(label.clusters) || (label.clusters[1] & is.character(z))){
    if(is.character(label.clusters)){
      clusters = label.clusters
      if(!is.null(spot.filter)){
        clusters = clusters[spot.filter]
      }
    }else
      clusters = z
    uclusters = unique(clusters)

    xl = sapply(split(xy[,1],clusters),mean,na.rm=T)[uclusters]
    yl = sapply(split(xy[,2],clusters),mean,na.rm=T)[uclusters]
    text(xl,yl,uclusters,adj = cluster.lab.adj,cex=cluster.lab.cex)
  }
  #cat(nrow(xy),',',ncol(xy),', ',length(z),',',length(col),'\n')
  invisible(data.frame(x=xy[,1],y=xy[,2],z=z,col=col))
}

#' Plot Visium sample on top of H&E
#'
#' Normally this function shouldn't be called directly. Use plotVisium instead
#'
#' @param xy seu@images$slice1@coordinates
#' @param img image (3D array)
#' @param scale.factor seu@images$slice1@scale.factors$lowres
#' @param cex spot size
#' @param col spot color
#' @param border color of spot border
#' @param spot.dist distance between spots (computed from xy by default)
#' @param img.alpha alpha level for H&E image
#' @param xlim,ylim,xlab,ylab parameters of \code{plot} function
#' @param symmetric.lims logical, specifies whether image crop should be square
#' @param pie.fraqs if specified plots pies instead of simple cycles. Matrix with number of rows equal to the number of spots, and number of columns equal to pie pieces.
#' @param pie.cols colors to be used for pie pieces (ncol(pie.fraqs) should be equal to length(pie.cols))
#' @param pie.min.fraq all pieces with relative size less than \code{pie.min.fraq} will be discared
#' @param he.img.width integer, defines width (in pixels) the H&E figure should be resized to. No resizing if NULL (default)
#' @param he.grayscale logical, specifies whether H&E image should be converted to grayscale
#' @param ... other parameters to be passed to \code{plot} function
#'
#' @return
#' @export
plotVisiumImg = function(xy,img,scale.factor,cex=1,col='red',border=NA,spot.dist=NULL,img.alpha=1,xlim=NULL,
                         ylim=NULL,symmetric.lims=TRUE,xlab='',ylab='',pie.fraqs=NULL,pie.cols=NULL,pie.min.fraq=0.05,
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
    spot.dist = min(dist(xy[,4:5]))*0.5
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
  if(any(cex>0) & is.null(pie.fraqs))
    symbols(xy$imagecol[f]*scale.factor,
            nrow(img)-xy$imagerow[f]*scale.factor,
            circles=cex[f]*spot.dist*scale.factor,
            bg=col[f],inches = FALSE,fg=border[f],add = T)

  if(any(cex>0) & !is.null(pie.fraqs)){
    pie.fraqs = sweep(pie.fraqs,1,apply(pie.fraqs,1,sum),'/')
    pie.fraqs[pie.fraqs<pie.min.fraq] = 0
    symbols.pie(x=xy$imagecol[f]*scale.factor,
                y=nrow(img)-xy$imagerow[f]*scale.factor,
                r=cex[f]*spot.dist*scale.factor,
                d=pie.fraqs[f,],
                cols= pie.cols,border = border[f])
  }
  invisible(data.frame(y=xy$imagecol*scale.factor,y=nrow(img)-xy$imagerow*scale.factor))
}

plotVisiumHex = function(xy,cex=1,col='red',border=NA,xlab='Cols',ylab='Rows',xlim=c(min(xy$col)-1,max(xy$col)+1),ylim=c(max(xy$row)+1,min(xy$row)-1),hexstep=0.25,transpose=TRUE,...){
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

#' Download Visium dataset from 10x website
#'
#' @param url url of the dataset (whithot dataset name, see example)
#' @param sample.name name of the dataset
#' @param outdir folder to save dataset (function will create subfolder named by \code{sample.name})
#'
#' @return
#' @export
#'
#' @examples
#' loadVisiumFrom10x("https://cf.10xgenomics.com/samples/spatial-exp/1.3.0","Visium_Mouse_Olfactory_Bulb",outdir="data")
loadVisiumFrom10x = function(url,sample.name,outdir){
  print(sample.name)
  url = paste0(url,'/',sample.name,'/',sample.name,'_')
  cmd = paste0(
    "rm -rf ",outdir,"/",sample.name,";
        mkdir ",outdir,"/",sample.name,";
        cd ",outdir,"/",sample.name,";
        curl -O ",url,"filtered_feature_bc_matrix.h5;
        curl -O ",url,"filtered_feature_bc_matrix.tar.gz;
        curl -O ",url,"raw_feature_bc_matrix.h5;
        curl -O ",url,"raw_feature_bc_matrix.tar.gz;
        curl -O ",url,"spatial.tar.gz;
        curl -O ",url,"analysis.tar.gz
        curl -O ",url,"web_summary.html;
        tar -xzvf ",sample.name,"_filtered_feature_bc_matrix.tar.gz;
        tar -xzvf ",sample.name,"_raw_feature_bc_matrix.tar.gz;
        tar -xzvf ",sample.name,"_spatial.tar.gz;
        tar -xzvf ",sample.name,"_analysis.tar.gz;
        mv ",sample.name,"_filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5;
        mv ",sample.name,"_raw_feature_bc_matrix.h5 raw_feature_bc_matrix.h5;
        rm ",sample.name,"_filtered_feature_bc_matrix.tar.gz;
        rm ",sample.name,"_raw_feature_bc_matrix.tar.gz;
        rm ",sample.name,"_spatial.tar.gz;
        rm ",sample.name,"_analysis.tar.gz;")
  #print(cmd)
  system(cmd)
}

# swaps and rotations

#' Extracts mean color of given spot
#'
#'
#' @param i image (3d array)
#' @param x,y coordinates of center of the spot
#' @param r radius of the spot
#'
#' @return numeric vector with three values (red, green, and blue)
#' @export
getMeanSpotColor1 = function(i,x,y,r){
  xr = round(x)
  yr = round(y)
  rn = ceiling(r)
  res = c(r=0,g=0,b=0)
  n = 0
  r2 = r^2
  for(k in (xr-rn):(xr+rn))
    for(m in (yr-rn):(yr+rn)){
      if((k-x)^2+(m-y)^2 <= r){
        n = n+1
        res = res + i[k,m,]
      }
    }
  res/n
}

#' Extracts mean color of spots
#'
#' @param vis Seurat object
#' @param scalefactors list of scale factors (see Details)
#'
#' @details Seurat version 4.1.0 for stores wrong value for spot size in .
#' So one needs to load scale factors themself using jsonlite::read_json('/spatial/scalefactors_json.json')).
#' Output of read_json can be used as scalefactors parameter.
#'
#' @return matrix with three coumns (red, green, blue) and number of rows equal to the number of spots in vis object
#' @export
getMeanSpotColor = function(vis,scalefactors){
  i = vis@images$slice1@image
  c = vis@images$slice1@coordinates
  r = scalefactors$spot_diameter_fullres/2*scalefactors$tissue_lowres_scalef
  c = c[,4:5]*scalefactors$tissue_lowres_scalef
  r = t(sapply(1:nrow(c),function(j)getMeanSpotColor1(i,c$imagerow[j],c$imagecol[j],r)))
  rownames(r) = colnames(vis)
  r
}



#' Title
#'
#' @param m
#'
#' @return
#' @export
getMaxInxByFirstDD = function(m){
  r = do.call(rbind,lapply(1:dim(m)[1],function(i){
    bout = which(m[i,,]==max(m[i,,]),arr.ind = T)[1,]
    bin  = which.max(m[i,i,])

    # if diag isnot worse then change to diag
    if(m[i,i,bin] == max(m[i,bout[1],bout[2]])){
      bout = c(i,bin)
    }
    data.frame(inx=i,
               best.match.inx=bout[1],
               best.match.tr = dimnames(m)[[3]][bout[2]],
               is.diag = i == bout[1],
               best.diag.tr = dimnames(m)[[3]][bin],
               ndiag.max=max(m[i,-i,]),
               diag.max =max(m[i, i,]))
  }))
  rownames(r) = dimnames(m)[[1]]
  r
}

#' Return list of 8 possible mirroring/rotations function
#'
#' @return list of functions
#' @export
getRotations = function(){
  t_ = function(i) aperm(i,c(2,1,3))
  c_ = function(i) i[,dim(i)[2]:1,]
  r_ = function(i) i[dim(i)[1]:1,,]

  transforms = list(o0 = list(identity),
                    o1 = list(r_,t_),
                    o2 = list(r_,c_),
                    o3 = list(t_,r_),
                    m0 = list(r_),
                    m1 = list(t_),
                    m2 = list(c_),
                    m3 = list(r_,c_,t_))
}

#' Apply list of transformations to image
#'
#' @param i image (3D array)
#' @param trs single transformation or list of transformations (see \code{link{getRotations}})
#' @param simplify logical, should image be returned in case if only one transformation was supplied
#'
#' @return list of transformed images
#' @export
applyTransforms = function(i,trs=getRotations(),simplify=TRUE){
  if(is.function(trs[[1]]))
    trs = list(trs)
  r = lapply(trs,function(t){
    for(f in t)
      i = f(i)
    i
  })
  if(simplify & length(r) == 1)
    r = r[[1]]
  r
}


#' Classifies visium spots according to its position relative to tissue slice border
#'
#' @param rc either seurat object or seu@images$slice1@coordinates dataframe
#'
#' @return list with two elements:
#' 1. augmented rc dataframe with spot coordinates. Following columns added:
#'  tissue.piece - number of tissue piece
#'  is.border - specifies whether spot is tissue border
#'  border.inx - consecutive number of border spots
#' 2. nj - list of spot neighbors
#'
#' @export
#'
#' @examples
findTissueBorder = function(rc){
  require(igraph)
  if(class(rc) == 'Seurat')
    rc = rc@images$slice1@coordinates
  rc$name = paste0(rc$row,'_',rc$col)
  njshifts = rbind(c(-1,-1),
                   c(-1, 1),
                   c( 1,-1),
                   c( 1, 1),
                   c( 0,-2),
                   c( 0, 2))
  nj = do.call(rbind,lapply(1:nrow(rc),function(i)data.frame(inx=i,nj.names=paste0(rc$row[i]+njshifts[,1],'_',rc$col[i]+njshifts[,2]))))
  nj = nj[nj$nj.names %in% rc$name,]
  nj$nj.inx = match(nj$nj.names,rc$name)
  nj$nj.is.tissue = rc$tissue[nj$nj.inx]
  intis = nj$nj.is.tissue==1 & rc$tissue[nj$inx]==1

  rc$tissue.piece = NA
  t = components(graph(t(as.matrix(nj[intis,c('inx','nj.inx')]))))$membership
  rc$tissue.piece[1:length(t)] = t
  rc$tissue.piece[rc$tissue==0] = NA

  nj = split(nj,rownames(rc)[nj$inx])
  rc$is.border = rc$tissue == 1 & sapply(nj,function(x)sum(x$nj.is.tissue==0)>0)
  rc$nnj = sapply(nj,nrow)[rownames(rc)]
  rc$nnj[is.na(rc$nnj)] = 0
  # number border
  rc$inx = 1:nrow(rc)
  rc$border.inx = NA
  no = 1
  for(p in as.numeric(names(sort(table(rc$tissue.piece),decreasing = TRUE)))){
    repeat{
      binx = which(!is.na(rc$tissue.piece) & rc$tissue.piece==p & is.na(rc$border.inx) & rc$is.border)
      if(length(binx)==0) break
      chain = binx[order(rc$nnj[binx])[1]]
      repeat{
        cinx = chain[length(chain)]
        if(is.na(rc$border.inx[cinx])){
          rc$border.inx[cinx] = no
          no = no + 1
        }
        cnj = rc[nj[[rownames(rc)[cinx]]]$nj.inx,]
        cenj = cnj[is.na(cnj$tissue.piece),] # empty njs
        cnj = cnj[cnj$is.border & is.na(cnj$border.inx),]
        if(nrow(cnj) == 0){
          if(length(chain)==1)
            break
          chain = chain[-length(chain)]
        }else{
          cnj$cmn.enj = sapply(rownames(cnj),function(i)sum(nj[[i]]$nj.inx %in% cenj$inx))
          cinx = cnj$inx[order(cnj$cmn.enj,decreasing = TRUE)[1]]
          chain = c(chain,cinx)
        }
      }
    }
  }
  list(rc=rc,nj=nj)
}

#' Classifies visium spots according to its position relative to tissue slice border
#'
#' @param rc,nj output of findTissueBorder
#'
#' @return augmented rc dataframe with spot coordinates. Following columns added:
#' nearest.border.inxs.graph - list of indexes (as specified by border.inx) of all nearestest border spots
#' dist2border.graph - distance to tissue border (in spots, by hex-graph of spots)
#' nearest.border.pos.graph - mean of nearest.border.inxs.graph
#'
#' please ignore other coumns added
#' @export
#'
#' @examples
calcDistance2border = function(rc,nj){
  # each tissue point should know its closest border and distance to it
  # by phisical distance
  rc$nearest.border.inx = NA
  rc$dist2border = NA
  dist = as.matrix(dist(rc[,c('imagerow','imagecol')],m='manh'))
  dist = dist/min(dist[dist>0])
  for(p in unique(na.omit(rc$tissue.piece))){
    f = !is.na(rc$tissue.piece) & rc$tissue.piece == p
    if(sum(rc$is.border & f)==0) next
    d = dist[f,rc$is.border & f,drop=F]
    d2b = do.call(rbind,apply(d,1,function(x){i=which.min(x);data.frame(i=i,dist=x[i])}))
    rc$nearest.border.inx[f] = rc[colnames(d)[d2b$i],'border.inx']
    rc$dist2border[f] = d2b$dist
  }
  # graph distance
  rc$nearest.border.inxs.graph = as.character(rc$border.inx)
  rc$dist2border.graph = NA

  spots = rownames(rc)[rc$is.border]
  d = 1
  rc$dist2border.graph[rc$is.border] = 0
  repeat{
    cnj = do.call(rbind,nj[spots])
    crc = rc[cnj$nj.inx,]
    f = !is.na(crc$tissue.piece) & is.na(crc$nearest.border.inxs.graph)
    if(sum(f)==0) break
    cnj = cnj[f,]
    crc = crc[f,]
    border.inxs = sapply(split(rc$nearest.border.inxs.graph[cnj$inx],cnj$nj.inx),function(x){paste(sort(unique(unlist(strsplit(x,',',TRUE)))),collapse = ',')})
    inxs = as.numeric(names(border.inxs))
    if(length(inxs)==0) break
    rc$dist2border.graph[inxs] = d
    rc$nearest.border.inxs.graph[inxs] =  border.inxs
    spots = rownames(rc)[inxs]
    d = d+1
  }
  rc$nearest.border.pos.graph = sapply(strsplit(rc$nearest.border.inxs.graph,',',TRUE),function(x)mean(as.numeric(x)))
  rc
}

#' Resizes H&E image
#'
#' @param v Seurat visium image
#' @param wpx desired image width (pixels)
#'
#' @return modified Seurat object
#' @export
#'
#' @examples
scaleVisiumImage = function(v,wpx=500){
  require(EBImage)
  coef = wpx/dim(v@images$slice1@image)[1]
  v@images$slice1@image = EBImage::resize(v@images$slice1@image,w=wpx)
  v@images$slice1@coordinates$imagerow = v@images$slice1@coordinates$imagerow*coef
  v@images$slice1@coordinates$imagecol = v@images$slice1@coordinates$imagecol*coef
  v
}

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
#' @return
#' @export
#'
#' @examples
plotNMFCons = function(coefs,cons,clcols=NULL,max.cex = 4.5/14*8,
                       colfun=function(x)num2col(x,c('blue','gray','orange','violet','black'),minx = 0,maxx = 1),
                       title='',ylab.cex=1,xlab.cex=1){
  cons = cons[colnames(coefs),colnames(coefs)]
  cls = apply(coefs,2,which.max)
  if(is.null(clcols))
    clcols = RColorBrewer::brewer.pal(nrow(coefs),'Set3')

  o = names(cls)[order(cls)]
  layout(matrix(1:3,ncol=3),widths = c(0.9,0.3,1.8))
  parl = par()
  par(mar=c(3,9,2.5,0),bty='n',oma=c(0,0,1,0),cex=1,tcl=-0.2,mgp=c(1.2,0.3,0),las=1)

  dotPlot(t(coefs[,o]),rowColours = cbind(clcols[cls[o]]),colColours = clcols,max.cex = max.cex,ylab.cex = ylab.cex,xlab.cex = xlab.cex,scaleWM=F,colfun=colfun)

  slh = cluster::silhouette(cls,dmatrix=1-cons)
  rownames(slh) = names(cls)


  par(mar=c(3,0,2.5,0.3))
  b=barplot(slh[rev(o),3],horiz = T,width = 1,space = 0, yaxs = "i",ylim=c(-1.5,length(o)+0.5),names.arg = '',border=NA,col=clcols[cls[rev(o)]],xlab='silhouette',xlim=c(-1,1),xaxt='n')
  par(cex=0.6)
  axis(1,)
  par(cex=1,xaxt='n',yaxt='n')
  dotPlot(cons[o,o],max.cex = 1.3,rowColours = cbind(clcols[cls[o]]),colColours = clcols[cls[o]],main='Consensus clustering matrix')

  mtext(title,3,outer = TRUE,line = 0)
  w = options(warn = -1)
  par(parl)
  options(w)
  invisible(list(clusters=cls,cols=clcols,slh=slh))
}

