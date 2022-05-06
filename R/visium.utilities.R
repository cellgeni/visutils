#' Crop image in Seuart Visium object
#'
#' Retain square image that contains all spots (listed in the object). Spot coordinates adjasted accordingly.
#'
#' @param v Seurat object
#'
#' @return
#' @export
#'
#' @examples
cropVisiumImage = function(v){
  c = v@images$slice1@coordinates
  scalefactors = v@images$slice1@scale.factors
  img = v@images$slice1@image

  c[,4:5] = c[,4:5]*scalefactors$lowres
  r = min(dist(c[,4:5]))*0.5

  rr = range(c$imagerow)
  rc = range(c$imagecol)

  cc = round(mean(rc))
  cr = round(mean(rr))

  #ir = ceiling(max(rr[2]/2 - rr[1]/2,rc[2]/2 - rc[1]/2)) + ceiling(r)
  irr = ceiling(rr[2]/2 - rr[1]/2 + r)
  irc = ceiling(rc[2]/2 - rc[1]/2 + r)

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
#' @export
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
        polygon(c(x[i],x[i]+r[i]*cos(as),x[i]),c(y[i],y[i]+r[i]*sin(as),y[i]),col = cols[j],border=border)
        a = a+d[i,j]
      }
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
#'
#' @examples
myLoad10X_Spatial = function(data.dir,filter.matrix=TRUE,ens_id=TRUE,...){
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

#' Adjast image brightnes and contrast
#'
#' @param p image (3d numeric array)
#' @param wb logical, specifies whether output image should be transformed to grayscale
#' @param pow
#' @param qs quantiles to trim. Numerical vector with two items. Trims all values outside of specified quantile range.
#'
#' @return image (3d numeric array)
#' @export
#'
#' @examples
enhanceImage = function(p,wb=FALSE,pow=3,qs=c(0.01,0.99)){
  p = (1-p)^pow
  pm = apply(p,1:2,max)
  qs = quantile(pm,probs = qs)
  pm = (pm-qs[1])/(qs[2]-qs[1])
  pm[pm>1] = 1
  pm[pm<0] = 0
  p = sweep(p,1:2, apply(p,1:2,max)/pm,'/')
  p[is.na(p)] = 1
  p = 1-p
  if(wb){
    z=apply(p,1:2,mean)
    for(j in 1:3)
      p[,,j] = z
  }
  p
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
#'
#' @examples
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
#' @param xy
#' @param img
#' @param scale.factor
#' @param cex
#' @param col
#' @param border
#' @param spot.dist
#' @param img.alpha
#' @param xlim
#' @param ylim
#' @param symmetric.lims
#' @param xlab
#' @param ylab
#' @param pie.fraqs if specified plots pies instead of simple cycles. Matrix with number of rows equal to the number of spots, and number of columns equal to pie pieces.
#' @param pie.cols colors to be used for pie pieces (ncol(pie.fraqs) should be equal to length(pie.cols))
#' @param pie.min.fraq all pieces with relative size less than \code{pie.min.fraq} will be discared
#' @param ... other parameters to be passed to \code{plot} function
#'
#' @return
#' @export
#'
#' @examples
plotVisiumImg = function(xy,img,scale.factor,cex=1,col='red',border=NA,spot.dist=NULL,img.alpha=1,xlim=NULL,ylim=NULL,symmetric.lims=TRUE,xlab='',ylab='',pie.fraqs=NULL,pie.cols=NULL,pie.min.fraq=0.05,...){
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
                d=pie.fraqs,
                cols= pie.cols)
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

applyTrans = function(i,trs=transforms){
  lapply(trs,function(t){
    for(f in t)
      i = f(i)
    i
  })
}
