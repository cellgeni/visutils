cropVisiumImage = function(v,r){
  c = v@images$slice1@coordinates
  scalefactors = v@images$slice1@scale.factors
  img = v@images$slice1@image

  c[,4:5] = c[,4:5]*scalefactors$lowres

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
    xy=plotVisiumImg(xy,v@images$slice1@image,v@images$slice1@scale.factors$lowres,cex=cex,col=col,border=border,xaxt='n',yaxt='n',...)
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

symbols.pie = function(x,y,r,d,cols,bg=NA,alpha0=0){
  r = recycle(r,1:length(x))
  d = sweep(d,1,apply(d,1,sum),'/')*2*pi
  for(i in 1:nrow(d)){
    a = alpha0
    for(j in 1:ncol(d)){
      if(d[i,j]>0){
        as = c(seq(a,a+d[i,j],by=pi/18),a+d[i,j])
        polygon(c(x[i],x[i]+r[i]*cos(as),x[i]),c(y[i],y[i]+r[i]*sin(as),y[i]),col = cols[j],border=bg)
        a = a+d[i,j]
      }
    }
  }
  # plot(1,t='n',xlim=c(0,5),ylim=c(0,10))
  # d = matrix(runif(6*4),ncol=4)
  # symbols.pie(c(1,2,3,1,2,3),c(1,1,1,1,1,1),1/2,d,1:4)
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

myLoad10X_Spatial = function(data.dir,filter.matrix=TRUE,ens_id=TRUE,...){
  d = Load10X_Spatial(data.dir,ifelse(filter.matrix,'filtered_feature_bc_matrix.h5','raw_feature_bc_matrix.h5'),filter.matrix=filter.matrix,...)
  if(ens_id){
    gids = read.table(paste0(data.dir,'/raw_feature_bc_matrix/features.tsv.gz'))
    d@assays$Spatial@counts@Dimnames[[1]] = gids$V1
    d@assays$Spatial@data@Dimnames[[1]] = gids$V1
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
