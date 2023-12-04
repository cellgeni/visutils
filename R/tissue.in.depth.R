#' Annotate spots on border between two features
#'
#' find spots that belongs to one feature but contact with another
#'
#' @param v Seurat visium object
#' @param ann.column name of meta.data column that holds feature annotation
#' @param which name if feature border spots should belong to
#' @param contactTo name if feature border spots should contact with
#' @param image.name image to use, specify if object contains more than one sample
#'
#' @return data.frame with two columns logical 'junction' that specifies whether spot belongs to junction and numerical 'dist2junction' expressed in inter-spot distances
#' @export
defineJunction = function(v,ann.column,which,contactTo,image.name=NULL){
  dist = calcSpotDistance(v,image.name)
  junction = v@meta.data[,ann.column] == which &  apply(dist[v@meta.data[,ann.column] == contactTo,],2,min) <= 1.2
  dist2junction = apply(dist[junction,],2,min)
  dist2junction[v@meta.data[,ann.column] != which] = -dist2junction[v@meta.data[,ann.column] != which]
  r = data.frame(junction=junction,dist2junction=dist2junction)
  rownames(r) = colnames(v)
  r
}


#' Calculates distance from each spot to the spot set
#'
#' Finds distance to the closest spot from the set and returns it expressed in inter-spot distances
#'
#' @param v Seurat visium object
#' @param spots set of spots (character, intereger, or logical index of the spots in v)
#' @param image.name image to use, specify if object contains more than one sample
#'
#' @return numeric vector with distances to the nearest spot from set (in interspot distances)
#' @export
calcDistance2SpotSet = function(v,spots,image.name=NULL){
  dist = calcSpotDistance(v,image.name)
  apply(dist[spots,],2,function(x){ifelse(length(x)>0,min(x),Inf)})
}

#' Calculates matrix of spatial distances between spot centers
#'
#'and normolazes it by inter-spot distance
#'
#' @param v Seurat visium object
#' @param image.name image to use, specify if object contains more than one sample
#'
#' @return distance matrix
#' @export
calcSpotDistance = function(v,image.name=NULL){
  image.name = getImageName(v,image.name,stopIfMoreThanOne=TRUE)
  dist = as.matrix(dist(v@images[[image.name]]@coordinates[,c('imagerow','imagecol')]))
  spot.dist = min(dist[upper.tri(dist,diag = FALSE)])
  dist = dist/spot.dist
  dist
}

#' Summarises feature expression along specified factor and sample
#'
#' summary statistics will be calculated for each dist-sample-feature combination
#'
#' @param dist factor to summarise data along. Normally it is binnarised distance to something.
#' @param sample character vector specifying sample
#' @param data numeric matrix (probably sparse) with data to be summarised.
#' @param filter spot filter, the variable to be used to subset spots (will be applied to dist, sample, and data). No filtering applied if NULL (default)
#' @param FUN function to summarise the data (mean by default)
#' @param min.spots minimal namber of spots for statistics to be calculated (5 by default)
#' @param ncores number of cores to be used
#' @param per.spot.norm whether to perform per-spot normalization (makes sence for celltype abundancies)
#'
#' @return 3d (distance/feature/sample) numeric matrix with summarized data
#' @export
makeDistFeatureSampleTable = function(dist,sample,data,filter=NULL,FUN=mean,min.spots=5,ncores=1,per.spot.norm=TRUE){
  require(plyr)
  require(doMC)
  if(ncores>=1)
    doMC::registerDoMC(ncores)
  if(per.spot.norm)
    data = sweep(data,1,rowSums(data),'/')
  if(!is.null(filter)){
    sample = sample[filter]
    dist = dist[filter]
    data = data[filter,]
  }
  features = colnames(data)
  sids = unique(sample)
  dists = sort(unique(dist))

  dnames = list(dists,features,sids)

  ldctmtx = llply(dists,function(d){
    res = array(NA,dim=sapply(dnames[-1],length),dimnames=dnames[-1])
    flt = dist == d
    for(f in features){
      r = tapply(data[flt,f],sample[flt],function(x)ifelse(length(x)<min.spots,NA,FUN(x)))
      res[f,names(r)] = r
    }
    res
  },.parallel = ncores>1)
  dctmtx = array(NA,dim=sapply(dnames,length),dimnames=dnames)
  for(i in 1:length(ldctmtx)){
    dctmtx[i,,] = ldctmtx[[i]]
  }
  dctmtx
}

#' Compare features bewteeb two conditions
#'
#' per distance bean
#'
#' @param m 3d (distance/feature/sample) numeric matrix with summarized data (output of makeDistFeatureSampleTable)
#' @param f1,f2 indexes that specify samples belonging to first and second conditions
#' @param test.fun function to be used to test difference. Should return list of p.value item. t.test is default
#'
#' @return list of pv, fdr, m1, and m2 (condition means) distance/feature matricex
#' @export
testTDConditions = function(m,f1,f2,test.fun=t.test){
  pv = m[,,1]
  pv[,] = NA
  m1 = m2 = pv
  for(i in 1:dim(m)[1]){
    for(j in 1:dim(m)[2]){
      x1 = m[i,j,f1]
      x2 = m[i,j,f2]
      if(sum(!is.na(x1))<2 || sum(!is.na(x2))<2){
        next
      }
      tryCatch({
        r = test.fun(x1,x2)$p.value
        pv[i,j] = r
        m1[i,j] = mean(x1,na.rm=T)
        m2[i,j] = mean(x2,na.rm=T)
      },error=function(e){})
    }
  }
  r = list(pv=pv,m1=m1,m2=m2)
  r$fdr = r$pv
  r$fdr[,] = p.adjust(r$pv,'fdr')
  r
}

#' Plots results of testTDConditions
#'
#' as three heatmaps: two shows mean per-condition spatial feature distribution and one show log fold change
#'
#' @param comp list, output of testTDConditions
#' @param cond.titles condition titles to be plot on top of heatmaps
#' @param order feature order, NULL (default) for no reordering
#' @param l2fc.zlim range to trim log fold change
#' @param feature.class list (or matrix of columns) with feature annotations to be shown by color
#' @param cols cols to be used for condition heatmaps
#' @param fdr.thrs named numerical list used to show levels of significance on log fold change plot, default is c('*'=0.05,'.'=0.2)
#'
#' @export
plotTD.HM = function(comp,cond.titles=c('cond1','cond2'),order=NULL,l2fc.zlim=NULL,feature.class=NULL,log=TRUE,
                     col = hcl.colors(100, "YlOrRd", rev = TRUE),fdr.thrs=c('*'=0.05,'.'=0.2)){
  if(!is.null(order)){
    comp = lapply(comp,function(x)x[,order])
  }
  if(!is.null(feature.class) & !is.list(feature.class))
    feature.class = list(feature.class)
  if(log)
    l2fc = log2(comp$m2/comp$m1)
  else
    l2fc = comp$m2 - comp$m1
  fdr.thrs = sort(fdr.thrs,decreasing = TRUE)
  pvt = comp$pv
  pvt[] = ''
  for(i in 1:length(fdr.thrs)){
    pvt[comp$fdr<=fdr.thrs[i]] = names(fdr.thrs)[i]
  }

  # scale per feature
  m1 = sweep(comp$m1,2,apply(comp$m1,2,max,na.rm=TRUE),'/')
  m2 = sweep(comp$m2,2,apply(comp$m2,2,max,na.rm=TRUE),'/')

  xlim = c(0.5,0.5+nrow(l2fc))
  imageWithText(m1,'',main=cond.titles[1],xlab='Distance to epidermis-dermis interface (spots)',col=col,rowAnns = feature.class)
  imageWithText(m2,'',main=cond.titles[2],xlab='Distance to epidermis-dermis interface (spots)',col=col,rowAnns = feature.class)
  if(is.null(l2fc.zlim)){
    l2fc.zlim=range(l2fc,na.rm=T)
  }
  l2fc[l2fc<l2fc.zlim[1]] = l2fc.zlim[1]
  l2fc[l2fc>l2fc.zlim[2]] = l2fc.zlim[2]
  zmax=max(abs(l2fc.zlim))
  imageWithText(l2fc,pvt,zlim=c(-zmax,zmax),main=paste0('log2(',cond.titles[2],'/',cond.titles[1],')'),
                xlab='Distance to epidermis-dermis interface (spots)',rowAnns = feature.class,col=hcl.colors(100, "Blue-Red 3"))

  if(any(!is.na(l2fc))){
    zlim=range(l2fc,na.rm=T)
    plotColorLegend2(1,1.2,0.2,0.8,c(-zmax,zmax),zlim,identity,function(x)num2col(x),title='log2FC')
  }
  legend(grconvertX(1,'npc','user'),grconvertY(1,'npc','user'),xpd=NA,pch=names(fdr.thrs),legend=paste0('<',fdr.thrs),bty='n',title='FDR')
}


#' Plots features profiles
#'
#' along distance. Shows mean abundance as line and CI as shade.
#'
#' @param m 3d (distance/feature/sample) numeric matrix with summarized data (output of makeDistFeatureSampleTable)
#' @param features features (indexes along third dimension of m)
#' @param cols names (by feature names) colour vector
#' @param sd.mult width of CI shown by shade expressed in standart deviations of mean
#' @param legend. logical, whether logend should be plotted, or list with parameters to be passed to legend function.
#' @param ylim ylim to be set, set by data if NULL (default)
#' @param scaleY whether all features should be scaled by its max value (to bring all features to the same scale)
#' @param area.opacity opacity of CI shade
#' @param lwd line width
#' @param xlab,ylab,main plot annotation
#' @param ... other parameters to be passed to plot
#'
#' @export
plotFeatureProfiles = function(m,features,cols=NULL,sd.mult=2,legend.=TRUE,ylim=NULL,scaleY=TRUE,area.opacity = 0.2,lwd=2,xlab='Distance (spots)',
                               ylab='Relative abundance',main='',...){
  if(is.null(cols))
    cols = char2col(features)
  x = as.numeric(dimnames(m)[[1]])
  areas = lapply(features,function(ct){
    do.call(rbind,apply(m[,ct,],1,function(x){
      x = x[!is.na(x)]
      data.frame(mean=mean(x),sd=sd(x)/sqrt(length(x)),n = length(x))
    }))
  })
  names(areas) = features
  if(scaleY)
    areas = lapply(areas,function(a){a[,1:2]=a[,1:2]/max(a$mean);a})
  if(is.null(ylim))
    ylim = range(sapply(areas,function(a)c(min(a$mean-a$sd*sd.mult),max(a$mean+a$sd*sd.mult))))
  plot(1,t='n',xlim=range(x),ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
  for(n in names(areas))
    plotArea(x,areas[[n]][,1:2],sd.mult = sd.mult,col=cols[n],new = FALSE,lwd=lwd,area.transp=area.opacity)
  if(!is.null(legend.) && (!is.logical(legend.) || legend.)){
    if(!is.list(legend.))
      legend. = list()
    if(is.null(legend.$x)){
      legend.$x = grconvertX(1,'npc','user')
      legend.$y = grconvertY(1,'npc','user')
    }
    legend.$fill = cols
    legend.$legend = names(cols)
    legend.$bty=par('bty')
    legend.$border = NA
    legend.$xpd = NA
    do.call(legend,legend.)
  }
}

#' Plot feature profile along distance for multiple conditions
#'
#' @param m 3d (distance/feature/sample) numeric matrix with summarized data (output of makeDistFeatureSampleTable)
#' @param feature single character, name of feature to be plotted (from second dimension of m)
#' @param conditions vector specifying conditions of all samples (along third dimension of m)
#' @param cols named (by conditions) color vector. Only conditions mentioned here will be shown. Colours are aoutogenerated by char2col if null (default)
#' @param sd.mult width of CI shown by shade expressed in standart deviations of mean
#' @param legend. logical, whether logend should be plotted, or list with parameters to be passed to legend function.
#' @param ylim ylim to be set, set by data if NULL (default)
#' @param scaleY whether all features should be scaled by its max value (to bring all features to the same scale)
#' @param area.opacity opacity of CI shade
#' @param lwd line width
#' @param xlab,ylab,main plot annotation
#' @param ... other parameters to be passed to plot
#'
#' @export
plotConditionsProfiles = function(m,feature,conditions,cols=NULL,sd.mult=2,legend.=TRUE,ylim=NULL,scaleY=FALSE,
                                  area.opacity = 0.2,lwd=2,xlab='Distance (spots)',
                                  ylab='Abundance',main='',...){

  if(is.null(cols)){
    uniq.conds = unique(conditions)
    cols = char2col(uniq.conds)
  }else{
    uniq.conds = intersect(names(cols),conditions)
  }

  x = as.numeric(dimnames(m)[[1]])
  areas = lapply(uniq.conds,function(cnd){
    do.call(rbind,apply(m[,feature,conditions==cnd,drop=FALSE],1,function(x){
      x = x[!is.na(x)]
      data.frame(mean=mean(x),sd=sd(x)/sqrt(length(x)),n = length(x))
    }))
  })
  names(areas) = uniq.conds
  if(scaleY)
    areas = lapply(areas,function(a){a[,1:2]=a[,1:2]/max(a$mean);a})
  if(is.null(ylim))
    ylim = range(sapply(areas,function(a)c(min(a$mean-a$sd*sd.mult),max(a$mean+a$sd*sd.mult))))
  plot(1,t='n',xlim=range(x),ylim=ylim,xlab=xlab,ylab=ylab,main=main)
  for(n in names(areas))
    plotArea(x,areas[[n]][,1:2],sd.mult = sd.mult,col=cols[n],new = FALSE,lwd=lwd,area.transp=area.opacity)
  if(!is.null(legend.) && (!is.logical(legend.) || legend.)){
    if(!is.list(legend.))
      legend. = list()
    if(is.null(legend.$x)){
      legend.$x = grconvertX(1,'npc','user')
      legend.$y = grconvertY(1,'npc','user')
    }
    legend.$fill = cols
    legend.$legend = names(cols)
    legend.$bty=par('bty')
    legend.$border = NA
    legend.$xpd = NA
    do.call(legend,legend.)
  }
}
