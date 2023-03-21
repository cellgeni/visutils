#' Makes one step coarser spot mesh
#'
#' @param coors Seurat object or coors@images$slice1@coordinates
#' @param cstart column coordinate of first (in -2th row) cluster center. Numbers for 0 to 13 will give different groupings.
#' @param type orientation of mesh, 1 or -1
#' @param to.merge logical vector specifying spots to be merged, use v$nCount_Spatial<500 to merge spots with less than 500 UMIs. NULL is identicall to all TRUE.
#'
#' @return modified spot coordinate table with group column added
#' @export
#'
#' @details There are severals ways to do so in dependence on starting point, this is controlled by cstart and type parameters. Maybe some combinations of these parameters results in identicall groupings, I didn't check it properly. Anyway I would suggest to always use defaults since it should not be very important.
#'
#' @examples
#' c2 = getCenters(v2)
#' plotVisium(v2,c2$group)
getCenters = function(coors,cstart=0L,type=1L,to.merge=NULL){
  if(class(coors) == 'Seurat')
    coors = coors@images$slice1@coordinates
  if(!(type %in% c(-1,1))) stop("type should be either 1 or -1")
  cstep = 14L
  cstart = cstart %% cstep
  col.shift = 5L*type

  # row,col
  area = rbind(c(-1L, 1L),
               c(-1L,-1L),
               c( 0L,-2L),
               c( 0L, 2L),
               c( 1L,-1L),
               c( 1L, 1L))
  r = NULL
  for(row in (-2):(max(coors$row)+1)){
    for(cs in seq(cstart-cstep,max(coors$col)+cstep,by=cstep)){
      center = c(row,cs)
      r = rbind(r,cbind(rbind(as.data.frame(sweep(area,2,center,'+')),center),center.row=center[1],center.col=center[2]))
    }
    cstart = (cstart + col.shift) %% cstep
  }
  r = setNames(paste0(r[,3],'_',r[,4]),paste(r[,1],r[,2]))[paste(coors$row,coors$col)]
  names(r) = rownames(c)
  if(!is.null(to.merge))
    r = ifelse(to.merge,r,paste0(coors$row,'_',coors$col))
  coors$group = r
  r = strsplit(r,'_',fixed = T)
  coors$group.row = as.integer(sapply(r,'[',1))
  coors$group.col = as.integer(sapply(r,'[',2))
  coors
}

#' Merge grouped spots
#'
#' @param v Seurat object
#' @param gr output of getCenters
#'
#' @return Seurat object with merged spots
#' @export
#'
#' @examples
#' c2 = getCenters(v2,to.merge = v2$nCount_Spatial<500)
#' v2m = mergeSpots(v2,c2)
#' v2m@meta.data[1:20,]
#' cex = scaleTo(sqrt(v2m$nspots),1,sqrt(8),minx = 1,maxx = sqrt(7))*0.9
#' plotVisium(v2m,v2m$nCount_Spatial,zfun = log1p,cex=cex)
mergeSpots = function(v,gr){
  gr$nCount_Spatial = v$nCount_Spatial
  gr$nspots = as.numeric(table(gr$group)[gr$group])
  gr$barcode = rownames(gr)
  mtx_ = calcColSums(v@assays$Spatial@counts,gr$group,mean = FALSE,verbose = FALSE)

  # combine spot info
  # either center (if exists) or max covered spot
  gr_ = do.call(rbind,lapply(split(gr,gr$group),function(x){
    # if center is to be merged
    f = which(x$group == paste0(x$row,'_',x$col))
    if(length(f)==1){
      r = x[f,]
    }else{ # use max covered spot othervise
      r = x[which.max(x$nCount_Spatial),]
    }
    r$merged_spots = paste(x$barcode,collapse=';')
    r
  }))
  gr_ = gr_[colnames(mtx_),]
  rownames(gr_) = gr_$barcode
  colnames(mtx_) = gr_$barcode

  object <- CreateSeuratObject(counts = mtx_, assay = "Spatial")

  image = new(Class = "VisiumV1", image = v@images$slice1@image,
              scale.factors = v@images$slice1@scale.factors,
              coordinates = gr_[,1:6],
              spot.radius = v@images$slice1@spot.radius)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) = "Spatial"
  object[["slice1"]] = image
  object@assays$Spatial@meta.features = v@assays$Spatial@meta.features

  cols = setdiff(colnames(v@meta.data),c(colnames(object@meta.data)))
  object@meta.data = cbind(object@meta.data,v@meta.data[rownames(object@meta.data),cols])

  # add info about merged spots
  object@meta.data$nspots = gr_$group.size
  object@meta.data$merged_spots = gr_$merged_spots

  object
}
