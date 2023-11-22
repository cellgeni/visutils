#' Classifies visium spots according to its position relative to tissue slice border
#'
#' @param rc either seurat object or seu@images[[.]]@coordinates dataframe
#'
#' @return list with two elements:
#' 1. augmented rc dataframe with spot coordinates. Following columns added:
#'  tissue.piece - number of tissue piece
#'  is.border - specifies whether spot is tissue border
#'  border.inx - consecutive number of border spots
#' 2. nj - list of spot neighbors
#'
#' @export
findTissueBorder = function(rc,image.name=NULL){
  require(igraph)

  if(class(rc) == 'Seurat'){
    image.name = getImageName(rc,image.name,stopIfMoreThanOne=TRUE)
    rc = rc@images[[image.name]]@coordinates
  }

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

  # problem: isolated spots (no idea how, but spaceranger makes them) have no items in nj
  emptynj = nj[c(),]
  nj = split(nj,rownames(rc)[nj$inx])
  # lets add empty items for them
  for(n in setdiff(rownames(rc),names(nj)))
    nj[[n]] = emptynj
  nj = nj[rownames(rc)]

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
calcDistance2border = function(rc,nj){
  warning("Not sure this fnction is useful. One can just use Euclidian distance")
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
