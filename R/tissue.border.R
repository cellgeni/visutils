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
