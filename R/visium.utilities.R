#' Rotates Seurat visium slide
#'
#' @param v Seurat object
#' @param n number of counter-clockwise rotations
#' @param mirror logical, whether to slide should be flipped horisontally
#'
#' @return rotated Seurat object
#' @export
#'
#' @examples
#' vm = rotateVisium(v,n=0,m=T)
#' SpatialFeaturePlot(v,'lnCount_Spatial')+SpatialFeaturePlot(vm,'lnCount_Spatial')
rotateVisium = function(v,n=1,mirror=FALSE,image.name=NULL){
  image.name = getImageName(v,image.name,stopIfMoreThanOne=TRUE)
  image = v@images[[image.name]]
  ts = getRotations()
  if(mirror){
    image@image = applyTransforms(image@image,ts$m0)
    t = image@coordinates$imagerow
    image@coordinates$imagerow = dim(image@image)[1]/image@scale.factors$lowres - t + 1/image@scale.factors$lowres
    image@coordinates$row = max(image@coordinates$row) - image@coordinates$row
  }
  if(n>0){
    image@image = applyTransforms(image@image,ts$o3)
    t = image@coordinates$imagecol
    image@coordinates$imagecol = image@coordinates$imagerow
    image@coordinates$imagerow = dim(image@image)[1]/image@scale.factors$lowres - t + 1/image@scale.factors$lowres

    t = image@coordinates$col
    image@coordinates$col = image@coordinates$row
    image@coordinates$row = max(t) - t
  }
  v@images[[image.name]] = image
  if(n>1)
    v = rotateVisium(v,n=n-1,mirror = FALSE,image.name=image.name)
  v
}


getImageName = function(v,image.name,stopIfMoreThanOne = TRUE){
  if(is.null(image.name)){
    if(length(v@images)==1 || !stopIfMoreThanOne)
      image.name = names(v@images)[1]
    else if(stopIfMoreThanOne)
      stop("This object contains more than one slide, please select which one you need")
  }
  image.name
}


#' Crop image in Seuart Visium object
#'
#' Retain square image that contains all spots (listed in the object). Spot coordinates adjasted accordingly.
#'
#' @param v Seurat object
#'
#' @return visium object with cropped image
#' @export
cropVisiumImage = function(v,image.name=NULL){
  image.name = getImageName(v,image.name,stopIfMoreThanOne=TRUE)
  image = v@images[[image.name]]

  c = image@coordinates
  scalefactors = image@scale.factors
  img = image@image

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
  spot.radius = v@images[[image.name]]@spot.radius * max(dim(v@images[[image.name]]@image))
  v@images[[image.name]]@spot.radius = spot.radius / max(dim(img))
  v@images[[image.name]]@coordinates = c
  v@images[[image.name]]@image = img

  v
}

#' Resizes H&E image
#'
#' @param v Seurat visium image
#' @param wpx desired image width (pixels)
#'
#' @return modified Seurat object
#' @export
scaleVisiumImage = function(v,wpx=500,image.name=NULL){
  require(EBImage)
  image.name = getImageName(v,image.name,stopIfMoreThanOne=TRUE)
  image = v@images[[image.name]]

  coef = wpx/dim(image@image)[1]
  image@image = EBImage::resize(image@image,w=wpx)
  image@coordinates$imagerow = image@coordinates$imagerow*coef
  image@coordinates$imagecol = image@coordinates$imagecol*coef
  v@images[[image.name]] = image
  v
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



#' Adjust image brightnes and contrast
#'
#' @param p image (3d numeric array)
#' @param wb logical, specifies whether output image should be transformed to grayscale
#' @param qs quantiles to trim. Numerical vector with two items. Trims all values outside of specified quantile range.
#' @param trim01 logical, specifies whether pixels with zero (black) and maximal (white) intensity should be trimed ahead of quantile trimming.
#'
#' @return image (3d numeric array)
#' @export
enhanceImage = function(p,wb=FALSE,qs=NULL,trim01 = TRUE){
  pm0 = pm = pmax(p[,,1],p[,,2],p[,,3])
  if(!is.null(qs)){
    if(trim01){
      f = pm == 0 | pm == 1
      qs = quantile(pm[!f],probs = qs)
    }else
      qs = quantile(pm,probs = qs)
    pm = (pm-qs[1])/(qs[2]-qs[1])
    pm[pm>1] = 1
    pm[pm<0] = 0
  }
  f = pm0/pm
  f[is.na(f)] = 1
  p = sweep(p,1:2, f,'/')
  if(wb){
    z=(p[,,1]+p[,,2]+p[,,3])/3
    for(j in 1:3)
      p[,,j] = z
  }
  p
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
myLoad10X_Spatial = function(data.dir,filter.matrix=TRUE,ens_id=TRUE,slice='slice1',...){
  #TODO: rename)
  if(dir.exists(paste0(data.dir,'/outs'))){
    data.dir = paste0(data.dir,'/outs')
  }
  d = Seurat::Load10X_Spatial(data.dir = data.dir,
                      filename = ifelse(filter.matrix,'filtered_feature_bc_matrix.h5','raw_feature_bc_matrix.h5'),
                      filter.matrix=filter.matrix,
                      use.name=!ens_id,
                      slice=slice,...)
  d[['is.tissue']] = d@images[[slice]]@coordinates$tissue
  d
}


#' Download Visium dataset from 10x website
#'
#' @param url url of the dataset (whithot dataset name, see example)
#' @param sample.name name of the dataset
#' @param outdir folder to save dataset (function will create subfolder named by \code{sample.name})
#'
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
getMeanSpotColor = function(vis,scalefactors,image.name=NULL){
  image.name = getImageName(v,image.name,stopIfMoreThanOne=TRUE)
  image = v@images[[image.name]]

  i = image@image
  c = image@coordinates
  r = scalefactors$spot_diameter_fullres/2*scalefactors$tissue_lowres_scalef
  c = c[,4:5]*scalefactors$tissue_lowres_scalef
  r = t(sapply(1:nrow(c),function(j)getMeanSpotColor1(i,c$imagerow[j],c$imagecol[j],r)))
  rownames(r) = colnames(vis)
  r
}



#' Find best match between image and RNA counts
#' to deal with image swaps/rotations
#'
#' @param m
#'
#' @return data.frame with best mtches
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




