---
title: "Visualization"
output: 
  github_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{Tissue carving}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---






```r
# devtools::install_github("cellgeni/visutils",force = TRUE)
library(visutils)
library(Seurat)
```

# Load data

Lets take one random skin sample from <https://spatial-skin-atlas.cellgeni.sanger.ac.uk/>

```r
sid = 'WSSKNKCLsp10446623'
tmpfile = tempfile()
download.file(paste0('https://cellgeni.cog.sanger.ac.uk/spatial-skin-atlas/download/',sid,'.h5ad'),tmpfile,quiet = TRUE)
vis = schard::h5ad2seurat_spatial(tmpfile,use.raw = TRUE,img.res = 'hires')
file.remove(tmpfile)
#> [1] TRUE
```


```r
par(mar=c(0,0,1,5),bty='n')
plotVisium(vis,vis$nCount_Spatial,legend.args = list(title='UMI'))
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

# Plot multiple microenvironments


```r
c2l = as.matrix(vis@meta.data[,grep('c2l',colnames(vis@meta.data))])
colnames(c2l) = sub('c2l_','',colnames(c2l))
celltypes = char2col(c('Suprabasal keratinocytes','APOD+ fibroblasts','Basal keratinocytes','Melanocytes'))
par(mfrow=c(2,1),mar=c(0,0,1,20),bty='n')
plotVisiumMultyColours(vis,c2l[,names(celltypes)],cols = celltypes,img.alpha=0.5,legend.ncol = 2,min.opacity = 100)
# or as pie charts
c2l = sweep(c2l,2,apply(c2l,2,max),'/')
plotVisium(vis,pie.fracs = c2l[,names(celltypes)],pie.cols = celltypes,img.alpha=0.5)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

# Plot multiple genes


```r
cnts = vis@assays$Spatial@layers$counts
rownames(cnts) = vis@assays$Spatial@meta.data$`_index`

gids = char2col(c('KRT5','KRT10','COL1A2','PERP'))
cpm = t(as.matrix(cnts[names(gids),]))
cpm = sweep(cpm,1,vis$nCount_Spatial,'/')*1e4
par(mar=c(0,0,1,20),bty='n')
plotVisiumMultyColours(vis,cpm,cols = gids,img.alpha=0.5,legend.ncol = 2,min.opacity = 100)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)
