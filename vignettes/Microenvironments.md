---
title: "Microenvironments"
#output: rmarkdown::html_vignette
output: 
  github_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{Microenvironments}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
# install packages if necessary
# devtools::install_github("cellgeni/visutils")
library(visutils)
# package to load h5ad file as Seurat objects
library(schard) 
library(Seurat)
library(Matrix)
library(NMF)
library(plyr)
```

# Introduction
Visium data analyses frequently include deconvolution of gene expression matrix into celltype abundances using tools such as [cell2location](https://github.com/BayraktarLab/cell2location) or [RCTD](https://github.com/dmcable/spacexr). These methods allows to predict per spot abundance matrix using reference single cell dataset. The next step is to identify microenvironments: groups of co-locating cells. This can be achived using non-negative matrix factorization (NMF). This vignete will demonstrate to use NMF and visutils packages to do so, 

# Download the data
We will use data from <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13084>, lets first load metadata:

```r
meta = read.table('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/084/E-MTAB-13084/Files/E-MTAB-13084.sdrf.txt',header = TRUE,check.names = FALSE,sep='\t',quote = '')
# samples are duplicated for each fastq, lets collapse them
colnames(meta) = gsub('Characteristics\\[|]','',colnames(meta))
meta = unique(meta[,c('Source Name','age','sex','sampling site','disease','sample id')])
meta$`body part` = splitSub(meta$`sample id`,'_',1)
rownames(meta) = meta$`Source Name`
ord = c(body=1,face=2,bcc=3)
meta = meta[order(ord[meta$`body part`]),]
meta[1:3,]
#>                             Source Name age  sex            sampling site disease       sample id body part
#> WSSKNKCLsp104466211 WSSKNKCLsp104466211  60 male                     back  normal     body_back1a      body
#> WSSKNKCLsp10446623   WSSKNKCLsp10446623  55 male inguinal part of abdomen  normal body_inguinal1a      body
#> WSSKNKCLsp10767965   WSSKNKCLsp10767965  47 male                  abdomen  normal  body_abdomen1b      body
```


We'll take h5ad data from <https://spatial-skin-atlas.cellgeni.sanger.ac.uk/> and load them using hchard package. These h5ad contains cell2location celltype abundance predictions that we need as input for NMF. For sake of time we will use samples taken from temple:


```r
# download data to to temporary location and load as Seurat object
metat = meta[meta$`sampling site`=='temple',]
tmpfile = tempfile()
vs = list()
for(i in 1:nrow(metat)){
  tryCatch({
    download.file(paste0('https://cellgeni.cog.sanger.ac.uk/spatial-skin-atlas/download/',metat$`Source Name`[i],'.h5ad'),
                        tmpfile,quiet = TRUE)
    vs[[metat$`Source Name`[i]]] = schard::h5ad2seurat_spatial(tmpfile,use.raw = TRUE,img.res = 'hires')
    file.remove(tmpfile)
    cat('.')
  },warning=function(w){cat('!')})
}
#> ....
print('\n# spots:')
#> [1] "\n# spots:"
sapply(vs,ncol)
#> WSSKNKCLsp10446613 WSSKNKCLsp10446614 WSSKNKCLsp10446615 WSSKNKCLsp10446616 
#>                986                927               1547               1494
```

# Run NMF
lets first extract cell2location prediction from Seurat objects, they are in metadata prefixed with 'c2l_'

```r
celltypes = colnames(vs[[1]]@meta.data)
celltypes = celltypes[startsWith(celltypes,'c2l_')]
c2l = do.call(rbind,lapply(vs,function(v)v@meta.data[,celltypes]))
colnames(c2l) = sub('c2l_','',colnames(c2l))
c2l = as.matrix(c2l)
# we will use per-spot normalazed celltype abundancies
c2l = sweep(c2l,1,rowSums(c2l),'/')

# we will run nmf N times to assess results stubility
N = 50 # consider to increase it
rank = 6 # number of factors is almost arbitrary and to be set manually in dependence on desired granularity of microenvironments. here it is 6 assuming about 5 celltypes per ME
set.seed(1234)
doMC::registerDoMC(4)
nmf = llply(1:N,function(i){nmf(c2l, rank = rank)},.parallel = T)
```
Next we will calculate consensus clustering matrix across all nmf runs and choose the one that agees with consensus matrix most:

```r
best.nmf = nmfGetBest(nmf,getNMFNormFs('max'))
plotNMFCons(best.nmf$coefn,best.nmf$cons/N,ylab.cex=0.7,max.cex = 2)
```

![plot of chunk unnamed-chunk-5](Microenvironments_files/unnamed-chunk-5-1.png)
Lets visualize spatial ME distribution

```r
par(mfrow=c(4,6),mar=c(0.5,0.5,1,5),bty='n')
for(i in 1:nrow(metat)){
  for(me in 1:rank){
    x = best.nmf$basisn[paste0(rownames(metat)[i],'.',colnames(vs[[i]])),me]
    plotVisium(vs[[i]],x,main=paste0(metat$`sample id`[i],'; ME',me))
  }
}
```

![plot of chunk unnamed-chunk-6](Microenvironments_files/unnamed-chunk-6-1.png)
