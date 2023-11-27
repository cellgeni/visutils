library(visutils)
# devtools::load_all()
library(Seurat)

# load data ########
vl = schard::h5ad2seurat_spatial('../test.h5ad/vis.heart.h5ad',simplify = FALSE)
vm = merge(vl[[1]],vl[-1])
v1 = schard::h5ad2seurat_spatial('../test.h5ad/single.vis.fat.WSSKNKCLsp12887264.h5ad',use.raw = TRUE)

mbf = readRDS('../test.h5ad/slide_seq_ens.obs.f1.rds')
mbf.tiles = readRDS('../test.h5ad/slide_seq_tiles.f1.rds')
rects = data.frame(x=rep(1:100,each=100),y=rep(1:100,times=100))
rects$z = sin(abs(rects$x-50)^2/185 + abs(rects$y-50)^3/1137)

# test plotting #######
# _single sample ########

pdf('plot.test.pdf',w=8*3.5,h=4*3)
#par(mfrow=c(2,2),bty='n',mar=c(0,0,1,5))
par(mfrow=c(4,8),bty='n',mar=c(0,0,1,5))
plotVisium(v1,v1$total_counts,cex=0.7)
plotVisium(v1,v1$nCount_Spatial,type='hex')
plotVisium(v1,v1$nCount_Spatial,type='xy')
e = v1@meta.data[,c('c2l_Diff and undiff Keratinocytes','c2l_POSTN+ fibroblasts','c2l_Melanocytes')]
plotVisium(v1, pie.fracs=e,pie.col=c('red','blue','green'))
plotVisiumMultyColours(v1, e,col=c('red','blue','green'),min.opacity = 50,img.alpha=0.3,he.grayscale=TRUE)

# _from multiple ######
plotVisium(vm,vm$total_counts,cex=0.7,image.name = 'HCAHeartST11702009')
e = vm@meta.data[,c('prop_Adip1','prop_EC5_art','prop_ILC')]
plotVisium(vm, pie.fracs=e,pie.col=c('red','blue','green'),image.name = 'HCAHeartST11702009')
plotVisiumMultyColours(vm, e,col=c('red','blue','green'),min.opacity = 250,img.alpha=0.3,he.grayscale=TRUE,image.name = 'HCAHeartST11702009')

# _rect, tiles #####
plotVisium(mbf.tiles,mbf$n_counts,type='tiles',zfun=log1p)
plotVisium(rects,rects$z,type='rect')

# test crop and rotation
# _single sample ######
v1_1 = rotateVisium(v1,1)
v1_3 = rotateVisium(v1,3)
v1_m = rotateVisium(v1,0,mirror = TRUE)
plotVisium(v1,v1$nCount_Spatial,cex=0.6)
plotVisium(v1_1,v1_1$nCount_Spatial,cex=0.6)
plotVisium(v1_m,v1_m$nCount_Spatial,cex=0.6)
plotVisium(v1_3,v1_3$nCount_Spatial,cex=0.6)

# _mult samples #######
image.name = 'HCAHeartST11702010'
vm_1 = rotateVisium(vm,1,image.name = image.name)
vm_3 = rotateVisium(vm,3,image.name = image.name)
vm_m = rotateVisium(vm,0,mirror = TRUE,image.name = image.name)
plotVisium(vm,vm$nCount_Spatial,cex=0.6,image.name=image.name)
plotVisium(vm_1,vm_1$nCount_Spatial,cex=0.6,image.name=image.name)
plotVisium(vm_m,vm_m$nCount_Spatial,cex=0.6,image.name=image.name)
plotVisium(vm_3,vm_3$nCount_Spatial,cex=0.6,image.name=image.name)

# _crop
v1_c = cropVisiumImage(v1)
vm_c = cropVisiumImage(vm,image.name = image.name)
plotVisium(v1,z='red',cex=0.6)
plotVisium(v1_c,z='red',cex=0.6)
plotVisium(vm,z='red',cex=0.6,image.name = image.name)
plotVisium(vm_c,z='red',cex=0.6,image.name = image.name)

# scale
v1_s = scaleVisiumImage(v1,100)
vm_s = scaleVisiumImage(vm,100,image.name = image.name)
plotVisium(v1,z='red',cex=0.3)
plotVisium(v1_s,z='red',cex=0.3)
plotVisium(vm,z='red',cex=0.3,image.name = image.name)
plotVisium(vm_s,z='red',cex=0.3,image.name = image.name)


# border
v1_b = findTissueBorder(v1)
# vm (visium from cellxgene) is not working because the object has no mesh coordinates
v1_b$rc[1:2,]
plotVisium(v1,z=as.character(v1_b$rc$tissue.piece))
plotVisium(v1,z=v1_b$rc$is.border) # not working because there are no empty spots


# merge spots
plotVisium(v1,v1$nCount_Spatial,zfun=log1p)
c = getCenters(v1,to.merge = v1$total_counts<1000)
plotVisium(v1,c$group)
v1_mr = mergeSpots(v1,c)
range(colSums(v1_mr))
plotVisium(v1_mr,v1_mr$nCount_Spatial,cex=v1_mr$cex)
dev.off()
print('all done')
