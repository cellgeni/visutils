#devtools::create('~/nfs/rcode/visutils/')
devtools::document()
devtools::load_all()
?log10p1

# usethis::create_github_token()
# gitcreds::gitcreds_set()
# usethis::edit_r_environ()
options(timeout=400)
devtools::install_github("iaaka/visutils")

devtools::install_local()
# git add -u
# git commit -m "fix "
# git push -u origin main
?dotPlot

library(visutils)
library(Seurat)
v = myLoad10X_Spatial('~/nfs/projects/2302.fetal.skin/data.nfs/my/visium/face.body/spacerangers/spaceranger200_count_46862_HCA_rFSKI13460601_GRCh38-2020-A/',filter.matrix = T)
v1 = rotateVisium(v,1)
vm1 = rotateVisium(v,1,mirror = T)
par(mar=c(0,0,1,1),mfrow=c(2,3))
plotVisium(v,v$nCount_Spatial,cex=0.5,zfun=log10)
plotVisium(v1,v1$nCount_Spatial,cex=0.5,zfun=log10)
plotVisium(vm1,vm1$nCount_Spatial,cex=0.5,zfun=log10)

plotVisium(v,v$nCount_Spatial,cex=1,zfun=log10,t='hex')
plotVisium(v1,v1$nCount_Spatial,cex=1,zfun=log10,t='hex')
plotVisium(vm1,vm1$nCount_Spatial,cex=1,zfun=log10,t='hex')
