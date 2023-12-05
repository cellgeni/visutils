# visutils
The package provides set of functions to facilitate Visium data analyses, QC, and visualization. The package is functional (that is not OOP), it uses base graphic (that is not ggplot) and it is Seurat-friendly.

It provides functions to perform analyses of gene expression or celltype abundancies (predicted by methods such as cell2location) in dependence on distance to some histological feature defined as set of spots ("tissue in depth"), microenvironment analyses (via NMF), merging low-covered spots to achieve reasonable coverage, splirring slides by tissue pieces and defining tissue border. See tutorials below for more information
# Installation
```
devtools::install_github("iaaka/visutils")
```
# Tutorials
* [Tissue in depth](vignettes/TissueInDepth.md) shows how to analyse celltype abundance and gene expression in dependence on distance to dermsit-to-epidermis junction in public skin dataset
Coming soon
* __Spot merge__. Normally spots with low coverage (say below 500 UMI) are excluded from analyses. In some sample it could result in removal of whole areas since some tissues such as dermis, adipose, or cartilage are frequently have low read counts. Merging adjacent low-covered spot together is an alternative to filter then away.
* __Microenvironments__. Soft-clustering of celltype based on their collocalization using non-negative matrix factorization (NMF).
* __Tissue cutting__. Semi-automatic debris-removal, splitting sample by tissue pieces (when multiple pieces we analused on same slide), tissue border demarcation.

