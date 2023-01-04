------------------------------------------------------------------------

# CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation, and visualization

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/taoshengxu/CancerSubtypes?branch=master&svg=true)](https://ci.appveyor.com/project/taoshengxu/CancerSubtypes)
[![bioc](http://www.bioconductor.org/shields/downloads/CancerSubtypes.svg)](http://bioconductor.org/packages/stats/bioc/CancerSubtypes.html)
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/CancerSubtypes.svg)](http://bioconductor.org/packages/CancerSubtypes/)
[![bioc](http://bioconductor.org/shields/availability/devel/CancerSubtypes.svg)](http://bioconductor.org/packages/CancerSubtypes/)
[![bioc](http://www.bioconductor.org/shields/build/release/bioc/CancerSubtypes.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/CancerSubtypes.html)

# This is a fork repo, but you are welcome to use

The original R/Bioconductor package `CancerSubtypes` integrates the current common computational biology methods for cancer subtypes identification and provides a standardized framework for cancer subtype analysis based multi-omics data, such as gene expression, miRNA expression, DNA methylation and others. This is a forked repo that I added additional features such as random seed (for reproducibility) and dependency updates (iCluster is removed from CRAN and there is an updated iClusterPlus).

------------------------------------------------------------------------

## Installation

```{r,eval=FALSE,warning=FALSE,message=FALSE}
renv::install("ConsensusClusterPlus")
renv::install("hsiaoyi0504/CancerSubtypes")
```
------------------------------------------------------------------------

## Manual
Tutorial and examples can be found in original CancerSubtypes [here](https://bioconductor.org/packages/devel/bioc/vignettes/CancerSubtypes/inst/doc/CancerSubtypes-vignette.html).

<!--(http://htmlpreview.github.io/?https://github.com/taoshengxu/Documents/blob/master/CancerSubtypes-vignette.html)-->


------------------------------------------------------------------------

## Citation
Please cite the original CancerSubtypes article and this customized fork [https://github.com/hsiaoyi0504/CancerSubtypes](https://github.com/hsiaoyi0504/CancerSubtypes) when using this customized version:

[![doi](https://img.shields.io/badge/doi-10.1093/bioinformatics/btx378-green.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx378) [![citation](https://img.shields.io/badge/cited%20by-34-green.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx378) [![Altmetric](https://img.shields.io/badge/Altmetric-2-green.svg?style=flat)](https://www.altmetric.com/details/21038105)

Xu, T. et al. CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation, and visualization. Bioinformatics (2017) [doi:10.1093/bioinformatics/btx378](https://doi.org/10.1093/bioinformatics/btx378).
