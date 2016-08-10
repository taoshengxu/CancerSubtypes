#' Dataset: Gene expression
#'
#'A glioblastoma (GBM) gene expression dataset downloaded from TCGA. This is a small dataset with 1500 genes 
#'and  100 cancer samples extracted from gene expression data for examples.
#' 
#'
#'\itemize{
#'  \item Rows are genes
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name GeneExp
#'@examples
#' data(GeneExp)
NULL


#' Dataset: miRNA expression
#'
#'A glioblastoma (GBM) miRNA expression dataset downloaded from TCGA. This is a small miRNA expression dataset with 470 miRNAs and 
#' 100 cancer samples extracted from miRNA expression data for examples.
#'
#'\itemize{
#'  \item Rows are miRNAs
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name miRNAExp
#'@examples
#' data(miRNAExp)
NULL


#' Dataset: Survival time
#'
#'\itemize{
#'  \item A vector representing the right censored survival time (days) for GBM cancer patients matched with the "GeneExp" and "miRNAExp"
#' datasets.
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A numeric vector
#'@name time
#'@examples
#' data(time)
NULL

#' Dataset: Survival status
#'
#'\itemize{
#'  \item A vector representing the survival status for GBM cancer patients matched with the "GeneExp" and "miRNAExp"
#'. 0=alive or censored, 1=dead
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A numeric vector
#'@name status
#'@examples
#' data(status)
NULL
