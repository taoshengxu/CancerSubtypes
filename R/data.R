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

#' Dataset: A default ranking of features for the fuction ExecuteWSNF()
#' 
#' A dataframe represents the regulatory ranking for features(mRNA,TF,miRNA) caculated based on the miRNA-TF-miRNA regulatory network which was promoted in our published work: 
#' Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data(PLos One,2016).
#'
#'\itemize{
#'  \item  mRNA_TF_miRNA_ID : ENTREZID for genes(mRNA,TF)  and miRBase Accession ID for miRNAs.
#'  \item  mRNA_TF_miRNA.v21._SYMBOL: gene symbol and miRNA names(miRBase Version 21)
#'  \item  feature_ranking:  the numeric values represents regulatory ranking for each feature.
#'}
#'
#'@docType data
#'@keywords datasets
#'@format dataframe
#'@name Ranking
#'@references 
#' Xu, T., Le, T. D., Liu, L., Wang, R., Sun, B., & Li, J. (2016). Identifying cancer subtypes from mirna-tf-mrna regulatory networks and expression data. PloS one, 11(4), e0152792.
#'@examples
#' data(Ranking)
NULL