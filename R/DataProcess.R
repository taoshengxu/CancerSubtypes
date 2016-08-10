#' Data check distribution
#' 
#' @param Data A matrix representing the genomic data such as gene expression data, miRNA expression data.\cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.

#' @examples 
#' data(GeneExp)
#' data.checkDistribution(GeneExp)
#' @return 
#' A plot describes the mean, variance and Median Absolute Deviation (MAD) distribution of features.
#' 
#' @export
data.checkDistribution<-function(Data)
{
  mean_features=apply(Data,1,function(x) mean(x,na.rm = TRUE))
  #which(is.na(mean_features))
  variance_features=apply(Data,1, function(x) var(x,na.rm = TRUE))
  #which(is.na(variance_features))
  MAD_features=apply(Data,1, function(x) mad(x,na.rm = TRUE))
  #which(is.na(MAD_features))
  feature_num=nrow(Data)
  par(mfrow=c(3,1))
  ###1. check the expression distribution by average the absolut value of each feature.
  hist(mean_features, breaks=feature_num*0.1, col="red",
       main="Data (mean) distribution",
       xlab="The average value of features")
  ###2. check the data variance of each feature
  hist(variance_features, breaks=feature_num*0.1, col="red",
       main="Data (Variance) distribution",
       xlab="The variance of features")
  ###3. check the data MAD of each feature
  hist(MAD_features, breaks=feature_num*0.1, col="red",
       main="Data (MAD) distribution",
       xlab="The Median Absolute Deviation of features")
  
  par(mfrow=c(1,1))
}


#'Data normalization
#'
#'Conduct normalization for dataset.
#'
#'@param Data A matrix representing the genomic data such as gene expression data, miRNA expression data. \cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#'@param type A character value representing the normalization type. The optional values are shown below:
#'\itemize{
#' \item "feature_Median". The default value. Normalize dataset by sweeping the median values of each feature.
#' \item "feature_Mean". Normalize dataset by sweeping the mean values of each feature.    
#' \item "feature_zscore". Conduct z_score normalization for each feature.
#' \item "sample_zscore".  Conduct z_score normalization for each samples.
#' 
#'}
#'@param log2 A logical value. If TRUE, the data is transform as log2(x+1). This is commonly used for RNAseq data.
#'@return
#' The normalized data matrix.
#'@examples
#' data(GeneExp)
#' result=data.normalization(GeneExp,type="feature_Median",log2=FALSE)
#' 
#'@export
#'
data.normalization<-function(Data,type="feature_Median",log2=FALSE)
{
  if(log2)
    data=log2(Data+1)
  else
    data=Data
  
  if(type=="feature_Median")
  {
    result=sweep(data,1,apply(data,1,function(x) median(x, na.rm = TRUE)))
  }
  else if(type=="feature_Mean")
  {
    result=sweep(data,1,apply(data,1,function(x) mean(x, na.rm = TRUE)))
  }
  else if(type=="feature_zsocre")
  {
    result=t(scale(t(data)))
  }
  else if(type=="sample_zsocre")
  {
    result=scale(data)
  }
  result
}

#'Data imputation
#'
#'Data imputation for features with missing values
#'@importFrom impute impute.knn
#'@param Data A matrix representing the genomic data such as gene expression data, miRNA expression data. \cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#'@param fun A character value representing the imputation type. The optional values are shown below:
#'\itemize{
#' \item "median". The NAs will be replaced by the median of the existing values of this feature in all samples.
#' \item "mean". The NAs will be replaced by the mean of the existing values of this feature in all samples.
#' \item "microarray". It will apply the "impute" package to impute the missing values. This is a common way to process 
#' the missing observation for MicroArray dataset.
#' 
#'}
#'@return
#' The data matrix after imputation (without NAs).
#'@examples
#'Data=matrix(runif(1000),nrow = 50,ncol = 20)
#'geneName=paste("Gene", 1:50, sep = " ")
#'sampleName=paste("Sample", 1:20, sep = " ")
#'rownames(Data)=geneName
#'colnames(Data)=sampleName
#'index=sample(c(1:1000),60)
#'Data[index]=NA
#'result=data.imputation(Data,fun="median")
#'@export
#'
data.imputation<-function(Data,fun="median")
{
  if(fun=="median")
  {
    result=apply(Data,1,function(x){
      x<-as.numeric(x)
      x[is.na(x)] =median(x, na.rm=TRUE) 
      x
    })
    result=t(result)
  }
  else if(fun=="mean")
  {
    result=apply(Data,1,function(x){
      x<-as.numeric(x)
      x[is.na(x)] =mean(x, na.rm=TRUE)
      x
    })
    result=t(result)
  }
  else if(fun=="microarray")
  {
    result=impute.knn(Data)$data
  }
  result
}
  
  
  