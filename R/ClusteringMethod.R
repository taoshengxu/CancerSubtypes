#' Execute Consensus Clustering
#'
#' This function is based on the R package "ConsensusClusterPlus". 
#' We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The parameters are compatible to the original R package "ConsensusClusterPlus" function "ConsensusClusterPlus()".\cr
#' Please note: we add a new parameter "clusterNum" which represents the result with cancer subtypes group we want to return. 
#' 
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @param clusterNum A integer representing the return cluster number, this value should be less
#' than maxClusterNum(maxK). This is the only additional parameter in our function compared to the original
#' R package "ConsensusClusterPlus". All the other parameters are compatible to the function "ConsensusClusterPlus().
#' @param d data to be clustered; either a data matrix where columns=items/samples and rows are features. For example, a gene expression matrix of genes in rows and microarrays in columns, or ExpressionSet object, or a distance object (only for cases of no feature resampling)  
#' 
#' Please Note: We add a new data type (list) for this parameter. Please see details and examples.
#' @param maxK integer value. maximum cluster number  for Consensus Clustering Algorithm to evaluate.
#' @param reps  integer value. number of subsamples(in other words, The iteration number of each cluster number) 
#' @param clusterAlg character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
#' @param distance character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param  title character value for output directory. This title can be an absolute or relative path
#' @param pItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param pFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param plot Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param innerLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param finalLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param writeTable Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param verbose Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param corUse Please refer to the "ConsensusClusterPlus" package for detailed information.

#' @return A list with the following elements.
#'\itemize{
#'  \item group 
#'  
#'   A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.
#'  \item distanceMatrix
#'  
#'   It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'  \item originalResult
#'  
#'  The clustering result of the original function "ConsensusClusterPlus()"
#'  }
#'  
#' @details 
#'  If the data is a list containing the matched mutli-genomics  data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   we use "z-score" to normalize features for each data matrix first. Then all the normalized data matrices from the data list are concatenated
#'   according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.
#'   
#' @seealso \code{ConsensusClusterPlus}
#'   
#' @examples
#' ### The input dataset is a single gene expression matrix.
#' data(GeneExp)
#' data(miRNAExp)
#' result1=ExecuteCC(clusterNum=3,d=GeneExp,maxK=10,clusterAlg="hc",distance="pearson",title="GBM")
#' result1$group
#' 
#' ### The input dataset is multi-genomics data as a list
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result2=ExecuteCC(clusterNum=3,d=GBM,maxK=5,clusterAlg="hc",distance="pearson",title="GBM")
#' result2$group
#' 
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., Golub, T. (2003) Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Machine Learning, 52, 91-118.
#' @export
#'
ExecuteCC<-function(clusterNum,
                    d,maxK=10,clusterAlg="hc",
                    distance="pearson",title="ConsensusClusterResult",
                    reps=500, pItem=0.8, pFeature=1,plot="png",
                    innerLinkage="average", finalLinkage="average",
                    writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,
                    verbose=FALSE,corUse="everything")
{
  if(is.list(d))
  {
    temp=NULL
    for(i in 1: length(d))
    {
      temp=rbind(temp,d[[i]])
    }
    temp=t(scale(t(temp)))
  }
  else
   temp=d
  
  originalResult=ConsensusClusterPlus(
      temp, maxK=maxK,clusterAlg=clusterAlg,
      distance=distance,title=title,
      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
      innerLinkage=innerLinkage, finalLinkage=finalLinkage,
      writeTable=writeTable,weightsItem=weightsItem,weightsFeature=weightsFeature,
      verbose=verbose,corUse=corUse)
  
  group=originalResult[[clusterNum]][["consensusClass"]]
  
  distanceMatrix=originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix,'class')="Similarity"
    
  #icl=calcICL(result,title =fileName,plot="png" )
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=originalResult)
  result
}



#' Execute iCluster (Integrative clustering of multiple genomic data)
#'
#' Shen (2009) proposed a latent variable regression with a lasso constraint for joint modeling of multiple omics 
#' data types to identify common latent variables that can be used to cluster patient samples into biologically and clinically relevant disease subtypes.\cr
#' This function is based on the R package "iCluster". 
#' The R package "iCluster" should be installed. 
#' We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The parameters is compatible to the original R package "iCluster" function "iCluster2()".\cr
#' Please note: The data matrices are transposed in our function comparing to the original R package "iCluster" on the behalf of the unified input format with other functions.
#' We try to build a standardized flow for cancer subtypes analysis and validation.
#'
#' @importFrom iCluster iCluster2 plotiCluster
#' @param datasets A list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' In order to unify the input parameter with other clustering methods, the data matrices are transposed comparing to the definition in the original "iCluster" package.
#' @param k Number of subtypes for the samples
#' @param lambda Penalty term for the coefficient matrix of the iCluster model
#' @param scalar Logical value. If true, a degenerate version assuming scalar covariance matrix is used.
#' @param max.iter  maximum iteration for the EM algorithm
#' @param scale Logical value. If true, the genomic features in the matrix is centered.

#' @return A list with the following elements.
#'\itemize{
#'  \item group  
#'  
#'  A vector represent the group of cancer subtypes. The order is corresponding to the samples in the data matrix.
#'  \item originalResult
#'  
#'  The clustering result of the original function "iCluster2()"
#'  }
#'  
#'  
#' @details 
#'  For iCluster algorithm, it cannot process high-dimensional data, otherwise it is very very time-consuming or reports a mistake.  
#'  Based on test, it could smoothly run for the matrix with around 1500 features. Normally it need feature selection step first to reduce feature number.
#' @references
#' Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. Bioinformatics 25, 2906-2912.\cr
#' Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. PLoS ONE 7, e35236
#' 
#' @seealso \code{\link{iCluster2}}
#' 
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' data1=FSbyVar(GeneExp, cut.type="topk",value=500)
#' data2=FSbyVar(miRNAExp, cut.type="topk",value=100)
#' GBM=list(GeneExp=data1,miRNAExp=data2)
#' result=ExecuteiCluster(datasets=GBM, k=3, lambda=list(0.44,0.33,0.28))
#' result$group
#' @export
#'
ExecuteiCluster<-function(datasets, k, lambda=NULL, scale=TRUE, scalar=FALSE, max.iter=10)
{
  data1=list()
  for(i in 1:length(datasets))
  {
    data1[[i]]=t(datasets[[i]])
  }
  
  fit=iCluster2(datasets=data1, k=k, lambda=lambda, scale=scale, scalar=scalar, max.iter=10) 
  
  plotiCluster(fit=fit, label=rownames(data1[[1]]))
  group=fit$clusters
  result=list(group=group,originalResult=fit)
  result
}


#' Execute SNF(Similarity Network Fusion )
#'
#' SNF is a multi-omics data processing method that constructs a fusion patient similarity network
#' by integrating the patient similarity obtained from each of the genomic data types. 
#' SNF calculates the similarity between patients using each single data type separately. The similarities
#' between patients from different data types are then integrated by a cross-network diffusion process to construct the fusion patient similarity matrix. 
#' Finally, a clustering method is applied to the fusion patient similarity matrix to cluster patients into different groups, which imply different cancer subtypes.
#' This function is based on the R package "SNFtool". 
#' The R package "SNFtool" should be installed. 
#' We write a function to integrate the clustering process and unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation.\cr
#' Please note: The data matrices are transposed in our function comparing to the original R package "SNFtools".
#' We try to build a standardized flow for cancer subtypes analysis and validation.
#'
#' @importFrom SNFtool dist2 affinityMatrix SNF spectralClustering displayClusters
#' @param datasets A list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' @param clusterNum A integer representing the return cluster number
#' @param K Number of nearest neighbors
#' @param alpha Variance for local model
#' @param t Number of iterations for the diffusion process
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#' @return A list with the following elements.
#'\itemize{
#'  \item group
#'  
#'      A vector represents the group of cancer subtypes. The order is corresponding to the samples in the data matrix.
#'  \item distanceMatrix  
#'  
#'  It is a samplesimilarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'  }
#' @references
#' B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014
#' 
#' @seealso \code{\link{affinityMatrix}} \code{\link{SNF}}
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' result$group
#' @export
#'
ExecuteSNF<-function(datasets, clusterNum, K=20, alpha=0.5, t=20,plot=TRUE)
{
 # library("SNFtool")
  W_temp=list()
  for(i in 1:length(datasets))
  {
    distance=dist2(as.matrix(t(datasets[[i]])), as.matrix(t(datasets[[i]])))
    W_temp[[i]] = affinityMatrix(distance, K, alpha)
  }
  W = SNFtool::SNF(W_temp, K=K, t=t)
  group =spectralClustering(W,clusterNum)
  
  diag(W)=0
  diag(W)=max(W)
  distanceMatrix=W
  attr(distanceMatrix,'class')="Similarity"
  
  if(plot)
    displayClusters(W, group)
  result=list(group=group,distanceMatrix=distanceMatrix)
  result
}


#' Execute the combined SNF (Similarity Network Fusion) and Consensus clustering
#'
#' This function is a combined process of SNF and Consensus Clustering for cancer subtypes identification.
#' First it applied SNF to get the fusion patients similarity matrix. Then use this 
#' fusion patients similarity matrix as the sample distance for Consensus Clustering.
#'
#' @importFrom SNFtool dist2 affinityMatrix SNF 
#' @param datasets Same as ExecuteSNF
#' @param clusterNum Same as ExecuteSNF
#' @param K Same as ExecuteSNF
#' @param alpha Same as ExecuteSNF
#' @param t Same as ExecuteSNF
#' @param maxK Same as ExecuteCC
#' @param pItem Same as ExecuteCC
#' @param reps Same as ExecuteCC
#' @param title Same as ExecuteCC
#' @param plot Same as ExecuteCC
#' @param finalLinkage Same as ExecuteCC
#' 
#' @return Same as the ExecuteCC(). A list with the following elements.
#'\itemize{
#'  \item group
#'  
#'   A vector represent the group of cancer subtypes. The order is corresponding to the samples in the data matrix.
#'  \item distanceMatrix
#'  
#'   It is a sample dissimilarity matrix. The more large value between samples in the matrix, the more dissimilarity the samples are.
#'  \item originalResult
#'  
#'  The clustering result of the original function "ConsensusClusterPlus()"
#'  }
#'  
#' @seealso \code{\link{ExecuteSNF}} \code{\link{ExecuteCC}}
#' @examples
#' 
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp,miRNAExp)
#' result=ExecuteSNF.CC(GBM, clusterNum=3, K=20, alpha=0.5, t=20,
#'                     maxK = 5, pItem = 0.8,reps=500, 
#'                     title = "GBM", plot = "png", 
#'                     finalLinkage ="average")
#' result$group
#' @export
#'
ExecuteSNF.CC<-function(datasets, clusterNum, K=20, alpha=0.5, t=20,
                        maxK = 10, pItem = 0.8,reps=500, 
                        title = "ConsensusClusterResult", plot = "png", 
                        finalLinkage ="average")
{
  W_temp=list()
  for(i in 1:length(datasets))
  {
    distance=dist2(as.matrix(t(datasets[[i]])), as.matrix(t(datasets[[i]])))
    W_temp[[i]] = affinityMatrix(distance, K, alpha)
  }
  W = SNFtool::SNF(W_temp, K=K, t=t)
  W=as.dist(W)
  result=ExecuteCC(clusterNum=clusterNum,d=W,maxK=maxK,
                   clusterAlg="spectralAlg",title=title,reps=reps,
                   pItem=pItem,plot=plot,
                   finalLinkage=finalLinkage)
  result
}

#' Execute Consensus NMF (Nonnegative matrix factorization)
#'
#' Brunet applied nonnegative matrix factorization (NMF) to analyze the Gene MicroArray dataset in 2004. In the original paper, the author
#' proved that NMF is an efficient method for distinct molecular patterns identification and provides a powerful method 
#' for class discovery. This method was implemented in an R package "NMF". Here we applied the "NMF" package to 
#' conduct the cancer subtypes identification. We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The R package "NMF" should be installed. 
#' 
#' @importFrom NMF nmf predict
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' If the matrices have negative values, first the negative values will be set to zero to get a matrix 1;
#' all the positive values will be set to zero to get the matrix 2; then a new matrix with all positive values will be
#' get by concatenating matrix1 and -maxtrix2.
#'  
#' 
#' @param clusterNum Number of subtypes for the samples
#' @param nrun Number of runs to perform NMF. A default of 30 runs are performed, allowing the computation of a consensus matrix that is used in selecting the best result for cancer subtypes identification as Consensus Clustering method.
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item group  
#'  
#'  A vector represent the group of cancer subtypes. The order is corresponding to the samples in the data matrix.
#'  \item distanceMatrix
#'  
#'   It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'  \item originalResult

#'  A NMFfitX class from the result of function "nmf()".
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' [1] Brunet, Jean-Philippe, Pablo Tamayo, Todd R Golub, and Jill P Mesirov. "Metagenes and Molecular Pattern Discovery Using Matrix Factorization." Proceedings of the National Academy of Sciences 101, no. 12 (2004):4164-69.
#' 
#' [2] Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for Nonnegative Matrix Factorization." BMC Bioinformatics 11 (2010): 367. doi:10.1186/1471-2105-11-367.
#' @seealso \code{\link{nmf}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' #To save the  execution time, the nrun is set to 5, but the recommended value is 30.
#' result=ExecuteCNMF(GeneExp,clusterNum=3,nrun=5)
#' result$group
#' 
#' @export
#'
ExecuteCNMF<-function(datasets, clusterNum,nrun=30 )
{
  if(is.list(datasets))
  {
    temp=NULL
    for(i in 1: length(datasets))
    {
      temp=rbind(temp,datasets[[i]])
    }
  }
  else
    temp=datasets
  ## change all value to positive
  data1=rbind(pmax(temp,0),-pmin(temp,0))
  index=which(rowSums(data1)==0)
  data1=data1[-index,]
  res=nmf(data1,rank=clusterNum,nrun=nrun)
  
  distanceMatrix=slot(res,"consensus")
  attr(distanceMatrix,'class')="Similarity"
  
  group=as.numeric(as.vector(predict(res)))
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=res)
}


