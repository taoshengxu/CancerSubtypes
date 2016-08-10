#' This is an internal function but need to be exported for the function ExecuteSNF.CC() call.
#' 
#' This is  Spectral Clustering Algorithm extracted from SNFtools package spectralClustering() with
#' a tiny modification. 
#' @param affinity Similarity matrix
#' @param K Number of clusters
#' @param type The variants of spectral clustering to use.
#' @examples
#' ####see the spectralClustering() in SNFtool package for the detail example.
#' data(miRNAExp)
#' Dist1=SNFtool::dist2(t(miRNAExp),t(miRNAExp))
#' W1 = SNFtool::affinityMatrix(Dist1, 20, 0.5)
#' group=spectralAlg(W1,3, type = 3)
#' 
#' @return 
#' A vector consisting of cluster labels of each sample.
#' @export
spectralAlg <- function (affinity, K, type = 3) 
{
  affinity=as.matrix(affinity)
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  }
  else if (type == 2) {
    Di = diag(1/d)
    NL = Di %*% L
  }
  else if (type == 3) {
    Di = diag(1/sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values), index.return = TRUE)
  U = eig$vectors[, res$ix[1:K]]
  normalize <- function(x) x/sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U, 1, normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete, 1, which.max)
  return(labels)
}

.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  return(Y)
}