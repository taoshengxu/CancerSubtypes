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
#' #Dist1=dist2(t(miRNAExp),t(miRNAExp))
#' #W1 = affinityMatrix(Dist1, 20, 0.5)
#' #group=spectralAlg(W1,3, type = 3)
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

.distanceWeighted2<-function(X,weight)  ##X is the expression Matrix(Row is sample, column is feature) 
{
  if(length(weight)==ncol(X))
  {
    X_row = nrow(X)      
    weight_diag<-diag(weight)
    X2<-(X^2)%*%weight_diag 
    sumsqX = rowSums(X2)
    X1<-X%*%weight_diag
    XY = X1 %*% t(X)
    XX=matrix(rep(sumsqX, times = X_row), X_row, X_row)
    res=XX+t(XX)-2*XY
    res[res < 0] = 0
    diag(res)=0
    #res<-sqrt(res)
    return(res)
  }
  else
  {
    stop("The number of weights is not equal to the number of features in the data matix")
  }
}

## Followings are SNFtool functions, because SNFtool is not  available in CRAN for a long time

#' This is the distance function extracted from SNFtool package.
#' Computes the Euclidean distances between all pairs of data point given.
#' @param X A data matrix where each row is a different data point
#' @param C A data matrix where each row is a different data point. If this matrix is the same as X, pairwise distances for all data points are computed.
#' @examples
#' data(miRNAExp)
#' Dist1=dist2(t(miRNAExp),t(miRNAExp))
#' 
#' @return 
#' Returns an N x M matrix where N is the number of rows in X and M is the number of rows in M. element (n,m) is the squared Euclidean distance between nth data point in X and mth data point in C
#' @export
dist2 <- function(X,C) {
  ndata = nrow(X)
  ncentres = nrow(C)
  
  sumsqX = rowSums(X^2)
  sumsqC = rowSums(C^2)
  
  XC = 2 * (X %*% t(C))
  
  res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
  res[res < 0] = 0
  return(res)
}

#' This is the affinity Matrix function extracted from SNFtool package.
#' Computes affinity matrix from a generic distance matrix
#' @param Diff Distance matrix
#' @param K Number of nearest neighbors
#' @param sigma Variance for local model
#' @examples
#' data(miRNAExp)
#' Dist1=dist2(t(miRNAExp),t(miRNAExp))
#' W1 = affinityMatrix(Dist1, 20, 0.5)
#' 
#' @return 
#' Returns an affinity matrix that represents the neighborhood graph of the data points.
#' @export
#' 
affinityMatrix <- function(Diff,K=20,sigma=0.5) {
  ###This function constructs similarity networks.
  N = nrow(Diff)
  
  Diff = (Diff + t(Diff)) / 2
  diag(Diff) = 0;
  sortedColumns = as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }
  means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
  
  avg <- function(x,y) ((x+y)/2)
  Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
  densities = dnorm(Diff,0,sigma*Sig,log = FALSE)
  
  W = (densities + t(densities)) / 2
  return(W)
}

SNF <- function(Wall,K=20,t=20) {
  
  ###This function is the main function of our software. The inputs are as follows:
  # Wall : List of affinity matrices
  # K : number of neighbors
  # t : number of iterations for fusion
  
  ###The output is a unified similarity graph. It contains both complementary information and common structures from all individual network. 
  ###You can do various applications on this graph, such as clustering(subtyping), classification, prediction.
  
  LW = length(Wall)
  #normalize <- function(X) X / rowSums(X)
  
  #New normalization method
  normalize <- function(X){
    X <- X/(2*(rowSums(X) - diag(X)))
    diag(X) <- 0.5
    return(X)
  }
  # makes elements other than largest K zero
  
  
  newW <- vector("list", LW)
  nextW <- vector("list", LW)
  ###First, normalize different networks to avoid scale problems.
  for( i in 1: LW){
    Wall[[i]] = normalize(Wall[[i]]);
    Wall[[i]] = (Wall[[i]]+t(Wall[[i]]))/2;
  }
  
  ### Calculate the local transition matrix.
  for( i in 1: LW){
    newW[[i]] = (.dominateset(Wall[[i]],K))
  }
  
  # perform the diffusion for t iterations
  for (i in 1:t) {
    for(j in 1:LW){
      sumWJ = matrix(0,dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      for(k in 1:LW){
        if(k != j) {
          sumWJ = sumWJ + Wall[[k]]
        }
      }
      nextW[[j]] = newW[[j]] %*% (sumWJ/(LW-1)) %*% t(newW[[j]]);
    }
    ###Normalize each new obtained networks.
    for(j in 1 : LW){
      #Adding normalization after each iteration
      Wall[[j]] <- normalize(nextW[[j]])
      Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2;
    }
  }
  
  # construct the combined affinity matrix by summing diffused matrices
  W = matrix(0,nrow(Wall[[1]]), ncol(Wall[[1]]))
  for(i in 1:LW){
    W = W + Wall[[i]]
  }
  W = W/LW;
  W = normalize(W);
  # ensure affinity matrix is symmetrical
  
  #Removed addition of diag(W) before centering
  W = (W + t(W)) / 2;
  
  return(W)  
}

spectralClustering <- function(affinity, K, type=3) {
  
  ###This function implements the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
  ###THe inputs are as follows:
  
  #affinity: the similarity matrix;
  #K: the number of clusters
  # type: indicators of variants of spectral clustering 
  
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  
  
  
  return(labels)
}

displayClusters <- function(W, group) {
  normalize <- function(X) X / rowSums(X)
  ind = sort(as.vector(group),index.return=TRUE)
  ind = ind$ix
  diag(W) = 0
  W = normalize(W);
  W = W + t(W);
  image(1:ncol(W),1:nrow(W),W[ind,ind],col = grey(100:0 / 100),xlab = 'Patients',ylab='Patients');
}


.csPrediction <- function(W,Y0,method){
  ###This function implements the label propagation to predict the label(subtype) for new patients.	
  ### note method is an indicator of which semi-supervised method to use
  # method == 0 indicates to use the local and global consistency method
  # method >0 indicates to use label propagation method.
  
  alpha=0.9;
  P= W/rowSums(W)
  if(method==0){
    Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
  } else {
    NLabel=which(rowSums(Y0)==0)[1]-1;
    Y=Y0;
    for (i in 1:1000){
      Y=P%*%Y;
      Y[1:NLabel,]=Y0[1:NLabel,];
    }
  }
  return(Y);
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

.dominateset <- function(xx,KK=20) {
  ###This function outputs the top KK neighbors.	
  
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
    
  }
  
  
  return(normalize(A))
}

# Calculate the mutual information between vectors x and y.
.mutualInformation <- function(x, y) {
  classx <- unique(x)
  classy <- unique(y)
  nx <- length(x)
  ncx <- length(classx)
  ncy <- length(classy)
  
  probxy <- matrix(NA, ncx, ncy)
  for (i in 1:ncx) {
    for (j in 1:ncy) {
      probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx
    }
  }
  
  probx <- matrix(rowSums(probxy), ncx, ncy)
  proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
  result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
  return(result)
}

# Calculate the entropy of vector x.
.entropy <- function(x) {
  class <- unique(x)
  nx <- length(x)
  nc <- length(class)
  
  prob <- rep.int(NA, nc)
  for (i in 1:nc) {
    prob[i] <- sum(x == class[i])/nx
  }
  
  result <- -sum(prob * log(prob, 2))
  return(result)
}

.repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  if (is.null(dim(X))) {
    mx = length(X)
    nx = 1
  } else {
    mx = dim(X)[1]
    nx = dim(X)[2]
  }
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE)
}