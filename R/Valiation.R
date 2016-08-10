
#'Compute or Extract Silhouette Information from Clustering based on similarity matrix.
#'
#'Silhouette refers to a method of interpretation and validation of consistency within clusters of data. 
#'The technique provides a succinct graphical representation of how well each object lies within its cluster (From Wiki).\cr
#'Note that: This function is a rewriting version of the function "silhouette()" in R package cluster.
#'   The original function "silhouette()" is to compute the silhouette information based on a dissimilarity matrix.
#'   Here the silhouette_SimilarityMatrix() is to solve the computation based on the similarity matrix.
#'   The result of the silhouette_SimilarityMatrix() is compatible to the function "Silhouette()".


#'
#'@param group A vector represent the cluster label for a set of samples.
#'@param similarity_matrix A similarity matrix between samples
#'@return 
#' An object, sil, of class silhouette which is an [n x 3] matrix with
#' attributes. The colnames correspondingly are c("cluster", "neighbor", "sil_width").
#' @details
#' For each observation i, the return sil[i,] contains the cluster to which i belongs as well as the neighbor 
#' cluster of i (the cluster, not containing i, for which the average 
#' dissimilarity between its observations and i is minimal), 
#' and the silhouette width s(i) of the observation. 
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
#' plot(sil)
#' ###If use the silhouette(), the result is wrong because the input is a similarity matrix.
#' sil1=silhouette(result$group, result$distanceMatrix)
#' plot(sil1)  ##wrong result
#' 
#' @seealso \code{\link{silhouette}}
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'  
#'@references
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. J. Comput. Appl. Math., 20, 53-65.
#'@export
#'
silhouette_SimilarityMatrix<-function(group, similarity_matrix)
{
  similarity_matrix=as.matrix(similarity_matrix)
  similarity_matrix<-(similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix)=0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix<-normalize(similarity_matrix)
  
  n <- length(group)
  if(!all(group == round(group))) stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if(k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if(doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  
  wds <- matrix(NA, n,3, dimnames =list(names(group), c("cluster","neighbor","sil_width")))  
  for(j in 1:k)
  { 
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, group[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else 0
    wds[index,"sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}


#' Survival analysis(Survival curves, Log-rank test) and compute Silhouette information for cancer subtypes
#' 
#' Survival analysis is a very common tool to explain and validate the cancer subtype identification result. It provides the significance testing and 
#' graphical display for the verification of the survival patterns between the identified cancer subtypes.
#' 
#' @importFrom survival survfit survdiff
#' @importFrom NMF consensusmap
#' @importFrom cluster silhouette
#'
#' @param mainTitle A character will display in the result plot.
#' @param time A numeric vector representing the survival time (days) of a set of samples.
#' @param status A numeric vector representing the survival status of a set of samples. 0=alive/censored, 1=dead.
#' @param group A vector represent the cluster label for a set of samples.
#' @param distanceMatrix A data matrix represents the similarity matrix or dissimilarity matrix between samples.\cr
#' If NULL, it will not compute silhouette width and draw the plot.
#' @param similarity A logical value. If TRUE, the distanceMatrix is a similarity distance matrix between samples. Otherwise a dissimilarity distance matrix between samples
#'
#' @return
#' The log-rank test p-value
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'@examples
#'### SNF result analysis
#' data(GeneExp)
#' data(miRNAExp)
#' data(time)
#' data(status)
#' data1=FSbyCox(GeneExp,time,status,cutoff=0.05)
#' data2=FSbyCox(miRNAExp,time,status,cutoff=0.05)
#' GBM=list(GeneExp=data1,miRNAExp=data2)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' group=result$group
#' distanceMatrix=result$distanceMatrix
#' p_value=survAnalysis(mainTitle="GBM",time,status,group,
#'                      distanceMatrix=distanceMatrix,similarity=TRUE)
#' @export
#'
survAnalysis<-function(mainTitle="Survival Analysis",time,status,group,distanceMatrix=NULL,similarity=TRUE)
{
  
  clusterNum=length(unique(group))
  dataset=list(time,status,x=group)  
  surv=survfit(Surv(time, status) ~ x,dataset)
  if(clusterNum>1)
  {
    sdf=NULL
    sdf=survdiff(Surv(time, status) ~ group)##log-rank test
    cat("                                                     \n")
    cat("*****************************************************\n")
    cat(paste(mainTitle,"Cluster=",clusterNum,"  "))
    print(sdf)
    p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  }
  else
  {
    cat("There is only one cluster in the group")
    p_value=1
  }
  
  if(!is.null(distanceMatrix[1,1]))
  {
    layout(matrix(c(1,2,3,3), 2, 2, byrow = FALSE),widths=c(2.2,2), heights=c(2,2))
  }
  title=paste(mainTitle,"Cluster =",clusterNum)
  plot(surv, lty = 1,col=2:(clusterNum+1),lwd=2,xscale=30,xlab="Survival time (Months)", ylab="Survival probability",
       main= title,font.main=4,cex.main=0.9)
  
  legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=0.3, paste("Subtpye", 1:clusterNum),
         lty=1,lwd=3, cex=0.8,text.font=4,text.col=2:(clusterNum+1), bty="n",col=2:(clusterNum+1),
         seg.len = 0.3)
  
  digit=ceiling(-log10(p_value)+2)
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.9,paste("p-value=",round(p_value,digit)),col="red",font=3,cex=1)   
  
  if(!is.null(distanceMatrix[1,1]))
  {
    #####
    if(class(distanceMatrix)=="Similarity")
    {
      si=silhouette_SimilarityMatrix(group,distanceMatrix)
    }
    else
    {
      si=silhouette(group,distanceMatrix)
    }
     
    attr(distanceMatrix,'class')=NULL
    
    ind=order(group,-si[, "sil_width"])
    
    num=length(unique(group))
    annotation=data.frame(group=as.factor(group))
    Var1 = c(palette()[2:(num+1)])
    names(Var1) = sort(unique(group))
    ann_colors =  list(group=Var1)
    
    consensusmap(distanceMatrix,Rowv=ind,Colv=ind,main = "Clustering dispaly",
                 annCol = annotation,annColors=ann_colors,
                 labRow ="Sample", labCol = "Sample",scale="none")
    
    plot(si,col =2:(clusterNum+1))
  }
  par(mfrow=c(1,1))
  p_value
}


#' Generate heatmaps
#' 
#' Generate heatmap for datasets.
#' @importFrom NMF aheatmap
#' @importFrom grDevices bitmap colorRampPalette dev.copy dev.off palette pdf png postscript rainbow rgb
#' @param data A matrix representing the genomic data such as gene expression data, miRNA expression data.\cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#' @param group A vector representing the subtype on each sample.  The default is NULL. If it is not
#' NULL, the samples will be rearrangement according to the subtypes in the heatmap.
#' @param silhouette An object of class silhouette. It is a result from function 
#' silhouette() or silhouette_SimilarityMatrix(). The default is NULL. If it is not NULL,  an annotation will be drawn to show the silhouette width for each sample.
#' @param scale A string for data normalization type before heatmap drawing.
#' The optional values are shown below:
#' \itemize{
#' \item "no". No normalization. This is default.
#' \item "z_score". Normalize data by z_score of features. 
#' \item "max_min". Normalize each feature by (value-min)/(max-min).
#' }
#' @param labRow labels for the rows. Possible values are:
#' \itemize{
#' \item NULL. The default value. It will use the row names of the matrix for the heatmap labels.
#' \item NA. No row label will be shown. 
#' \item A list of labels.
#' }
#' @param labCol labels for the columns.  See labRow.
#' @param color color specification for the heatmap.
#' @param Title A string for the Main title of the heatmap.
#' @details 
#' We applied the R package "NMF" function "aheatmap()" as the heatmap drawer.
#' @author
#' Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#' @references 
#' Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 1.
#' @return 
#' A heatmap
#' @examples 
#' ### SNF result analysis
#' data(GeneExp)
#' data(miRNAExp)
#' data(time)
#' data(status)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' group=result$group
#' distanceMatrix=result$distanceMatrix
#' silhouette=silhouette_SimilarityMatrix(group, distanceMatrix)
#' drawHeatmap(GeneExp,group,silhouette=silhouette,scale="max_min",Title="GBM Gene Expression")
#' drawHeatmap(GeneExp,group,silhouette=silhouette,scale="max_min",
#'             color="-RdYlBu",Title="GBM Gene Expression")
#' @export
#' 
drawHeatmap<-function(data,group=NULL,silhouette=NULL,scale="no",
                      labRow = NULL, labCol = NULL,
                      color=colorRampPalette(c("green","black","red"))(300),
                      Title=NA)
{
  if(!is.matrix(data))
    stop("The input dataset is not a matrix")
  if(scale=="z_score")
  {
    data1=data.normalization(data, type = "feature_zsocre")
  }
  else if(scale=="max_min")
  {
    data1=t(apply(data,1,function(x) (x-min(x))/(max(x)-min(x))))
  }
  else
  {
    data1=data
  }
  
  num=length(unique(group))
  
  if(!is.null(group))
  {
    if(!is.null(silhouette))
    {
      ind=order(group,-silhouette[, "sil_width"])
      annotation=data.frame(gruop=as.factor(group),silhouette=silhouette[, "sil_width"])
      
      Var1 = c(palette()[2:(num+1)])
      names(Var1) = sort(unique(group))
      Var2 = c("deeppink", "yellow")
      
      ann_colors = list(Var1, Var2)
    }
    else
    {
      ind=order(group)
      annotation=annotation=data.frame(group=as.factor(group))
      
      Var1 = c(palette()[2:(num+1)])
      names(Var1) = sort(unique(group))
      ann_colors =  list(group=Var1)
    }
      
  }
  else
  {
    ind=NA
    annotation=NA
    ann_colors=NA
  }
    
  if(is.na(Title))
    Title="Heatmap"
  
  aheatmap(data1,Colv = ind,treeheight=0,labRow = labRow, labCol = labCol,color = color,
           annCol = annotation,annColors=ann_colors,main=Title)
  
  
}

#' This function save the figure in the current plot.
#' @param foldername Character values. It specifies the folder name which will be created in the present working path.
#' @param filename Character values. It specifies the saved file name.
#' @param image_width the figure width
#' @param image_height the figure height
#' @param image_res the figure resolution
#' 
#' @return A * .png file in the specified folder.
#' @author
#' Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}

#' @examples 
#' data(GeneExp)
#' data(miRNAExp)
#' data(time)
#' data(status)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' group=result$group
#' distanceMatrix=result$distanceMatrix
#' p_value=survAnalysis(mainTitle="GBM",time,status,group,
#'       distanceMatrix=distanceMatrix,similarity=TRUE)
#' saveFigure(foldername="GBM",filename="GBM",image_width=10,image_height=10,image_res=300)
#' @export
saveFigure<-function(foldername=NULL,filename="saveFig",image_width=10,image_height=10,image_res=300)
{
  mainDir<-getwd()
  filename1=paste(filename,".png",sep="")
  if(!is.null(foldername))
  {
    subDir=paste(foldername,"Figures")
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    path=file.path(mainDir, subDir,filename1)
  }
  else
    path=file.path(mainDir,filename1)  
  dev.copy(png,path,width=image_width,height=image_height,res=image_res,units="in")
  dev.off();
}

#' 
#'
#'
#' A statistical method for testing the significance of clustering results. 
#' 
#' SigClust (Statistical significance of clustering) is a statistical method for testing the significance of clustering results. SigClust can be applied to
#' assess the statistical significance of splitting a data set into two clusters. 
#' SigClust studies whether clusters are really there, using the 2-means (k = 2) clustering index as a statistic. It assesses the significance of clustering by 
#' simulation from a single null Gaussian distribution. Null Gaussian parameters are estimated from the data.
#' Here we apply the SigClust to assess the statistical significance of pairwise subtypes. "sigclust" package should be installed.
#' @importFrom sigclust sigclust
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param group The subtypes label of each sample
#' @param nsim This is a parameter inherited from sigclust() in "sigclust" Package. 
#' Number of simulated Gaussian samples to estimate the distribution of the clustering index for the main p-value computation.
#' @param nrep This is a parameter inherited from sigclust() in "sigclust" Package. 
#' Number of steps to use in 2-means clustering computations (default=1, chosen to optimize speed).
#' @param icovest This is a parameter inherited from sigclust() in "sigclust" Package.
#' Covariance estimation type: 1. Use a soft threshold method as constrained MLE (default); 
#' 2. Use sample covariance estimate (recommended when diagnostics fail); 3. Use original background noise threshold estimate (from Liu, et al, (2008)) ("hard thresholding").
#' @return 
#' A matrix indicates the p-value between pairwise subtypes.
#' @references 
#' Liu, Yufeng, Hayes, David Neil, Nobel, Andrew and Marron, J. S, 2008, Statistical Significance of Clustering for High-Dimension, Low-Sample Size Data, Journal of the American Statistical Association 103(483) 1281-1293.\cr
#' Huang, Hanwen, Yufeng Liu, Ming Yuan, and J. S. Marron. "Statistical Significance of Clustering Using Soft Thresholding." Journal of Computational and Graphical Statistics, no. just-accepted (2014): 00-00.
#' @seealso \code{\link{sigclust}}
#' @author
#' Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}

#' @examples 
#' data(GeneExp)
#' data(miRNAExp)
#' data(time)
#' data(status)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' group=result$group
#' sigclust1=sigclustTest(miRNAExp,group, nsim=500, nrep=1, icovest=3)
#' sigclust2=sigclustTest(miRNAExp,group, nsim=1000, nrep=1, icovest=1)
#' @export
sigclustTest<-function(Data,group, nsim=1000, nrep=1, icovest=1)
{
  groupN=sort(unique(group))
  len=length(groupN)
  pvalue=matrix(data = NA, nrow =len , ncol = len)
  name=paste("Subtype",groupN)
  rownames(pvalue)=name
  colnames(pvalue)=name
  
  group_temp=sort(group,index.return = TRUE)
  data=t(Data[,group_temp$ix])
  group_temp=group_temp$x
  
  for(i in 1: (len-1))
  {
    for(j in (i+1): len)
    {
      index=which(group_temp==groupN[i] | group_temp==groupN[j])
      data1=data[index,]
      label=group_temp[index]
      label[which(label == min(label))]=1
      label[which(label == max(label))]=2 
      
      res<-sigclust(x=data1, nsim=nsim, nrep=nrep, labflag=1, label=label, icovest=icovest)
      
      plot(res,arg="pvalue",sub=paste("subtype",groupN[i],"and","subtype",groupN[j]))
      pvalue[i,j]=res@pval
      pvalue[j,i]=res@pval
    }
  }
  diag(pvalue)=1
  pvalue
}


#' DiffExp.limma
#' 
#' Differently Expression Analysis for genomic data. We apply limma package to conduct the analysis.
#' @importFrom limma lmFit eBayes makeContrasts contrasts.fit voom topTable
#' @param Tumor_Data A matrix representing the genomic data of cancer samples such as gene expression data, miRNA expression data.\cr
#' For the matrix, the rows represent the genomic features, and the columns represent the cancer samples.
#' @param Normal_Data A matrix representing the genomic data of Normal samples.\cr
#' For the matrix, the rows represent the genomic features corresponding to the Tumor_Data, and the columns represent the normal samples.
#' @param group A vector representing the subtype of each tumor sample in the Tumor_Data. The length of group is equal to the column number of Tumor_Data.  
#' @param topk The top number of different expression features that we want to extract in the return result.
#' @param RNAseq A bool type representing the input datatype is a RNASeq or not. Default is FALSE for microarray data.
#' @return 
#' A list representing the differently expression for each subtype comparing to the Normal group.
#' @examples 
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' group=result$group
#' ######Fabricate a normal group by extracting some samples from the cancer dataset 
#' ######for demonstrating the examples.
#' Normal_Data=GeneExp[,sample(1:100,20)]
#' result=DiffExp.limma(Tumor_Data=GeneExp,Normal_Data=Normal_Data,group=group,topk=NULL,RNAseq=FALSE)

#' @author
#' Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#' @references 
#' Smyth, Gordon K. "Limma: linear models for microarray data." Bioinformatics and computational biology solutions using R and Bioconductor. 
#' Springer New York, 2005. 397-420.
#' @export

DiffExp.limma<-function(Tumor_Data,Normal_Data,group=NULL,topk=NULL,RNAseq=FALSE)
{ 
  if(is.null(group))
    group=rep(1,ncol(Tumor_Data))
  groupN=sort(unique(group))
  len=length(groupN)
  
  mylist.names <- paste("Subtype",order(groupN))
  result <- vector("list", length(mylist.names))
  names(result) <- mylist.names
  
  num=ncol(Normal_Data)
  for(i in 1:len)
  {
    index=which(group==groupN[i])
    #label=c(rep(0,num),rep(1,length(index)))
    #label=as.factor(label)
    #design <- model.matrix(~-1+label)
    #colnames(design)=c("Normal","Cancer")
    
    Normal=NULL
    Cancer=NULL
    design=cbind(Normal=c(rep(1,num), rep(0,length(index))), Cancer=c(rep(0,num), rep(1,length(index))))
    
    Data=cbind(Normal_Data,Tumor_Data[,index])
    if(RNAseq)
      mR <- voom(Data, design, plot=TRUE)
    else
      mR=Data
    
    ###In order to return the index of features, set the same name for two feautre
    ### then restore the name
    name1=rownames(mR)[1]
    name2=rownames(mR)[2]
    rownames(mR)[1]="repeat"
    rownames(mR)[2]="repeat" 
    ########## mR ############
    mRfit=lmFit(mR, design)
    mRfit=eBayes(mRfit)
    contrast.matrix=makeContrasts(CancervNormal=Cancer - Normal, levels=design)
    mRfit2=contrasts.fit(mRfit, contrast.matrix)
    mRfit2=eBayes(mRfit2)
    
    if(is.null(topk))
    {
      topk=nrow(mR)
    }
    mRresults=topTable(mRfit2, number= topk, sort.by="p", adjust.method="BH")
    ######restore feature name
    mRresults[which(rownames(mRresults)=="1"),1]=name1
    mRresults[which(rownames(mRresults)=="2"),1]=name2
    
    result[[i]]=mRresults
  }
  result
}

