## ----message = FALSE,warning=FALSE---------------------------------------
### Prepare a TCGA gene expression dataset for analysis. 
library(CancerSubtypes)
library("RTCGA.mRNA")
rm(list = ls())
data(BRCA.mRNA)
mRNA=t(as.matrix(BRCA.mRNA[,-1]))
colnames(mRNA)=BRCA.mRNA[,1]

###To observe the mean, variance and Median Absolute Deviation distribution of the dataset, it helps users to get the distribution characteristics of the data, e.g. To evaluate whether the dataset fits a normal distribution or not.
data.checkDistribution(mRNA)

## ----eval = FALSE--------------------------------------------------------
#  index=which(is.na(mRNA))
#  res1=data.imputation(mRNA,fun="median")
#  res2=data.imputation(mRNA,fun="mean")
#  res3=data.imputation(mRNA,fun="microarray")

## ------------------------------------------------------------------------
result1=data.normalization(mRNA,type="feature_Median",log2=FALSE)
result2=data.normalization(mRNA,type="feature_zscore",log2=FALSE)

## ----eval = FALSE--------------------------------------------------------
#  ###The top 1000 most variance features will be selected.
#  data1=FSbyVar(mRNA, cut.type="topk",value=1000)
#  ###The features with (variance>0.5) are selected.
#  data2=FSbyVar(mRNA, cut.type="cutoff",value=0.5)

## ----eval = FALSE--------------------------------------------------------
#  data1=FSbyMAD(mRNA, cut.type="topk",value=1000)
#  data2=FSbyMAD(mRNA, cut.type="cutoff",value=0.5)

## ----eval = FALSE--------------------------------------------------------
#  mRNA1=data.imputation(mRNA,fun="microarray")
#  data1=FSbyPCA(mRNA1, PC_percent=0.9,scale = TRUE)

## ------------------------------------------------------------------------
data(GeneExp)
data(time)
data(status)
data1=FSbyCox(GeneExp,time,status,cutoff=0.05)

## ----eval = FALSE--------------------------------------------------------
#  ### The input dataset is single gene expression matrix.
#  data(GeneExp)
#  result=ExecuteCC(clusterNum=3,d=GeneExp,maxK=10,clusterAlg="hc",distance="pearson",title="GBM")
#  
#  ### The input dataset is multi-genomics data as a list
#  data(GeneExp)
#  data(miRNAExp)
#  GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#  result=ExecuteCC(clusterNum=3,d=GBM,maxK=10,clusterAlg="hc",distance="pearson",title="GBM")

## ----eval = FALSE--------------------------------------------------------
#  ### The input dataset is single gene expression matrix.
#  data(GeneExp)
#  result=ExecuteCNMF(GeneExp,clusterNum=3,nrun=30)
#  
#  ### The input dataset is multi-genomics data as a list
#  data(GeneExp)
#  data(miRNAExp)
#  GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#  result=ExecuteCNMF(GBM,clusterNum=3,nrun=30)

## ----eval = FALSE--------------------------------------------------------
#  data(GeneExp)
#  data(miRNAExp)
#  data1=FSbyVar(GeneExp, cut.type="topk",value=1000)
#  data2=FSbyVar(miRNAExp, cut.type="topk",value=300)
#  GBM=list(GeneExp=data1,miRNAExp=data2)
#  lambda=alist()
#  lambda[[1]] = 30
#  lambda[[2]] = c(20,1)
#  lambda[[3]] = c(20,20)
#  lambda[[4]] = 30
#  lambda[[5]] = c(30,20)
#  method = c('lasso', 'enet', 'flasso', 'glasso', 'gflasso')
#  result=ExecuteiCluster(GBM, 3, lambda=lambda, method=method)

## ----eval = FALSE--------------------------------------------------------
#  data(GeneExp)
#  data(miRNAExp)
#  GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#  result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)

## ----eval = FALSE--------------------------------------------------------
#  data(GeneExp)
#  data(miRNAExp)
#  data(time)
#  data(status)
#  data1=FSbyCox(GeneExp,time,status,cutoff=0.05)
#  data2=FSbyCox(miRNAExp,time,status,cutoff=0.05)
#  GBM=list(GeneExp=data1,miRNAExp=data2)
#  result=ExecuteSNF.CC(GBM, clusterNum=3, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500,
#                       title = "GBM", plot = "png", finalLinkage ="average")

## ----eval = FALSE--------------------------------------------------------
#  data(GeneExp)
#  data(miRNAExp)
#  data(Ranking)
#  ####Retrieve there feature ranking for genes
#  gene_Name=rownames(GeneExp)
#  index1=match(gene_Name,Ranking$mRNA_TF_miRNA.v21._SYMBOL)
#  gene_ranking=data.frame(gene_Name,Ranking[index1,],stringsAsFactors=FALSE)
#  index2=which(is.na(gene_ranking$ranking_default))
#  gene_ranking$ranking_default[index2]=min(gene_ranking$ranking_default,na.rm =TRUE)
#  
#  ####Retrieve there feature ranking for genes
#  miRNA_ID=rownames(miRNAExp)
#  index3=match(miRNA_ID,Ranking$mRNA_TF_miRNA_ID)
#  miRNA_ranking=data.frame(miRNA_ID,Ranking[index3,],stringsAsFactors=FALSE)
#  index4=which(is.na(miRNA_ranking$ranking_default))
#  miRNA_ranking$ranking_default[index4]=min(miRNA_ranking$ranking_default,na.rm =TRUE)
#  ###Clustering
#  ranking1=list(gene_ranking$ranking_default ,miRNA_ranking$ranking_default)
#  GBM=list(GeneExp,miRNAExp)
#  result=ExecuteWSNF(datasets=GBM, feature_ranking=ranking1, beta = 0.8, clusterNum=3,
#                     K = 20,alpha = 0.5, t = 20, plot = TRUE)

## ------------------------------------------------------------------------
data(GeneExp)
data(miRNAExp)
GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20,plot = FALSE)

###Similarity smaple matrix
sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
plot(sil)

## ------------------------------------------------------------------------
sil1=silhouette(result$group, result$distanceMatrix)
plot(sil1)  ##wrong result

## ------------------------------------------------------------------------
data(GeneExp)
data(miRNAExp)
data(time)
data(status)
data1=FSbyCox(GeneExp,time,status,cutoff=0.05)
data2=FSbyCox(miRNAExp,time,status,cutoff=0.05)
GBM=list(GeneExp=data1,miRNAExp=data2)

#### 1.ExecuteSNF
result1=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20,plot = FALSE)
group1=result1$group
distanceMatrix1=result1$distanceMatrix
p_value=survAnalysis(mainTitle="GBM1",time,status,group1,
                     distanceMatrix1,similarity=TRUE)

## ----results="hide",message=FALSE----------------------------------------
#### 2.ExecuteSNF.CC
result2=ExecuteSNF.CC(GBM, clusterNum=3, K=20, alpha=0.5, t=20,
                      maxK = 5, pItem = 0.8,reps=500, 
                      title = "GBM2", plot = "png", 
                      finalLinkage ="average")

## ------------------------------------------------------------------------
group2=result2$group
distanceMatrix2=result2$distanceMatrix
p_value=survAnalysis(mainTitle="GBM2",time,status,group2,
                     distanceMatrix2,similarity=TRUE)


## ----warning=FALSE-------------------------------------------------------
data(GeneExp)
data(miRNAExp)
data(time)
data(status)
GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20,plot = FALSE)
group=result$group
sigclust=sigclustTest(miRNAExp,group, nsim=1000, nrep=1, icovest=1)
sigclust

## ----results="hide",message=FALSE,fig.show="hide"------------------------
library("RTCGA.mRNA")
#require(TCGAbiolinks)
rm(list = ls())
data(BRCA.mRNA)
mRNA=t(as.matrix(BRCA.mRNA[,-1]))
colnames(mRNA)=BRCA.mRNA[,1]
mRNA1=data.imputation(mRNA,fun="microarray")
mRNA1=FSbyMAD(mRNA1, cut.type="topk",value=5000)

###Split the normal and tumor samples
index=which(as.numeric(substr(colnames(mRNA1),14,15))>9)
mRNA_normal=mRNA1[,index]
mRNA_tumor=mRNA1[,-index]

### Remove the duplicate samples
index1=which(as.numeric(substr(colnames(mRNA_tumor),14,15))>1)
mRNA_tumor=mRNA_tumor[,-index1]

##### Identify cancer subtypes
result=ExecuteCC(clusterNum=5,d=mRNA_tumor,maxK=5,clusterAlg="hc",distance="pearson",title="BRCA")
group=result$group
res=DiffExp.limma(Tumor_Data=mRNA_tumor,Normal_Data=mRNA_normal,group=group,topk=NULL,RNAseq=FALSE)

## ------------------------------------------------------------------------
## Differently expression genes in subtype 1
head(res[[1]])

## Differently expression genes in subtype 2
head(res[[2]])

