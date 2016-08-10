

#'Biological feature (such as gene) selection based on Cox regression model.
#'
#' Cox model (Proportional hazard model) is a statistical approach for survival risk analysis. We applied the univariate Cox model
#' for feature selection. The proportional hazard assumption test is used to evaluate the significant level of each biological feature
#' related to the survival result for samples. 
#' Eventually, the most significant genes are selected for clustering analysis. 
#' 
#' @importFrom survival coxph Surv
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#' @param time A numeric vector representing the survival time (days) of a set of samples. Note that the
#' order of the time should map the samples in the Data matrix.
#' @param status A numeric vector representing the survival status of a set of samples. 0=alive/censored, 1=dead. Note that the
#' order of the time should map the samples in the Data matrix.
#' @param cutoff A numeric value in (0,1) representing whether the significant feature Xi is selected
#' according to the Proportional Hazards Assumption p-value of the feature Xi. If p-value(Xi)<cutoff, the
#' features Xi will be selected for downstream analysis. Normally the significant level is set to 0.05.
#' @return A data matrix, extracted a subset with significant features from the input data matrix.
#' The rows represent the significant features, and the columns represents the samples.
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'  
#'  
#' @examples
#' data(GeneExp)
#' data(time)
#' data(status)
#' data1=FSbyCox(GeneExp,time,status,cutoff=0.05)
#' @references
#' Andersen, P. and Gill, R. (1982). Cox's regression model for counting processes, a large sample study. Annals of Statistics 10, 1100-1120.\cr
#' Therneau, T., Grambsch, P., Modeling Survival Data: Extending the Cox Model. Springer-Verlag, 2000.
#' @export
FSbyCox<-function(Data,time,status,cutoff=0.05)
{
  feature_cox_test=NULL;
  for(i in 1:nrow(Data))
  {
    dataset=list(time,status,x<-Data[i,])
    temp<-coxph(Surv(time, status) ~ x, dataset)
    feature_cox_test=rbind(feature_cox_test,summary(temp)$coefficients)
  }
  indexNum=which(feature_cox_test[,5]<cutoff)
  selectData=Data[indexNum,]
  selectData
}



#' Biological feature (such as gene) selection based on the most variant Median Absolute Deviation (MAD).
#'
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param cut.type A character value representing the selection type. The optional values are shown below:
#' \itemize{
#' \item "topk"
#' \item "cutoff"
#' }
#' @param value A numeric value.\cr
#' If the cut.type="topk", the top number of value features are selected.\cr
#' If the cut.type="cutoff", the features with (MAD>value) are selected.
#' @return An extracted subset data matrix with the most variant MAD features from the input data matrix.
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'  
#'  
#' @examples
#' data(GeneExp)
#' data1=FSbyMAD(GeneExp, cut.type="topk",value=1000)
#' 
#' @export
FSbyMAD<-function(Data, cut.type="topk", value)
{
  mads=apply(Data,1,mad)
  feature_num=length(mads)
  hist(mads, breaks=feature_num*0.1, col="red",
       main="Expression (MAD) distribution",
       xlab="The MAD of feature")
  if(cut.type=="topk")
  {
    index= sort(mads,decreasing = TRUE,index.return=TRUE)
    if(value>nrow(Data))
    {
      value=nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff=index$x[value]
    abline(v=cutoff,col = "blue",lty = 5,lwd=1.5)
    index=index$ix[1:value]
    selectData=Data[index,]
  }
  if(cut.type=="cutoff")
  {
    abline(v=value,col = "blue",lty = 5,lwd=1.5)
    index=which(mads>value)
    selectData=Data[index,]
  }
  selectData
  
  
}


#' Biological feature (such as gene) selection based on the most variance.
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param cut.type A character value representing the selection type. The optional values are shown below:
#' \itemize{
#' \item "topk"
#' \item "cutoff"
#' }
#' @param value A numeric value.\cr
#' If the cut.type="topk", the top number of value features are selected.\cr
#' If the cut.type="cutoff", the features with (var>value) are selected.
#' @return An extracted subset data matrix with most variance features from the input data matrix.
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'  
#'  
#' @examples
#' data(GeneExp)
#' data1=FSbyVar(GeneExp, cut.type="topk",value=1000) 
#' 
#' @export
FSbyVar<-function(Data, cut.type="topk", value)
{
  vars=apply(Data,1,var)
  feature_num=length(vars)
  hist(vars, breaks=feature_num*0.1, col="red",
       main="Expression (Variance) distribution",
       xlab="The Variance of feature")
  if(cut.type=="topk")
  {
    index= sort(vars,decreasing = TRUE,index.return=TRUE)
    if(value>nrow(Data))
    {
      value=nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff=index$x[value]
    abline(v=cutoff,col = "blue",lty = 5,lwd=1.5)
    index=index$ix[1:value]
    selectData=Data[index,]
  }
  if(cut.type=="cutoff")
  {
    abline(v=value,col = "blue",lty = 5,lwd=1.5)
    index=which(vars>value)
    selectData=Data[index,]
  }
  selectData
}


#' Biological feature (such as gene) dimension reduction and extraction based on Principal Component Analysis.
#' 
#' This function is based on the prcomp(), we write a shell for it and make it easy to use on genomic data.
#'
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param PC_percent A numeric values in [0,1]  representing the ratio  of principal component is seclected.
#' @param scale A bool variable, If true, the Data is normalized before PCA.
#'
#' @return
#' A new matrix with full or part Principal Component in new projection space.
#' @author
#'  Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'  
#'  
#' @examples
#' data(GeneExp)
#' data1=FSbyPCA(GeneExp, PC_percent=0.9,scale = TRUE)
#'  
#' @export
#'
FSbyPCA<-function(Data,PC_percent=1,scale = TRUE)
{
  newData=t(Data)
  aa=prcomp(newData,scale = scale)
  vars <- apply(aa$x, 2, var)
  props <- vars / sum(vars)
  xx=as.vector(cumsum(props))
  num=which(xx>PC_percent)[1]
  coeff=aa$ro
  score=(newData%*%coeff)[,1:num]
  result=t(score)
  result
}
