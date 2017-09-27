### merge methods
# mergeMeth <- function(a,b){
#   aMean <- mean(a)
#   bMean <- mean(b)
#   num <- length(a)
#   c <- array(NA,num)
# 
#   for (i in 1:num){
#     if(a[i] > aMean & b[i] > bMean){
#       c[i] <- mean(c(a[i],b[i]))
#     }else if(a[i] < aMean & b[i] < bMean){
#       c[i] <- mean(c(a[i],b[i]))
#     }else{
#       if (abs(a[i]-aMean) < abs(b[i]-aMean)){
#         c[i] <- a[i]
#       }else{
#         c[i] <- b[i]
#       }
#     }
#   }
#   c
# }

################## transform the result predlist to a matrix ###################
#' @title Transform Prediction Result List to Matrix
#' @description This function provides a way to get prediction matix from prediction list.
#' @param predList  prediction list.
#' @return
#' the prediction result matrix.
#' @author Chuang Ma  , Qian Cheng , Zhixu Qiu , Jie Song
#' @keywords transform
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#' 
#' ## cross validation ##
#' predlist <- G2PCrossvalidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
#'                 pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP", "RFC"),
#'                 outputModel = FALSE)
#' resultMat <- resultsMerge(predlist)
#' }
#' 
resultsMerge <- function(predList){
  total_i <- dim.data.frame(predList)[2]
  predMatrix <- rbind(predList[[1]],predList[[2]])
  
  if (total_i > 2){ 
    for(i in 3:total_i){
      predMatrix <- rbind(predMatrix,predList[[i]])
    }
  }else{
    predMatrix <- predMatrix
  }
  predMatrix
}

############################# merge function ###########################
#' @title Prediction Integration from Two Methods 
#' @description This function provides a strategy to integrate the results of two or more algorithms.
#' @param predResMat  (numeric, matrix)the prediction results of algorithms which you want to merge, the first column is the real value of trait.
#' @param ratio  (numeric,array) the weights of every algorithms.
#' @param autoOptimize  (logical)if auto select two method results from multi-results and then compute the mean of two methods results (1:1), default FALSE.
#' @details 
#' The predResMat must including real value in first clumn, and if you set "autoOptimize = T", the count of algorithms must more than 2.
#' 
#' In this function, if autoOptimize = T, the final two algorithms merge are selected from multi methods by following strategy:
#' Firstly, compute the pearson's correlation of predResMat, choose the best correlation between real value and prediction scores, named method 1.
#' Secondly, choose the best correlation between method 1 and other methods, named method 2.
#' Finally, merge method 1 and method 2 with 1:1(mean).
#' 
#' This function auto merge only provide integrate base \bold{pearson's correlation} evaluation.
#' @seealso 
#' \link{GSEnsemble}
#' @return
#' a matrix:  
#' involve real value of trait, merge algorithms and the merge result.
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu, Jie Song
#' @keywords integration, merge
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#' 
#' ## cross validation ##
#' predlist <- G2PCrossvalidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
#'                 pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP", "RFR"),
#'                 outputModel = FALSE)
#' resultMat <- resultsMerge(predlist)
#' 
#' ## merge ##
#' inter <- GSIntegrate(predResMat = resultMat[,1:6],
#'                      ratio = c(2,3,4,4,5), autoOptimize = F)
#' interAuto <- GSIntegrate(predResMat = resultMat, autoOptimize = T,
#'                          allMethodPredResMat = resultMat)
#' }

GSIntegrate <- function(predResMat, ratio, autoOptimize = F ){
  methodCount <- ncol(predResMat) - 1
  ## normalization
  rowCount <- nrow(predResMat)
  
  ## Auto optimize
  if(autoOptimize == T){
    
    # probMat <- transValMat2ProbMat(evalMat = predResMat,BestIndividuals = BestIndividuals)
    corMat <- cor(predResMat[,-1])
    evalresult <- evaluateGS(realScores = predResMat[,1],predScores = predResMat[,2:ncol(predResMat)],evalMethod = "pearson",topAlpha = 1:90)
    a <- evalresult$corMethosds[1,]
    setMax <- which(a == max(a))
    b <- rank(a) - rank(corMat[,setMax])
    setPair <- which(b == max(b))
    
    c <- (predResMat[,setMax + 1]+predResMat[,setPair + 1])/2
    mergePredRes <- as.matrix(c)
    mergePredRes <- cbind(predResMat[,c(1,setMax + 1,setPair + 1)],mergePredRes)
    colnames(mergePredRes)[4] <- "merge"
  }else{
    ## User optimize 
    # probMat <- transValMat2predResMat(evalMat = probMat,BestIndividuals = BestIndividuals)
    mergePredRes <- as.matrix(rep(ratio, each = rowCount) * predResMat[,-1])
    mergePredRes <- as.matrix(apply(mergePredRes,1,sum))/sum(ratio)
    colnames(mergePredRes) <- "merge"
    mergePredRes <- cbind(predResMat,mergePredRes)
  }
  mergePredRes
}


