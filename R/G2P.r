################################################## Part1 GSDataQC ###################################################
#' @title Genomic Selection Data Quality control
#' @description  This function can examine and summary the quality of GS data. And can be used for imputation.  
#' @param markers  (numeric, matrix)row is sample well column is SNP information (feature).Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param phenotype  (numeric)the phenotype value of each individual.
#' @param impute  (logical)if TRUE, imputation, default F.
#' @param filter (logical)if TRUE, filter the markers data with the setting of MAF and missing value.
#' @param NABound (numeric, [0,1])max missing value percentage.
#' @param MAFBound (numerix, [0,1])min MAF percentage in each marker local sit.
#' @param imputeMethod  (character)the method of imputation, "mean", "median" or "KNN", default "mean".
#' @param round  (numeric)rounds the values in its first argument to the specified number of decimal places, default 4.
#' @param k,maxp,rowmax,colmax,rng.seed (numeric)KNN method parameters.
#' @seealso
#' \pkg{impute}
#' @return
#' A list of the data quality information.
#' 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords feature selction, Gini, Accuracy, RR-BLUP
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ## generate missing value
#' misIndex <- sample(1:242000,100000)
#' Markers[misIndex] <- NA
#' 
#' ## GSDataQC, not impute ##
#' QCRes <- GSDataQC(markers = Markers, phenotype = phenotype, impute = F)
#' 
#' ## GSDataQC, not impute ##
#' QCResImpute <- GSDataQC(markers = Markers, phenotype = phenotype, impute = T, 
#'                                     imputeMethod = "mean")
#' }
GSDataQC <- function(markers, phenotype, impute = F, filter = F, NABound = 0.8, MAFBound = 0.005, imputeMethod = "mean", round = 4 ,k = 10, maxp = 1500, rowmax = 0.5, colmax = 0.8, rng.seed = 362436069){

  if (! is.matrix(markers)){
    markerFormat <- warning(paste0("The formation of markers is wrong! Please check!","\n"))
  }else{
    markerFormat <- cat(paste0("The formation of markers is right!","\n"))
  }
  if (! is.numeric(markers)){
    markersClass <- warning(paste0("The data in matrix is not numeric! Please transform!","\n"))
  }else{
    markersClass <- cat(paste0("The data in matrix is right!","\n"))
  }
  if(! is.numeric(phenotype)){
    phenotypeFormat <- warning(paste0("The phenotype data is not numeric!Please transform!","\n"))
  }
  if (is.matrix(markers) & is.numeric(markers) & is.numeric(phenotype)) {
    MAF <- round(MAFSummary(markers),2)
    if (length(which(is.na(markers)) == TRUE) == 0){
      NACount <- 0 
      NAPercentTotal <- 0
      NACountRow <- 0
      NACountCol <- 0
      NACountRowPercent <- 0
      NACountColPercent <- 0
      resImpute <- markers
      NASummary <- cat(paste0("The marker matrix have no missing value!","\n"))
      NAIdx <- NA
      tableRes <- round(table(markers))
    }else{
      NACount <- round(length(which(is.na(markers)) == TRUE))
      NAPercentTotal <- round(NACount/(ncol(markers)*nrow(markers)),round)
      
      NACountRow <- round(apply(markers,1,function(x){length(which(is.na(x)) == TRUE)}))
      NACountCol <- round(apply(markers,2,function(x){length(which(is.na(x)) == TRUE)}))
      
      NACountRowPercent <- round(NACountRow/ncol(markers),2)
      NACountColPercent <- round(NACountCol/nrow(markers),2)
      
      NAIdx <- which(is.na(markers) == TRUE)
      
      tableRes <- round(table(markers))
      
      if(filter){
        NAFliIdx <- which(NACountColPercent >= NABound)
        MAFFliIdx <- which(MAF <= MAFBound)
        fliterIdx <- unique(c(NAFliIdx,MAFFliIdx))
        markers <- markers[,-fliterIdx]
      }
      
      if(impute){
        resImpute <- GSImputation(x = markers, imputeMethod = imputeMethod , k = k,rowmax = rowmax,colmax = colmax, maxp = maxp, rng.seed = rng.seed)
      }else{
        resImpute <- markers
      }
    }
    
    phenotypeNACount <- round(length(which(is.na(phenotype) == TRUE)))
    phenotypeNAIdx <- which(is.na(phenotype) == TRUE)
    
    NASum <- c(NACount,NAPercentTotal)
    names(NASum) <- c("missValueCount","missValuePercent") 
    
    markersTable <- c(tableRes,NASum)
    
    ## summarize
    #     summarize <- paste0(markerFormat,markersClass,"The data have ",NACount," missing value!","\n","Missing value percent:",NAPercentTotal,"\n",
    #                             "The phenotype data have ",phenotypeNACount," missing value!","\n",
    #                             "The data have ",length(tableRes)," element.","\n","The detail in markerReport!",
    #                             "\n")
    markerReport <- list(Global = markersTable, Imputation = resImpute,MAF = MAF, NACountRow = NACountRow,NACountRowPercent = NACountRowPercent,NACountCol = NACountCol,
                         NACountColPercent = NACountColPercent, NAIdx = NAIdx,penotypeNAIdx = phenotypeNAIdx)
    cat(paste0(markerFormat,markersClass,"The data have ",NACount," missing value!","\n","Missing value percent:",NAPercentTotal,"\n",
               "The phenotype data have ",phenotypeNACount," missing value!","\n",
               "The data have ",length(tableRes)," element.","\n","The detail in markerReport!",
               "\n"))
    markerReport
  }else{
    stop(paste0("The format of data is wrong! There may be the following above cases:","\n","The format of maker is not a matrix;","\n",
                "The class of markers or phenotype is not numeric."))
  }
}
GSImputation <- function(x, imputeMethod = "mean", k = 10,rowmax = 0.5,colmax = 0.8,maxp = 1500, rng.seed = 362436069){
  require("impute")
  if(is.numeric(x)){
    if(imputeMethod == "KNN"){
      requireNamespace("impute", quietly = TRUE)
      x <- impute.knn(data = t(x), k = k, rowmax = rowmax, colmax = colmax, maxp = maxp, rng.seed = rng.seed)
      x <- round(t(x$data))
    }else{
      imputeMethod <- imputeMethod
    }
    if(imputeMethod == "mean"){
      x[which(is.na(x))] <- mean(x,na.rm=TRUE)
    }else if(imputeMethod=="median"){
      x[which(is.na(x))] <- median(x,na.rm=TRUE)
    }
  }else{
    if(imputeMethod =="mean"){
      stop("Method 'mean' is not available for non-numeric vectors.",call. = FALSE)
    }else if(imputeMethod=="median"){
      tt <- table(x)
      x[which(is.na(x))] <-  names(tt)[which(tt== max(tt))]
    }else{
      x[which(is.na(x))] <-  names(tt)[which(tt==max(tt))]
    }
  }
  return(x)
}
MAFSummary <- function(x){
  MAF <- apply(x,2,function(x){
    level <- table(x)
    if (length(level) == 2){
      index <- order(level,decreasing = T)
      MAF <- level[index[2]]/sum(level)
    } else if(length(level) == 3){
      index <- order(level[-2],decreasing = T)
      MAF <- (level[-2][index[2]]*2 + level[2])/(sum(level)*2)
    }
  })
  MAF
}

##################################################################################################333
############################## Feature Selection ############################
#' @title Feature Selection
#' @description  This function scores each marker,so that reduce the data dimension and perform feature selection.(Methods including Gini,Accuracy and rrBLUP).
#' @param markers  (numeric, matrix)row is sample well column is SNP information (feature).Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param phenotype  (numeric)the phenotype value of each individual.
#' @param method  (character)the method of feature selction including "Gini" "Accuracy" "rrBLUP", default "RR-BLUP"
#' @param ntree  (numeric)the number of random forest decision tree, default 500
#' @param importance  (logical)whether the results of variable importance,default TRUE
#' @param posPercentage  (numeric,[0,1])phenotype of extreme individuals which expected, default 0.4
#' @param BestIndividuals (character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",defult "top".
#' @return
#' A numeric mode score of each position of SNPs
#' 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords feature selction, Gini, Accuracy, RR-BLUP
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ## feature selection with Gini ##
#' Gini_selection <- feature_assess(markers = Markers, phenotype = phenotype, method = "Gini", 
#'                                  ntree = 500, importance = TRUE, posPercentage = 0.40,
#'                                  BestIndividuals = "top")
#'
#' ## feature selection with Acc ##
#' Acc_selection <- feature_assess(markers = Markers, phenotype = phenotype, method = "Accuracy", 
#'                                 ntree = 500, importance = TRUE, posPercentage = 0.40,
#'                                 BestIndividuals = "top")
#'
#' ## feature selection with rrBLUP ##
#' rrBLUP_selection <- feature_assess(markers = Markers, phenotype = phenotype, method = "rrBLUP", 
#'                                    posPercentage = 0.40, BestIndividuals = "top")
#' }
##SNP feature order Gini Accuracy rrBLUP
feature_assess <- function(markers, phenotype, method = c("rrBLUP", "Gini", "Accuracy"), ntree = 500, importance = TRUE, posPercentage = 0.40, BestIndividuals = c("top", "middle", "buttom")){
  require("rrBLUP")
  require("randomForest")

  if(length(method) > 1){
    method = "rrBLUP"
  }
  if(length(BestIndividuals > 1)){
    BestIndividuals <- "top"
  }
  posNegSampleList <- sampleClassify(phenotype = phenotype, posPercentage = posPercentage, BestIndividuals = BestIndividuals )
  Xtrain <- markers[c(posNegSampleList$posSampleIndex,posNegSampleList$negSampleIndex),]
  YClasif <-  as.factor( c( rep("1", length(posNegSampleList$posSampleIndex)), rep("0", length(posNegSampleList$negSampleIndex)) ) )
  switch(method,
         rrBLUP = trainModel_RRBLUP( markers , phenotype  )$phD,
         Accuracy =  randomForest( x = Xtrain, y = YClasif, ntree = ntree, importance = importance)$importance[,"MeanDecreaseAccuracy"],
         Gini = randomForest( x = Xtrain, y = YClasif, ntree = ntree, importance = importance)$importance[,"MeanDecreaseGini"]
  )
}

########################################################################################################3
# ######################### get system time for seed and then generate random index ###########
#' @title Generate Random Seed
#' @description This funcation is appplied for generating random seed with current system time
#' @export
#' @return  
#' (numeric) A random seed 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}
########################## load example data ###########################
#' @name GYSS
#' @title Example Data for G2P
#' @description
#' The data of maize yield under drought stressed.(SNP genotypes informations) 
#'  \itemize{
#'  \item{Markers} {A numeric matrix, each row is the each individual's SNP genotypes informations.}
#'  \item{phenotype} {The real phenotype Value of each individual.}
#'  }
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords datasets
#' @docType data
#' @examples
#' \dontrun{
#' ## load maize yield data sets
#' data(GYSS)
#' }
NULL
########################## marker data ###########################
#' @name Markers
#' @title Example Data for G2P
#' @description
#' A numeric matrix.
#'  \itemize{
#'  \item{each row represents each sample}
#'  \item{each column represents each SNP locus}
#'  }
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords datasets
#' @docType data
#' @format matrix
NULL
########################## phenotype data ###########################
#' @name phenotype
#' @title Example Data for G2P
#' @description
#' The real phenotype Value of each individual
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords data , phenotype value
#' @docType data
#' @format vector
NULL
########### generated positive and negative samples for training ####################
#' @title Generate Positive and Negative Samples for Training 
#' @description This function can be use to generate positive and negative samples for training.The positive samples represent the excellent individuals which's breeding values we expect to obtain in your research.And the negative samples represent the lower breeding values of individuals. 
#' @param phenotype  (numeric)the breeding values of each individual.
#' @param posPercentage  (numeric,[0,1])phenotype of extreme individuals which expected, default 0.4
#' @param BestIndividuals (character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".
#' @return  
#' A list of row number of positive and negative samples
#' $posSampleIndex Index of positive samples
#' $negSampleIndex  Index of negative samples
#' 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords positive , negative 
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ## percentage of positive samples is 0.4 ##
#' sampleCly <- sampleClassify(phenotype, posPercentage = 0.4, BestIndividuals = "top")
#' } 
sampleClassify <- function(phenotype, posPercentage = 0.4, BestIndividuals = c("top", "middle", "buttom")){
  trainSampleSize <- length(phenotype)
  posSampleSize <- round( posPercentage*trainSampleSize )
  
  if( BestIndividuals == "top" ){
    trainPhenOrder <- order( phenotype, decreasing = TRUE )
  }else if( BestIndividuals == "buttom" ){
    trainPhenOrder <- order( phenotype, decreasing = FALSE )   
  }else if( BestIndividuals == "middle"){
    trainPhenOrder <- order( abs(phenotype), decreasing = FALSE ) 
  }
  posIdx <- trainPhenOrder[1:posSampleSize]
  negIdx <- setdiff( c(1:trainSampleSize), posIdx )
  
  res <- list( posSampleIndex = posIdx, negSampleIndex = negIdx )
  res
}
########################generate train idx and test idx ##########################
#' @title Generate Sample Indices for Training Sets and Testing Sets
#' @description  This function generates indices for samples in training and testing sets for performing the N-fold cross validation experiment.
#' @param sampleNum  The number of samples needed to be partitioned into training and testing sets.
#' @param cross  The fold of cross validation.
#' @param seed  An integer used as the seed for data partition. The default value is 1.
#' @param randomSeed  Logical variable, default FALSE.
#' @return 
#' A list and each element including $trainIdx $testIdx and $cvIdx
#' 
#' $trainIdx  The index of training samples.
#' 
#' $testIdx   The index of testing samples.
#' 
#' $cvIdx     The cross validation index.
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#' ## leave-one out cross validation
#' a <- cvSampleIndex(sampleNum = nrow(Markers), cross = nrow(Markers), seed = 1)
#' 
#' ## 5-fold cross validation
#' b <- cvSampleIndex(sampleNum = nrow(Markers), cross = 5, seed = 1)
#' }
##get sample idx for training and testing
cvSampleIndex <- function( sampleNum, cross = 5, seed = 1,randomSeed = FALSE ) {
  if(randomSeed == TRUE){
    seed <- randomSeed()
  }
  cv <- cross
  resList <- list()
  
  # leave-one-out
  if( cv == sampleNum ){
    vec <- 1:sampleNum
    for( i in 1:sampleNum ){
      resList[[i]] <- list( trainIdx = vec[-i], testIdx = i, cvIdx = i)
    }
  }else {
    #random samples 
    set.seed(seed)
    index <- sample(1:sampleNum, sampleNum, replace = FALSE )
    step = floor( sampleNum/cv )
    
    start <- NULL
    end <- NULL
    train_sampleNums <- rep(0, cv)
    for( i in c(1:cv) ) {
      start <- step*(i-1) + 1
      end <- start + step - 1
      if( i == cv ) 
        end <- sampleNum
      
      testIdx <- index[start:end]
      trainIdx <- index[-c(start:end)]
      resList[[i]] <- list( trainIdx = trainIdx, testIdx = testIdx, cvIdx = i)
    }
  }
  names(resList) <- paste0("cv",1:cross)
  resList
}
################################################################################################################################
####################################################### Fit regression model ###################################################
### models from BGLR package
## model: BayesA (Scaled-t prior), BayesB (point of mass at zero + scaled-t slab), BayesC, Bayesian LASSO (BL), Bayesian Ridge Regression(BR; (Gaussian prior), equivalent to G-BLUP), RKHS
#' @title Fit Regression Model  
#' @description This function can fit several regression models of genomic selection such as BayesA, BayesB, BayesC, BRR(BayesBayesian Ridge Regression) and BL(Bayesian LASSO).
#' @param trainedMarkerMat  (numeric, matrix)each row is the each training sets individual's SNP genotypes informations.Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedphenotype  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric, matrix)each row is the each testing sets individual's SNP genotypes informations.Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param modelMethods  (character)the model to fit. "BayesA", "BayesB", "BayesC", "BL", "BRR".
#' @param nIter,burnIn,thin  (integer)the number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.
#' @param saveAt  (string)this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default "".
#' @param S0,df0  (numeric)the scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Default S0 NULL,df0 5.
#' @param R2  (numeric, (0,1))the proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Defult 0.5
#' @param weights  (numeric, n)a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.
#' @param verbose  (logical)if TRUE the iteration history is printed, default FALSE.
#' @param rmExistingFiles  (logical)if TRUE, removes existing output files from previous runs, default TRUE.
#' @param groups  (factor)a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.
#' @param outputModel  (logical)if TRUE, return the list of training model and prediction result, default FALSE.
#' @seealso 
#' \pkg{BGLR}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords regression model, BayesA, BayesB, BayesC, BRR, BL
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## Fit Bayes A model and prediction ##
#' BayesA_model <- fit.BGLR(trainedMarkerMat = Markers[1:200,], predictMarkerMat = Markers[201:242,],
#'                          trainedPheVal = phenotype[1:200], modelMethods = "BayesA", outputModel = TRUE)
#' }
fit.BGLR <- function( trainedMarkerMat, trainedPheVal, predictMarkerMat, modelMethods,  
                      outputModel = FALSE, nIter = 1500, burnIn = 500, thin = 5,
                      saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                      verbose = FALSE, rmExistingFiles = TRUE, groups=NULL){
  require("BGLR")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  numT <- nrow(predictMarkerMat)
  MarkerMat <- rbind(predictMarkerMat,trainedMarkerMat)
  pheVal <- c(rep(NA,numT),trainedPheVal)
  
  checkModel <- modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")
  if( ! checkModel ) {
    stop("Error: not defined model!")
  }
  
  BGLRModel.fit <- NULL
  tryCatch({
    if( modelMethods != "RKHS") {
      BGLRModel.fit <- BGLR( y = pheVal, ETA=list(list(X=MarkerMat, model= modelMethods)),  nIter = nIter, burnIn = burnIn, thin = thin, 
                             saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                             verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups )
      
      
    }else{
      # Computing the genomic relationship matrix
      X <-scale(x = MarkerMat,center=TRUE,scale=TRUE)
      G <-tcrossprod(X)/ncol(X)
      
      BGLRModel.fit <- BGLR( y=pheVal, ETA=list(list(K=G, model= modelMethods)),  nIter = nIter, burnIn = burnIn, thin = thin, 
                             saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                             verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups )
      
    }
  }, error=function(e){
    gc(verbose=F)
    BGLRModel.fit <- NULL
  }, warning=function(w){
    gc(verbose=F)
    BGLRModel.fit <- NULL
  }); gc(verbose=F)
  
  Res <- BGLRModel.fit$yHat[1:numT]
  if(outputModel){
    Res <- list(model = BGLRModel.fit,predictRes = Res)
  }
  Res
}
################# 
trainModel_RRBLUP <- function( markerMat, phenVec ){
  phen_answer<-mixed.solve(phenVec, Z = markerMat, K=NULL, SE = FALSE, return.Hinv=FALSE)
  beta <- phen_answer$beta
  phD <- phen_answer$u
  e <- as.matrix(phD)
  return( list(beta = beta, e = e, phD = phD) ) 
}

#' @title Fit Regression Model  
#' @description This function can fit several regression models of genomic selection including SPLS, rrBLUP, LASSO and bigRR.
#' @param markers  (numeric)a matrix, each row is the each individual's SNP genotypes informations.Genotypes should be coded as {0,1,2}or{-1,0,1};0(-1) represent AA(homozygote),2(1) represent BB(homozygote) and 1(0) represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param pheVal  (numeric)the phenotype value of each individual.
#' @param K  (integer)SPLS model parameter: number of hidden components.
#' @param eta	 (numeric)SPLS model parameter: thresholding parameter. eta should be between 0 and 1.
#' @param select  (character)SPLS model parameter: PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2".
#' @param fit	 (character)SPLS model parameter: PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".
#' @param scale.x	 (character)SPLS model parameter: scale predictors by dividing each predictor variable by its sample standard deviation?
#' @param scale.y  (character)SPLS model parameter: scale responses by dividing each response variable by its sample standard deviation?
#' @param eps  (character)SPLS model parameter: an effective zero. Default is 1e-4.
#' @param maxstep  (integer)SPLS model parameter: maximum number of iterations when fitting direction vectors. Default is 100.
#' @param trace  (logical)SPLS model parameter: print out the progress of variable selection?
#' @param alpha  (numeric)LASSO model parameter: the elasticnet mixing parameter.Detail in glmnet.
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param family the distribution family of y, see help('family') for more details. 
#' @param lambda the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.
#' @param tol.err internal tolerance level for extremely small values; default value is 1e-6.
#' @param tol.conv tolerance level in convergence; default value is 1e-8.
#' @param ... arguments passed to or from other package.
#' @seealso
#' \pkg{rrBLUP}
#' \pkg{hglm}
#' \pkg{glmnet}
#' \pkg{spls}
#' @return  
#' A regression model which is enable to predict.
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords regression model, rrBLUP, LASSO, SPLS, bigRR
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## Fit rrBLUP model ##
#' rrBLUP_model <- GSReModel(markers = Markers,pheVal = phenotype,modelMethods = "rrBLUP")
#' }

GSReModel <- function(markers, pheVal, modelMethods,
                      K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,
                      alpha = 1,
                      X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8, weights = NULL,
                      ...){  
  
  checkModel <- modelMethods %in% c("rrBLUP","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GSReModel!")
  }
  #   if (modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")){
  #     BGLRmethods <- modelMethods
  #     modelMethods <- "BGLRModel"
  #   }
  
  
  switch(modelMethods,
         #        BGLRModel  = trainedPredictModel_BGLR(trainedMarkerMat = markers, trainedPhenVec = pheVal, modelMethods = BGLRmethods ,nIter = nIter, burnIn = burnIn, thin = thin, 
         #                                                saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
         rrBLUP  = trainModel_RRBLUP(markerMat = markers, phenVec = pheVal),
         LASSO = trainModel_LASSO(markers,pheVal,alpha = alpha, ...),
         SPLS = trainModel_spls(markers,pheVal,K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,maxstep = maxstep, ...),
         bigRR = trainModel_bigRR(markers = markers, X = X, pheVal = pheVal, weight = weights,family = family, lambda = lambda, tol.err = tol.err,tol.conv = tol.conv, ...)
  )
  
}


################################  LASSO ##############################
trainModel_LASSO <- function(markers,pheVal,alpha = 1, ...){
  require("glmnet")
  #glmnet fits a lasso model when we specify that alpha=1
  LASSO.fit <- glmnet(y=pheVal,x=markers,alpha=1, ...)
  #cv.glmnet finds the lambda value such that the the cvm value is the minimum
  cv <- cv.glmnet(y = pheVal, x=markers)
  LASSO_Res <- list(LASSO.fit = LASSO.fit,cv = cv)
  LASSO_Res
}

############################### spls #################################
trainModel_spls <- function(markers,pheVal,K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,trace = FALSE, maxstep = 100, ...){
  require("spls")
  f <- spls(markers,pheVal,K = K,eta = eta,select = select,fit = fit,scale.x =scale.x,scale.y = scale.y,eps = eps,trace = trace, ...)
  f
}

############################## bigRR #################################
trainModel_bigRR <- function (markers, X = NULL, pheVal, weight = NULL,
                              family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, 
                              tol.conv = 1e-8, ...){
  require("hglm.data")
  require("hglm")
  if(is.null(X)){
    X <- as.matrix(rep(1,length(pheVal)))
  }
  
  y <- (pheVal - min(pheVal))/(max(pheVal) - min(pheVal)) 
  Z <- markers
  
  Call <- match.call()
  if (!(is.matrix(X))) 
    stop("X should be a matrix.")
  if (!(is.matrix(Z))) 
    stop("Z should be a matrix.")
  if (!(is.vector(y))) 
    stop("y should be a vector.")
  
  if (any(is.na(y))) {
    naidx <- which(is.na(y))
    y <- y[-naidx]
    X <- X[-naidx,]
    Z <- Z[-naidx,]
  }
  
  N <- n <- nrow(X)
  p <- ncol(X)
  k <- ncol(Z)
  if (N != nrow(Z) | N != length(y)) 
    stop("Sizes of y, X, and Z are not all equal.")
  if (is.null(weight)) w <- rep(1, k) else w <- weight
  
  #G <- crossprod(sqrt(w)*t(Z)) ## bug fixed 111201 -- Xia
  wZt <- sqrt(w)*t(Z)
  G <- crossprod(wZt)
  ############ Bending to allow for p<n problems -- Lars (Xia added SVD)
  if (k < n) {
    eigen.values <- eigen(G)$values
    min.eigen <- min(eigen.values)
    if (min.eigen < tol.err) G <- G + diag(N)*(abs(min.eigen) + tol.err) 
  }
  ############
  #invG <- solve(G)
  #L <- t(chol(G))
  svdG <- svd(G)
  L <- svdG$u %*% diag(sqrt(svdG$d))
  invG <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
  phi0 <- sa0 <- 1
  if (is.null(lambda)) {
    hm <- hglm(y = y, X = X, Z = L, family = family, conv = tol.conv, bigRR = TRUE,...) ## checked with old emme code, conv = 1e-6 removed -- Xia
  }
  else {
    start.beta = c(rep(0, p))
    start.v = c(rep(0.01, n))
    start.lambda = lambda
    start.sigma2e = 1
    cat("Only 1 iteration applied for fixed lambda")
    hm <- hglm(y = y, X = X, Z = L, family = family, startval = c(start.beta, 
                                                                  start.v, start.lambda, start.sigma2e), maxit = 1,bigRR = TRUE, ...)
  }
  phi <- as.numeric(hm$varFix)
  sa <- as.numeric(hm$varRanef)
  a <- L%*%hm$ranef
  tZinvG <- crossprod(Z, invG)
  u <- (w*tZinvG)%*%a
  qu <- GCV <- NULL
  result <- list(phi = phi, lambda = hm$varRanef, beta = hm$fixef, hglm = hm,
                 u = u, leverage = qu, GCV = GCV, Call = Call, y = y, X = X)
  class(result) <- "bigRR"
  return(result)
}

#########################################################################################################################
############################### Fit a machine learning  model ############################
#' @title Fit machine learning model  
#' @description This function can fit several machine learning models of genomic selection such as "SVR", "SVC", "RFR" and "RFC".
#' @param markers  (numeric)a matrix, each row is the each individual's SNP genotypes informations.Genotypes should be coded as {0,1,2};0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param pheVal  (numeric)the phenotype value of each individual.
#' @param posPercentage  (numeric,[0,1])phenotype of extreme individuals which expected, default 0.4.
#' @param BestIndividuals (character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".
#' @param modelMethods  (character)alternative machine learning models. "SVR" and "SVC" from SVM, "RFR" and "RFC" from RF.
#' @param ntree  (integer)ramdomforest parameter. Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times,default 500.
#' @param nodesize	(integer)ramdomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param kernel  (numeric)svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial)Default "linear".
#' @param gamma  (numeric)svm parameter parameter needed for all kernels except linear, default 1.
#' @param cost  (numeric)svm parameter cost of constraints violation, default: 2^(-9), it is the 'C'-constant of the regularization term in the Lagrange formulation.
#' @param ... other parameters.
#' @details 
#' SVM (support vector machine) and RF (random forest) models can be fitted in this function, including 2 classification (RFC, SVC) and 2 regression (SVR, SVC) models.
#' @seealso 
#' \pkg{e1017}
#' \pkg{randomForest}
#' @return  
#' a machine learning model
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords model, randomforest, svm
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## Fit RFR model ##
#' machine_model <- GSmachine(markers = Markers, pheVal = phenotype, modelMethods = "RFR")
#' 
#' ## Fit classification model(RFC) ##
#' machine_model <- GSmachine(markers = Markers, pheVal = phenotype, modelMethods = "RFC",
#'                            posPercentage = 0.4, ntree = 500)
#' }


GSmachine <- function(markers, pheVal, modelMethods ="SVC", posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500,
                      nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9), ...){
  
  if(is.null(nodesize)){
   if(modelMethods == "SVC"){
   nodesize <- 1
   }else if(modelMethods == "SVR"){
   nodesize <- 5
   }
  }else{
  nodesize <- nodesize
  }
  require("e1071")
  require("randomForest") 
  if( !modelMethods%in% c("SVR","SVC","RFR","RFC") ) {
    stop("Error: not defined category")
  }
  if (modelMethods %in%  c("SVC","RFC")){
    posNegSampleList <- sampleClassify(phenotype = pheVal ,posPercentage = posPercentage ,BestIndividuals = BestIndividuals )
    markers <- markers[c(posNegSampleList$posSampleIndex,posNegSampleList$negSampleIndex),]
    pheVal <-  as.factor( c( rep("1", length(posNegSampleList$posSampleIndex)), rep("0", length(posNegSampleList$negSampleIndex)) ) )
    
  }
  
  if(modelMethods %in%  c("SVR","SVC")){
    modelMethods <- "svm"
  }
  
  if(modelMethods %in% c("RFR","RFC")){
    modelMethods <- "randomforest"
  }
  switch(modelMethods,
         svm = svm(x= markers, y = pheVal,kernel = kernel,cost=cost,gamma = gamma,probability = TRUE, ...),
         randomforest = randomForest( x = markers, y = pheVal, ntree = ntree, importance = F,nodesize = nodesize), ...)
}

##############################################################################################################################3
########################### the prediction of genomic selection ###########################
#' @title Prediction with Trained Model from Geomic Selection Model
#' @description This function can give the prediction score of a new GS data by using already model.
#' @param testMat  (numeric)a matrix, each row is the each testing sets or new GS data individual's SNP genotypes informations.Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param modelMethods  (character)the type name of training model including "bigRR","rrBLUP","LASSO","SPLS","SVC","SVR","RFC","RFR".
#' @param trainModel  (model)the trained model.It's type must be similar whith modelMethods.
#' @param type  (character)SPLS parameter,detail see package spls.\pkg{spls}.
#' @return
#' a list:  
#' The prediction result of testing sets which predicted through already models
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords model , predicttion
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## Fit rrBLUP model ##
#' rrBLUP_model <- GSReModel(markers = Markers, pheVal = phenotype, modelMethods = "rrBLUP")
#'
#' ## Predict 1-20 subset of all example data with already rrBLUP model ## 
#' res <- predictGS(testMat = Markers[1:20,], trainModel = rrBLUP_model, modelMethods = "rrBLUP")
#' }

predictGS <- function(testMat, trainModel, modelMethods = "SVC",type = "fit"){
  ########## check the methods of GS
  checkModel <- modelMethods %in% c("rrBLUP","SVR","SVC","RFR","RFC","LASSO","SPLS","bigRR")
  if( ! checkModel ) {
    stop("Error: not defined models for implementing GS Model")
  }
  ####### 
  #   if (modelMethods %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")){
  #     Methods <- "BGLRModel"
  #   }else{
  Methods <- modelMethods
  # }
  
  ####### check testset 
  if(!is.matrix(testMat)) {
    testMat <- rbind(testMat,testMat)
    testnum <- 1 
  }else{
    testnum <- nrow(testMat)
  }
  
  
  #############
  if (modelMethods %in% c("rrBLUP","SVR","RFR","LASSO","SPLS","bigRR")){
    predresult <- switch(Methods,
                         #BGLRModel = {mkr.effs <- as.numeric(trainModel$ETA[[1]]$b); testMat %*% mkr.effs},
                         rrBLUP = {pred_phenVec <-  testMat %*% trainModel$e; pred_phenVec[,1] + trainModel$beta },
                         SVR = predict( trainModel, testMat),
                         RFR = predict( object = trainModel,testMat ),
                         LASSO = pred_LASSO(trainModel,testMat),
                         SPLS = pred_SPLS(trainModel,testMat,type= type),
                         bigRR = as.numeric(trainModel$beta + testMat %*% trainModel$u)
    )
  }
  else if (modelMethods %in% c("SVC","RFC")){
    predresult <- switch(modelMethods,
                         SVC = {obj_pred <- predict(trainModel,testMat, probability = TRUE); as.matrix(attr(obj_pred, "probabilities")[,"1"])},
                         RFC = predict(trainModel, testMat, type= "vote" )[,"1"])
  }
  predresult[1:testnum]
}


################### Modeling and predicting using RR ########################
#' @title Modeling and Predicting using ridge regression  
#' 
#' @description This function can fit ridge regression model and export the prediction value of testing sets. 
#' @param trainedMarkerMat  (numeric,matrix)each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric,matrix)each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param cpus  (integer)number of cpu cores to use for calculations (only available in UNIX-like operating systems), default 1.
#' @param outputModel (logical)if true, return the list of training model.
#' @seealso 
#' \pkg{rrBLUP}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords RR model predict
#' @export
#' @examples
#' \dontrun{
#' ## Not run: 
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' rr_Res <- fit.RR(trainedMarkerMat = Markers, trainedPheVal = phenotype,
#'                  predictMarkerMat = Markers[1:10,], cpus = 1 )
#' }
###predicted model using ridge regression

fit.RR <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat, cpus = 1,outputModel = FALSE){
  require("rrBLUP")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  yield_answer <- kinship.BLUP( y = trainedPheVal, G.train = trainedMarkerMat, G.pred= predictMarkerMat, K.method="RR", n.core = cpus)
  Res <-  as.numeric(yield_answer$g.pred) + as.numeric(yield_answer$beta)
  if(outputModel){
    Res <- list(model = yield_answer,predictRes = Res)
  }
  Res
}

#' @title Modeling and predicting using Bayesian Regularization Neural Networks(BRNN)
#' 
#' @description This function can fit BRNN model and export the prediction value of testing sets. 
#' @param trainedMarkerMat  (numeric)a matrix, each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric)a matrix, each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote),2 represent BB(homozygote) and 1 represent AB(heterozygote);missing (NA) alleles are not allowed.
#' @param cpus  (integer)number of cpu cores to be used for calculations (only available in UNIX-like operating systems), default 1.
#' @param neurons  (integer)indicates the number of neurons,defult 4.
#' @param epochs  (integer)maximum number of epochs(iterations) to train, default 30.
#' @param verbose (logical)if TRUE, will print iteration history.
#' @param outputModel (logical)if true, return the list of training model.
#' @param ... other parameters, details see package brnn
#' @seealso 
#' \pkg{brnn}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma  , Qian Cheng , Zhixu Qiu
#' @keywords BRNN, model predict
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' BRNN_Res <- fit.BRNN(trainedMarkerMat = Markers,trainedPheVal = phenotype,
#'                      predictMarkerMat = Markers[1:10,],cpus = 1 )
#' }

################### Modeling and predicting using BRNN#############################
fit.BRNN <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat, outputModel = FALSE,verbose=TRUE, neurons=4, epochs=30, cpus = 1, ...){
  require("brnn")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  X <- rbind(trainedMarkerMat,predictMarkerMat)
  trainNum <- nrow(trainedMarkerMat)
  yTRN <- trainedPheVal
  n<-nrow(X) 
  p<-ncol(X)
  
  for(i in 1:ncol(X)){ (X[,i]<-X[,i]-mean(X[,i]))/sd(X[,i])}
  G<-tcrossprod(X)/ncol(X)
  trainedMarkerMat <- G[1:trainNum,]
  predictMarkerMat <- G[-(1:trainNum),]
  
  NN<-brnn(y=yTRN,x=trainedMarkerMat,neurons=neurons, epochs=epochs,verbose=verbose,cores = cpus,...)  
  Res <- predict(NN, newdata = predictMarkerMat)
  if(outputModel){
    Res <- list(model = NN,predictRes = Res)
  }
  Res
}
####################### LASSO pred ######################### 
pred_LASSO <- function(trainModel,testMat){
  predict(object = trainModel$LASSO.fit,testMat,s = trainModel$cv$lambda.min)
}

###################### pls and spls###########################
pred_SPLS <- function(trainModel,testMat,type = "fit"){
  predict( trainModel,testMat,type=type )
}

#' @title Modeling and predicting using Reproducing Kernel Hilbert Space(RKHS).
#' 
#' @description This function can fit RKHS model and export the prediction value of testing sets. 
#' @param trainedMarkerMat  (numeric)a matrix, each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric)a matrix, each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param nIter,burnIn,thin  (integer)the number of iterations, burn-in and thinning,default nIter 1500, burnIn 500, thin 5.
#' @param saveAt  (string)this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default "".
#' @param S0,df0  (numeric)the scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf. Default S0 NULL,df0 5.
#' @param R2  (numeric, (0,1))the proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf. Defult 0.5
#' @param weights  (numeric, n)a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.
#' @param verbose  (logical)if TRUE, the iteration history is printed, default FALSE.
#' @param rmExistingFiles  (logical)if TRUE, removes existing output files from previous runs, default TRUE.
#' @param groups  (factor)a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.
#' @param outputModel (logical)if true, return the list of training model.
#' @seealso 
#' \pkg{BGLR}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma  , Qian Cheng , Zhixu Qiu
#' @keywords RKHS, model predict
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' RKHS_Res <- fit.RKHS(trainedMarkerMat = Markers, trainedPheVal = phenotype, 
#'                      predictMarkerMat = Markers[1:10,],nIter = 1500, burnIn = 500)
#' }

########################################### modeling and predicting RKHS ###########################################
fit.RKHS <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat, outputModel = FALSE,nIter = 1500, burnIn = 500, thin = 5,
                     saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                     verbose = FALSE, rmExistingFiles = TRUE, groups=NULL){
  require("BGLR")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  numT <- nrow(predictMarkerMat)
  MarkerMat <- rbind(predictMarkerMat,trainedMarkerMat)
  pheVal <- c(rep(NA,numT),trainedPheVal)
  X <-scale(x = MarkerMat,center=TRUE,scale=TRUE)
  G <-tcrossprod(X)/ncol(X)
  
  RKHS.Fit <- BGLR(y=pheVal, ETA=list(list(K=G, model= "RKHS")),  nIter = nIter, burnIn = burnIn, thin = thin, 
                   saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,
                   verbose = verbose, rmExistingFiles = rmExistingFiles, groups = groups)
  Res <- RKHS.Fit$yHat[1:numT]
  if(outputModel){
    Res <- list(model = RKHS.Fit,predictRes = Res)
  }
  Res
}

########################################## modeling and predicting mmer models #######################################
#' @title Modeling and predicting using methods AI, NR , EM , EMMA.
#' @description This function can fit AI, NR , EM , EMMA model and export the prediction values of testing sets. 
#' @param trainedMarkerMat  (numeric)a matrix, each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric)a matrix, each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param Z  A 2-level list,incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. This needs to be provided as a 2-level list structure, defult NULL. 
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param method  (string)select the algorithm to fit the model including NR, AI, EM, EMMA."AI and EM method has been discontinued because of its instability. Try 'EMMA' and 'NR'.See details in the G2P help page.If you insist AI & EM, please use function fit.EM and fit.AI to debug.
#' @param effectConcider  (string)if Z = NULL, random effects are auto generated.
#' @param iters  (numeric)a scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.
#' @param cpus  (numeric)number of cores to use when the user passes multiple responses in the model for a faster execution of univariate models. The default is 1
#' @param verbose  (logical)if TRUE, the iteration history is printed, default FALSE.
#' @param rmExistingFiles  (logical)if TRUE removes existing output files from previous runs, default TRUE.
#' @param outputModel (logical)if true, return the list of training model.
#' @param ... other parameters, details see package sommer.
#' @seealso 
#' \pkg{sommer}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords AI, NR, EM, EMMA, model predict
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' EM_Res <- fit.mmer(trainedMarkerMat = Markers[1:200,], trainedPheVal = phenotype, 
#'                    predictMarkerMat = Markers[201:242,], method = "EM", 
#'                    effectConcider = "A", iters = 20,
#'                    )
#' }
fit.mmer <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat, Z = NULL, X = NULL, method = "NR", effectConcider = "A", outputModel = FALSE, iters = 20){
require("sommer")
if (!is.matrix(predictMarkerMat)) {
  predictMarkerMat <- t(as.matrix(predictMarkerMat))
}
x <- rbind(trainedMarkerMat, predictMarkerMat)
Y <- c(trainedPheVal, rep(NA,nrow(predictMarkerMat)))
idx <- which(is.na(Y) == T)

if (is.null(Z)){
  effect <- c("A","D","E","AD","AE","DE","ADE")
  if (effectConcider %in% effect){
    Z = diag(dim(x)[1])
    switch (effectConcider,
            A = {A <- A.mat(x);ETA = list(A = list(Z = Z,K = A));ETA <- ETA},
            D =  {D <- D.mat(x);ETA = list(A = list(Z = Z,K = D));ETA <- ETA},
            E = {E <- E.mat(x);ETA = list(A = list(Z = Z,K = E));ETA <- ETA},
            AD = {A <- A.mat(x);D <- D.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D));ETA <- ETA},
            AE = {A <- A.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),E = list(Z = Z,K = E));ETA <- ETA},
            DE = {D <- D.mat(x);E <- E.mat(x);ETA = list(D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA},
            ADE = {A <- A.mat(x);D <- D.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA}
    )
  }else{
    stop("The setting of random effects is wrong!")
  }
} else {
  ETA = Z
}

 if(method %in% c("EMMA","NR")){
  mmerModel <- mmer(X = X, Y = Y, Z = ETA, method = method, iters = iters, silent = T)
 }else if(method %in% c("EM","AI")){
   # mmerModel <- New_EM(X = X,y = Y,ZETA = ETA,R = NULL,iters = iters,silent = TRUE,draw = F)
   stop("AI and EM method has been discontinued because of its instability. Try 'NR' and 'AI'.\nSee details in the G2P help page. \n If you insist AI & EM, please use function fit.EM and fit.AI to debug.")
 }
Res <- mmerModel$fitted.y[idx]
if(outputModel){
  Res <- list(model = mmerModel,predictRes = Res)
}
Res
}

############################################################# New EM ##################################################
########################################## modeling and predicting EM models #######################################
#' @title Modeling and predicting using methods EM model.
#' @description This function can fit EM model and export the prediction values of testing sets. 
#' @param trainedMarkerMat  (numeric)a matrix, each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric)a matrix, each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param Z  A 2-level list,incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. This needs to be provided as a 2-level list structure, defult NULL. 
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param effectConcider  (string)if Z = NULL, random effects are auto generated.
#' @param iters  (numeric)a scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.
#' @param draw  (logical)indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE, default FALSE.
#' @param silent (logical)indicating if the function should draw the progress bar or iterations performed while working or should not be displayed, default FALSE
#' @param constraint (logical)indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way, default FALSE.
#' @param init vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.
#' @param forced a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.
#' @param tolpar tolerance parameter for convergence in the models.
#' @param tolparinv tolerance parameter for matrix inversion in the models.
#' @param outputModel (logical)if true, return the list of training model.
#' @param ... other parameters, details see package sommer.
#' @seealso 
#' \pkg{sommer}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords AI, NR, EM, EMMA, model predict
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' EM_Res <- fit.EM(trainedMarkerMat = Markers[1:200,], trainedPheVal = phenotype, 
#'                    predictMarkerMat = Markers[201:242,], method = "EM", 
#'                    effectConcider = "A", iters = 20,
#'                    )
#' }
fit.EM <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat,effectConcider = "A",outputModel = F,Z = NULL, X = NULL, iters = 20, draw = FALSE, 
                  silent = TRUE, constraint = TRUE, init = NULL, forced = NULL, 
                  tolpar = 1e-04, tolparinv = 1e-06, ...) 
{
  require("sommer")
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  x <- rbind(trainedMarkerMat, predictMarkerMat)
  Y <- c(trainedPheVal, rep(NA,nrow(predictMarkerMat)))
  idx <- which(is.na(Y) == T)
  
  if (is.null(Z)){
    effect <- c("A","D","E","AD","AE","DE","ADE")
    if (effectConcider %in% effect){
      Z = diag(dim(x)[1])
      switch (effectConcider,
              A = {A <- A.mat(x);ETA = list(A = list(Z = Z,K = A));ETA <- ETA},
              D =  {D <- D.mat(x);ETA = list(A = list(Z = Z,K = D));ETA <- ETA},
              E = {E <- E.mat(x);ETA = list(A = list(Z = Z,K = E));ETA <- ETA},
              AD = {A <- A.mat(x);D <- D.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D));ETA <- ETA},
              AE = {A <- A.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),E = list(Z = Z,K = E));ETA <- ETA},
              DE = {D <- D.mat(x);E <- E.mat(x);ETA = list(D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA},
              ADE = {A <- A.mat(x);D <- D.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA}
      )
    }else{
      stop("The setting of random effects is wrong!")
    }
  } else {
    ETA = Z
  }

ZETA = ETA
R = NULL
y <- c(trainedPheVal, rep(NA,nrow(predictMarkerMat)))
convergence <- FALSE
if (is.null(X)) {
  X <- matrix(1, nrow = length(y))
}
if (is.null(R)) {
  R <- list(units = diag(length(y)))
}
y.or <- y
X.or <- X
ZETA.or <- ZETA
R.or <- R
if (is.null(names(ZETA))) {
  varosss <- c(paste("u.", 1:length(ZETA), sep = ""))
}else {
  varosss <- c(names(ZETA))
}
varosssZ <- varosss
if (is.null(names(R))) {
  varosss <- c(varosss, paste("Res.", 1:length(R), sep = ""))
}else {
  varosss <- c(varosss, names(R))
}
good <- which(!is.na(y))
y <- y[good]
X <- as.matrix(X[good, ])
ZETA <- lapply(ZETA, function(x, good) {
  x[[1]] <- x[[1]][good, ]
  x[[2]] <- x[[2]]
  return(x)
}, good = good)
R <- lapply(R, function(x, good) {
  x <- x[good, good]
  return(x)
}, good = good)
qr <- qr(X)
ranx <- length(y) - qr$rank
nz <- length(ZETA)
nr <- length(R)
nx <- dim(X)[2]
dimzs <- lapply(ZETA, function(x) {
  dim(x$Z)[2]
})
dimrs <- lapply(R, function(x) {
  dim(x)[1]
})
N <- length(y)
sta <- 1
ind1 <- numeric()
for (u in 1:length(dimzs)) {
  sta <- dimzs[[u]] + sta
  ind1[u] <- sta
}
ind1 <- c(1, ind1[-c(length(ind1))])
ind2 <- (c(ind1 - 1, sum(unlist(dimzs))))[-1]
ind1g <- ind1 + nx
ind2g <- ind2 + nx
ind1 <- c(1, ind1 + nx)
ind2 <- c(nx, ind2 + nx)
if (is.null(init)) {
  varcom <- rep(var(y, na.rm = TRUE)/(nz + nr), nz + nr)
  names(varcom) <- c(rep("z", nz), rep("r", nr))
  varcomz <- varcom[which(names(varcom) == "z")]
  varcomr <- varcom[which(names(varcom) == "r")]
}else {
  varcom <- init
  if (length(init) != (nr + nz)) {
    stop("Please provide initial values for all variance components", 
         call. = FALSE)
  }else {
    names(varcom) <- c(rep("z", nz), rep("r", nr))
    varcomz <- varcom[which(names(varcom) == "z")]
    varcomr <- varcom[which(names(varcom) == "r")]
  }
}
ll2 = -1e+07
ll.stored <- ll2
conv = 0
wi = 0
taper <- rep(1, iters)
Riw <- lapply(R, function(x) {
  solve(x)
})
Ki <- lapply(ZETA, function(x) {
  findrank <- qr(x$K)$rank
  if (findrank < dim(x$K)[1]) {
    return(solve(x$K + diag(tolparinv, dim(x$K)[1])))
  }else {
    return(solve(x$K))
  }
})
varso <- as.matrix(varcom)
if (!silent) {
  count <- 0
  tot <- iters + 1
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
}
if (!is.null(forced)) {
  varcom <- forced
  if (length(init) != (nr + nz)) {
    stop("Please provide initial values for all variance components", 
         call. = FALSE)
  }
  else {
    names(varcom) <- c(rep("z", nz), rep("r", nr))
    varcomz <- varcom[which(names(varcom) == "z")]
    varcomr <- varcom[which(names(varcom) == "r")]
  }
  iters = 1
}
while (conv == 0) {
  wi = wi + 1
  if (!silent) {
    count <- count + 1
  }
  Rse <- R[[1]] * 0
  for (u in 1:length(R)) {
    Rse <- Rse + R[[u]]
  }
  Rsei <- solve(Rse)
  varoz <- as.list(varcomz)
  Gi <- lapply(as.list(c(1:nz)), function(x, Ko, v) {
    oo = Ko[[x]] * as.numeric(varcomr[1]/(v[[x]]))
    return(oo)
  }, Ko = Ki, v = varoz)
  Zbind <- do.call(cbind, lapply(ZETA, function(x) {
    x$Z
  }))
  XZbind <- as(cbind(X, Zbind), Class = "sparseMatrix")
  CM <- t(XZbind) %*% Rsei %*% XZbind
  for (u in 1:length(Gi)) {
    rox <- ind1g[u]
    cox <- ind2g[u]
    CM[rox:cox, rox:cox] <- CM[rox:cox, rox:cox] + Gi[[u]]
  }
  RHS <- t(XZbind) %*% Rsei %*% y
  ge <- solve(CM, RHS)
  CMi <- solve(CM)
  na <- sum(unlist(dimzs))
  na
  se <- varcomr
  left <- t(XZbind) %*% Rsei %*% y
  top <- t(left)
  corner <- t(y) %*% Rsei %*% y
  M <- rbind(cbind(corner, top), cbind(left, CM))
  M[1:4, 1:4]
  vb <- try(chol(as(M, Class = "sparseMatrix"), pivot = FALSE), 
            silent = TRUE)
  if (class(vb) == "try-error") {
    vb <- try(chol(as(M + diag(tolparinv, dim(M)[1]), 
                      Class = "sparseMatrix")), silent = TRUE)
    if (class(vb) == "try-error") {
      ypy <- 1
      logdc <- 2
    }else {
      L <- t(vb)
      L[1:4, 1:4]
      ypy <- (L[1, 1]^2)
      ypy
      logdc <- 2 * sum(log(diag(L)[-1]))
      logdc
    }
  }else {
    L <- t(vb)
    L[1:4, 1:4]
    ypy <- (L[1, 1]^2)
    ypy
    logdc <- 2 * sum(log(diag(L)[-1]))
    logdc
  }
  logda <- sum(unlist(lapply(ZETA, function(x) {
    determinant(x$K, logarithm = TRUE)$modulus[[1]]
  })))
  nlogsu <- log(varcomz) * unlist(dimzs)
  ll <- -0.5 * (((N - ranx - na) * log(se)) + logdc + logda + 
                  sum(nlogsu) + ypy)
  ll
  if (wi == iters) {
    conv = 1
    #       if (abs(ll - ll2) < tolpar) {
    #         convergence <- TRUE
    #       }
  }
  ll2 = (ll)
  ll.stored <- c(ll.stored, ll2)
  if (!silent) {
    setTxtProgressBar(pb, (count/tot))
  }
  ulist <- list()
  for (k in 1:length(ind1g)) {
    ulist[[k]] <- as.matrix(ge@x[ind1g[k]:ind2g[k]])
  }
  b <- as.matrix(ge@x[1:nx])
  Xb <- (X %*% b)
  Zulist <- list()
  Welist <- list()
  for (k in 1:length(ind1g)) {
    Zulist[[k]] <- ZETA[[k]]$Z %*% ulist[[k]]
    Welist[[k]] <- Zulist[[k]]/varcomz[k]
  }
  zu <- do.call(cbind, Zulist)
  zu <- rowSums(zu)
  now <- 0
  for (f in 1:length(ind1g)) {
    now <- now + crossprod(as.matrix(Zulist[[f]]), y)
  }
  b <- as.matrix(ge@x[ind1[1]:ind2[1]])
  now <- now + crossprod(X %*% b, y)
  varcomr[1] <- ((t(y) %*% as.matrix(Rsei) %*% y) - now)/(length(y) - 
                                                            nx)
  for (k in 1:length(ind1g)) {
    Kinv <- Ki[[k]]
    nk <- nrow(Kinv)
    rox <- ind1g[k]
    cox <- ind2g[k]
    varcomz[k] <- ((t(ulist[[k]]) %*% Kinv %*% ulist[[k]]) + 
                     (as.numeric(varcomr[1]) * sum(diag(Kinv %*% CMi[rox:cox, 
                                                                     rox:cox]))))/nk
  }
  varcom <- c(varcomz, varcomr)
  zeros <- which(varcom <= 0)
  if (length(zeros) > 0 & constraint) {
    varcom[zeros] <- varcomr[1] * (1.011929e-07^2)
  }
  names(varcom) <- c(rep("z", nz), rep("r", nr))
  varcomz <- varcom[which(names(varcom) == "z")]
  varcomr <- varcom[which(names(varcom) == "r")]
  varso <- cbind(varso, as.matrix(varcom))
}
e <- y - (Xb + zu)
We <- (e)/varcomr[1]
B <- cbind(do.call(cbind, Welist), We)
left <- t(XZbind) %*% Rsei %*% B
top <- t(left)
corner <- t(B) %*% Rsei %*% B
M <- rbind(cbind(corner, top), cbind(left, CM))
vb <- try(chol(as(M, Class = "sparseMatrix")), silent = TRUE)
if (class(vb) == "try-error") {
  aii <- diag(nz + nr)
}else {
  L <- t(vb)
  ai <- L[1:(nz + nr), 1:(nz + nr)]^2
  ai
  ai <- as.matrix(ai)
  ai[upper.tri(ai)] <- t(ai)[upper.tri(ai)]
  aii <- solve(ai)
}
if (draw) {
  layout(matrix(1:2, 2, 1))
  plot(y = ll.stored[-1], x = 1:length(ll.stored[-1]), 
       type = "l", main = "Log-likelihood", xlab = "iteration", 
       ylab = "Log-likelihood")
  for (l in 1:dim(varso)[1]) {
    ub <- max(unlist(varso), na.rm = TRUE)
    if (l == 1) {
      plot(varso[1, ], ylim = c(0, ub), type = "l", 
           main = "MME-EM results", ylab = "varcomp", 
           xlab = "iteration")
    }else {
      lines(varso[l, ], x = 1:dim(varso)[2], col = l)
    }
  }
}
out1 <- as.matrix(varcom)
rownames(out1) <- varosss
colnames(out1) <- "component"
Vinv <- XZbind %*% CMi %*% t(XZbind)
names(ulist) <- varosssZ
for (f in 1:length(ZETA)) {
  rownames(ulist[[f]]) <- colnames(ZETA[[f]]$Z)
}
Var.u <- list()
for (f in 1:length(ind1g)) {
  Var.u[[f]] <- diag(CMi[ind1g[f]:ind2g[f], ind1g[f]:ind2g[f]])
  names(Var.u[[f]]) <- colnames(ZETA[[f]]$Z)
}
rownames(b) <- colnames(X)
xvxi <- CMi[ind1[1]:ind2[1], ind1[1]:ind2[1]]
rownames(xvxi) <- colnames(xvxi) <- colnames(X)
ee <- y - (Xb)
ll <- as.vector(ll)
AIC = as.vector((-2 * ll) + (2 * nx))
BIC = as.vector((-2 * ll) + (log(length(y)) * nx))
Ksp <- do.call(adiag1, lapply(ZETA, function(x) {
  x$K
}))
fit0 <- Xb + zu
Xor.b <- X.or %*% b
Zor.u <- lapply(as.list(1:nz), function(x, z, u) {
  z[[x]]$Z %*% ulist[[x]]
}, z = ZETA.or, u = ulist)
Zor.u2 <- do.call(cbind, Zor.u)
Zor.u2 <- rowSums(Zor.u2)
fit1 <- Xor.b + Zor.u2
if (!silent) {
  setTxtProgressBar(pb, (tot/tot))
}
out1 <- as.data.frame(out1)
out1[, "constraint"] <- "Positive"
if (length(zeros) > 0) {
  out1[zeros, "constraint"] <- "Boundary"
  out1[zeros, 1] <- 0
}
sigma.scaled <- out1[, 1]/var(y, na.rm = TRUE)
res <- list(var.comp = out1, V.inv = Vinv, u.hat = ulist, 
            Var.u.hat = Var.u, beta.hat = b, Var.beta.hat = xvxi, 
            residuals = ee, cond.residuals = e, LL = ll, sigma.scaled = sigma.scaled, 
            AIC = AIC, BIC = BIC, fish.inv = aii, fitted.y.good = fit0, 
            X = X.or, Z = Zbind, K = Ksp, ZETA = ZETA, fitted.y = fit1, 
            fitted.u = Zor.u2, forced = forced, convergence = convergence)
layout(matrix(1, 1, 1))
Res <- res$fitted.y[idx]
if(outputModel){
  Res <- list(model = res,predictRes = Res)
}
Res
}

########################################## modeling and predicting AI models #######################################
#' @title Modeling and predicting using methods AI model.
#' @description This function can fit AI model and export the prediction values of testing sets. 
#' @param trainedMarkerMat  (numeric)a matrix, each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainedPheVal  (numeric)the phenotype Value of each individual.
#' @param predictMarkerMat  (numeric)a matrix, each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param Z  A 2-level list,incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. This needs to be provided as a 2-level list structure, defult NULL. 
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param effectConcider  (string)if Z = NULL, random effects are auto generated.
#' @param iters  (numeric)a scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.
#' @param draw  (logical)indicating if a plot of updated values for the variance components and the likelihood should be drawn or not. The default is TRUE. COMPUTATION TIME IS SMALLER IF YOU DON'T PLOT SETTING draw=FALSE, default FALSE.
#' @param silent (logical)indicating if the function should draw the progress bar or iterations performed while working or should not be displayed, default FALSE
#' @param constraint (logical)indicating if the program should use the boundary constraint when one or more variance component is close to the zero boundary. The default is TRUE but needs to be used carefully. It works ideally when few variance components are close to the boundary but when there are too many variance components close to zero we highly recommend setting this parameter to FALSE since is more likely to get the right value of the variance components in this way, default FALSE.
#' @param init vector of initial values for the variance components. By default this is NULL and variance components are estimated by the method selected, but in case the user want to provide initial values this argument is functional.
#' @param forced a vector of numeric values for variance components including error if the user wants to force the values of the variance components. On the meantime only works for forcing all of them and not a subset of them. The default is NULL, meaning that variance components will be estimated by REML/ML.
#' @param tolpar tolerance parameter for convergence in the models.
#' @param tolparinv tolerance parameter for matrix inversion in the models.
#' @param outputModel (logical)if true, return the list of training model.
#' @param ... other parameters, details see package sommer.
#' @seealso 
#' \pkg{sommer}
#' @return  
#' a list including model and prediction result (outputModel = TRUE)
#' a array of prediction result (outputModel = FALSE)
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords AI, NR, EM, EMMA, model predict
#' @export
#' @examples
#' \dontrun{
#' ## Load example data ##
#' data(GYSS)
#'
#' ## use RR model to modeling and predict ##
#' EM_Res <- fit.AI(trainedMarkerMat = Markers[1:200,], trainedPheVal = phenotype, 
#'                 predictMarkerMat = Markers[201:242,], method = "EM", 
#'                 effectConcider = "A", iters = 20,
#'                 )
#' }

fit.AI <- function(trainedMarkerMat, trainedPheVal, predictMarkerMat,effectConcider = "A",outputModel = F,
               X=NULL,Z=NULL,iters=20,draw=FALSE,silent=FALSE, constraint=TRUE, init=NULL, forced=NULL,
               tolpar = 1e-04, tolparinv = 1e-06){
  require("sommer")
  trainedMarkerMat <- Markers[1:200,]
  predictMarkerMat <- Markers[201:242,]
  trainedPheVal <- phenotype[1:200]
  
  if (!is.matrix(predictMarkerMat)) {
    predictMarkerMat <- t(as.matrix(predictMarkerMat))
  }
  x <- rbind(trainedMarkerMat, predictMarkerMat)
  Y <- c(trainedPheVal, rep(NA,nrow(predictMarkerMat)))
  idx <- which(is.na(Y) == T)
  
  Z = NULL
  effectConcider = "A"
  if (is.null(Z)){
    effect <- c("A","D","E","AD","AE","DE","ADE")
    if (effectConcider %in% effect){
      Z = diag(dim(x)[1])
      switch (effectConcider,
              A = {A <- A.mat(x);ETA = list(A = list(Z = Z,K = A));ETA <- ETA},
              D =  {D <- D.mat(x);ETA = list(A = list(Z = Z,K = D));ETA <- ETA},
              E = {E <- E.mat(x);ETA = list(A = list(Z = Z,K = E));ETA <- ETA},
              AD = {A <- A.mat(x);D <- D.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D));ETA <- ETA},
              AE = {A <- A.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),E = list(Z = Z,K = E));ETA <- ETA},
              DE = {D <- D.mat(x);E <- E.mat(x);ETA = list(D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA},
              ADE = {A <- A.mat(x);D <- D.mat(x);E <- E.mat(x);ETA = list(A = list(Z = Z,K = A),D = list(Z = Z,K = D),E = list(Z = Z,K = E));ETA <- ETA}
      )
    }else{
      stop("The setting of random effects is wrong!")
    }
  } else {
    ETA = Z
  }
  
  ZETA = ETA
  R = NULL
  y <- c(trainedPheVal, rep(NA,nrow(predictMarkerMat)))
  convergence <- FALSE
  
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### fix that some initially were not specified
  if(is.null(X)){
    X <- matrix(1,nrow=length(y))
  }
  if(is.null(R)){
    R <- list(units=diag(length(y)))
  }
  y.or <- y
  X.or <- X
  ZETA.or <- ZETA
  R.or <- R
  if(is.null(names(ZETA))){
    varosss <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA))
  }; varosssZ <- varosss
  if(is.null(names(R))){
    varosss <- c(varosss,paste("Res.",1:length(R),sep=""))
  }else{
    varosss <- c(varosss,names(R))
  }
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  
  #### now reduce the original inputs
  #y <- scale(y)
  good <- which(!is.na(y))
  y <- y[good]
  X <- as.matrix(X[good,])
  ZETA <- lapply(ZETA, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
  R <- lapply(R, function(x,good){x <- x[good,good]; return(x)}, good=good)
  #### get sizes of reduced
  qr <- qr(X)
  ranx <- qr$rank # length(good)
  nz <- length(ZETA)
  nr <- length(R)
  nx <- dim(X)[2]
  dimzs <- lapply(ZETA,function(x){dim(x$Z)[2]})
  dimrs <- lapply(R,function(x){dim(x)[1]})
  N <- length(y)
  
  #### get the indexes for all effects
  sta <- 1
  ind1 <- numeric()
  for(u in 1:length(dimzs)){
    sta <- dimzs[[u]]+sta
    ind1[u] <- sta
  }; ind1<- c(1,ind1[-c(length(ind1))]);
  ind2 <- (c(ind1-1, sum(unlist(dimzs))))[-1]
  
  ind1g <- ind1 + nx ## indexes for Gi
  ind2g <- ind2 + nx ## indexes for Gi
  ind1 <- c(1,ind1 + nx) ## indexes for all including x
  ind2 <- c(nx,ind2 + nx) ## indexes for all including x
  
  #### initialize var components
  if(is.null(init)){
    varcom <- rep(var(y,na.rm=TRUE)/(nz+nr),nz+nr)
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
  }else{
    varcom <- init 
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
  }
  
  ll2=-10000000 # initial log-likelihood
  ll.stored <- ll2
  conv=0 # convergence
  wi=0 # counter
  taper <- rep(1, iters) # weighting parameter for updates
  #taper[1:2] <- c(.5,.7)#c(0.5, 0.7) # c(0.5, 0.7)
  
  #### get inverse of covariance and residual matrices
  Riw <-  lapply(R,function(x){solve(x)})
  Ki <- lapply(ZETA,function(x){
    findrank <- qr(x$K)$rank
    if(findrank < dim(x$K)[1]){
      return(solve(x$K + diag(tolparinv,dim(x$K)[1])))
    }else{
      return(solve(x$K))
    }
  })
  
  varso <- as.matrix(varcom)
  #### start
  if(!silent){
    count <- 0
    tot <- iters+1
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!is.null(forced)){
    varcom <- forced
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
    iters=1
  }
  
  ups <- as.matrix(varcom)*0
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  while (conv==0) { # ==================== START AI ALGORITHM =========================
    wi=wi+1
    if(!silent){
      count <- count + 1
    }
    ##############################
    ### for R * se as a direct sum
    ##############################
    ### for R * se as a direct sum
    Rse <- R[[1]]*0
    for(u in 1:length(R)){
      Rse <- Rse + R[[u]]#*(1/varcomr[u]))
    }
    ##############################
    ## R inverse
    ## R inverse
    varor <- as.list(varcomr) 
    Rsei <- lapply(as.list(c(1:nr)),function(x,K,v){
      oo=K[[x]]*as.numeric(1/(v[[x]])); return(oo)
    }, K=Riw, v=varor) ## K*v(u)
    
    Rsei <- Reduce('+', Rsei)
    #Rsei <- solve(Rse) *(1/varcomr[u])
    ##############################
    ## G inverse
    se <- varcomr[1]
    varoz <- as.list(varcomz)  # variance components as list, no error included
    Gi <- lapply(as.list(c(1:nz)),function(x,Ko,v){
      oo=Ko[[x]]*as.numeric(1/(v[[x]])); return(oo)
    }, Ko=Ki, v=varoz) ## K*v(u)
    
    Zbind <- do.call(cbind,lapply(ZETA, function(x){x$Z}))
    XZbind <- as(cbind(X,Zbind), Class="sparseMatrix")
    CM <- t(XZbind) %*% Rsei %*% XZbind
    
    ## add Gi's to the C matrix
    for(u in 1:length(Gi)){
      rox <- ind1g[u]; cox <- ind2g[u]
      CM[rox:cox,rox:cox] <- CM[rox:cox,rox:cox] + Gi[[u]]
    }
    
    ## do gaussian eliminations (absorption)
    RHS <- t(XZbind) %*% Rsei %*% y
    ge <- solve(CM,RHS)
    
    #### Calculate ypy by building M
    left <- t(XZbind) %*% Rsei %*% y
    top <- t(left)
    corner <- t(y)%*%Rsei%*%y
    
    M <- rbind(cbind(corner,top),cbind(left,CM))
    M[1:4,1:4]
    vb <- chol(as(M, Class="sparseMatrix"))
    L <- t(vb); L[1:4,1:4]
    ypy <- (L[1,1]^2); ypy
    logdc <- 2* sum(log(diag(L)[-1])); logdc
    logda <- sum(unlist(lapply(ZETA,function(x){determinant(x$K, logarithm = TRUE)$modulus[[1]]})))
    nlogsu <- log(varcomz)*unlist(dimzs) # n*log(su)
    ##%%%%%%%%%%%%%%
    ### likelihood
    ##%%%%%%%%%%%%%%
    na <- sum(unlist(dimzs));na
    se <- varcomr
    ll <- - 0.5 * ( ((N - ranx - na)*log(se)) + logdc + logda + sum(nlogsu) +  ypy);ll
    if (abs(ll-ll2) < tolpar | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
      conv=1
      if (abs(ll-ll2) < tolpar){
        convergence <- TRUE
      }
    }
    ll2=(ll) # replace the initial logLik value for the value reached
    ll.stored <- c(ll.stored,ll2)
    ##%%%%%%%%%%%%%%
    ##%%%%%%%%%%%%%%
    M <- NULL; vb <- NULL
    
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    
    ##%%%%%%%%%%%%%%
    ### fill AImatrix
    ##%%%%%%%%%%%%%%
    ai <- matrix(NA,nz+nr,nz+nr)
    
    ## extract BLUPs and BLUEs
    ulist <-list()
    for(k in 1:length(ind1g)){
      ulist[[k]] <- as.matrix(ge@x[ind1g[k]:ind2g[k]])
    }
    b <- as.matrix(ge@x[1:nx])
    Xb <- (X%*%b)
    ## get fitted Xb and Zu, 
    ## and working variates HiPy
    Zulist <- list()
    Welist <- list()
    for(k in 1:length(ind1g)){
      Zulist[[k]] <- ZETA[[k]]$Z %*% ulist[[k]]
      Welist[[k]] <- Zulist[[k]]/varcomz[k]
    }
    ## add all Zu's
    zu <- do.call(cbind,Zulist); zu <- rowSums(zu)
    ## get residuals
    e <- y - (Xb+zu)
    ## get working variate for residuals
    We <- (e)/varcomr[1]; #length(ye)
    ## matrix of working variates
    B <- cbind(do.call(cbind,Welist),We)
    
    left <- t(XZbind) %*% Rsei %*% B #left of the M matrix
    top <- t(left) # top of M
    corner <- t(B)%*%Rsei%*%B # left corner of M
    
    M <- rbind(cbind(corner,top),cbind(left,CM))
    #M[1:5,1:5]
    
    vb <- chol(as(M, Class="sparseMatrix"))
    L <- t(vb); #L[1:4,1:4]
    ai <- L[1:(nz+nr),1:(nz+nr)];ai
    
    ai <- as.matrix(ai)#*.5 #originally should not be multiplied
    ai[upper.tri(ai)] <- t(ai)[upper.tri(ai)]
    ai <- ai^2
    
    #print(ai)
    aii <- solve(ai)
    
    ##%%%%%%%%%%%%%%
    ### first derivatives dL/ds
    ##%%%%%%%%%%%%%%
    
    CMi <- solve(CM)
    
    fdlist <- list()
    delist <- list() # store the pieces to substract in the dL/dse
    for(i in 1:(nz)){
      
      caa <- CMi[ind1g[i]:ind2g[i],ind1g[i]:ind2g[i]]
      fac <- sum(diag(Ki[[i]]%*%caa)) # to use in dL/dse as well
      fdlist[[i]] <- -0.5*( (dimzs[[i]]/varcomz[i]) - 
                              (fac/(varcomz[i]^2)) - 
                              (((t(e)/varcomr[1])%*%(Zulist[[i]]/varcomz[i])))  ) # pay attention to this "se"
      # for first derivatives
      delist[[i]] <- dimzs[[i]] - (fac/(varcomz[i])) # Na - tr(A%*%Caa)/su
    }
    
    xprov <- lapply(delist, function(x){x*(1/varcomr[1])}) # (Na - tr(A%*%Caa)/su) * (1/se)
    
    fd <- ( ((N-ranx)/varcomr[1]) - (((t(e)/varcomr[1])%*%(e/varcomr[1])))  )
    for(o in 1:length(xprov)){
      fd <- fd - xprov[[o]]  
    }# this is the derivative for error
    fd <- -0.5* (fd)
    
    FD <- c(unlist(fdlist),fd) # all first derivatives
    
    up <- (aii%*%FD)
    
    ups <- cbind(ups,up)
    
    ######## reduced updates issue
    #     if(wi>0){
    #       reduced <- which(scale(up)>1)
    #       if(length(reduced)>0){
    #         if(wi==1){
    #           up[reduced] <- 0
    #         }else{
    #           up[reduced] <- ups[reduced,(wi-1)]*.01
    #         }
    #         
    #       }
    #     }
    ##########
    #print(up)
    
    varcom <- as.vector(c(varcomz,varcomr) + (taper[wi]*as.matrix(up)))
    
    zeros <- which(varcom <= 0)
    
    if(length(zeros)>0 & constraint){
      varcom[zeros] <- varcomr[1] * (1.011929e-07^2)
    }
    ## move back to the z and r separation
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
    #y <- Xb+zu
    varso <- cbind(varso,as.matrix(varcom))# just to keep track
  } ################# ======================== END OF ALGORITHM =======================
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  
  
  ## plot likelihood
  if(draw){
    layout(matrix(1:2,2,1))
    plot(y=ll.stored[-1], x=1:length(ll.stored[-1]), type="l", main="Log-likelihood", xlab="iteration",ylab="Log-likelihood")
    ## plot varcomp
    for(l in 1:dim(varso)[1]){
      #lb <- min(unlist(varso),na.rm = TRUE)
      ub <- max(unlist(varso),na.rm = TRUE)
      if(l==1){plot(varso[1,],ylim=c(0,ub),type="l", main="MME-Average Information results",ylab="varcomp", xlab="iteration")}else{
        lines(varso[l,],x=1:dim(varso)[2], col=l)
      }
    }
  }
  ## ======================== ##
  ## ======================== ##
  ## PROVIDE EXTRA PARAMETERS
  ## ======================== ##
  ## ======================== ##
  
  ## variance components
  out1 <- as.matrix(varcom); rownames(out1) <- varosss; colnames(out1) <- "component"
  ## inverse of the phenotypic variance (ZKZ'+R)-
  Vinv <- XZbind %*% CMi %*% t(XZbind)
  ## BLUPs
  names(ulist) <- varosssZ
  for(f in 1:length(ZETA)){
    rownames(ulist[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  ## VarBLUPs
  Var.u <- list()
  for(f in 1:length(ind1g)){
    #Var.u[[f]] <- diag(CMi[ind1g[f]:ind2g[f],ind1g[f]:ind2g[f]])
    Var.u[[f]] <- ZETA[[f]]$K * varcom[f]#diag(CMi[ind1g[f]:ind2g[f],ind1g[f]:ind2g[f]])
    names(Var.u[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  names(Var.u) <- varosss[1:length(ZETA)]
  ## PEV.BLUPs
  Pev.u <- list()
  for(f in 1:length(ind1g)){
    Pev.u[[f]] <- diag(CMi[ind1g[f]:ind2g[f],ind1g[f]:ind2g[f]])
    names(Pev.u[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  names(Pev.u) <- varosss[1:length(ZETA)]
  ## betas
  rownames(b) <- colnames(X)
  ## var.betas
  xvxi <- CMi[ind1[1]:ind2[1],ind1[1]:ind2[1]]
  rownames(xvxi) <- colnames(xvxi) <- colnames(X)
  ## cond. residuals
  #e
  ## residuals
  ee <- y - (Xb)
  ## log likelihood
  ll <- as.vector(ll)
  ## AIC and BIC
  AIC = as.vector((-2 * ll ) + ( 2 * nx))
  BIC = as.vector((-2 * ll ) + ( log(length(y)) * nx))
  ## fish.inv
  # aii
  ## Ksp
  Ksp <- do.call(adiag1,lapply(ZETA,function(x){x$K}))
  ## fitted values only for non-missing data
  fit0 <- Xb + zu
  
  ## ======================== ##
  ## ======================== ##
  ## 2. PROVIDE EXTRA PARAMETERS
  ## using original data
  ## ======================== ##
  ## ======================== ##
  Xor.b <- X.or%*%b 
  Zor.u <-  lapply(as.list(1:nz),function(x,z,u){z[[x]]$Z %*% ulist[[x]]},z=ZETA.or,u=ulist)
  Zor.u2 <- do.call(cbind,Zor.u); Zor.u2 <- rowSums(Zor.u2)
  fit1 <- Xor.b + Zor.u2 # fitted values
  
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  
  ### let user knows the constraints
  out1 <- as.data.frame(out1)
  out1[,"constraint"] <- "Positive"
  if(length(zeros)>0){
    out1[zeros,"constraint"] <- "Boundary"
    out1[zeros,1] <- 0
  }
  res <- list(var.comp=out1, V.inv = Vinv, u.hat=ulist, Var.u.hat=Var.u, 
              PEV.u.hat=Pev.u, 
              beta.hat=b, Var.beta.hat=xvxi, residuals=ee, cond.residuals=e,
              LL=ll, AIC=AIC, BIC=BIC, fish.inv=aii,fitted.y.good=fit0, 
              X=X.or, Z=Zbind, K=Ksp, ZETA=ZETA,
              ## go back to original data
              fitted.y=fit1, fitted.u=Zor.u2, 
              forced=forced, convergence=convergence)
  
  layout(matrix(1,1,1))
  Res <- res$fitted.y[idx]
  if(outputModel){
    Res <- list(model = res,predictRes = Res)
  }
  Res
}
####################################################################################################################################
############################################################# parallel #############################################################
#' @title Genotype to Phenotype
#' @description  this function is to predict phenotype from genotype.
#' @param trainMarker  (numeric,matrix)each row is the each training sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param trainPheno  (numeric)the phenotype value of each individual.
#' @param testMarker  (numeric,matrix)each row is the each testing sets individual's SNP genotypes informations. Genotypes should be coded as {0,1,2}; 0 represent AA(homozygote), 2 represent BB(homozygote) and 1 represent AB(heterozygote); missing (NA) alleles are not allowed.
#' @param testPheno  (numeric)the phenotype value of test population individual, default NULL.
#' @param modelMethods  the model to fit."BayesA", "BayesB", "BayesC", "BL", "BRR","RKHS","rrBLUP","LASSO","SPLS","SVC","SVR","RFR","RFC", "RR", "RKHS", "BRNN", "EM", "EMMA", "AI", "NR".
#' @param nIter,burnIn,thin  (integer) the number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.
#' @param saveAt  (string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default ""
#' @param S0,df0  (numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Default S0 NULL,df0 5.
#' @param R2  (numeric, (0,1)) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Defult 0.5
#' @param weights  (numeric, n) a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.
#' @param verbose  (logical) if TRUE the iteration history is printed, default FALSE
#' @param rmExistingFiles  (logical) if TRUE removes existing output files from previous runs, default TRUE.
#' @param groups  (factor) a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.
#' @param ntree  RandomForest parameter:Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.Defualt 500
#' @param nodesize Randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param importance  RandomForest parameter:Should importance of predictors be assessed?Defualt FALSE
#' @param posPercentage  (numeric)the percentage positive samples in all samples.1 > posPercentage > 0.
#' @param BestIndividuals  It is a position that the best individuals (positive samples) in a training group, according to the breeding values of a training group's trait.
#'                          if the trait was yield,flowering or disease resistance,and  male flowering time to female flowering time,it is "top"(default), "buttom",and "middle" of the breeding values, respectively.
#' @param kernel  svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial)Default "linear".
#' @param gamma  svm parameter parameter needed for all kernels except linear (default: 1/(data dimension))
#' @param cost  svm cost,default 2^(-9)
#' @param outputModel  (logical)if true, return the list of training model.
#' @param K  Number of hidden components
#' @param eta	 Thresholding parameter. eta should be between 0 and 1.
#' @param select  PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2"
#' @param fit	 PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".
#' @param scale.x	 Scale predictors by dividing each predictor variable by its sample standard deviation?
#' @param scale.y  Scale responses by dividing each response variable by its sample standard deviation?
#' @param eps  An effective zero. Default is 1e-4
#' @param maxstep  Maximum number of iterations when fitting direction vectors. Default is 100.
#' @param trace  Print out the progress of variable selection?
#' @param alpha  The elasticnet mixing parameter.Detail in glmnet.
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param family the distribution family of y, see help('family') for more details. 
#' @param lambda the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.
#' @param tol.err internal tolerance level for extremely small values; default value is 1e-6.
#' @param tol.conv tolerance level in convergence; default value is 1e-8.
#' @param neurons  (integer)indicates the number of neurons,defult 4.
#' @param epochs  (integer)maximum number of epochs(iterations) to train, default 30.
#' @param Z  (2-level list)incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. This needs to be provided as a 2-level list structure, defult NULL. 
#' @param effectConcider  (string)if Z = NULL, random effects are auto generated.
#' @param iters  (numeric)a scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.
#' @param ...  arguments passed to or from other methods.
#' @seealso 
#' \link[G2P]{GSmachine}
#' \link[G2P]{GSReModel}
#' \link[G2P]{fit.RR}
#' \link[G2P]{fit.BRNN}
#' \link[G2P]{fit.RKHS}
#' \link[G2P]{fit.mmer}
#' \link[G2P]{predictGS}
#' @return 
#' a matrix:
#' The prediction results of multi-methods
#' 
#' @author Chuang Ma ,Qian Cheng
#' @keywords genomic prediction
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ########## predicting breeding value
#' predRes <- G2P(Markers[1:200,],phenotype[1:200],Markers[201:242,],
#'                 phenotype[201:242],modelMethods = c("rrBLUP", "RFC"),
#'                 outputModel = FALSE)
#' }
G2P <- function(trainMarker,trainPheno,testMarker,testPheno = NULL,modelMethods ="BayesA",outputModel = FALSE,  # main parameters
                nIter = 1500, burnIn = 500, thin = 5, 
                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                posPercentage = 0.3,BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9),  # machine learing parameters
                K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,trace = FALSE,maxstep = 100, # SPLS parameters
                alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 20,
                ...){   # LASSO parameters
  if(is.null(testPheno)){
    ntestSample <- nrow(testMarker)
    if (is.null(ntestSample)) {
      testPheno <- rep(0,1)
    }else{
      testPheno <- rep(0,ntestSample)
    }
  }
  
  if(length(intersect(modelMethods,c("AI","EM"))) != 0){
   modelMethods <- setdiff(modelMethods,c("AI","EM"))
   cat("AI and EM method has been discontinued because of its instability. Try 'NR' and 'AI'.\nSee details in the G2P help page. \n If you insist 'AI' & 'EM', please use function fit.EM and fit.AI to debug.")
  }
  GSReMethods <- c("rrBLUP","LASSO","SPLS","bigRR") 
  MachineMethods <- c("SVC","RFC","SVR","RFR")
  # BGLRMethods <- c("BayesA", "BayesB", "BayesC", "BL", "BRR")
  aloneMethods <- c("BayesA", "BayesB", "BayesC", "BL", "BRR","RR","BRNN","RKHS","AI","NR","EM","EMMA")
  NamesMethods <- NULL
  
  ########################## building regessoion model for GS
  if (length(intersect(modelMethods,GSReMethods)) != 0){
    GSReMethods  <- intersect(modelMethods,GSReMethods)
    NamesMethods <- GSReMethods 
    trainedModelList <- lapply(GSReMethods,function(ii){GSReModel(modelMethods = ii,markers = trainMarker,pheVal = trainPheno,
                                                                  K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
                                                                  alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep
                                                                  , ...)})
  }else{
    NamesMethods <- NULL
    trainedModelList <- NULL
  }
  
  ###############################3  building  machine model for GS
  if (length(intersect(modelMethods,MachineMethods)) != 0){
    MachineMethods <- intersect(modelMethods,MachineMethods)
    NamesMethods <- c(NamesMethods,MachineMethods)
    machinemodel <- lapply(MachineMethods,function(ii){GSmachine(markers = trainMarker,pheVal = trainPheno,posPercentage = posPercentage ,BestIndividuals = BestIndividuals,
                                                                 modelMethods = ii ,ntree = ntree ,nodesize = nodesize,
                                                                 kernel = kernel,gamma = gamma, cost = cost)})
    trainedModelList <- c(trainedModelList,machinemodel)
    
  } else {
    NamesMethods <- NamesMethods
    trainedModelList <- trainedModelList
  } 
  
  if(length(intersect(modelMethods,aloneMethods)) != 0){
    aloneMethods <- intersect(modelMethods,aloneMethods)
    #    NamesMethods <- c(NamesMethods,aloneMethods)
    aloneMethodsResList <- lapply(aloneMethods, function(ii){switch(ii,
                                                                    RR = fit.RR(trainMarker,trainPheno,testMarker,outputModel = TRUE),
                                                                    BRNN = fit.BRNN(trainMarker,trainPheno,testMarker,outputModel = TRUE,verbose = verbose, neurons = neurons, epochs = epochs, ...),
                                                                    RKHS = fit.RKHS(trainMarker,trainPheno,testMarker,outputModel = TRUE, nIter = nIter, burnIn = burnIn,thin = thin,
                                                                                    saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
                                                                    #AI = fit.mmer(trainMarker,trainPheno,testMarker,outputModel = TRUE,method = "AI",Z = Z, effectConcider = effectConcider,iters = mmerIters),
                                                                    NR = fit.mmer(trainMarker,trainPheno,testMarker,outputModel = TRUE,method = "NR",Z = Z, effectConcider = effectConcider,iters = mmerIters),
                                                                    #EM = fit.EM(trainMarker,trainPheno,testMarker,outputModel = TRUE,method = "EM",Z = Z, effectConcider = effectConcider,iters = mmerIters),
                                                                    EMMA = fit.mmer(trainMarker,trainPheno,testMarker,outputModel = TRUE,method = "EMMA",Z = Z, effectConcider = effectConcider,iters = mmerIters),
                                                                    BayesA = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesA",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                                                                                      saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
                                                                    BayesB = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesB",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                                                                                      saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
                                                                    BayesC = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BayesC",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                                                                                      saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
                                                                    BL = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BL",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                                                                                  saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups),
                                                                    BRR = fit.BGLR(trainMarker, trainPheno, testMarker, modelMethods = "BRR",outputModel = TRUE,nIter = nIter, burnIn = burnIn, thin = thin,
                                                                                   saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups)
    )
    })
    #     aloneMethodsList <- list()
    #     for (i in 1:length(aloneMethodsResList)) {
    #       aloneMethodsList <- c(aloneMethodsList,)
    #       
    #     }
    mmerModelList <- lapply(aloneMethodsResList, function(x){a <- x[[1]]})
    alonePredResList <- lapply(aloneMethodsResList, function(x){b <- x[[2]]})
    lengtList <- length(alonePredResList)
    alonepredResMat <- as.matrix(testPheno)
    for(i in 1:lengtList){
      alonepredResMat <- cbind(alonepredResMat,as.matrix(alonePredResList[[i]]))
    }
    colnames(alonepredResMat) <- c("realPhenScore",aloneMethods)
    names(mmerModelList) <- aloneMethods
    aloneRes <- list(model = mmerModelList, predScore = alonepredResMat)
  }
  ######### predict breeding value of test sample
  
  if(!is.null(NamesMethods)){
    names(trainedModelList) <- NamesMethods
    predscores <- sapply(NamesMethods,function(ii){predictGS(testMat = testMarker,trainModel = trainedModelList[[ii]],modelMethods = ii )})
    if(!is.matrix(testMarker)){
      predscores <- matrix(c(testPheno,predscores),nrow = 1)
      colnames(predscores) <-  c("realPhenScore",NamesMethods)
      rownames(predscores) <- NULL
    }else{
      predscores <- cbind(testPheno,predscores)
      colnames(predscores) <-  c("realPhenScore",NamesMethods)
    }
    if (length(intersect(modelMethods,aloneMethods)) != 0) {
      predscores <- cbind(predscores,aloneRes[[2]][,-1])
      colnames(predscores) <- c("realPhenScore",NamesMethods,aloneMethods)
      rownames(predscores) <- rownames(testMarker)
      trainedModelList <- c(trainedModelList,aloneRes[[1]])
    }
  }else{
    predscores <- aloneRes[[2]]
    trainedModelList <- aloneRes[[1]]
  }
  
  predscores <- predscores[,c("realPhenScore",modelMethods)]
  result <- list()
  
  ####### output result
  if(outputModel){
    result[["trainModel"]] <- trainedModelList
    result[["predscores"]] <- predscores 
    rm(trainMarker,testMarker)
    return(result)
  }
  else{
    return(predscores)
  }
}
################################################ eval ###############################################
meanNDCG <- function( realScores, predScores, topK = 10 ){
  resVec <- rep(0, topK )
  for( idx in 1:topK ){
    resVec[idx] <- NDCG( realScores, predScores, topK = idx )
  }
  meanNDCG <- mean(resVec)
  names(meanNDCG ) <- paste0("meanNDCG","_top",topK)
  return (meanNDCG)
}
##from plos one, 2015, A Ranking Approach to Genomic Selection
NDCG <- function( realScores, predScores, topK = 10){
  
  if( length(realScores) != length(predScores)){
    stop("Error: different length between realScores and predScores")
  }
  if( length(realScores) < topK ) {
    stop("Error: too large topK")
  }
  
  scoreMat <- cbind(realScores,predScores)
  scoreMatSortbyPred <- scoreMat[order(scoreMat[,2],decreasing = TRUE),]
  scoreMatSortByReal <- scoreMat[order(scoreMat[,1],decreasing = TRUE),]
  
  DCG <- rep(0, topK)
  IDCG <- rep(0, topK)
  for(idx in 1:topK){
    DCG[idx] <-  scoreMatSortbyPred[idx,1]/log(idx+1,2)
    IDCG[idx] <- scoreMatSortByReal[idx,1]/log(idx+1,2)
  }
  
  NDCG <- sum(DCG)/sum(IDCG) 
  names(NDCG) <- paste0("NDCG","_top",topK)
  return(NDCG)
}
##evaluation method: pearson, spearman, kendall, MSE
corEvaluation <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),BestIndividuals,Probability = FALSE){
  # Probability handle
  if (Probability) {
    if(BestIndividuals == "top"){
      realScores <- realScores
      predScores <- predScores
    }else if(BestIndividuals == "buttom"){
      realScores <- realScores
      predScores <- 1 - predScores
    }else if(BestIndividuals == "middle"){
      realScores <- abs(realScores)
      predScores <- 1 - predScores
    }
  }else{
    realScores <- realScores
    predScores <- predScores
  }
  
  if( length(method) > 1 ){
    method <- method[0]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  res <- ""
  if( (method == "pearson") | (method == "kendall") | (method == "spearman") ){
    res <- cor( realScores, predScores,  use="complete", method = method  )
  }else if(method == "MSE") {
    res <- apply(predScores,2,function(ii){
      deltaVec <- abs( realScores - ii)
      deltaVec <- deltaVec^2
      mean(deltaVec)})
  }else if(method == "R2"){
    res <-apply(predScores,2,function(ii){
      R2 <- summary(lm(realScores ~ ii))$r.squared })
  }
  
  res <- matrix(res,nrow = 1,dimnames = list(method,colnames(predScores)))
  res
}
####################
## evaluation method alpha : pearson, spearman, kendall, MSE
corEvaluationAlpha <- function( realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha,BestIndividuals){
  topNum <- 1: round(length(realScores)*(topAlpha/100))
  if( length(method) > 1 ){
    method <- method[0]
  }
  checkMethodType <- method %in% c("pearson", "kendall", "spearman", "MSE","R2")
  if( !checkMethodType ){
    stop("Error: undefined method in corEvaluation")
  }
  realScores <- as.matrix(realScores)
  predScores <- as.matrix(predScores)
  mat <- cbind(realScores,predScores)
  
  if(BestIndividuals == "top"){
    mat <- mat[order(realScores,decreasing = T),]
  }else if(BestIndividuals == "buttom"){
    mat <- mat[order(realScores,decreasing = F),]
  }else if(BestIndividuals == "middle"){
    mat <- mat[order(abs(realScores),decreasing = F),]
  }
  
  mat <- mat[topNum,]
  
  realScores <- as.matrix(mat[,1])
  predScores <- as.matrix(mat[,2])
  
  res <- ""
  if( (method == "pearson") | (method == "kendall") | (method == "spearman") ){
    res <- cor( realScores, predScores,  use="complete", method = method  )
  }else if(method == "MSE") {
    res <- apply(predScores,2,function(ii){
      deltaVec <- abs( realScores - ii)
      deltaVec <- deltaVec^2
      mean(deltaVec)})
  }else if(method == "R2"){
    res <-apply(predScores,2,function(ii){
      R2 <- summary(lm(realScores ~ ii))$r.squared })
  }
  
  res <- matrix(res,nrow = 1,dimnames = list(method,colnames(predScores)))
  res
}
#a <- A(realScores = Maize_LOOCV_GYSSRes_NC[,1],predScores = Maize_LOOCV_GYSSRes_NC[,-1],method = "pearson",topAlpha = 1:90)
multiAlphaCor <- function(realScores, predScores, method = c("pearson", "kendall", "spearman", "MSE","R2"),topAlpha,BestIndividuals){
  res <- lapply(method,function(meth){
    cormat <- apply(predScores,2,function(x){
      sapply(topAlpha,function(ii){
        corEvaluationAlpha(realScores = realScores,predScores = x,method = meth,topAlpha = ii,BestIndividuals = BestIndividuals)})})
    rownames(cormat) <- paste0("top",topAlpha)
    cormat
  })
  names(res) <- method
  res
}
####
## function: RE and kappa from Heredity paper
## type:regression
## rank: top:(TRUE,TRUE); middle and down :(FALSE,FALSE)
## type: classfication
## rank: top:(TRUE,TRUE): middle and down :top:(FALSE,TRUE)
classEvaluation <- function(realScores, predScores, topAlpha = 15, Beta = 1,Probability = TRUE,evalMethod = "RE",BestIndividuals = c("top", "middle", "buttom") ) {
  
  #if( is.null( names(realScores)) ){
  #names(realScores) <- paste("sample", 1:length(realScores), sep = "")
  #}
  if( length(realScores) != length(predScores) ){
    stop("Error: different length between realScores and predScores")
  }
  
  if( length(BestIndividuals) > 1 ) {
    BestIndividuals <- BestIndividuals[1]
  }
  
  realScores <- as.numeric(realScores)
  predScores <- as.numeric(predScores)
  total <- length(realScores)
  topK <- round( total*topAlpha/100 )
  classVec <- c( rep(1, topK), rep(0, total-topK))
  ##
  if(BestIndividuals == "top"){
    decreaseReal <- TRUE
    decreasePred <- TRUE
  }else if(BestIndividuals == "buttom"){
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }else if(BestIndividuals == "middle"){
    realScores <- abs(realScores)
    predScores <- abs(predScores)
    if(Probability){
      decreaseReal <- FALSE
      decreasePred <- TRUE
    }else{
      decreaseReal <- FALSE
      decreasePred <- FALSE
    }
  }
  
  # if(Probability == TRUE){
  #   if(BestIndividuals != "top")  {
  #     decreaseReal <- FALSE
  #     decreasePred <- FALSE
  #   }
  # }else if (Probability == FALSE)  {
  #   if(BestIndividuals != "top") {
  #     decreaseReal <- FALSE
  #     decreasePred <- TRUE
  #     if(BestIndividuals == "middle") {
  #       realScores <- abs(realScores)
  #       predScores <- abs(predScores)
  #     }
  #   }
  # }else {
  #   stop("Error: undefined model type")
  # }
  
  scoreMat <- cbind( realScores, predScores )
  newScoreMat <- scoreMat[order(scoreMat[,1], decreasing = decreaseReal),] 
  newScoreMat <- cbind( newScoreMat, classVec )
  topRealMean <- mean( newScoreMat[1:topK,1] )
  #
  threshold <- newScoreMat[topK,1]
  #
  
  ##### RE ,kappa
  newScoreMat <- newScoreMat[order(newScoreMat[,2], decreasing = decreasePred),]
  # classVec <- c(rep(1,length(which(newScoreMat[,2] <= threshold))),rep(0, total-(length(which(newScoreMat[,2] <= threshold)))))
  newScoreMat <- cbind( newScoreMat, classVec )
  colnames(newScoreMat) <- c("real", "predicted", "realClass", "predClass")
  PredRealMean <- mean( newScoreMat[1:topK,1] ) 
  
  TP <- sum(newScoreMat[1:topK, 3])
  FP <- topK - TP
  FN <- topK - TP
  TN <- total - topK - FN
  Po <- (TP+TN)/total
  Pe <- ((FP+TN)/total)*((FN+TN)/total) + ((TP+FN)/total)*((TP+FP)/total)
  allRealMean <- mean( scoreMat[,1] )
  precision = TP/(TP + FP)
  recall = TP/(TP +FN)
  ###  the area under the receiver operating characteristics curve(AUC),and  the area under the precision-recall  curve (AUCpr)
  result <- sapply(evalMethod,function(one_method){
    switch(one_method,
           Kappa =  (Po-Pe)/(1-Pe),
           RE = ( PredRealMean - allRealMean)/(topRealMean - allRealMean),
           #AUC = roc.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE)$AUC ,
           AUC = roc(newScoreMat[,3],newScoreMat[,2])$auc[[1]],
           AUCpr = pr.curve(scores.class0 = newScoreMat[newScoreMat[,3] == 1,2],scores.class1 = newScoreMat[newScoreMat[,3] == 0,2],curve = TRUE,sorted = TRUE)$auc.integral,
           accuracy = (TP + TN)/total,
           F1 = (1 + Beta^2)*precision *recall/(precision + recall )
    )
  })
  names(result) <- paste(evalMethod,"_top", topAlpha, sep = "" )
  return(result) 
}
##################################### to evaluate for multiple parameters set
multiParameters <- function(realScores,predScores, Probability = TRUE,evalMethod = c("RE"),topAlpha ,Beta = 1,BestIndividuals = c("top"),probIndex = NULL){
  
  ############ calculate one or multiple topAlpha 
  multiTopAlpha <- function(realScores,predScores, Probability ,evalMethod ,topAlpha,Beta,BestIndividuals){
    ######################### one or multiple topAlpha for classifiction evaluation
    if (length(intersect(evalMethod,c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" ))) != 0){
      result <-  sapply(topAlpha,function(ii)
        classEvaluation(realScores = realScores,predScores = predScores,Probability = Probability,
                        evalMethod = evalMethod,topAlpha = ii,Beta = Beta,BestIndividuals = BestIndividuals)
      ) }
    
    ######################### one or multiple topAlpha for NDCG evaluation
    ################ set format of output #############
    if(length(evalMethod) > 1){
      result <- t(result)
    }else{
      result <- as.matrix(result)
    }
    dimnames(result) <- list(paste0("top",topAlpha),evalMethod )
    return(result)
  }
  
  
  predScores <- as.matrix(predScores)
  evalNum <- 1:ncol(predScores)
  
  ## class ##
  if(!is.null(probIndex)){
    classMethSite <- probIndex
    evalNum1 <- evalNum[-classMethSite]
    multiPredresultS <- lapply(evalNum1,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    classPredresult <- lapply(classMethSite,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = TRUE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
    multiPredresult <- list()
    length(multiPredresult) <- ncol(predScores)
    multiPredresult[evalNum1] <- multiPredresultS
    multiPredresult[classMethSite] <- classPredresult
  }else{
    ############ evaluate one or multiple prediction result
    multiPredresult <- lapply(evalNum,function(ii){
      multiTopAlpha(realScores = realScores,predScores = predScores[,ii],Probability = FALSE,evalMethod = evalMethod ,topAlpha = topAlpha,Beta = Beta,BestIndividuals = BestIndividuals )
    })
  }
  ######### set format of output
  result <- lapply(evalMethod,function(ii){
    one_evalMethod <- sapply(evalNum,function(jj) multiPredresult[[jj]][,ii])
    if(length(topAlpha) > 1){
      colnames(one_evalMethod) <- colnames(predScores)
    }else{
      one_evalMethod <- matrix(one_evalMethod,nrow = 1,dimnames = list(paste0("top",topAlpha),colnames(predScores)))
    }
    one_evalMethod
  })
  names(result) <- evalMethod
  return(result)
}
################################## evaluateGS ###################################
#' @title evaluateGS
#' @description  this function is used to evaluete the accuracy of predicted by genomic selection model.
#' @param realScores  (numeric,vector)vector is the real breeding values of the validation individual for a trait.
#' @param predScores  (numeric,vector or matrix)the prediction breeding value predicted by genomic selection model of the individuals.
#' @param Probability  (logical)whether the predScores is probability? Default FALSE.
#' @param evalMethod  (character)the evaluation methods selected to evaluate, which include "pearson", "kendall", "spearman", "MSE","R2"
#' "RE", "Kappa", "AUC","AUCpr","accuracy","F1","meanNDCG", "NDCG".
#' @param Beta (numeric)the parameter of "F1".
#' @param BestIndividuals (character)the position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".
#' @param topAlpha  (numeric,(0,100])a vector is the proportion of excellent individuals,default 1:90.
#' @param globalAlpha (logical)indicates if evaluate global methods(pearson, kendall, spearman, MSE and R2) by alpha,default FALSE.
#' @param probIndex (integer)indicates the column index which prediction result is probability.
#' @return 
#' a list inculding  evaluation results with methods which user selected.
#' @author Chuang Ma, Zhixu Qiu, Qian Cheng
#' @keywords evaluation, pearson, kendall, spearman, MSE, R2, RE, Kappa, AUC, AUCpr, accuracy, F1, meanNDCG, NDCG
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ########## predicting breeding value
#' predlist <-  G2PCrossValidation(cross = 10,seed = 1 ,cpus = 3,markers  = Markers,pheVal  = phenotype,
#'                  modelMethods = c("rrBLUP","RFC"),outputModel = FALSE)
#' predMartix <- NULL
#' for(ii in 1:10){predMartix <- rbind(predMartix,predlist[[ii]])}
#' ######## evaluate the accuracy of the prediction result
#'
#' evaluareTest <- evaluateGS(realScores = predMartix[,1], predScores = predMartix[,2:3], 
#'                            evalMethod = c("pearson", "kendall","spearman","RE","Kappa",
#'                                          "AUC","AUCpr","NDCG","meanNDCG",
#'                                          "MSE","R2","F1","accuracy"),topAlpha = 1:90, probIndex = 2)
#' }
evaluateGS <- function(realScores, predScores, Probability = FALSE, evalMethod = "RE", Beta = 1, BestIndividuals = "top", topAlpha = 1:90, allNDCG = F, globalAlpha = F, probIndex = NULL){
  require(pROC)
  require(PRROC)
  ## data process
  evalMat <- cbind(as.matrix(realScores),as.matrix(predScores))
  # probMat <- transValMat2ProbMat(evalMat = evalMat,BestIndividuals = BestIndividuals)
  # realScores <- probMat[,1]
  # predScores <- probMat[,-1]
  ## Evaluation
  selectFun <- NULL
  predScores <- as.matrix(predScores)
  globalMethods <-  c("pearson", "kendall", "spearman", "MSE","R2")
  thresholdMethods <- c("RE", "Kappa", "AUC","AUCpr","accuracy" ,"precision","recall","F1" )
  NDCGMethods <- c("meanNDCG", "NDCG")
  ###################### correlation methods
  if (length(intersect(evalMethod,globalMethods)) != 0){
    if(globalAlpha){
      selectFun <- c(selectFun,"globalEvaluation","globalEvaluationAlpha")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }else{
      selectFun <- c(selectFun,"globalEvaluation")
      globalMethods  <- intersect(evalMethod,globalMethods)
    }
  }
  
  ################# classifiction evaluation
  if (length(intersect(evalMethod,thresholdMethods)) != 0){
    selectFun <- c(selectFun,"thresholdEvaluation")
    thresholdMethods <- intersect(evalMethod,thresholdMethods)
  }
  
  result <- lapply(selectFun,function(one_fun){
    switch(one_fun,
           globalEvaluation = {corResult <- sapply(globalMethods,function(one_methods){corEvaluation( realScores, predScores, method = one_methods,BestIndividuals = BestIndividuals,Probability = FALSE)});
           matrix(t(corResult),ncol= ncol(predScores),dimnames = list(globalMethods,colnames(predScores)))},
           globalEvaluationAlpha = multiAlphaCor(realScores,predScores,topAlpha = topAlpha,method = globalMethods,BestIndividuals = BestIndividuals),
           thresholdEvaluation = multiParameters(realScores, predScores, topAlpha = topAlpha, Probability = Probability,evalMethod = thresholdMethods,Beta =Beta , BestIndividuals = BestIndividuals, probIndex = probIndex)
    )})
  ############ the format of output
  finalresult <- list()
  nList <- length(result)
  nList <- c("one","two","three")[nList]
  if(length(result) != 0){
    #     id <- 1
    #     if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];id <- 2};finalresult <- c(finalresult,result[[id]])
    switch(nList,
           one = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]]}else{finalresult <- c(finalresult,result[[1]])}},
           two = {if(!is.list(result[[1]])){finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]])}else{finalresult <- c(finalresult,result[[1]],result[[2]])}},
           three = {finalresult[["corMethods"]] <- result[[1]];finalresult <- c(finalresult,result[[2]],result[[3]])})
  }else{finalresult <- list()}
  #################### NDCG evaluation
  if (length(intersect(evalMethod,NDCGMethods)) != 0){
    selectFun <- intersect(evalMethod,NDCGMethods)
    result2 <- lapply(selectFun,function(one_fun){
      switch (one_fun,
              NDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "NDCG",topAlpha = topAlpha,realScores = realScores,predScores = x,allNDCG = allNDCG)}),
              meanNDCG = apply(predScores,2,function(x){NDCGEvaluation(method = "meanNDCG",topAlpha = topAlpha,realScores = realScores,predScores = x,allNDCG = allNDCG)})
      )})
    names(result2) <- intersect(evalMethod,NDCGMethods)
    finalresult <- c(finalresult,result2)
  }
  ### RFC SVC
  if ((!is.null(probIndex)) & ("corMethods" %in% names(finalresult))){
    probColIndex <- probIndex
    predScores <- as.matrix(predScores[,probColIndex])
    corResult <- sapply(globalMethods,function(one_methods){corEvaluation( realScores, predScores, method = one_methods,BestIndividuals = BestIndividuals,Probability = Probability)})
    corResult <- t(corResult)
    finalresult$corMethods[,probColIndex] <- corResult
  }
  finalresult
}
NDCGEvaluation <- function(method = NULL,topAlpha,realScores,predScores, allNDCG = T){
  
  if(allNDCG == T){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
  }else{
    ndcg <- round(length(realScores)*(topAlpha/100))
  }
  
  if(method == "NDCG"){
    result <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
      
    })
    if(allNDCG == F){
      names(result) <- paste0("NDCG_top",topAlpha)
    }
  }else if(method == "meanNDCG"){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
    NDCGEval <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
    })
    result <- sapply(topAlpha,function(ii){
      iii <- round(ii*length(realScores)/100)
      sum(NDCGEval[1:iii])/iii
    })
    names(result) <- paste0("meanNDCG_top",topAlpha)
  }
  result
}


one_cross_validation_GS <- function(cvSampleList,markers,pheVal,modelMethods ="SVC",outputModel = FALSE,
                                    nIter = 7000, burnIn = 500, thin = 5, 
                                    saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                                    verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,
                                    posPercentage = 0.4,BestIndividuals = c("top"),ntree = 500,nodesize = NULL,kernel = c("linear"),gamma = 1, cost = 2^(-9),
                                    K = 8,eta = 0.7,select = "pls2",fit = "simpls",scale.x = FALSE,scale.y = FALSE,eps = 1e-4,trace = FALSE, # SPLS parameters
                                    alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,maxstep = 100,
                                    epochs = 20, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 8,
                                    ...){
  ########### set training data and testing data
  trainIdx <- cvSampleList$trainIdx
  testIdx <-  cvSampleList$testIdx
  trainMarker <- markers[trainIdx,]
  trainPheno <- pheVal[trainIdx]
  testMarker <- markers[testIdx,]
  testPheno <- pheVal[testIdx]
  ####
  ##################  add of feature selection
  ############## 
  
  ################  run GS
  G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker,testPheno =testPheno ,modelMethods = modelMethods ,BestIndividuals = BestIndividuals,
      nIter = nIter, burnIn = burnIn, thin = thin,  saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,posPercentage = posPercentage,
      verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,ntree = ntree  ,nodesize = nodesize ,importance = importance,
      kernel = kernel,gamma = gamma, cost = cost,outputModel = outputModel,
      K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
      alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep,
      epochs = epochs, neurons = neurons, Z = Z,  effectConcider = effectConcider, mmerIters = mmerIters,
      ...)
}


##################################### cross validation for genomic selection  ####################################
#' @title G2PCrossValidation
#' @description  this function is apply cross validation to test Genomic Selection model trained by different methods and datas.
#' @param cross  (numeric)the fold number of cross validation.
#' @param seed  (numeric)random number options,defult 1.
#' @param cpus  (numeric)number of cpu cores to be used for calculations.
#' @param markers  (numeric) a matrix, each row is the each individual's SNP genotypes informations.Genotypes should be coded as {0,1,2};0 represents AA(homozygote),2 represents BB(homozygote) and 1 represents AB(heterozygote);missing (NA) alleles are not allowed
#' @param pheVal  (numeric)the phenotype Value of each individual.
#' @param modelMethods  the model to fit."BayesA", "BayesB", "BayesC", "BL", "BRR","RKHS","rrBLUP","LASSO","SPLS","SVC","SVR","RFR","RFC", "RR", "RKHS", "BRNN", "EM", "EMMA", "AI", "NR".
#' @param nIter,burnIn,thin  (integer) the number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.
#' @param saveAt  (string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs,default ""
#' @param S0,df0  (numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Default S0 NULL,df0 5.
#' @param R2  (numeric, 0<R2<1) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or http://genomics.cimmyt.org/BGLR-extdoc.pdf.Defult 0.5
#' @param weights  (numeric, n) a vector of weights, may be NULL. If weights is not NULL, the residual variance of each data-point is set to be proportional to the square of the weight. Only used with Gaussian outcomes.
#' @param verbose  (logical) if TRUE the iteration history is printed, default FALSE
#' @param rmExistingFiles  (logical) if TRUE removes existing output files from previous runs, default TRUE.
#' @param groups  (factor) a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.
#' @param ntree  RandomForest parameter:Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.Defualt 500
#' @param nodesize Randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param importance  RandomForest parameter:Should importance of predictors be assessed?Defualt FALSE
#' @param posPercentage  (numeric)the percentage positive samples in all samples.1 > posPercentage > 0.
#' @param BestIndividuals  It is a position that the best individuals (positive samples) in a training group, according to the breeding values of a training group's trait.
#'                          if the trait was yield,flowering or disease resistance,and  male flowering time to female flowering time,it is "top"(default), "buttom",and "middle" of the breeding values, respectively.
#' @param kernel  svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial)Default "linear".
#' @param gamma  svm parameter parameter needed for all kernels except linear (default: 1/(data dimension))
#' @param cost  svm cost,default 2^(-9)
#' @param outputModel  if true return the list of training model.
#' @param K  Number of hidden components
#' @param eta	 Thresholding parameter. eta should be between 0 and 1.
#' @param select  PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2"
#' @param fit	 PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".
#' @param scale.x	 Scale predictors by dividing each predictor variable by its sample standard deviation?
#' @param scale.y  Scale responses by dividing each response variable by its sample standard deviation?
#' @param eps  An effective zero. Default is 1e-4
#' @param maxstep  Maximum number of iterations when fitting direction vectors. Default is 100.
#' @param trace  Print out the progress of variable selection?
#' @param alpha  The elasticnet mixing parameter.Detail in glmnet.
#' @param X (numeric, matrix) design matrix related to the parameters not to be shrunk (i.e. fixed effects in the mixed model framework),defult no shrink.
#' @param family the distribution family of y, see help('family') for more details. 
#' @param lambda the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.
#' @param tol.err internal tolerance level for extremely small values; default value is 1e-6.
#' @param tol.conv tolerance level in convergence; default value is 1e-8.
#' @param neurons  (integer)indicates the number of neurons,defult 4.
#' @param epochs  (integer)maximum number of epochs(iterations) to train, default 30.
#' @param Z  (2-level list)incidence matrices and var-cov matrices for random effects. This works for ONE OR MORE random effects. This needs to be provided as a 2-level list structure, defult NULL. 
#' @param effectConcider  (string)if Z = NULL, random effects are auto generated.
#' @param iters  (numeric)a scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.
#' @param ...  arguments passed to or from other methods.
#' @seealso 
#' \link[G2P]{GSmachine}
#' \link[G2P]{GSReModel}
#' \link[G2P]{fit.RR}
#' \link[G2P]{fit.BRNN}
#' \link[G2P]{fit.RKHS}
#' \link[G2P]{fit.mmer}
#' \link[G2P]{predictGS}
#' \link[G2P]{G2P}
#' @return 
#' a list:
#' The prediction results of input GS method with cross validation.
#' 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu
#' @keywords cross validation,Genomic Selection
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ########## predicting breeding value
#' predlist <- G2PCrossValidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
#'                 pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
#'                 outputModel = FALSE)
#' }

# G2PCrossValidation <- function(cross = 5, seed = 1, cpus = 1, markers, pheVal, modelMethods ="SVC", outputModel = FALSE,
#                                nIter = 1500, burnIn = 500, thin = 5, 
#                                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
#                                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
#                                posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500, nodesize = 1, kernel = c("linear"), gamma = 1, cost = 2^(-9),
#                                K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,  # SPLS parameters
#                                alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
#                                epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 20,
#                                ...){
#   
#   cvSampleList <- cvSampleIndex(sampleNum = length(pheVal),cross = cross,seed = seed)
#   
#   sfInit(parallel = TRUE, cpus = cpus) 
#   sfLibrary( "brnn", character.only=TRUE)
#   sfLibrary( "glmnet", character.only=TRUE)
#   sfLibrary( "spls", character.only=TRUE)
#   sfLibrary( "pls", character.only=TRUE)
#   sfLibrary( "e1071", character.only=TRUE)
#   sfLibrary( "BGLR", character.only=TRUE)
#   sfLibrary( "rrBLUP", character.only=TRUE)
#   sfLibrary( "randomForest", character.only=TRUE)
#   sfLibrary( "G2P", character.only=TRUE)
#   sfLibrary( "sommer", character.only=TRUE)
#   sfLibrary( "hglm", character.only=TRUE)
#   #sfExport("RKHS.F2P")
#   
#   #sfExport("one_cross_validation_GS")
#   #sfExport("GSReModel")
#   #sfExport("trainedPredictModel_BGLR")
#   #sfExport("trainModel_RRBLUP")
#   #sfExport("GSmachine")
#   #sfExport("sampleClassify")
#   #sfExport("predictGS")
#   cvresult <- sfClusterApplyLB(cvSampleList,one_cross_validation_GS, markers = markers,pheVal = pheVal,modelMethods = modelMethods ,BestIndividuals = BestIndividuals,
#                                nIter = nIter, burnIn = burnIn, thin = thin,  saveAt = saveAt, S0 = S0, df0 = df0, R2 = R2, weights = weights,posPercentage = posPercentage,
#                                verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,ntree = ntree ,nodesize = nodesize ,importance = importance,
#                                kernel = kernel,gamma = gamma, cost = cost,outputModel = outputModel,
#                                K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,maxstep = maxstep,
#                                alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,
#                                epochs = epochs, neurons = neurons, Z = Z,  effectConcider = effectConcider, mmerIters = mmerIters, ...)
#   
#   sfStop()
#   
#   task_names <- paste0("cv",1:cross)
#   names(cvresult) <- task_names 
#   if(outputModel){
#     modelList <- lapply(task_names,function(cv) cvresult[[cv]]$trainModel)
#     predscores <- lapply(task_names,function(cv) cvresult[[cv]]$predscores)
#     names(modelList) <- task_names 
#     names(predscores ) <- task_names
#     cvresult <- list(modelList,predscores)
#     names(cvresult) <- c("modelList","predscores")
#     return(cvresult)
#     
#   }else{
#     predscores <- lapply(task_names,function(cv) cvresult[[cv]])
#     names(predscores ) <- task_names
#     return(predscores)
#   }
# }
# 

G2PCrossValidation <-function(cross = 5, seed = 1, cpus = 1, markers, pheVal, modelMethods ="SVC", outputModel = FALSE,
                                nIter = 1500, burnIn = 500, thin = 5, 
                                saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
                                posPercentage = 0.4, BestIndividuals = c("top"), ntree = 500, nodesize = NULL, kernel = c("linear"), gamma = 1, cost = 2^(-9),
                                K = 8, eta = 0.7, select = "pls2", fit = "simpls", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100,  # SPLS parameters
                                alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                                epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 20,
                                ...){
  require("parallel")
  sampleNum <- nrow(markers)
  phenotypeNum <- length(pheVal)
  cross = cross; seed = seed; cpus = cpus; markers = markers;pheVal = pheVal;
  modelMethods = modelMethods; BestIndividuals = BestIndividuals;
  nIter = nIter; burnIn = burnIn; thin = thin; saveAt = saveAt; S0 = S0; df0 =df0; R2 = R2; weights = weights;posPercentage = posPercentage;
  verbose = verbose; rmExistingFiles = rmExistingFiles; groups=groups;ntree = ntree  ;nodesize = nodesize ;importance = importance;
  kernel = kernel;gamma = gamma; cost = cost;outputModel = outputModel;
  K = K;eta = eta;select = select;fit = fit;scale.x = scale.x;scale.y = scale.y;eps = eps;trace = trace;
  alpha = alpha;X = X;family = family; lambda = lambda; tol.err = tol.err; tol.conv = tol.conv;maxstep = maxstep;
  epochs = epochs; neurons = neurons; Z = Z;  effectConcider = effectConcider; mmerIters = mmerIters;
  if(sampleNum != phenotypeNum) {
    stop("Marker count is not equal to phenotype count!")
  }
  
  if(!is.numeric(markers) | !is.numeric(pheVal)){
    stop("Marker or phenotype is not numeric, please check it!")  
  }
  
  cl <- makeCluster(cpus)
  cat(cpus," cores were used for cross validation ... \n")
  cat("Start cross validation ... \n")
  results <- parLapply(cl,1:cross,function(x){
    library("BGLR")
    library("G2P")
    library("brnn")
    library("glmnet")
    library("spls")
    library("pls")
    library("e1071")
    library("BGLR")
    library("rrBLUP")
    library("randomForest")
    library("sommer")
    library("hglm")
    cat("All needed package have loaded in all cores! \n")
    
    cvSampleList <- cvSampleIndex(sampleNum = sampleNum,cross = cross,seed = seed)
    trainIdx <- cvSampleList[[x]]$trainIdx
    testIdx <-  cvSampleList[[x]]$testIdx
    trainMarker <- markers[trainIdx,]
    trainPheno <- pheVal[trainIdx]
    testMarker <- markers[testIdx,]
    testPheno <- pheVal[testIdx]
    G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker,testPheno =testPheno ,modelMethods = modelMethods ,BestIndividuals = BestIndividuals,
        nIter = nIter, burnIn = burnIn, thin = thin,  saveAt = saveAt, S0 = S0, df0 =df0, R2 = R2, weights = weights,posPercentage = posPercentage,
        verbose = verbose, rmExistingFiles = rmExistingFiles, groups=groups,ntree = ntree  ,nodesize = nodesize ,importance = importance,
        kernel = kernel,gamma = gamma, cost = cost,outputModel = outputModel,
        K = K,eta = eta,select = select,fit = fit,scale.x = scale.x,scale.y = scale.y,eps = eps,trace = trace,
        alpha = alpha,X = X,family = family, lambda = lambda, tol.err = tol.err, tol.conv = tol.conv,maxstep = maxstep,
        epochs = epochs, neurons = neurons, Z = Z,  effectConcider = effectConcider, mmerIters = mmerIters,
        ...)
  }) # lapply的并行版本
  stopCluster(cl)
  cat(cross,"fold cross validation is done! \n")
  names(results) <- paste0("CV",1:cross)
  results
}

######################################### heatMapDataProcess ######################################
heatMapDataProcess <- function(x,highBound,lowBound,alpha,basedMethod = "best"){
  x <- rbind(x[["corMethosds"]],t(sapply(x[-1],function(x){x[alpha,]})))
  for(i in 1:nrow(x)){
    if(length(which(x[i,] == "NaN")) >= 1){
      x <- x[-i,]
    }
  }
  if(basedMethod == "best"){
    highBound = 0
    lowBound = lowBound
    best <- apply(x,1,function(x){max(x)})
    rowNames <- rownames(x)
    colNames <- colnames(x)
    finalMat <- matrix(NA,1,ncol(x))
    for (i in 1:length(best)) {
      if (rowNames[i] == "MSE") { 
        changeRes <- -((x[i,] - min(x[4,]))/min(x[4,]))*100
      }else{
        changeRes <- ((x[i,] - best[i])/best[i])*100  
      }
      finalMat <- rbind(finalMat, changeRes)
    }
    finalMat <- finalMat[-1,]
    colnames(finalMat) <- colNames
    rownames(finalMat) <- rowNames
    
    finalMat[which(finalMat >= highBound)] <- highBound
    finalMat[which(finalMat <= lowBound)] <- lowBound
  }else{
    heat_increase <- ((x - x[,basedMethod])/x[,basedMethod])*100
    heat_increase["MSE",] <- -heat_increase["MSE",]
    heat_increase[which(heat_increase >= highBound)] <- highBound
    heat_increase[which(heat_increase <= lowBound)] <- lowBound
    finalMat <- heat_increase
  }
  finalMat
}
