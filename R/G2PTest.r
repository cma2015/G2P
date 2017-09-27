##################################### G2P CV test for genomic selection  ####################################
#' @title G2PTest
#' @description  Testing the program and get the time of one cross validation.
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
#' $runOneFoldTime  The prediction results of input GS method with cross validation.
#' $res  Test program results
#' 
#' @author Chuang Ma ,Qian Cheng , Zhixu Qiu ,Jie Song
#' @keywords cross validation, genomic selection, test 
#' @export
#' @examples
#' \dontrun{
#' data(GYSS)
#' ########## predicting breeding value
#' test <-  G2PTest(cross = 10, seed = 1, cpus = 3, markers = Markers,
#'                  pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
#'                  outputModel = FALSE)
#' }
G2PTest <- function(cross = 10,ncross = NULL,seed = 1,cpus = 1, markers, pheVal, modelMethods ="SVC",
                    nIter = 1500, burnIn = 500, thin = 5, 
                    saveAt = "", S0 = NULL, df0 =5, R2 = 0.5, weights = NULL,
                    verbose = FALSE, rmExistingFiles = TRUE, groups=NULL, importance = FALSE,
                    posPercentage = 0.4, BestIndividuals = c("top"),
                    ntree = 500,nodesize = 1,
                    kernel = c("radial"), gamma = 1, cost = 2^(-9), outputModel = FALSE,
                    K = 8, eta = 0.7, select = "pls2", fit = "simples", scale.x = FALSE, scale.y = FALSE, eps = 1e-4, trace = FALSE, maxstep = 100, # SPLS parameters
                    alpha = 1,X = NULL,family = gaussian(link = identity), lambda = NULL, tol.err = 1e-6, tol.conv = 1e-8,
                    epochs = 30, neurons = 4, Z = NULL,  effectConcider = "A", mmerIters = 8,
                    ...){
  require(parallel)
  startTime <- Sys.time()
  cross = cross; seed = seed; cpus = cpus; markers = markers;pheVal = pheVal;
  modelMethods = modelMethods; BestIndividuals = BestIndividuals;
  nIter = nIter; burnIn = burnIn; thin = thin; saveAt = saveAt; S0 = S0; df0 =df0; R2 = R2; weights = weights;posPercentage = posPercentage;
  verbose = verbose; rmExistingFiles = rmExistingFiles; groups=groups;ntree = ntree  ;nodesize = nodesize ;importance = importance;
  kernel = kernel;gamma = gamma; cost = cost;outputModel = outputModel;
  K = K;eta = eta;select = select;fit = fit;scale.x = scale.x;scale.y = scale.y;eps = eps;trace = trace;
  alpha = alpha;X = X;family = family; lambda = lambda; tol.err = tol.err; tol.conv = tol.conv;maxstep = maxstep;
  epochs = epochs; neurons = neurons; Z = Z;  effectConcider = effectConcider; mmerIters = mmerIters;
  cvSampleList <- cvSampleIndex(sampleNum = length(pheVal),cross = cross,seed = seed)
  if(is.null(ncross)){
    cvSampleList <- cvSampleList[1]
  }else{
    cvSampleList <- cvSampleList[1:ncross]
  }

  cl <- makeCluster(cpus)
  cat(cpus," cores were used for cross validation ... \n")
  cat("Start cross validation ... \n")
  results <- parLapply(cl,1:ncross,function(x){
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
    markers <- markers
    cvSampleList <- cvSampleList
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
  names(results) <- paste0("CV",1:ncross)
  results
  
  totalTime <- .f(startTime)
  allTime <- (as.numeric(substr(totalTime,1,2))*3600 + as.numeric(substr(totalTime,4,5))*60 + as.numeric(substr(totalTime,7,8))) * cross/ncross
  testRes <- list(runOneFoldTime = totalTime,CVtime = allTime,res = results)
  
  cat(paste0("The time of run one fold is: ",totalTime,"\n","The total time of this system will spare is :",allTime,"\n"))
  testRes
}


.f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
