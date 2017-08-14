###### G2P test code
## 0. load example data
data(GYSS)
# Genotypic data
Markers[1:10,1:10]
# Phenotypic data
phenotype[1:10]

## 1. Data and check(GSDataQC)
## GSDataQC, not impute ##
QCRes <- GSDataQC(markers = Markers, phenotype = phenotype, impute = F)

## GSDataQC, not impute ##
misIndex <- sample(1:242000,100000)
misMarkers <- Markers
misMarkers[misIndex] <- NA
QCResImpute <- GSDataQC(markers = Markers, phenotype = phenotype, impute = T, 
                        imputeMethod = "mean")

## 2.Feature selection with rrBLUP
rrBLUP_selection <- feature_assess(markers = Markers, phenotype = phenotype, method = "rrBLUP", 
                                   posPercentage = 0.40, BestIndividuals = "top")
#This function return a numeric array indicates the score of each position of SNPs

## 3.Modeling

# Fit a regression model (modelMethods including "BayesA", "BayesB", "BayesC", "BL", "BRR",  # "rrBLUP",    
# "LASSO", "SPLS", "bigRR". ) 
rrBLUP_model <- GSReModel(markers = Markers, pheVal = phenotype, modelMethods = "rrBLUP")
# Fit a machine learning model (modelMethods including "SVR" , "SVC", "RFR", "RFC")
# Fit RFR model
machine_model <- GSmachine(markers = Markers, pheVal = phenotype, modelMethods = "RFR")
# Fit classification model(RFC)
machine_model <- GSmachine(markers = Markers, pheVal = phenotype, modelMethods = "RFC",
                           posPercentage = 0.4, ntree = 500)
# Fit other models ("BRNN", "RKHS", "RR", "AI", "NR", "EM", "EMMA"), set parameter "outputModel = TRUE" to get a list including prediction results and model, otherwise, only output prediction results.  
BRNN_Res <- fit.BRNN(trainedMarkerMat = Markers, trainedPheVal = phenotype,
                     predictMarkerMat = Markers[1:10,], outputModel = TRUE,verbose = F)
RKHS_Res <- fit.RKHS(trainedMarkerMat = Markers, trainedPheVal = phenotype, 
                     predictMarkerMat = Markers[1:10,],nIter = 1500, burnIn = 500, outputModel = TRUE)
RR_Res <- fit.RR(trainedMarkerMat = Markers, trainedPheVal = phenotype,
                 predictMarkerMat = Markers[1:10,], outputModel = TRUE )
# Fit mmer models (method including "AI", "NR", "EM", "EMMA" )
mmer_Res <- fit.mmer(trainedMarkerMat = Markers, trainedPheVal = phenotype, 
                     predictMarkerMat = Markers[1:10,], method = "NR", effectConcider = "A", outputModel = TRUE)

## 4. prediction 
# testMat is a new marker matrix which need to predict
# trainModel is the already model in 3 and the modelMethods indicates the name of method.
rrBLUP_Res <- predictGS(testMat = Markers[1:10,], trainModel = rrBLUP_model, modelMethods = "rrBLUP")

## 5.G2P
G2P(trainMarker = Markers, trainPheno = phenotype, testMarker = Markers[1:10,], testPheno = phenotype[1:10], modelMethods =c("BayesA", "BayesB", "BayesC", "BL", "BRR", "rrBLUP","RFC"), outputModel =FALSE)

## 6. G2P cross validation 
predlist <- G2PCrossValidation(cross = 10, seed = 1 , cpus = 3, markers = Markers,
                               pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
                               outputModel = FALSE)

## 7. Evaluation 
# Merge the results
predMat <- resultsMerge(predList = predlist)
# Evaluation
evalTest <- evaluateGS(realScores = predMat [,1], predScores = predMat [,2:3], 
                       evalMethod = c( "pearson", "kendall","spearman", "RE", "Kappa",
                                      "AUC", "AUCpr", "NDCG", "meanNDCG",
                                      "MSE", "R2", "F1", "accuracy"), topAlpha = 1:90)
# This function return a list and each element indicates evaluation results of one assessment strategy.
## 8 Visualization
## 
G2PCVRes <-  G2PCrossValidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
                                pheVal = phenotype, modelMethods = c("BayesA", "BayesB", "BayesC", "BL", "BRR","RR","RKHS","rrBLUP","LASSO","SPLS","bigRR","SVC","RFC","SVR","RFR"),
                                outputModel = FALSE)
CVres <- resultsMerge(predList = G2PCVRes)
evalTest <- evaluateGS(realScores = CVres[,1], predScores = CVres[,2:20], 
                       evalMethod = c( "pearson", "kendall","spearman", "RE", "Kappa",
                                       "AUC", "AUCpr", "NDCG", "meanNDCG",
                                       "MSE", "R2", "F1", "accuracy"), topAlpha = 1:90)
### 8.1 row data visulization
## phenotype distribution  plot
rowDataPlot(y = phenotype,show.line = T,barCol = "blue",lineCol = "red")
## PCA 3-D plot 
htmlwidgets::saveWidget(as_widget(rowDataPlot(markers = Markers,y = phenotype,plot.type = "PCA")), file="3-D_PCA.html",selfcontained=T)

### 8.2 scatter plot 
scatterPlot(CVres,x1 = "BayesA",x2 = "RFC",show.line = F,color_szie = T,make.plotly = F,sizeRange = c(4,6))
### 8.3 lines plot 
linePlot(evalMat = evalTest$RE,size = 1)
### 8.4 bar plot 
barPlot(evalTest$corMethosds,other = "sector")
### 8.5 heat map 
#### pred res heatmap 
heatmapPlot(predmat = CVres,make.plotly = F,col.low = "green",col.high = "red")
#### eval res heatmap 
heatmapEval <- heatMapDataProcess(x = evalTest,highBound = 0,lowBound = -30,alpha = 15,basedMethod = "best")
heatmapPlot(predmat = heatmapEval,make.plotly = F,col.low = "green",col.high = "red")
## 9 integrate 
# Preparation of Multi-methods results 
predlist <- G2PCrossValidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
                pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP"),
                outputModel = FALSE)
# Merge results
resultMat <- resultsMerge(predlist)
inte <- GSIntegrate(predResMat = resultMat[,1:5], ratio = c(2,3,4,4), autoOptimize = F)