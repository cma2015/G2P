# ___G2P : Genotypes to Phenotypes___<br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
![](https://encrypted-tbn3.gstatic.com/images?q=tbn:ANd9GcS3RzhXKSfXpWhWhvClckwi1Llj1j3HvjKpjvU8CQv4cje23TwS "windows logo")
<br>
The R package "G2P" is an integrated _Genomic Selection_ (GS) package for predicting phenotypes from genotypes,
which includes 20 GS algorithms and 13 evaluation measurements. G2P provide a comprehensive but easy-to-use platform  
for _Genomic Selection_ researchers.Besides, G2P also provide a interactive UI based on ___shiny___ which can easily operated and learned.
<br>
## Version and download <br>
* [Version 1.0](https://github.com/cma2015/G2P/blob/master/G2P_1.0.tar.gz) -First version released on Feb, 28th, 2017<br>
* [Version 1.1](https://github.com/cma2015/G2P/blob/master/G2P_1.1.0.tar.gz) -Second version released on August, 14th, 2017<br>
## Depends <br>
* [R](https://www.r-project.org/) (>= 3.3.1)
* BGLR(>= 2.10)
* pROC(>= 1.8)
* PRROC(>= 1.1)
* e1071(>= 1.6-7)
* glmnet(>= 2.0)
* pls(>= 2.5-0)
* randomForest(>= 4.6-12)
* rrBLUP(>= 4.4)
* snowfall(>= 1.84-6.1)
* spls(>= 2.2-1)
* brnn(>= 0.6)
* sommer(>= 2.9)
* hglm(>=2.1-1)
* hglm.data(>=1.0-0)
* snow(>=0.4-1)
* snowfall(>=1.84-6.1)
<br>
Please assure that you have installed all depends packages before you use G2P <br>

## Suggests
* rgl(>= 0.97.0)
* pheatmap(>= 1.0.8)
* shiny(>= 1.0.3)
* plotly(>= 4.7.0)
* shinythemes(>= 1.1.1)
* ggplot2(>= 2.2.1.9000)
* impute(>=1.46.0)
<br>
Thesew packages for results display and G2P.app, if you want use interactive UI and plot in G2P,please assure you have
install these packages.Besides,___impute___ for "GSDataQC" function to impute with ___knn___ methods. <br>

## Installation <br>
### Install dependency packages
```R
install.packages(c("BGLR", "PRROC","e1071","glmnet","spls","randomForest","rrBLUP","snowfall","pls","brnn","sommer","hglm","hglm.data","snow","snowfall"), dependencies=TRUE)
```
### Install G2P
```R
install.packages("path/G2P_1.1.0.tar.gz", repos = NULL, type = "source")
# The path cant including space.
```
### Install suggested packages
```R
install.packages(c("rgl", "pheatmap","shiny","plotly","shinythemes","ggplot2","impute"), dependencies=TRUE)
```
## Contents

#### Main functions
* GS Data summary <br>
* Example data <br>
* Trainning model  <br>
* Performance assement <br>
* Cross validation <br>
* Feature selection <br>
* Results display <br>
* Package help<br>
* G2P.app<br>

#### GS algorithms (20)
* __Statistics based methods__<br>
BayesA, BayesB, BayesC, BRR, BL, RKHS, RR, rrBLUP, SPLS, LASSO, BRNN, AI, NR, EM, EMM, bigRR
* __Machine-learning based methods__<br>
RFC, RFR, SVC, SVR  <br>

#### Evaluation measures (13)
* __Global measures__ <br>
 Pearson correlation, Kendall rank correlation, Spearman Correlation, Mean squared error (MSE), R2
* __Threshold-based measures__ <br>
Normalized discounted cumulative gain (NDCG), meanNDCG, AUC, AUCpr, F1, Kappa, Relative efficiency (RE), Accuracy <br>
#### Example data
The G2P package have built-in example dataset ___GYSS___, but it only a subset of grain yield under drought stressed(GYSS) dataset with 242 samples and 1000 SNPs. GYSS dataset ware from International Centre for the Improvement of Maize and Wheat (CIMMYT). Complete dataset could be [downloaded](https://github.com/cma2015/G2P) with 242 samples and 46373 SNPs. We use subset not complete dataset in order to shorten the compute time.
## Quick start
More details please install [G2P](https://github.com/cma2015/G2P/blob/master/G2P_1.1.tar.gz) and see help in R.<br>

## G2P Tutorial
### Command line
#### 0 Setting up the R session
Before starting, the user should choose a working directory, preferably a directory devoted exclusively for this tutorial. After starting an R session, change working directory, load the requisite packages and set standard options:
```R
# Display the current working directory
getwd();
# Set working directory
workingDir = ".";
setwd(workingDir)
# Library G2P package
Library(G2P)
```
#### 1 Load example datasets and quick look at the format of example dataset.
Load datasets and check the input of G2P. This function provide two methods to impute the missing value and finally giving a summary list including multifarious information about input data. For example the count and percent of miss value, the minor allele frequency (MAF) etc.
data(GYSS)
```R
# Genotypic data
Markers[1:10,1:10]
# Phenotypic data
phenotype[1:10]
## GSDataQC, not impute ##
QCRes <- GSDataQC(markers = Markers, phenotype = phenotype, impute = F)

## GSDataQC, impute ##
misIndex <- sample(1:242000,100000)
Markers[misIndex] <- NA
QCResImpute <- GSDataQC(markers = Markers, phenotype = phenotype, impute = T, 
                        imputeMethod = "mean")
```
#### 2 Feature selection
To score each SNP set, you can screen high grade of SNP for subsequent modeling, in order to simplify the operation and improve the precision of feature selection.Parameter "method" including "Gini","rrBLUP" and "Accuracy".
```
# Feature selection with rrBLUP
rrBLUP_selection <- feature_assess(markers = Markers, phenotype = phenotype, method = "rrBLUP")
# This function return a numeric array indicates the score of each position of SNPs
```
#### 3 Modeling
Details for parameter setting and illustration please refer to G2P package’s help in R or reference manual in pdf.
```R
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
```
#### 4 Prediction
Genomic selection methods including "BRNN", "RKHS", "RR", "AI", "NR", "EM", "EMMA" have got prediction results above. The following prediction function were apply to methods "BayesA", "BayesB", "BayesC", "BL", "BRR", "rrBLUP", "LASSO", "SPLS", "bigRR", "SVR" , "SVC", "RFR" and "RFC". Details for parameter setting and illustration please refer to G2P package’s help in R or reference manual in pdf.
```R
# testMat is a new marker matrix which need to predict
# trainModel is the already model in 3 and the modelMethods indicates the name of method.
rrBLUP_Res <- predictGS(testMat = Markers[1:10,], trainModel = rrBLUP_model, modelMethods = "rrBLUP")
This function return a numeric array indicates the prediction results of each predicted sample.
```
#### 5 G2P
Multi-methods genotype to phenotype.
```R
> G2P(trainMarker = Markers, trainPheno = phenotype, testMarker = Markers[1:10,], testPheno = phenotype[1:10], modelMethods =c("BayesA", "BayesB", "BayesC", "BL", "BRR", "rrBLUP","RFC"), outputModel =FALSE)
      realPhenScore    BayesA    BayesB    BayesC        BL       BRR    rrBLUP   RFC
DT10      0.7045494 0.6638141 0.6645500 0.6512245 0.6496041 0.6553746 0.6574948 0.820
DT100     0.3886150 0.4270637 0.4152724 0.4201070 0.4275309 0.4197220 0.4077610 0.096
DT101     0.1165665 0.3105247 0.3044984 0.2987718 0.3117646 0.2924395 0.2930300 0.104
DT102     0.3747609 0.4727063 0.4783084 0.4745835 0.4849857 0.4759522 0.4694439 0.162
DT103     0.4136186 0.4770212 0.4841845 0.4843115 0.4933730 0.4942586 0.4833140 0.166
DT104     0.6908267 0.6321259 0.6272183 0.6269430 0.6128531 0.6268135 0.6300071 0.778
DT105     0.5520021 0.5536663 0.5411970 0.5555117 0.5449907 0.5263780 0.5400619 0.144
DT109     0.5161542 0.4841527 0.4866321 0.4982482 0.4744184 0.4835680 0.4826495 0.110
DT110     0.5486291 0.6010441 0.6076803 0.5893248 0.6028943 0.6064955 0.5979550 0.124
DT111     0.7276410 0.5877352 0.5966961 0.6118398 0.6075789 0.6014777 0.6036835 0.786
```
#### 6 Cross validation
Run cross validation (CV) of GS methods. Parameter “cross” indicates the folds of CV, and “seed” sets random seed. “cpus” sets the core numbers of parallel.
```R
predlist <- G2PCrossValidation(cross = 10, seed = 1 , cpus = 3, markers = Markers,
                               pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
                               outputModel = FALSE)
# This function return a list and each element indicates one fold results of CV.
> predlist$cv1
       testPheno    rrBLUP   RFC
DT187 0.27187247 0.4862278 0.316
DT21  0.08312729 0.4208843 0.254
DT257 0.25201410 0.4822461 0.302
DT72  0.47947204 0.6243357 0.798
DT17  0.50980385 0.4601366 0.154
DT68  0.45270891 0.4154981 0.176
DT78  0.80182295 0.4193034 0.332
DT273 0.61428571 0.5145006 0.344
DT266 0.40762308 0.7258433 0.528
DT116 0.82464237 0.5425566 0.362
DT169 0.29187668 0.4181612 0.284
DT158 1.01129286 0.6474776 0.428
DT276 0.50811050 0.5782443 0.396
DT208 0.58764420 0.6148920 0.402
DT294 0.93534274 0.6184581 0.522
DT233 1.23127979 0.5756115 0.394
DT282 0.51583669 0.5852271 0.374
DT79  0.30790897 0.5193164 0.342
DT206 0.38559096 0.5057513 0.332
DT292 0.49675427 0.5991949 0.704
DT63  0.84989371 0.6700435 0.820
DT167 0.25885476 0.6123705 0.452
DT262 0.56417625 0.5828048 0.208
DT134 0.39868010 0.4842858 0.348
```
#### 7 Evaluation
Evaluation function for GS, 13 assessment methods can be set including "pearson", "kendall", "spearman", "MSE","R2" "RE", "Kappa", "AUC","AUCpr","accuracy","F1","meanNDCG", "NDCG". Parameter “topAlpha” indicates the proportion of excellent individuals.
```R
# Merge the results
predMat <- resultsMerge(predList = predlist)
# Evaluation
evalTest <- evaluateGS(realScores = predMat [,1], predScores = predMat [,2:3], 
                       evalMethod = c( "pearson", "kendall","spearman", "RE", "Kappa",
                                      "AUC", "AUCpr", "NDCG", "meanNDCG",
                                      "MSE", "R2", "F1", "accuracy"), topAlpha = 1:90)
This function return a list and each element indicates evaluation results of one assessment strategy.
> evalTest$corMethosds
             rrBLUP        RFC
pearson  0.35861146 0.34324085
kendall  0.22355200 0.21771634
spearman 0.32433095 0.32423002
MSE      0.07325809 0.09292603
R2       0.12860218 0.11781428
> evalTest$RE
         rrBLUP        RFC
top1  0.4623412 0.06032638
top2  0.3552597 0.18740391
top3  0.5868788 0.29001966
top4  0.5445002 0.32067173
top5  0.5474425 0.38223681
top6  0.5146146 0.35913813
top7  0.5288127 0.36279797
top8  0.4777331 0.32782280
top9  0.4415385 0.35465968
top10 0.3841515 0.37457313
top11 0.3953131 0.44471255
top12 0.4269692 0.41753296
top13 0.3919942 0.46025826
top14 0.3324953 0.50017893
top15 0.3418756 0.47511276
top16 0.3586926 0.46829151
top17 0.3621604 0.41855918
top18 0.3478598 0.42316310
top19 0.3563579 0.38201784
top20 0.3588952 0.37693601
top21 0.3095465 0.37088407
top22 0.3259495 0.35949169
top23 0.3308104 0.38578747
top24 0.3305805 0.37949081
top25 0.3430151 0.42111577
top26 0.3597613 0.43218671
top27 0.3730551 0.40364655
top28 0.3608020 0.39152421
top29 0.3769560 0.39211294
top30 0.3618355 0.39188584
```
#### 8 Visualization
```R
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
```
![](https://github.com/cma2015/G2P/blob/master/vignettes/density%20plot.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/scatter%20plot.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/bar%20plot.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/lines%20plot.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/heatmap1.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/heatmap2.png)
#### 9 Integration
Parameter “ratio” set the weight of each methods which need to be integrated. Auto optimize will select two methods which are the most suitable, then assign weights 1:1.
```R
# Preparation of Multi-methods results 
predlist <- G2PCrossValidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
                pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP"),
                outputModel = FALSE)
# Merge results
resultMat <- resultsMerge(predlist)
inte <- GSIntegrate(predResMat = resultMat[,1:5], ratio = c(2,3,4,4), autoOptimize = F)
> inte[1:10,]
      realPhenScore    BayesA    BayesB    BayesC    rrBLUP     merge
DT187    0.27187247 0.4854131 0.5022911 0.4905175 0.4862278 0.4911293
DT21     0.08312729 0.4062116 0.4054856 0.4023935 0.4208843 0.4093839
DT257    0.25201410 0.4416190 0.4693871 0.4514048 0.4822461 0.4635387
DT72     0.47947204 0.6265084 0.6269129 0.6308500 0.6243357 0.6272691
DT17     0.50980385 0.4457797 0.4529249 0.4495355 0.4601366 0.4530017
DT68     0.45270891 0.3990259 0.4056420 0.4045821 0.4154981 0.4073307
DT78     0.80182295 0.4042674 0.4313791 0.4183998 0.4193034 0.4194989
DT273    0.61428571 0.5154637 0.5050100 0.5040182 0.5145006 0.5092333
DT266    0.40762308 0.7369941 0.7325658 0.7329628 0.7258433 0.7313008
DT116    0.82464237 0.5539202 0.5514976 0.5504473 0.5425566 0.5487961
```
#### 10 Other functions
Other practical functions in G2P: “resultsMerge”, “G2PTest”.
resultsMerge : merge CV results list and create a matrix.
G2PTest: test run time of G2P cross validation
```R.
# Preparation of Multi-methods results 
predlist <- G2P(cross = 10,seed = 1 , cpus = 3, markers = Markers,
            pheVal = phenotype, modelMethods = c("BayesA","BayesB","BayesC","rrBLUP", "RFC"),
            outputModel = FALSE)
# Merge results
resultMat <- resultsMerge(predlist)
# G2PTest
test <-  G2PTest(cross = 10, seed = 1, cpus = 3, markers = Markers,
                 pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
                 outputModel = FALSE)
> test$runOneFoldTime
[1] "00:00:05"
> test$CVtime
[1] 50
```
### G2P.app (Interaction under shiny)
#### 0 Launching G2P.app
To make the implementation more user friendly, TSIS analysis is integrated into a [Shiny App](https://shiny.rstudio.com/) ([Chang, et al., 2016](https://shiny.rstudio.com/)). By typing
```R
G2P.app()
```
The App is opened in the default web browser. Users can upload input datasets, set parameters for switch analysis, visualize and save the results easily. The [G2P](https://github.com/cma2015/G2P/) App includes four tab panels (see __Figure 1__ to __Figure 5__ ).
#### Tab panel 1
The first panel including G2P data input and quality control
![](https://github.com/cma2015/G2P/blob/master/vignettes/%E5%9B%BE%E7%89%871.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/%E5%9B%BE%E7%89%873.png)
Input data example
##### Markers structure
![](https://github.com/cma2015/G2P/blob/master/vignettes/markers.png)
##### Phenotype structure
![](https://github.com/cma2015/G2P/blob/master/vignettes/phenotypes.png)
#### Tab panel 2
The second panel including G2P, G2P cross validation and evalation parameter setting and runing.
![](https://github.com/cma2015/G2P/blob/master/vignettes/panel2.png)
#### Tab panel 3
The third panel including visulization of G2P cross validation and evaluation results.
![](https://github.com/cma2015/G2P/blob/master/vignettes/%E5%9B%BE%E7%89%874.png)
Result table
![](https://github.com/cma2015/G2P/blob/master/vignettes/table.png)
Plots:
![](https://github.com/cma2015/G2P/blob/master/vignettes/visualScatter.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/visualLine.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/visualBar.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/visualHeat1.png)
![](https://github.com/cma2015/G2P/blob/master/vignettes/visualHeat2.png)
#### Tab panel 4 
The last tab panel includes this user manual.
## Ask questions
Please use [G2P/issues](https://github.com/cma2015/G2P/issues) for how to use DeepGS and reporting bugs.
