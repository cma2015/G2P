# ___G2P : Genotypes to phenotypes___<br>
![](https://halobi.com/wp-content/uploads/2016/08/r_logo.png "R logo")
![](https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcSvCvZWbl922EJkjahQ5gmTpcvsYr3ujQBpMdyX-YG99vGWfTAmfw "linux logo")
![](https://encrypted-tbn3.gstatic.com/images?q=tbn:ANd9GcS3RzhXKSfXpWhWhvClckwi1Llj1j3HvjKpjvU8CQv4cje23TwS "windows logo")
<br>
The R package "G2P" is an integrated _Genomic Selection_ (GS) package for predicting phenotypes from genotypes,
which includes 15 GS algorithms and 13 evaluation measurements. G2P provide a comprehensive but easy-to-use platform  
for _Genomic Selection_ researchers.
<br>
## Version and download <br>
* [Version 1.0](https://github.com/cma2015/G2P/blob/master/G2P_1.0.tar.gz) -First version released on Feb, 28th, 2017<br>

## Installation <br>
```R
install.package("Download path/G2P_1.0.tar.gz")
```
<br>
## Depends
* [R](https://www.r-project.org/) (>= 3.3.1)
* BGLR(>= 2.10)
* PRROC(>= 1.1)
* e1071(>= 1.6-7)
* glmnet(>= 2.0)
* pls(>= 2.5-0)
* randomForest(>= 4.6-12)
* rrBLUP(>= 4.4)
* snowfall(>= 1.84-6.1)
* spls(>= 2.2-1)
<br>
Please assure that you have installed all depends packages before you use G2P <br>

## Suggests
* rgl(>= 0.97.0)
* pheatmap(>= 1.0.8)
<br>
This two packages for result display, if you want use function "result_diplay()" in G2P,please assure you have
install _"rgl"_ and _"pheatmap"_. <br>

## Contents

#### Main functions
* GS Data summary <br>
* Example data <br>
* Trainning model  <br>
* Performance assement <br>
* Cross validation <br>
* Feature selection <br>
* Results display <br>
* [User manual](https://github.com/cma2015/DeepGS/blob/master/DeepGS.pdf)<br>

#### GS algorithms (15)
BayesA, BayesB, BayesC, BRR, BL, RKHS, RR, rrBLUP, GBLUP, RFC, RFR, SVC, SVR, SPLS, LASSO, BRNN <br>

#### Evaluation measurements (13)
Pearson correlation, Kendall rank correlation, Spearman Correlation, Mean squared error (MSE), R2,
Normalized discounted cumulative gain (NDCG), meanNDCG, AUC, AUCpr, F1, Kappa, Relative efficiency (RE), Accuracy <br>

## Quick start
More details please see download  [G2P]() see help in R.<br>

#### Do cross validation by using _G2P_ function 
```R
predlist <- G2P(cross = 10,seed = 1 , cpus = 3, markers = Markers,
                pheVal = phenotype, modelMethods = c("rrBLUP", "RFC"),
                outputModel = FALSE)
```
#### Modeling 
```R
## Fit regression model ##
rrBLUP_model <- GSReModel(markers = Markers,pheVal = phenotype,modelMethods = "rrBLUP")
## Fit classification model(RFC) ##
machine_model <- GSmachine(markers = Markers, pheVal = phenotype, modelMethods = "RFC",
                           posPercentage = 0.4, ntree = 500)
```
#### predict 
```R
res <- predictGS(testMat = Markers[1:20,], trainModel = rrBLUP_model, modelMethods = "rrBLUP")
```
#### Performance assement
```R
########## predicting breeding value
predlist <-  G2P(cross = 10,seed = 1 ,cpus = 3,markers  = Markers,pheVal  = phenotype,
                 modelMethods = c("rrBLUP","RFC"),outputModel = FALSE)
predMartix <- NULL
for(ii in 1:10){predMartix <- rbind(predMartix,predlist[[ii]])}
######## evaluate the accuracy of the prediction result

evaluareTest <- evaluateGS(realScores = predMartix[,1], predScores = predMartix[,2:3], 
                           evaMethod = c("pearson", "kendall","spearman","RE","Kappa",
                                         "auc","AUCpr","NDCGEvaluation","meanNDCGEvaluation",
                                         "MSE","R2","F1","accuracy"),topAlpha = 1:90)
```

## Ask questions
Please use [DeepGS/issues](https://github.com/cma2015/G2P/issues) for how to use DeepGS and reporting bugs
