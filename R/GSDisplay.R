###############
color <- function(){
  color = colorRampPalette(rev(palette()))(8000)
  color
}
########################## exhibit evaluation result and data structure###########################
#' @title The Visualization of Evaluation Result and Data Structure 
#' @description This function is designed for visualization of multiple results of G2P.
#' @usage rowDataPlot(y,show.line = T,barCol = "blue",lineCol = "red")
#'        rowDataPlot(markers,y,plot.type = "PCA")
#'        scatterPlot(predmat,x1 ,x2 = ,show.line = F,color_szie = T,make.plotly = F,sizeRange = c(4,6))
#'        linePlot(evalMat,size = 1)
#'        barPlot(data,other = "sector")
#'        heatmapEval <- heatMapDataProcess(x,highBound = 0,lowBound = -30,alpha = 15, basedMethod = "best")
#' @return  
#' plotd of result
#' @export 
#' @author Chuang Ma, Qian Cheng, Zhixu Qiu, Jie Song
#' @keywords plot,result visualization
#' @examples 
#' \dontrun{
#' ############# PCA analysis ############
#' data(GYSS)
#' G2PCVRes <-  G2PCrossValidation(cross = 10,seed = 1 , cpus = 3, markers = Markers,
#' pheVal = phenotype, modelMethods = c("BayesA", "BayesB", "BayesC", "BL", "BRR","RR",
#'                      "RKHS","rrBLUP","LASSO","SPLS","bigRR","SVC","RFC","SVR","RFR"),
#' outputModel = FALSE)
#' CVres <- resultsMerge(predList = G2PCVRes)
#' evalTest <- evaluateGS(realScores = CVres[,1], predScores = CVres[,2:20], 
#'                       evalMethod = c( "pearson", "kendall","spearman", "RE", "Kappa",
#'                                       "AUC", "AUCpr", "NDCG", "meanNDCG",
#'                                       "MSE", "R2", "F1", "accuracy"), topAlpha = 1:90)
#' ### row data visulization
#' ## phenotype distribution  plot
#' rowDataPlot(y = phenotype,show.line = T,barCol = "blue",lineCol = "red")
#' ## PCA 3-D plot 
#' htmlwidgets::saveWidget(as_widget(rowDataPlot(markers = Markers,y = phenotype,
#'                       plot.type = "PCA")), file="3-D_PCA.html",selfcontained=T)
#'
#' ### scatter plot 
#' scatterPlot(CVres,x1 = "BayesA",x2 = "RFC",show.line = F,color_szie = T,make.plotly = F,
#'            sizeRange = c(4,6))
#' ### lines plot 
#' linePlot(evalMat = evalTest$RE,size = 1)
#' ### bar plot 
#' barPlot(evalTest$corMethosds,other = "sector")
#' ### heat map 
#' #### pred res heatmap 
#' heatmapPlot(predmat = CVres,make.plotly = F,col.low = "green",col.high = "red")
#' #### eval res heatmap 
#' heatmapEval <- heatMapDataProcess(x = evalTest,highBound = 0,lowBound = -30,alpha = 15,
#'                                  basedMethod = "best")
#' heatmapPlot(predmat = heatmapEval,make.plotly = F,col.low = "green",col.high = "red")
#' }

rowDataPlot <- function (y, markers, make.plotly = T, plot.type = "density",pointSize = 5, 
                         show.line = T, title = "Density of switch points",bins = 30,color_low ='#FFE1A1',color_high= '#683531',barCol = "blue",lineCol = "red", ...) 
{
  if(plot.type %in% c("density", "frequency")){
    plot.type <- match.arg(plot.type, c("density", "frequency","PCA-3-D"))
    data2plot <- data.frame(phenotype = y,a = rep("a",length(y)))
    if (plot.type == "density") {
      g <- ggplot(data2plot, aes(y, fill = a)) + geom_histogram(aes(y = ..density..), 
                                                                bins = bins, closed = "left") + theme_bw() + 
        theme(legend.position = "none") + labs(y = "Density") +
        scale_fill_manual(values = barCol)
      
      if (show.line) 
        g <- g + stat_density(geom = "line", size = 1, color = lineCol)
    }
    else {
      g <- ggplot(data2plot, aes(y, fill = a)) + geom_histogram(aes(y = ..count..), 
                                                                bins = bins, closed = "left") + theme_bw() + 
        theme(legend.position = "none")  + labs(y = "Count")+
        scale_fill_manual(values = barCol)
      if (show.line) 
        g <- g + geom_freqpoly(binwidth = (range(y)[2]-range(y)[1])/bins, color = lineCol, 
                               size = 1)
    }
    g <- g + coord_cartesian(xlim = c(range(phenotype)[1],range(phenotype)[2])) + 
      labs(x = "Phenotypic Values ", title = title)
    if (make.plotly) 
      plotly::ggplotly(g)
    else g
  }else{
    ## PCA compute:
    PCA_Marker <- prcomp(markers)
    ############ kmean for cluster
    PCA_Marker <- PCA_Marker$x[,1:3]
    PCA_Marker <- as.data.frame(PCA_Marker)
    PCA_Marker <- cbind(PCA_Marker,y)
    PCA.ana <- PCA_Marker
    
    g <- plot_ly(PCA.ana, x = ~PC1, y = ~PC2, z = ~PC3, 
                 marker = list(size = pointSize,color = ~phenotype, colors = c(color_low,color_high), showscale = TRUE)) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3')),
             annotations = list(
               x = 1.13,
               y = 1.05,
               text = 'Phenotype',
               xref = 'paper',
               yref = 'paper',
               showarrow = FALSE
             ))
  }
  g
}

scatterPlot <- function (predmat,x1,x2,col.low = "blue",col.high = "red",col.mid = NULL,col = "blue",show.line = FALSE,
                         color_szie = FALSE,alpha = 0.8,make.plotly = TRUE,sizeRange) 
{
  plot.dataframe <- as.data.frame(predmat[,c(x1,x2)])
  namesData <- colnames(plot.dataframe)
  colnames(plot.dataframe) <- c("x1","x2")
  if (color_szie) {
    g <- ggplot(plot.dataframe, aes(x = x1, y = x2,color = x1,size = x2)) + 
      geom_point(alpha = alpha) + 
      scale_colour_gradientn(colours = c(col.low,col.mid,col.high))+
      labs(x = namesData[1],y = namesData[2],color = namesData[1], size = namesData[2])+
      scale_size(name = namesData[2],range = sizeRange)
  }else{
    g <- ggplot(plot.dataframe, aes(x = x1, y = x2)) + 
      geom_point(alpha = alpha,color = col,size = sizeRange[1])+
      labs(x = namesData[1],y = namesData[2],color = namesData[1], size = namesData[2])
  }
  if(show.line){
    g <- g + geom_smooth(method = "lm",se = F)
  }
  g
  if (make.plotly){ 
    plotly::ggplotly(g)}else{g}
}

heatmapPlot <- function (predmat,col.low = "blue",col.high = "red",col.mid = NULL,
                         x_cluster = "euclidean",y_cluster = "euclidean",make.plotly = TRUE,
                         xlab = "x",ylab = "y",legend.title ="value",boundSize = 0
) 
{
  require(reshape)
  xOrder <- hclust(dist(predmat,method = x_cluster))$order
  yOrder <- hclust(dist(t(predmat),method = y_cluster))$order
  plot.dataframe <- predmat[xOrder,yOrder]
  plot.dataframe <- melt(plot.dataframe)
  g <- ggplot(plot.dataframe, aes(X1,X2))+ 
    geom_tile(aes(fill = value),colour = "white",size=boundSize)+ 
    scale_fill_gradientn(colours = c(col.low,col.mid,col.high))+
    labs(x = xlab,y = ylab,fill = legend.title)
  if (make.plotly){ 
    plotly::ggplotly(g)}else{g}
}

linePlot <- function(evalMat,xlab = "Top alpha (%)",ylab = "Value",legendTittle = "GS Method",
                     size = 2,alpha = 0.8,make.plotly = FALSE){
  color <- rainbow(ncol(evalMat))
  times <- nrow(evalMat)
  plot.dataframe <- melt(evalMat)
  g <- ggplot(plot.dataframe, aes(x = rep(1:times,ncol(evalMat)), y = value, colour=X2,group=X2)) + 
    geom_line(size=size,alpha = alpha)+
    labs(x = xlab,y = ylab,colour = legendTittle)+
    scale_colour_manual(values = color)
  if (make.plotly){ 
    plotly::ggplotly(g)}else{g}
}

barPlot <- function(data,
                    xlab = "GS Methods",ylab = "Value",legendTittle = "Measures",
                    make.plotly = FALSE,other = "normal"){
  color <- rainbow(nrow(data))
  plot.dataframe <- melt(data)
  g <- ggplot(plot.dataframe, aes(x = X2, y = value, fill=X1)) + 
    geom_bar(stat="identity", position="dodge")+
    labs(x = xlab,y = ylab,fill = legendTittle)
  scale_colour_manual(values= color)
  switch (other,
          normal = {g <- g},
          parallel = {g <- g + coord_flip()},
          sector = {g <- g + coord_polar(theta = "x")},
          pie = {g <- g + coord_polar(theta = "y")}
  )
  
  if (make.plotly){ 
    plotly::ggplotly(g)}else{g}
}

