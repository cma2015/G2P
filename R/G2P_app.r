################################################### G2P #########################################################
#' @title Genotypes to phenotypes with Shiny App
#'
#' @return the Shiny App.
#'
#' @examples TSIS.app()
#'
#' @seealso \code{\link{shiny}}
#'
#' @export
G2P.app <- function() {
  
  if(!file.exists('tutorial'))
    dir.create('tutorial')
  if(file.exists('tutorial/tutorial.md'))
    file.remove('tutorial/tutorial.md')
  
  download.file('https://github.com/cma2015/G2P/blob/master/G2P.Rmd',
                destfile = 'tutorial/tutorial.md',quiet = T)
#   # unzip('tutorial/tutorial-shiny.zip',exdir = 'tutorial')
#   invisible(file.remove('tutorial/tutorial-shiny.zip'))
  
  # require(tools)
  require(shiny)
  # require(shinyFiles)
  library(shinythemes)
  # library(TSIS)
  require(plotly)
  
  message('Starting G2P App...')
  
  shinyApp(options = list(launch.browser=T),
           ui = navbarPage("Genotype to Phenotype[G2P]",
                           ##Page 2
                           tabPanel("G2P data input and quality control", #tabPanel the 1st module and title
                                    fluidPage(
                                      fluidRow(
                                        ##input gene expression data part
                                        column(3,                         # structure of window
                                               titlePanel('Input GS data files'),  # 2st module
                                               
                                               
                                               wellPanel( #inside panel
                                                 
                                                 ##input isoform expression data
                                                 h4('GS markers data:'),
                                                 fileInput('markersData','Select the GS markers data file',
                                                           accept = c(
                                                             'text/csv',
                                                             'text/comma-separated-values',
                                                             'text/tab-separated-values',
                                                             'text/plain',
                                                             '.csv',
                                                             '.tsv'
                                                           )),
                                                 p('The GS markers input is a data matrix with columns of markers(SNPs or DArTs) and rows of samples.'),
                                                 br(),
                                                 
                                                 ##phenotypic data
                                                 h4('GS phenotypic data:'),
                                                 fileInput('phenotypeData','Select the phenotypic data file',
                                                           accept = c(
                                                             'text/csv',
                                                             'text/comma-separated-values',
                                                             'text/tab-separated-values',
                                                             'text/plain',
                                                             '.csv',
                                                             '.tsv'
                                                           )),
                                                 p('The mapping input is a data array with phenotype.'),
                                                 br(),
                                                 
                                                 ##input subset of isoforms for investigation
                                                 h4('Other incidence matrices:'),
                                                 fileInput('file.Z','Select the fixed effects matrix data file',
                                                           accept = c(
                                                             'text/csv',
                                                             'text/comma-separated-values',
                                                             'text/tab-separated-values',
                                                             'text/plain',
                                                             '.csv',
                                                             '.tsv'
                                                           )),
                                                 HTML('The fixed effects incidence matrices for model.')
                                                 
                                               )
                                        ),
                                        column(9,
                                               titlePanel('Visualization and summary of input phenotype data'),
                                               wellPanel(
                                                 sidebarLayout(
                                                   sidebarPanel(
                                                     fluidRow(
                                                       column(6,
                                                              radioButtons("plotInputType", label = h5("Select plot type:"),
                                                                           choices = list("Frequency" = 'frequency',"Density" = 'density',"PCA-3-D" = 'PCA-3-D'),
                                                                           selected = 'frequency',inline = F)
                                                       ),
                                                       column(6,
                                                              radioButtons("showDensityLine", label = h5("Density/frequency line:"),
                                                                           choices = list("FALSE" = 'FALSE',"TRUE" = 'TRUE'),
                                                                           selected = 'FALSE',inline = F)
                                                       )
                                                     ),
                                                     
                                                     fluidRow(
                                                       
                                                       column(6,
                                                              radioButtons("save.input.type", label = h5("Select format to save:"),
                                                                           choices = list("html" = 'html', "png" = "png", "pdf" = 'pdf'),
                                                                           selected = 'html',inline = T)
                                                       ),
                                                       column(6,
                                                              br(),
                                                              downloadButton('download.input.plot', 'Save',class="btn btn-primary")
                                                       ),
                                                       column(12,
                                                              h4('Selection:'),
                                                              sliderInput("binsN", "Number of bins to view:", min = 0, max = 100,value = 20),
                                                              sliderInput("D3_size", "Points size:", min = 0, max = 20,value = 5,sep = 0.5),
                                                              conditionalPanel(
                                                                condition = "input.plotInputType == 'PCA-3-D'",
                                                                sliderInput('D3.color.low', label = 'Color low: ',min = 0,max = 8000,value = 3000),
                                                                sliderInput('D3.color.high', label = 'Color high: ',min = 0,max = 8000,value = 7000)
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.plotInputType == 'frequency'|input.plotInputType == 'density'",
                                                                sliderInput('barPlotCol', label = 'Color of bars',min = 0,max = 8000,value = 3000),
                                                                conditionalPanel(
                                                                  condition = "input.showDensityLine == 'TRUE'",
                                                                  sliderInput('barPlotLineCol', label = 'Color of lines',min = 0,max = 8000,value = 3000)
                                                                )
                                                              )
                                                       )
                                                     )
                                                   ),
                                                   mainPanel(
                                                     wellPanel(
                                                       plotlyOutput('density',height = "80%",width = "100%"),
                                                       
                                                       HTML('<b>Figure:</b> Density plot of phenotype or PCA analysis plot. The plot is made based
                                                      on the disbution of phenoypes.')),
                                                     wellPanel(
                                                       # h4('Phenotype summary:'),
                                                       verbatimTextOutput("summary"),
                                                       HTML('<b>Summary:</b> Basic information of phenotypes.')
                                                     )
                                                   )
                                                   
                                                 )))
                                        ,
                                        column(12,
                                               titlePanel('Data quality control'),
                                               wellPanel(
                                                 sidebarLayout(
                                                   sidebarPanel(width = 4,
                                                                fluidRow(
                                                                  column(12,
                                                                         #                                                                      h4('Global:'),
                                                                         #                                                                      verbatimTextOutput("Global.summary"),
                                                                         wellPanel(
                                                                           h4('parameter setting'),
                                                                           fluidRow(
                                                                             br(),
                                                                             column(12,
                                                                                    actionButton('inputQC','Parameter',width = "100%")
                                                                             ),
                                                                             column(12,
                                                                                    HTML('<b>Press GSDataQC button to set the parameters.</b>')      
                                                                             ),
                                                                             br(),
                                                                             wellPanel(" ",
                                                                                       fluidRow(
                                                                                         htmlOutput("QCparaSelect")
                                                                                       )
                                                                             ),
                                                                             column(12,
                                                                                    actionButton('inpute.GSDataQC','GSDataQC',icon("send outline icon"),class="btn btn-primary",width = "100%")
                                                                             ),
                                                                             column(12,
                                                                                    HTML('<b>Press GSDataQC button to control the raw data.</b>')
                                                                             )
                                                                           )
                                                                         ),
                                                                         wellPanel(
                                                                           fluidRow(
                                                                             h4('Resluts display'),
                                                                             column(12,
                                                                                    selectInput('QCres.output',"QC imformation:",c("Imputation matrix","MAF","NA in column","NA in row","NA index"))         
                                                                             ),
                                                                             column(12,
                                                                                    h4('selection'),
                                                                                    numericInput("QCres.obsR", "Number of rows to view:", 10),
                                                                                    numericInput("QCres.obsC", "Number of columns to view:", 10)
                                                                             )
                                                                           )
                                                                         )
                                                                         ,
                                                                         wellPanel(
                                                                           fluidRow(
                                                                             h4('Download'),
                                                                             column(6,
                                                                                    div(align='left',downloadButton('download.QC.inputation', 'Download inputation matrix',class="btn btn-primary")),
                                                                                    br()
                                                                             ),
                                                                             
                                                                             column(6,
                                                                                    div(align='left',downloadButton('download.QC.MAF', 'Download MAF',class="btn btn-primary")),
                                                                                    br()
                                                                             ),
                                                                             
                                                                             column(6,
                                                                                    div(align='left',downloadButton('download.QC.NACol', 'Download NA column',class="btn btn-primary")),
                                                                                    br()
                                                                             ),
                                                                             column(6,
                                                                                    div(align='left',downloadButton('download.QC.NARow', 'Download NA row',class="btn btn-primary")),
                                                                                    br()
                                                                             ),
                                                                             column(6,
                                                                                    div(align='left',downloadButton('download.QC.NAIndex', 'Download NA index',class="btn btn-primary")),
                                                                                    br()
                                                                             ) 
                                                                           )
                                                                         )
                                                                  )
                                                                )
                                                                #                                                               )
                                                                #                                                             )
                                                   ),
                                                   mainPanel(width = 8,
                                                             wellPanel(
                                                               tabsetPanel(type = "tabs",
                                                                           tabPanel("Markers data",
                                                                                    HTML('<b>Table:</b> The raw data of GS markers'),
                                                                                    tableOutput("markerView")
                                                                           ),
                                                                           tabPanel("Visualization after data QC",
                                                                                    HTML('<b>Table:</b> The GS data quliaty analysis results'),
                                                                                    br(),
                                                                                    tableOutput('GSDataQC.table')
                                                                           )
                                                               )
                                                             )
                                                   )
                                                 )
                                               )
                                        )
                                      )
                                    )
                                    # )
                                    
                           ),
                           ## page3
                           tabPanel("Genotype to phenotype[modeling,predict and cross validation]",
                                    fluidPage(
                                      fluidRow(
                                        #                   tabsetPanel(type = "tabs",
                                        #                               tabPanel("Parameter setting",
                                        column(12,
                                               titlePanel("G2P and Cross validattion"),
                                               wellPanel(
                                                 sidebarLayout(
                                                   sidebarPanel(width = 3,
                                                                fluidRow(
                                                                  column(12,
                                                                         h4('Input training sets boundary:'),
                                                                         numericInput("left.T", "TS start boundary: ", 0),
                                                                         numericInput("right.T", "TS end boundary:", 0),
                                                                         h4('Input testing sets boundary:'),
                                                                         numericInput("left.V", "VS start boundary: ", 0),
                                                                         numericInput("right.V", "VS end boundary:", 0)
                                                                         
                                                                  )),
                                                                fluidRow(
                                                                  column(12,
                                                                         br(),
                                                                         column(6,
                                                                                div(align = "left",actionButton('G2P.run','G2P',icon("send outline icon"),class="btn btn-primary")),
                                                                                br()
                                                                         ),
                                                                         column(6,
                                                                                div(align = "left",actionButton('G2PCV.run','G2P-CV',icon("send outline icon"),class="btn btn-primary")),
                                                                                br()
                                                                         ),
                                                                         br(),
                                                                         column(6,
                                                                                div(align = "left",actionButton('eval.run.G2P','EvaluationG2P',icon("send outline icon"),class="btn btn-primary"))
                                                                         ),
                                                                         column(6,
                                                                                div(align = "left",actionButton('eval.run.CV','EvaluationCV',icon("send outline icon"),class="btn btn-primary"))
                                                                         )
                                                                  )
                                                                )
                                                   )
                                                   ,
                                                   mainPanel(width = 9,
                                                             fluidRow(
                                                               h3('Parameters setting:'),
                                                               column(8,
                                                               column(12,
                                                                      div(align = "middle",actionButton('inputG2P','G2P and Cross Validation Parameters',width = "50%"))
                                                               ),
                                                               column(12,
                                                                      HTML('<b>Press Parameters button to set the parameters.</b>')      
                                                               ),
                                                               br(),
                                                               wellPanel(" ",
                                                                         fluidRow(
                                                                           column(12,
                                                                                  htmlOutput("G2PparaSelect")  
                                                                           ),
                                                                           column(6,
                                                                                  htmlOutput("G2PparaSelect1")  
                                                                                  ),
                                                                           column(6,
                                                                                  htmlOutput("G2PparaSelect2")
                                                                                  )
                                                                         )
                                                               )
                                                               ),
                                                               column(4,
                                                               column(12,
                                                                      actionButton('inputEval','Evaluation Parameters',width = "100%")
                                                               ),
                                                               column(12,
                                                                      HTML('<b>Press Parameters button to set the parameters.</b>')      
                                                               ),
                                                               br(),
                                                               wellPanel(" ",
                                                                         fluidRow(
                                                                           column(12,
                                                                           htmlOutput("EvalparaSelect")
                                                                           )
                                                                         )
                                                               )
                                                               ),
                                                               column(8,
                                                                      htmlOutput("G2PCVprocess")
                                                                      )
                                                             )
                                                   )
                                                 )
                                               )
                                        )
                                        # ),
                                        # tabPanel("Results display",
                                        #                               )
                                        #                   )
                                      )
                                      # )
                                    )
                           ),
                           ## page 4
                           tabPanel("Visualization",
                                    fluidPage(
                                      fluidRow(
                                        ##input gene expression data part
                                        tabsetPanel(type = "tabs",
                                                              tabPanel("Plot",
                                                                       column(4,                         # structure of window
                                                                              titlePanel('Switch plots and charts'),  # 2st module
                                                                              
                                                                              wellPanel( #inside panel
                                                                                fluidRow(
                                                                                  ##input isoform expression data
                                                                                  h3('Select plot data:'),
                                                                                  selectInput('plotData','Plot data',c("Prediction results","Global evaluation results","Threshold evaluation results")),
                                                                                  conditionalPanel(
                                                                                    condition = "input.plotData == 'Prediction results'",
                                                                                    selectInput("predPlotType", "Chart selection:",c("Scatter","Heatmap","Density")),
                                                                                    conditionalPanel(condition = "input.predPlotType == 'Scatter'",
                                                                                                     wellPanel(
                                                                                                       h4('Plot parameter setting:'),
                                                                                                       fluidRow(
                                                                                                         column(6,
                                                                                                                textInput('method1', label = 'Method 1: ', value = ' '),
                                                                                                                selectInput('scatter.showline',"Show the fitting line?",c("FALSE","TRUE")),
                                                                                                                sliderInput("scatter.alpha", "Diaphaneity", min = 0, max = 1,value = 0.8,sep = 0.01)
                                                                                                         ),
                                                                                                         column(6,
                                                                                                                textInput('method2', label = 'Method 2: ', value = ' '),
                                                                                                                selectInput('scatterColor_size',"Change of gradient?",c("TRUE","FALSE")),
                                                                                                                sliderInput("scatter.sizeRange", "Points size", sep = 0.05,min = 0, max = 20,value = c(2,5))
                                                                                                         ),
                                                                                                         conditionalPanel(
                                                                                                           condition = "input.scatterColor_size == 'TRUE'",
                                                                                                           column(12,
                                                                                                                  sliderInput('scatter.col.low', label = 'Color for min value:', min = 0,max = 8000,value = 5000),
                                                                                                                  sliderInput('scatter.col.mid', label = 'Color for median value:', min = 0,max = 8000,value = 0),
                                                                                                                  sliderInput('scatter.col.high', label = 'Color for max value:', min = 0,max = 8000,value = 7000)
                                                                                                           )
                                                                                                         ),
                                                                                                         conditionalPanel(
                                                                                                           condition = "input.scatterColor_size == 'FALSE'",
                                                                                                           column(12,
                                                                                                                  sliderInput('scatter.color',label = 'Color for points:',min = 0,max = 8000,value = 5000)
                                                                                                           )
                                                                                                         )
                                                                                                       ) 
                                                                                                     )
                                                                                                     
                                                                                    ),
                                                                                    # Only show this panel if Custom is selected
                                                                                    conditionalPanel(
                                                                                      condition = "input.predPlotType == 'Desity'",
                                                                                      
                                                                                      selectInput("predres.select.method", "Select the methods", c("BayesA","BayesB","RFR","SVC"))
                                                                                    ),
                                                                                    conditionalPanel(
                                                                                      condition = "input.predPlotType == 'Heatmap'",
                                                                                      tabsetPanel(type = "tabs",
                                                                                                  tabPanel("Basic",
                                                                                                           column(6,
                                                                                                                  selectInput('predHeatmap.x_cluster',"X cluster method:",c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  selectInput('predHeatmap.y_cluster',"Y cluster method:",c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) 
                                                                                                           ),
                                                                                                           column(12,
                                                                                                                  sliderInput('predHeatmap.col.low', label = 'Color for min value:', min = 0,max = 8000,value = 5000),
                                                                                                                  sliderInput('predHeatmap.col.mid', label = 'Color for median value:', min = 0,max = 8000,value = 0),
                                                                                                                  sliderInput('predHeatmap.col.high', label = 'Color for max value:', min = 0,max = 8000,value = 7000),
                                                                                                                  sliderInput('predHeatmap.boundSize', label = 'Bound line for heatmap:', min = 0,max = 1,value = 0,sep = 0.01)
                                                                                                           )
                                                                                                  ),
                                                                                                  tabPanel("Theme",
                                                                                                           column(6,
                                                                                                                  textInput('predHeatmap.xlab', label = 'X axis label ', value = 'X')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('predHeatmap.ylab', label = 'Y axis label ', value = 'Y')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('predHeatmap.legend.title', label = 'Legend title', value = 'Value')
                                                                                                           )
                                                                                                  )
                                                                                      )
                                                                                    )
                                                                                  ),
                                                                                  conditionalPanel(
                                                                                    condition = "input.plotData == 'Global evaluation results'",
                                                                                    selectInput("golbalPlotType", "Chart selection:",c("Barplot")),
                                                                                    conditionalPanel(condition = "input.golbalPlotType == 'Barplot'",
                                                                                                     tabsetPanel(type = "tabs",
                                                                                                                 tabPanel("Basic",
                                                                                                                          column(12,
                                                                                                                                 selectInput("globalBar.type","The type of bar plot",c("normal","parallel","sector","pie"))
                                                                                                                          )
                                                                                                                 ),
                                                                                                                 tabPanel("Theme",
                                                                                                                          column(6,
                                                                                                                                 textInput('globalBar.xlab', label = 'X axis label ', value = 'X')
                                                                                                                          ),
                                                                                                                          column(6,
                                                                                                                                 textInput('globalBar.ylab', label = 'Y axis label ', value = 'Y')
                                                                                                                          ),
                                                                                                                          column(6,
                                                                                                                                 textInput('globalBar.legend.title', label = 'Legend title', value = 'Value')
                                                                                                                          )
                                                                                                                 )
                                                                                                     )
                                                                                    )
                                                                                  ),
                                                                                  conditionalPanel(
                                                                                    condition = "input.plotData == 'Threshold evaluation results'",
                                                                                    selectInput("thresholdPlotType", "Chart selection:",c("Line","Heatmap","AUC","AUCpr")),
                                                                                    conditionalPanel(
                                                                                      condition = "input.thresholdPlotType == 'Line'",
                                                                                      tabsetPanel(type = "tabs",
                                                                                                  tabPanel("Basic",
                                                                                                           column(12,
                                                                                                                  selectInput("thresholdLine.measure","Which measure be plot?",c("RE","Kappa","AUC","AUCpr","F1","accuracy","NDCG","meanNDCG")),
                                                                                                                  sliderInput("thresholdLine.width","Size of lines ",min = 0.5,max = 10,value = 2),
                                                                                                                  sliderInput("thresholdLine.alpha","Diaphaneity ",min = 0,max = 1,value = 0.8,sep = 0.01)
                                                                                                           )
                                                                                                  ),
                                                                                                  tabPanel("Theme",
                                                                                                           column(6,
                                                                                                                  textInput('thresholdLine.xlab', label = 'X axis label ', value = 'Top alpha (%)')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('thresholdLine.ylab', label = 'Y axis label ', value = 'Measure value')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('thresholdLine.legend.title', label = 'Legend title', value = 'GS methods')
                                                                                                           )
                                                                                                  )
                                                                                      )
                                                                                    ),
                                                                                    conditionalPanel(
                                                                                      condition = "input.thresholdPlotType == 'Heatmap'",
                                                                                      tabsetPanel(type = "tabs",
                                                                                                  tabPanel("Basic",
                                                                                                           column(6,
                                                                                                                  selectInput('thresholdHeatmap.x_cluster',"X cluster method:",c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  selectInput('thresholdHeatmap.y_cluster',"Y cluster method:",c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) 
                                                                                                           ),
                                                                                                           column(12,
                                                                                                                  sliderInput('thresholdHeatmap.col.low', label = 'Color for min value:', min = 0,max = 8000,value = 5000),
                                                                                                                  sliderInput('thresholdHeatmap.col.mid', label = 'Color for median value:', min = 0,max = 8000,value = 0),
                                                                                                                  sliderInput('thresholdHeatmap.col.high', label = 'Color for max value:', min = 0,max = 8000,value = 7000),
                                                                                                                  sliderInput('thresholdHeatmap.boundSize', label = 'Bound line for heatmap:', min = 0,max = 1,value = 0,sep = 0.01)
                                                                                                           )
                                                                                                  ),
                                                                                                  tabPanel("Theme",
                                                                                                           column(6,
                                                                                                                  textInput('thresholdHeatmap.xlab', label = 'X axis label ', value = 'X')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('thresholdHeatmap.ylab', label = 'Y axis label ', value = 'Y')
                                                                                                           ),
                                                                                                           column(6,
                                                                                                                  textInput('thresholdHeatmap.legend.title', label = 'Legend title', value = 'Value')
                                                                                                           )
                                                                                                  ),
                                                                                                  tabPanel("Data Process",
                                                                                                           column(12,
                                                                                                                  selectInput('thresholdHeatmap.basedMethod', "Benchmark method",c("best","BayesA","BayesB","BayesC","rrBLUP","RFC","LASSO",
                                                                                                                                                                                   "BRR","SPLS","RFR","SVC","SVR","RR","NR","AI","BL","bigRR",
                                                                                                                                                                                   "BRNN","EM","EMMA","RKHS")),
                                                                                                                  sliderInput('thresholdHeatmap.highBound',"Upper bound (%)",min = 0,max = 100,value = 0),
                                                                                                                  sliderInput('thresholdHeatmap.lowBound',"Lower bound (%)",min = -100,max = 0,value = -15),
                                                                                                                  sliderInput('thresholdHeatmap.alpha',"Threshod top individuals (%)",min = 1,max = 99,value = 15)
                                                                                                           )
                                                                                                  )
                                                                                      )
                                                                                    )
                                                                                  )
                                                                                  #                                                   textInput('method1', label = 'Method 1: ', value = ''),
                                                                                  #                                                   textInput('method2', label = 'Method 2: ', value = ''),
                                                                                  #                                                   h4('Plot type:'),
                                                                                  #                                                   selectInput('visual.plot.type',label = "Select the plot type:",choices = c("Scatter plot","Bar plot","Heatmap","Line plot","Box plot","Pie chart")),
                                                                                  #                                                   h4('Plot type:'),
                                                                                  #                                                   selectInput('eval.data.type',label = "Select the plot type:",choices = c("Scatter plot","Bar plot","Heatmap","Line plot","Box plot","Pie chart")),
                                                                                  #                                                   
                                                                                  #                                                   h4('Results Selection'),
                                                                                  #                                                   selectInput('res.select',label = " ",choices = c("G2P results","G2P CV results"))
                                                                                  #                                                   fileInput('markersData','Select the GS markers data file',
                                                                                  #                                                             accept = c(
                                                                                  #                                                               'text/csv',
                                                                                  #                                                               'text/comma-separated-values',
                                                                                  #                                                               'text/tab-separated-values',
                                                                                  #                                                               'text/plain',
                                                                                  #                                                               '.csv',
                                                                                  #                                                               '.tsv'
                                                                                  #                                                             )),
                                                                                  #                                                   p('The GS markers input is a data matrix with columns of markers(SNPs or DArTs) and rows of samples.'),
                                                                                  #                                                   br(),
                                                                                  
                                                                                  ##phenotypic data
                                                                                  #                                                   h4('GS phenotypic data:'),
                                                                                  #                                                   fileInput('phenotypeData','Select the phenotypic data file',
                                                                                  #                                                             accept = c(
                                                                                  #                                                               'text/csv',
                                                                                  #                                                               'text/comma-separated-values',
                                                                                  #                                                               'text/tab-separated-values',
                                                                                  #                                                               'text/plain',
                                                                                  #                                                               '.csv',
                                                                                  #                                                               '.tsv'
                                                                                  #                                                             )),
                                                                                  # p('The mapping input is a data array with phenotype.')
                                                                                )
                                                                                
                                                                              )                                                
                                                                              #                                                 wellPanel(
                                                                              #                                                   ##input subset of isoforms for investigation
                                                                              #                                                   h4('Other incidence matrices:'),
                                                                              #                                                   fileInput('file.Z','Select the fixed effects matrix data file',
                                                                              #                                                             accept = c(
                                                                              #                                                               'text/csv',
                                                                              #                                                               'text/comma-separated-values',
                                                                              #                                                               'text/tab-separated-values',
                                                                              #                                                               'text/plain',
                                                                              #                                                               '.csv',
                                                                              #                                                               '.tsv'
                                                                              #                                                             )),
                                                                              #                                                   HTML('The fixed effects incidence matrices for model.')
                                                                              #                                                   
                                                                              #                                                 )
                                                                       ),                     
                                                                       column(8,
                                                                              titlePanel('Plot show window:'),
                                                                              wellPanel(
                                                                                fluidRow(
                                                                                  tabsetPanel(type = "pills",
                                                                                              tabPanel("Scatter plot",
                                                                                                       plotlyOutput('visualization.scatter',height = "100%",width = "1000px")
                                                                                              ),
                                                                                              tabPanel("Lines plot",
                                                                                                       plotlyOutput('visualization.line',height = "100%",width = "auto")
                                                                                              ),
                                                                                              tabPanel("Bar plot",
                                                                                                       plotlyOutput('visualization.bar',height = "100%",width = "100%")
                                                                                              ),  
                                                                                              tabPanel("Heatmap",
                                                                                                       plotlyOutput('visualization.heatmap',height = "100%",width = "1000px")
                                                                                              )
                                                                                  )
                                                                                )
                                                                              )
                                                                       )
                                                              ),
                                                    tabPanel("Data table",
                                                             column(12,
                                                                    titlePanel("Resluts display"),
                                                                    wellPanel(
                                                                      sidebarLayout(
                                                                        sidebarPanel(width = 3,
                                                                                     fluidRow(
                                                                                       column(12,
                                                                                              tabsetPanel(
                                                                                                type = "tabs",
                                                                                                tabPanel("Results Selection",
                                                                                                         column(12,
                                                                                                                h4('Results Selection'),
                                                                                                                selectInput('res.select',label = " ",choices = c("G2P results","G2P CV results")) 
                                                                                                         ),
                                                                                                         #                                                                          column(12,
                                                                                                         #                                                                                 selectInput('QCres.output',"QC imformation:",c("Imputation matrix","MAF","NA in column","NA in row","NA index"))         
                                                                                                         #                                                                          ),
                                                                                                         column(12,
                                                                                                                h4('Selection:'),
                                                                                                                numericInput("G2P.obsR", "Number of rows to view:", 10),
                                                                                                                numericInput("G2P.obsC", "Number of columns to view:", 5)
                                                                                                                
                                                                                                         ),
                                                                                                         column(12,
                                                                                                                h4('Threshold measures selection:'),
                                                                                                                selectInput('eval.res.select',"Measures:",c("RE","Kappa","AUC","AUCpr","NDCG","meanNDCG","F1","accuracy")) 
                                                                                                         ), 
                                                                                                         column(12,
                                                                                                                h4('Selection:'),
                                                                                                                numericInput("evalres.obsR", "Number of rows to view:", 10),
                                                                                                                numericInput("evalres.obsC", "Number of columns to view:", 5)
                                                                                                                
                                                                                                         )
                                                                                                ),
                                                                                                tabPanel("Download results",
                                                                                                         column(12,
                                                                                                                fluidRow(
                                                                                                                  column(6,
                                                                                                                         br(),
                                                                                                                         div(align='left',downloadButton('download.predres.mat', 'Prediction results',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         br(),
                                                                                                                         div(align='left',downloadButton('download.evalres.global', 'Golbal evaluation',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.RE', 'RE',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.Kappa', 'Kappa',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.AUC', 'AUC',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.AUCpr', 'AUCpr',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.NDCG', 'NDCG',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.meanNDCG', 'meanNDCG',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.F1', 'F1',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  ),
                                                                                                                  column(6,
                                                                                                                         div(align='left',downloadButton('download.evalres.accuracy', 'accuracy',class="btn btn-primary")),
                                                                                                                         br()
                                                                                                                  )
                                                                                                                )    
                                                                                                         )
                                                                                                )
                                                                                              )
                                                                                              
                                                                                              
                                                                                       )
                                                                                     )
                                                                        ),
                                                                        mainPanel(width = 9,
                                                                                  wellPanel(
                                                                                    HTML('<b>Table:</b> The G2P results'),
                                                                                    br(),
                                                                                    tableOutput('G2P.res.table')
                                                                                  ),
                                                                                  wellPanel(
                                                                                    h4('Table: The evaluation results'),
                                                                                    wellPanel(
                                                                                      HTML('<b>Table:</b> The glpbal results'),
                                                                                      br(),
                                                                                      tableOutput('eval.res.global.table')
                                                                                    ),
                                                                                    wellPanel(
                                                                                      HTML('<b>Table:</b> The threshold results'),
                                                                                      br(),
                                                                                      tableOutput('eval.res.table')
                                                                                    )
                                                                                  )
                                                                                  
                                                                        )
                                                                        
                                                                      )
                                                                    ))
                                                    )
                                        )
                                        
                                        #                                          column(8,
                                        #                                                 titlePanel('Visualization and summary of input phenotype data'),
                                        #                                                 wellPanel(
                                        #                                                   sidebarLayout(
                                        #                                                     sidebarPanel(
                                        #                                                       fluidRow(
                                        #                                                         column(6,
                                        #                                                                radioButtons("plotInputType", label = h5("Select plot type:"),
                                        #                                                                             choices = list("Frequency" = 'frequency',"Density" = 'density',"PCA-3-D" = 'PCA-3-D'),
                                        #                                                                             selected = 'frequency',inline = F)
                                        #                                                         ),
                                        #                                                         column(6,
                                        #                                                                radioButtons("showDensityLine", label = h5("Density/frequency line:"),
                                        #                                                                             choices = list("FALSE" = 'FALSE',"TRUE" = 'TRUE'),
                                        #                                                                             selected = 'FALSE',inline = F)
                                        #                                                         )
                                        #                                                       ),
                                        #                                                       
                                        #                                                       fluidRow(
                                        #                                                         
                                        #                                                         column(6,
                                        #                                                                radioButtons("save.input.type", label = h5("Select format to save:"),
                                        #                                                                             choices = list("html" = 'html', "png" = "png", "pdf" = 'pdf'),
                                        #                                                                             selected = 'html',inline = T)
                                        #                                                         ),
                                        #                                                         column(6,
                                        #                                                                br(),
                                        #                                                                downloadButton('download.input.plot', 'Save',class="btn btn-primary")
                                        #                                                         ),
                                        #                                                         column(12,
                                        #                                                                h4('selection:'),
                                        #                                                                sliderInput("binsN", "Number of bins to view:", min = 0, max = 100,value = 20),
                                        #                                                                h4('Phenotype summary:'),
                                        #                                                                verbatimTextOutput("summary")
                                        #                                                         ),
                                        #                                                         column(12,
                                        #                                                                br(),
                                        #                                                                br(),
                                        #                                                                br(),
                                        #                                                                br(),
                                        #                                                                br(),
                                        #                                                                br()
                                        #                                                         )
                                        #                                                         
                                        #                                                       )
                                        #                                                     ),
                                        #                                                     mainPanel(
                                        #                                                       wellPanel(
                                        #                                                         plotlyOutput('density',height = "80%",width = "100%"),
                                        #                                                         
                                        #                                                         HTML('<b>Figure:</b> Density plot of phenotype or PCA analysis plot. The plot is made based
                                        #                                                              on the disbution of phenoypes.')))
                                        #                                                         )))
                                        
                                      )
                                    )
                                    
                           ),
                           #Page 4
                           tabPanel("Manual", 
                                    htmlOutput("tutorial")
                           )
           ),
           ################################# server function ########################################################
           server = function(input, output,session) {
             
             ######################################## windows #################################
             dataModal <- function() {
               modalDialog(easyClose = TRUE,
                           tabsetPanel(type = "tabs",
                                       tabPanel("Parameter setting",
                                                wellPanel(
                                                  fluidRow(
                                                    column(6,
                                                           selectInput('if.inpute',"Inpute the marker data",c("FALSE","TRUE"))      
                                                    ),
                                                    column(6,
                                                           selectInput('inpute.method',"Inpute method",c("mean","median","KNN"))         
                                                    ),
                                                    column(6,
                                                           numericInput('inpute.round',label = 'Round',value = 4)
                                                    ),
                                                    column(6,
                                                           numericInput('inpute.k',label = 'K',value = 10)
                                                    ),
                                                    column(6,
                                                           numericInput('inpute.maxp',label = "Maxp",value = 1500)
                                                    ),
                                                    column(6,
                                                           numericInput('inpute.rowmax',label = "Rowmax",value = 0.5)
                                                    ),
                                                    column(6,
                                                           div(align = "left", numericInput('inpute.colmax',label = "Colmax",value = 0.8))
                                                    ),
                                                    column(6,
                                                           div(align = "left", numericInput('inpute.rng.seed',label = "RNG.seed",value = 362436069))
                                                    )
                                                  )
                                                )
                                       ),
                                       tabPanel("Parameter arguments",
                                                HTML('The details of parameters:
                                          <ul><li><b>Inpute the marker data:</b>(logical)if TRUE, imputation, default F.</li>
                                          <li><b>Inpute method:</b>(character)the method of imputation, "mean", "median" or "KNN", default "mean".</li>
                                          <li><b>Round:</b>(numeric)rounds the values in its first argument to the specified number of decimal places, default 4 </li>
                                          <li><b>Min time in interval:</b> The minimum time points for both columns "before.t.points" and "after.t.points" in the output table.</li>
                                          <li><b>K, Maxp, Rowmax, Colmax, RNG.seed:</b>(numeric)KNN method parameters, details see <a href="http://www.bioconductor.org/packages/release/bioc/manuals/impute/man/impute.pdf">KNN</a>.</li>
                                          </ul>')
                                       )
                           ),
                           footer = tagList(
                             modalButton("Cancel"),
                             actionButton("okQC", "OK")
                           )
                           # modalButton("OK"),
                           
                           #           if (failed)
                           #             div(tags$b("Invalid name of data object", style = "color: red;")),
                           
               )
             }
             
             observeEvent(input$inputQC, {
               showModal(dataModal())
             })
             
               observeEvent(input$okQC, {
                 removeModal()
               })
             
             output$QCparaSelect <- renderPrint({
               HTML('<h5><b>The setting of parameters:</h5></b>',
                    '<ul><li><b>Inpute the marker data:</b>',input$if.inpute,'</li>',
                    '<li><b>Inpute method:</b>',input$inpute.method,'</li>',
                    '<li><b>Round:</b>',input$inpute.round,'</li>',
                    '<li><b>K:</b>',input$inpute.k,'</li>',
                    '<li><b>Maxp:</b>',input$inpute.maxp,'</li>',
                    '<li><b>Rowmax:</b>',input$inpute.rowmax,'</li>',
                    '<li><b>Colmax:</b>',input$inpute.colmax,'</li>',
                    '<li><b>RNG.seed:</b>',input$inpute.rng.seed,'</li>',
                    '</ul>')
               
             })
## G2P and CV 
             G2PModal <- function() {
               modalDialog(easyClose = TRUE,size = "l",
                           tabsetPanel(type = "tabs",
                                       tabPanel("Parameter setting",
                                                wellPanel(
                                                  fluidRow(
                                                    column(12,
                                                           checkboxGroupInput("model.group.R",label = h4("Statistics models"),inline = T,
                                                                              choices = list(
                                                                                "Bayes A" = "BayesA",
                                                                                "Bayes B" = "BayesB",
                                                                                "Bayes C" = "BayesC",
                                                                                "RR-BLUP" = "rrBLUP",
                                                                                "LASSO  " = "LASSO",
                                                                                "BRR    " = "BRR",
                                                                                "SPLS   "= "SPLS",
                                                                                "RR     " = "RR",
                                                                                "RKHS   " = "RKHS",
                                                                                "BRNN   " = "BRNN",
                                                                                "NR     " = "NR",
                                                                                "AI     " = "AI",
                                                                                "EM     " = "EM",
                                                                                "EMMA   " = "EMMA",
                                                                                "BL     " = "BL",
                                                                                "BigRR  " = "bigRR"
                                                                                
                                                                              ))
                                                    ),
                                                    column(12,
                                                           checkboxGroupInput("model.group.C",label = h4("Machine-learning models"),inline = T,
                                                                              choices = list(
                                                                                "SVC    " = "SVC",
                                                                                "RFC    " = "RFC",
                                                                                "RFR    " = "RFR",
                                                                                "SVR    " = "SVR"
                                                                              )
                                                                              
                                                           )
                                                    ),
                                                    column(12,
                                                           column(12,
                                                                  h4("Statistics models parameters"),
                                                                  h5('Numberic parameters'),
                                                                  column(3,
                                                                         numericInput('G2P.nIter',label = "nIter",value = 1500)      
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.burnIn',label = "burnIn",value = 500)         
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.thin',label = 'thin',value = 5)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.S0',label = 'S0',value = NULL)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.df0',label = "df0",value = 5)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.R2',label = "R2",value = 0.5)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.K',label = 'K',value = 8)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.eta',label = 'eta',value = 0.7)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.eps',label = 'eps',value = 1e-4)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.alpha',label = 'alpha',value = 1)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.lambda',label = 'lambda',value = NULL)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.tol.err',label = 'tol.err',value = 1e-6)  
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.tol.conv',label = 'tol.conv',value = 1e-8)        
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.epochs',label = 'epochs',value = 30)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.neurons',label = 'neurons',value = 4)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.mmerIters',label = "mmerIters",value = 20)
                                                                  ),
                                                                  column(3,
                                                                         numericInput('G2P.maxstep',label = "maxstep",value = 100)
                                                                  )
                                                                  ),
                                                           column(12,
                                                                  h5('Logical parameters'),
                                                                  column(3,
                                                                         selectInput('G2P.outputModel',"Outpute the models",c("FALSE","TRUE"))
                                                                  ),
                                                                  column(3,
                                                                         selectInput('G2P.scale.x',"scale.x",c("FALSE","TRUE")) 
                                                                  ),
                                                                  column(3,
                                                                         selectInput('G2P.scale.y',"scale.y",c("FALSE","TRUE"))
                                                                  )),
                                                           column(12,
                                                                  h5('Character parameters'),
                                                                  column(3,
                                                                         selectInput('G2P.select',"SPLS variable selection",c("pls2","simpls"))
                                                                  ),
                                                                  column(3,
                                                                         selectInput('G2P.fit',"SPLS model fitting",c("simpls","kernelpls","kernelpls","oscorespls"))   
                                                                  ),
                                                                  column(3,
                                                                         selectInput('G2P.effectConcider',"effectConcider ",c("A","D","E","AD","AE","DE","ADE"))   
                                                                  )
                                                           )
                                                    ),
                                                    column(12,
                                                           h4('Machine-learning parameters'),
                                                           h5('Numberic parameters'),
                                                           column(3,
                                                                  numericInput('G2P.posPercentage',label = "posPercentage",value = 0.4)
                                                           ),
                                                           column(3,
                                                                  numericInput('G2P.ntree',label = "ntree",value = 500)
                                                           ),
                                                           column(3,
                                                                  numericInput('G2P.nodesize',label = "nodesize",value = 1)
                                                           ),
                                                           column(3,
                                                                  numericInput('G2P.gamma',label = "gamma",value = 1)       
                                                           ),
                                                           column(3,
                                                                  numericInput('G2P.cost',label = 'cost',value = 2^(-9))
                                                           )
                                                    ),
                                                    column(12,
                                                           h5('Character parameters'),
                                                           column(3,
                                                                  selectInput('G2P.BestIndividuals',"Except group",c("top","middle","buttom"))
                                                           ),
                                                           column(3,
                                                                  selectInput('G2P.kernel',"kernel",c("linear","polynomial","sigmoid","radial"))
                                                           )
                                                    ),
                                                    column(12,
                                                           h4('Cross validation parameters:'),
                                                           column(3,
                                                                  numericInput('G2PCV.cross',label = 'Cross',value = 10)
                                                           ),
                                                           column(3,
                                                                  numericInput('G2PCV.seed',label = 'Seed',value = 1)
                                                           ),
                                                           column(3,
                                                                  numericInput('G2PCV.cpus',label = 'CPUs',value = 1)
                                                           )
                                                           
                                                    )
                                                  )
                                                )
                                       ),
                                       tabPanel("Parameter arguments",
                                                HTML('The details of parameters:
                                                     <ul><li><b>nIter, burnIn, thin:</b>(integer)The number of iterations, burn-in and thinning,default nIter 7000,burnIn 500,thin 5.</li>
                                                     <li><b>S0, df0:</b>(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. In the parameterization of the scaled-inverse chi square in BGLR the expected values is S0/(df0-2). The default value for the df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*R2 (see below). For further details see the vignettes in the package or <a href="http://genomics.cimmyt.org/BGLR-extdoc.pdf">BGLR</a>.Default S0 NULL,df0 5.</li>
                                                     <li><b>R2:</b>(numeric, (0,1)) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode. For further details see the vignettes in the package or <a href="http://genomics.cimmyt.org/BGLR-extdoc.pdf">BGLR</a>.Defult 0.5.</li>
                                                     <li><b>K:</b>Number of hidden components.</li>
                                                     <li><b>eta:</b>(numeric)An effective zero. Default is 1e-4.</li>
                                                     <li><b>eps:</b>(numeric)Thresholding parameter. eta should be between 0 and 1.</li>
                                                     <li><b>maxstep:</b>(numeric)Maximum number of iterations when fitting direction vectors. Default is 100.</li>
                                                     <li><b>alpha:</b>The elasticnet mixing parameter.Detail in glmnet.</li>
                                                     <li><b>lambda:</b>The shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters.</li>
                                                     <li><b>tol.err:</b>Internal tolerance level for extremely small values; default value is 1e-6.</li>
                                                     <li><b>tol.conv:</b>Tolerance level in convergence; default value is 1e-8.</li>
                                                     <li><b>epochs:</b>(integer)Maximum number of epochs(iterations) to train, default 30.</li>
                                                     <li><b>neurons:</b>(integer)Indicates the number of neurons,defult 4.</li>
                                                     <li><b>mmerIters:</b>(numeric)A scalar value indicating how many iterations have to be performed if the optimization algorithms. There is no rule of tumb for the number of iterations but less than 8 is usually enough, default 20.</li>
                                                     <li><b>scale.x:</b>Scale predictors by dividing each predictor variable by its sample standard deviation?</li>
                                                     <li><b>scale.y:</b>Scale responses by dividing each response variable by its sample standard deviation?</li>
                                                     <li><b>SPLS variable selection:</b>PLS algorithm for variable selection. Alternatives are "pls2" or "simpls". Default is "pls2".</li>
                                                     <li><b>SPLS model fitting:</b>PLS algorithm for model fitting. Alternatives are "kernelpls", "widekernelpls", "simpls", or "oscorespls". Default is "simpls".</li>
                                                     <li><b>effectConcider:</b>(string)if Z = NULL, random effects are auto generated.</li>
                                                     <li><b>posPercentage:</b>(numeric)the percentage positive samples in all samples.1 > posPercentage > 0.</li>
                                                     <li><b>ntree:</b>RandomForest parameter:Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.Defualt 500.</li>
                                                     <li><b>nodesize:</b>Randomforest parameter Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).</li>
                                                     <li><b>kernel:</b>svm parameter the kernel used in training and predicting. You might consider changing some of the following parameters, depending on the kernel type.(linear,polynomial,sigmoid,radial)Default "linear".</li>                                                     
                                                     <li><b>gamma:</b>svm parameter parameter needed for all kernels except linear (default: 1/(data dimension)).</li>
                                                     <li><b>cost:</b>svm cost,default 2^(-9).</li>
                                                     <li><b>Cross:</b>(numeric)the fold number of cross validation.</li>
                                                     <li><b>Seed:</b>(numeric)random number options,defult 1.</li>
                                                     <li><b>CPUs:</b>(numeric)number of cpu cores to be used for calculations.</li>
                                                     <li><b>Outpute the models:</b>if true return the list of training model.</li>
                                                     </ul>')
                                                )
                                                ),
                           footer = tagList(
                             modalButton("Cancel"),
                             actionButton("okG2P", "OK")
                           )
                           # modalButton("OK"),
                           
                           #           if (failed)
                           #             div(tags$b("Invalid name of data object", style = "color: red;")),
                           
                                                )
             }
             
             observeEvent(input$inputG2P, {
               showModal(G2PModal())
             })
             
             observeEvent(input$okG2P, {
               removeModal()
             })
             
             output$G2PparaSelect <- renderPrint({
               HTML('<h5><b>The setting of parameters:</h5></b>',
                    '<li><b>Models:</b>',c(input$model.group.R,input$model.group.C),'</li>')
             })
             
             output$G2PparaSelect1 <- renderPrint({
               HTML('<li><b>nIter:</b>',input$G2P.nIter,'</li><li><b>burnIn:</b>',input$G2P.burnIn,'</li>',
                    '<li><b>thin:</b>',input$G2P.thin,'</li><li><b>S0:</b>',input$G2P.S0,'</li>',
                    '<li><b>df0:</b>',input$G2P.df0,'</li><li><b>R2:</b>',input$G2P.R2,'</li>',
                    '<li><b>K:</b>',input$G2P.K,'</li><li><b>eta:</b>',input$G2P.eta,'</li>',
                    '<li><b>eps:</b>',input$G2P.eps,'</li><li><b>alpha:</b>',input$G2P.alpha,'</li>',
                    '<li><b>lambda:</b>',input$G2P.lambda,'</li><li><b>tol.err:</b>',input$G2P.tol.err,'</li>',
                    '<li><b>tol.conv:</b>',input$G2P.tol.conv,'</li><li><b>epochs:</b>',input$G2P.epochs,'</li>',
                    '<li><b>neurons:</b>',input$G2P.neurons,'</li><li><b>mmerIters:</b>',input$G2P.mmerIters,'</li>',
                    '<li><b>maxstep:</b>',input$G2P.maxstep,'</li>'
                    )
                     
             })
             
             output$G2PparaSelect2<- renderPrint({
               HTML('<li><b>Outpute the models:</b>',input$G2P.outputModel,'</li><li><b>scale.x:</b>',input$G2P.scale.x,'</li>',
                    '<li><b>scale.y:</b>',input$G2P.scale.y,'</li><li><b>SPLS variable selection:</b>',input$G2P.select,'</li>',
                    '<li><b>effectConcider:</b>',input$G2P.effectConcider,'</li><li><b>posPercentage:</b>',input$G2P.posPercentage,'</li>',
                    '<li><b>ntree:</b>',input$G2P.ntree,'</li><li><b>nodesize:</b>',input$G2P.nodesize,'</li>',
                    '<li><b>gamma:</b>',input$G2P.gamma,'</li><li><b>cost:</b>',input$G2P.cost,'</li>',
                    '<li><b>Except group:</b>',input$G2P.BestIndividuals,'</li><li><b>kernel:</b>',input$G2P.kernel,'</li>',
                    '<li><b>Cross:</b>',input$G2PCV.cross,'</li><li><b>Seed:</b>',input$G2PCV.seed,'</li>',
                    '<li><b>CPUs:</b>',input$G2PCV.cpus,'</li><li><b>Inpute the marker data:</b>',input$if.inpute,'</li>')
             })

## Evaluation 
             evalModal <- function() {
               modalDialog(easyClose = TRUE,size = "l",
                           tabsetPanel(type = "tabs",
                                       tabPanel("Parameters setting",
                                                wellPanel(
                                                  fluidRow(
                                                    column(12,
                                                           h4('Global measures'),
                                                           checkboxGroupInput("measures.group.global",label = NULL,inline = T,
                                                                              choices = list(
                                                                                "PCC" = "pearson",
                                                                                "Kendall" = "kendall",
                                                                                "spearman" = "spearman",
                                                                                "MSE" = "MSE",
                                                                                "R2" = "R2")
                                                                              
                                                           )
                                                    ),
                                                    column(12,
                                                           h4('Threshold-based  measures'),
                                                           checkboxGroupInput("measures.group.threshold",label = NULL,inline = T,
                                                                              choices = list(
                                                                                "RE" = "RE",
                                                                                "Kappa" = "Kappa",
                                                                                "AUC" = "AUC",
                                                                                "AUCpr"= "AUCpr",
                                                                                "NDCG" = "NDCG",
                                                                                "MeanNDCG" = "meanNDCG",
                                                                                "F-score" = "F1",
                                                                                "Accuracy" = "accuracy")
                                                                              
                                                           )
                                                    ),
                                                    column(12,
                                                           h5('Evaluation parameters'),
                                                           column(3,
                                                                  sliderInput("eval.topAlpha", "Threshold: top alpha(%) ", min = 1, max = 100,value = c(1,90))
                                                           ),
                                                           column(3,
                                                                  numericInput('eval.Beta',label = 'Bata',value = 1)
                                                           ),
                                                           column(3,
                                                                  selectInput('eval.BestIndividuals',"BestIndividuals",c("top","middle","buttom"))
                                                           ),
                                                           column(3,
                                                                  selectInput('eval.Probability',"Probability",c("FALSE","TRUE"))
                                                           ),
                                                           column(3,
                                                                  selectInput('eval.allNDCG',"AllNDCG",c("FALSE","TRUE"))
                                                           ),
                                                           column(3,
                                                                  numericInput('eval.probIndex',label = "probIndex",value = 10 )
                                                           )
                                                    )
                                                  )
                                                )
                                       ),
                                       tabPanel("Parameter arguments",
                                                HTML('The details of parameters:
                                                     <ul><li><b>Threshold: top alpha(%):</b>(numeric,(0,100])A vector is the proportion of excellent individuals,default 1:90.</li>
                                                     <li><b>Bata:</b>(numeric)The parameter of "F1".</li>
                                                     <li><b>BestIndividuals:</b>(character)The position of expected phenotype in whole phenotypic data set."top","buttom" or "middle",default "top".</li>
                                                     <li><b>Probability:</b>(logical)Whether the predScores is probability? Default FALSE.</li>
                                                     <li><b>AllNDCG:</b>(logical)Compute all NDCG results one by one?default FALSE.</li>
                                                     <li><b>probIndex:</b>(integer)indicates the column index which prediction result is probability.</li>
                                                     </ul>')
                                                )
                                                ),
                           footer = tagList(
                             modalButton("Cancel"),
                             actionButton("okEval", "OK")
                           )
                           # modalButton("OK"),
                           
                           #           if (failed)
                           #             div(tags$b("Invalid name of data object", style = "color: red;")),
                           
                                                )
             }
             
             observeEvent(input$inputEval, {
               showModal(evalModal())
             })
             
             observeEvent(input$okEval, {
               removeModal()
             })
             
             output$EvalparaSelect <- renderPrint({
               HTML('<h5><b>The setting of parameters:</h5></b>',
                    '<li><b>Models:</b>',c(input$measures.group.global,input$measures.group.threshold),'</li>',
                    '<li><b>Threshold: top alpha(%):</b>',input$eval.topAlpha,'</li>',
                    '<li><b>Bata:</b>',input$eval.Beta,'</li>',
                    '<li><b>BestIndividuals:</b>',input$eval.BestIndividuals,'</li>',
                    '<li><b>Probability:</b>',input$eval.Probability,'</li>',
                    '<li><b>AllNDCG:</b>',input$eval.allNDCG,'</li>',
                    '<li><b>probIndex:</b>',input$eval.probIndex,'</li>'
#                   '</ul>'
                )
             })
             
             
             ################################
             # Output manual
             output$tutorial<-renderUI({
               withMathJax(includeMarkdown("tutorial/tutorial.md"))
             })
             
             
             # options(shiny.maxRequestSize = data.size.max*1024^2)
             
             
             data.exp0 <- reactive({
               infile.data <- input$markersData
               if (is.null(infile.data))
                 return(NULL)
               
               data.exp0 <- read.csv(file=infile.data$datapath,header=T)
               return(data.exp0)
               
             })
             
             ## markers data
             D.markers <- reactive({
               if (is.null(data.exp0))
                 return(NULL)
               
               D.markers <- data.exp0()
               
               D.markers <- as.matrix(D.markers)
               rownames(D.markers) <- D.markers[,1]
               D.markers <- D.markers[,-1]
               class(D.markers) <- "numeric"
               
               return(D.markers)
               
             })
             
             
             ## phenotype data
             D.phenotype <-reactive({
               infile.target<-input$phenotypeData
               if (is.null(infile.target))
                 return(NULL)
               
               D.phenotype<-read.csv(file=infile.target$datapath,header=T)
               D.phenotype<-D.phenotype[,2]
               return(D.phenotype)
             })
             
             ## extra matrix 
             D.Z<-reactive({
               infile.file.Z<-input$file.Z
               
               if (is.null(infile.target))
                 return(NULL)
               
               D.Z<-read.csv(file=infile.target$datapath,header=T)
               D.Z<- as.matrix(D.Z)[,-1]
               
               return(D.Z)
             })
             
             ## PCA analysis
             PCA.ana  <- reactive({
               PCA_Marker <- prcomp(D.markers())
               ############ kmean for cluster
               PCA_Marker <- PCA_Marker$x[,1:3]
               PCA_Marker <- as.data.frame(PCA_Marker)
               PCA_Marker <- cbind(PCA_Marker,D.phenotype())
               PCA.ana <- PCA_Marker
             })
             #######################################################################################################################             
             output$summary <- renderPrint({
               dataset <- D.phenotype()
               summary(dataset)
             })
             
             output$markerView <- renderTable({
               if (is.null(data.exp0())) {
                 
               }else{
                 x <- D.markers()[,(1:input$QCres.obsC)]
                 head(x,n = input$QCres.obsR)
               }
             },rownames = T
             )
             
             output$density <- renderPlotly({
               phenotype = D.phenotype()
               if(is.null(phenotype)){
                 return(h4("No phenotype data input! Please inpute phenotype data!"))
               }else{
               phenotype = phenotype
               barCol <- color()[input$barPlotCol]
               lineCol <- color()[input$barPlotLineCol]
               color_low <- color()[input$D3.color.low]
               color_high <- color()[input$D3.color.high]
               rowDataPlot(y =  phenotype,markers = D.markers(),plot.type = input$plotInputType, make.plotly = T,
                          show.line = input$showDensityLine,bins = input$binsN,title = '',pointSize = input$D3_size,color_low = color_low,color_high = color_high,
                          barCol = barCol,lineCol = lineCol
               )
               }
             })
             
             GSDataQC.res <- eventReactive(input$inpute.GSDataQC,{
               if(is.null(D.markers()) | is.null(D.phenotype())){
                 return()
               }else{
                 impute <- input$if.inpute
                 imputeMethod <- input$inpute.method
                 round <- input$inpute.round
                 k <- input$inpute.k
                 maxp <- input$inpute.maxp
                 rowmax <- input$inpute.rowmax
                 colmax <- input$inpute.colmax
                 rng.seed <- input$inpute.rng.seed
                 GSDataQC.res <- GSDataQC(markers = D.markers(), phenotype = D.phenotype(), impute = impute,
                                          imputeMethod = imputeMethod, round = round ,k = k,
                                          maxp = maxp, rowmax = rowmax, colmax = colmax, rng.seed = rng.seed)
                 
               }
               
             })
             
             
             output$Global.summary <- renderPrint({
               if(is.null(GSDataQC.res()[["Imputation"]])){
                 return()}else{
                   Globalres <- t(as.matrix(GSDataQC.res()$Global))
                   rownames(Globalres) <- "Statistic"
                   Globalres
                 }
             })
             
             output$GSDataQC.table <- renderTable({
               QCRes <- GSDataQC.res()
               statisticLabel <- input$QCres.output
               if(statisticLabel == "Imputation matrix"){
                 statisticRes <- QCRes[["Imputation"]][1:input$QCres.obsR,1:input$QCres.obsC]
               }else if(statisticLabel == "NA in column"){
                 statisticRes <- cbind(QCRes[["NACountCol"]],QCRes[["NACountColPercent"]])
                 colnames(statisticRes) <- c("Count","Percent")
                 statisticRes <- statisticRes[1:input$QCres.obsR,]
               }else if(statisticLabel == "NA in row"){
                 statisticRes <- cbind(QCRes[["NACountRow"]],QCRes[["NACountRowPercent"]])
                 colnames(statisticRes) <- c("Count","Percent")
                 statisticRes <- statisticRes[1:input$QCres.obsR,]
               }else if(statisticLabel == "MAF"){
                 statisticRes <- as.matrix(QCRes[["MAF"]])
                 colnames(statisticRes) <- c("Percent")
                 statisticRes <- statisticRes[1:input$QCres.obsR,]
               }else if(statisticLabel == "NA index"){
                 statisticRes <- as.matrix(QCRes[["NAIdx"]])
                 colnames(statisticRes) <- "Index"
                 statisticRes <- statisticRes[1:input$QCres.obsR,]
               }
               statisticRes
             },rownames = T
             )
             
             output$download.QC.inputation <- downloadHandler(
               filename=function(){
                 'inputation.csv'
               },
               content=function(file){
                 statisticRes <- GSDataQC.res()[["Imputation"]]
                 write.csv(statisticRes,file,row.names = T)
               }
             )
             
             output$download.QC.NACol <- downloadHandler(
               filename=function(){
                 'NAinCol.csv'
               },
               content=function(file){
                 statisticRes <- cbind(GSDataQC.res()[["NACountCol"]],GSDataQC.res()[["NACountColPercent"]])
                 colnames(statisticRes) <- c("Count","Percent")
                 write.csv(statisticRes,file,row.names = T)
               }
             )
             
             output$download.QC.NARow <- downloadHandler(
               filename=function(){
                 'NAinRow.csv'
               },
               content=function(file){
                 statisticRes <- cbind(GSDataQC.res()[["NACountRow"]],GSDataQC.res()[["NACountRowPercent"]])
                 colnames(statisticRes) <- c("Count","Percent")
                 write.csv(statisticRes,file,row.names = T)
               }
             )
             
             output$download.QC.MAF <- downloadHandler(
               filename=function(){
                 'MAF.csv'
               },
               content=function(file){
                 statisticRes <- as.matrix(GSDataQC.res()[["MAF"]])
                 colnames(statisticRes) <- c("Percent")
                 write.csv(statisticRes,file,row.names = T)
               }
             )
             
             output$download.QC.NAIndex <- downloadHandler(
               filename=function(){
                 'NAIdex.csv'
               },
               content=function(file){
                 statisticRes <- as.matrix(GSDataQC.res()[["NAIdx"]])
                 colnames(statisticRes) <- "Index"
                 write.csv(statisticRes,file,row.names = T)
               }
             )
             
             output$download.input.plot <- downloadHandler(
               
               filename = function() {
                 paste0('Density plot.',input$save.input.type)
               },
               content = function(file,format=input$save.input.type) {
                 if(format=='html')
                   suppressWarnings(htmlwidgets::saveWidget(as.widget(
                     rowDataPlot(y =  D.phenotype(),plot.type = input$plotInputType, make.plotly = T,
                                show.line = input$showDensityLine,bins = input$binsN,title = '')
                   ), file=file,selfcontained=T))
                 else ggsave(file,
                             rowDataPlot(y =  D.phenotype(),plot.type = input$plotInputType, make.plotly = F,
                                        show.line = input$showDensityLine,bins = input$binsN,title = ''),
                             width = 16,height = 12,units = "cm")
               })
             
             #### G2P 
             G2P.res <- eventReactive(input$G2P.run,{
#                if(is.null(D.markers()) | is.null(D.phenotype())){
#                  return()
#                }else{
                 # output$G2PCVprocess <- renderPrint(HTML('<h5><b>G2P is running ......</h5></b>'))
                 trainMarker <- D.markers()[input$left.T:input$right.T,]
                 trainPheno <- D.phenotype()[input$left.T:input$right.T]
                 testMarker <- D.markers()[input$left.V:input$right.V,]
                 testPheno <- D.phenotype()[input$left.V:input$right.V]
                 modelMethods <- c(input$model.group.R,input$model.group.C)
                 if(is.na(input$G2P.lambda)){lambda <- NULL}else{lambda <- input$G2P.lambda}
                 if(is.na(input$G2P.S0)){S0 <- NULL}else{S0 <- input$G2P.S0}
                 
                 G2P.res <- G2P(trainMarker = trainMarker,trainPheno = trainPheno,testMarker = testMarker, testPheno = testPheno,
                                modelMethods = modelMethods,outputModel = input$G2P.outputModel,  # main parameters
                                nIter = input$G2P.nIter, burnIn = input$G2P.burnIn, thin = input$G2P.thin, 
                                saveAt = "", S0 = S0, df0 = input$G2P.df0, R2 = input$G2P.R2, weights = NULL,
                                verbose = FALSE, rmExistingFiles = TRUE, groups=NULL,importance = FALSE,    # # # BGLR method parameters
                                posPercentage = input$G2P.posPercentage,BestIndividuals = input$G2P.BestIndividuals,ntree = input$G2P.ntree,
                                nodesize = input$G2P.nodesize,kernel = input$G2P.kernel,gamma = input$G2P.gamma, cost = input$G2P.cost,  # machine learing parameters
                                K = input$G2P.K,eta = input$G2P.eta,select = input$G2P.select,fit = input$G2P.fit,scale.x = input$G2P.scale.x,scale.y = input$G2P.scale.y,eps = input$G2P.eps,trace = FALSE, maxstep = input$maxstep,# SPLS parameters
                                alpha = input$G2P.alpha,X = NULL,family = gaussian(link = identity), lambda = lambda, tol.err = input$G2P.tol.err, tol.conv =input$G2P.tol.conv,
                                epochs = input$G2P.epochs, neurons =input$G2P.neurons, Z = NULL,  effectConcider = input$G2P.effectConcider, mmerIters = input$G2P.mmerIters)
                 if(is.null(rownames(G2P.res))){
                   rownames(G2P.res) <- 1:nrow(testMarker)
                 }else{
                   rownames(G2P.res) <- rownames(testMarker)
                 }
                 # output$G2PCVprocess <- renderPrint(HTML('<h5><b>G2P is Done !!!</h5></b>'))
                 G2P.res
#                }
             })
             
             ### G2P CV 
             G2PCV.res <- eventReactive(input$G2PCV.run,{
#                if(is.null(D.markers()) | is.null(D.phenotype())){
#                  return()
#                }else{
                 # output$G2PCVprocess <- renderPrint(HTML('<h5><b>G2P Cross Validation is running ......</h5></b>'))
                 markers <- D.markers()
                 pheVal <- D.phenotype()
                 cross <- input$G2PCV.cross
                 seed <- input$G2PCV.seed
                 cpus <- input$G2PCV.cpus
                 modelMethods <- c(input$model.group.R,input$model.group.C)
                 if(is.na(input$G2P.lambda)){lambda <- NULL}else{lambda <- input$G2P.lambda}
                 if(is.na(input$G2P.S0)){S0 <- NULL}else{S0 <- input$G2P.S0}
                 
                 G2PCV.res <- G2PCrossValidation(cross = cross,seed = seed,cpus = cpus,markers = markers, pheVal = pheVal,
                                                 modelMethods = modelMethods,outputModel = input$G2P.outputModel,  # main parameters
                                                 nIter = input$G2P.nIter, burnIn = input$G2P.burnIn, thin = input$G2P.thin, 
                                                 saveAt = "", S0 = S0, df0 = input$G2P.df0, R2 = input$G2P.R2, weights = NULL,
                                                 verbose = FALSE, rmExistingFiles = TRUE, groups = NULL,importance = FALSE,    # # # BGLR method parameters
                                                 posPercentage = input$G2P.posPercentage,BestIndividuals = input$G2P.BestIndividuals,ntree = input$G2P.ntree,
                                                 nodesize = input$G2P.nodesize,kernel = input$G2P.kernel,gamma = input$G2P.gamma, cost = input$G2P.cost,  # machine learing parameters
                                                 K = input$G2P.K,eta = input$G2P.eta,select = input$G2P.select,fit = input$G2P.fit,scale.x = input$G2P.scale.x,scale.y = input$G2P.scale.y,eps = input$G2P.eps,trace = FALSE,maxstep = input$maxstep, # SPLS parameters
                                                 alpha = input$G2P.alpha,X = NULL,family = gaussian(link = identity), lambda = lambda, tol.err = input$G2P.tol.err, tol.conv =input$G2P.tol.conv,
                                                 epochs = input$G2P.epochs, neurons =input$G2P.neurons, Z = NULL,  effectConcider = input$G2P.effectConcider, mmerIters = input$G2P.mmerIters)
                 G2PCV.res <- resultsMerge(G2PCV.res)
                 output$G2PCVprocess <- renderPrint(HTML('<h5><b> G2P Cross Validation have done !!!</h5></b>'))
                 G2PCV.res
               # }
             })
             
#              G2PCV.signal <- reactive({
#                cv.colnum <- ncol(G2PCV.res()) - 1
#                cv.cross  <- input$G2PCV.cross
#                output$G2PCVprocess <- renderPrint(HTML('<h5><b>',cv.cross,'fold G2P Cross Validation with',cv.colnum,'Method(s) have done!!!</b></h5>'))
#              })
             
             output$G2PCVprocess <- renderPrint({
               cv.colnum <- ncol(G2PCV.res()) - 1
               cv.cross  <- input$G2PCV.cross
               HTML('<h5><b>',cv.cross,'fold G2P Cross Validation with',cv.colnum,'Method(s) have done!!!</b></h5>')
             })
             
             ##### display 
             output$G2P.res.table <- renderTable({
               if(input$res.select == "G2P results"){
                 G2P.res()[1:input$G2P.obsR,1:input$G2P.obsC]
               }else if(input$res.select == "G2P CV results"){
                 G2PCV.res()[1:input$G2P.obsR,1:input$G2P.obsC]
               }
             },rownames = T
             )
             
             eval.res.G2P <- eventReactive(input$eval.run.G2P,{
               if(is.null(G2P.res())){
                 return()
               }else{
                 realScores <- G2P.res()[,1]
                 predScores = G2P.res()[,2:(ncol(G2P.res()))]
                 evalMethod <- c(input$measures.group.global,input$measures.group.threshold)
                 BestIndividuals <- input$eval.BestIndividuals
                 Probability <- input$eval.Probability
                 probIndex <- input$eval.probIndex
                 topAlpha <- input$eval.topAlpha[1]:input$eval.topAlpha[2]
                 Beta <- input$eval.Beta  
                 allNDCG <- input$eval.allNDCG
                 eval.res.G2P <- evaluateGS(realScores = realScores, predScores = predScores, 
                                            evalMethod = evalMethod,topAlpha = topAlpha,
                                            BestIndividuals = BestIndividuals,Probability = Probability,probIndex = probIndex,Beta = Beta,allNDCG = allNDCG )
                 eval.res.G2P
               }
             })
             
             eval.res.CV <- eventReactive(input$eval.run.CV,{
               if(is.null(G2PCV.res())){
                 return()
               }else{
                 realScores <- G2PCV.res()[,1]
                 predScores = G2PCV.res()[,2:(ncol(G2PCV.res()))]
                 evalMethod <- c(input$measures.group.global,input$measures.group.threshold)
                 BestIndividuals <- input$eval.BestIndividuals
                 Probability <- input$eval.Probability
                 probIndex <- input$eval.probIndex
                 topAlpha <- input$eval.topAlpha[1]:input$eval.topAlpha[2]
                 Beta <- input$eval.Beta  
                 allNDCG <- input$eval.allNDCG
                 eval.res.CV <- evaluateGS(realScores = realScores, predScores = predScores, 
                                           evalMethod = evalMethod,topAlpha = topAlpha,
                                           BestIndividuals = BestIndividuals,Probability = Probability,probIndex = probIndex,Beta = Beta,allNDCG = allNDCG )
                 eval.res.CV
               }
             })
             
             output$eval.res.global.table <- renderTable({
               #                  if(is.null(eval.res.G2P()[["corMethosds"]]) & is.null(eval.res.CV()[["corMethosds"]])){
               #                    return(NULL)
               #                  }else{
               if(input$res.select == "G2P results"){
                 eval.res.G2P()[["corMethosds"]]
               }else if(input$res.select == "G2P CV results"){
                 eval.res.CV()[["corMethosds"]]
               }
               # }
             },rownames = T
             )
             
             output$eval.res.table <- renderTable({
               if(input$res.select == "G2P results"){
                 evalRes <- eval.res.G2P()
               }else if(input$res.select == "G2P CV results"){
                 evalRes <- eval.res.CV()
               }
               evalMethodLabel <- input$eval.res.select
               if(evalMethodLabel == "RE"){
                 evalMethosRes <- evalRes[["RE"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "Kappa"){
                 evalMethosRes <- evalRes[["Kappa"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "AUC"){
                 evalMethosRes <- evalRes[["AUC"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "AUCpr"){
                 evalMethosRes <- evalRes[["AUCpr"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "NDCG"){
                 evalMethosRes <- evalRes[["NDCG"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "meanNDCG"){
                 evalMethosRes <- evalRes[["meanNDCG"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "F1"){
                 evalMethosRes <- evalRes[["F1"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }else if(evalMethodLabel == "accuracy"){
                 evalMethosRes <- evalRes[["accuracy"]][1:input$evalres.obsR,1:input$evalres.obsC]
               }
               evalMethosRes
             },rownames = T
             )
             
             output$download.predres.mat <- downloadHandler(
               filename=function(){
                 'predresMat.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   predres <- G2P.res()
                 }else if(input$res.select == "G2P CV results"){
                   predres <- G2PCV.res()
                 }
                 write.csv(predres,file,row.names = T)
               }
             )
             
             output$download.evalres.global <- downloadHandler(
               filename=function(){
                 'evalresGlobal.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["corMethosds"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["corMethosds"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.RE <- downloadHandler(
               filename=function(){
                 'evalresRE.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["RE"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["RE"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.Kappa <- downloadHandler(
               filename=function(){
                 'evalresKappa.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["Kappa"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["Kappa"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.AUC <- downloadHandler(
               filename=function(){
                 'evalresAUC.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["AUC"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["AUC"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.AUCpr <- downloadHandler(
               filename=function(){
                 'evalresAUCpr.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["AUCpr"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["AUCpr"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.NDCG <- downloadHandler(
               filename=function(){
                 'evalresNDCG.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["NDCG"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["NDCG"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.meanNDCG <- downloadHandler(
               filename=function(){
                 'evalresMeanNDCG.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["meanNDCG"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["meanNDCG"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.F1 <- downloadHandler(
               filename=function(){
                 'evalresF1.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["F1"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["F1"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             
             output$download.evalres.accuracy <- downloadHandler(
               filename=function(){
                 'evalresAcc.csv'
               },
               content=function(file){
                 if(input$res.select == "G2P results"){
                   evalRes <- eval.res.G2P()[["accuracy"]]
                 }else if(input$res.select == "G2P CV results"){
                   evalRes <- eval.res.CV()[["accuracy"]]
                 }
                 write.csv(evalRes,file,row.names = T)
               }
             )
             ################################################ visualiztion ##############################################
             output$visualization.scatter <- renderPlotly({
               if(is.null(G2PCV.res()) | input$method1 == " " | input$method2 == " "){
                 return(NULL)
               }else{
                 if(input$scatter.col.mid == 0){
                   col.mid = NULL}else{
                     col.mid = color()[input$scatter.col.mid]
                   }
                 if(input$scatter.col.low == 0){col.low = NULL}else{col.low = color()[input$scatter.col.low]}
                 if(input$scatter.col.high == 0){col.high = NULL}else{col.high = color()[input$scatter.col.high]}
                 scatterPlot(predmat = G2PCV.res(),x1 = input$method1,x2= input$method2,col.low = col.low,
                             col.high = col.high,col.mid = col.mid,col = input$scatter.color,show.line = input$scatter.showline,color_szie = input$scatterColor_size,
                             alpha = input$scatter.alpha,sizeRange = input$scatter.sizeRange,make.plotly = T
                 )
               }
               
             })
             
             output$visualization.line <- renderPlotly({
               if(is.null(eval.res.CV()[[input$thresholdLine.measure]])){
                 return(NULL)
               }else{
                 evalMat <- eval.res.CV()[[input$thresholdLine.measure]]
               }
               linePlot(evalMat = evalMat,xlab = input$thresholdLine.xlab,ylab = input$thresholdLine.ylab,
                        legendTittle = input$thresholdLine.legend.title,size = input$thresholdLine.width,
                        alpha = input$thresholdLine.alpha,
                        make.plotly = T
               )
             })
             
             output$visualization.bar <- renderPlotly({
               if(is.null(eval.res.CV()[["corMethosds"]])){
                 return(NULL)
               }else{
                 data <- eval.res.CV()[["corMethosds"]]
               }
               barPlot(data = data,xlab = input$globalBar.xlab,ylab = input$globalBar.ylab,
                       legendTittle = input$globalBar.legend.title,
                       other = input$globalBar.type,
                       make.plotly = T
               )
             })
             
             output$visualization.heatmap <- renderPlotly({
               if(is.null(G2PCV.res())){
                 return(NULL)
               }else{
                 if(input$plotData == "Prediction results"){
                   heatmapData <- G2PCV.res()
                   
                   if(input$predHeatmap.col.mid == 0){
                     col.mid = NULL}else{
                       col.mid = color()[input$predHeatmap.col.mid]
                     }
                   if(input$predHeatmap.col.low == 0){col.low = NULL}else{col.low = color()[input$predHeatmap.col.low]}
                   if(input$predHeatmap.col.high == 0){col.high = NULL}else{col.high = color()[input$predHeatmap.col.high]}
                   heatmapPlot(predmat = heatmapData,col.low = col.low,col.high = col.high,col.mid = col.mid,boundSize = input$predHeatmap.boundSize,
                               x_cluster = input$predHeatmap.x_cluster,y_cluster = input$predHeatmap.y_cluster,
                               xlab = input$predHeatmap.xlab,ylab = input$predHeatmap.ylab,legend.title = input$predHeatmap.legend.title,
                               make.plotly = T
                   )
                   
                 }else if(input$plotData == "Threshold evaluation results"){
                   heatmapData <- heatMapDataProcess(x = eval.res.CV(),highBound = input$thresholdHeatmap.highBound,
                                                     lowBound = input$thresholdHeatmap.lowBound,alpha = input$thresholdHeatmap.alpha,
                                                     basedMethod = input$thresholdHeatmap.basedMethod)
                   if(input$thresholdHeatmap.col.mid == 0){
                     col.mid = NULL}else{
                       col.mid = color()[input$thresholdHeatmap.col.mid]
                     }
                   if(input$thresholdHeatmap.col.low == 0){col.low = NULL}else{col.low = color()[input$thresholdHeatmap.col.low]}
                   if(input$thresholdHeatmap.col.high == 0){col.high = NULL}else{col.high = color()[input$thresholdHeatmap.col.high]}
                   heatmapPlot(predmat = heatmapData,col.low = col.low,col.high = col.high,col.mid = col.mid,boundSize = input$thresholdHeatmap.boundSize,
                               x_cluster = input$thresholdHeatmap.x_cluster,y_cluster = input$thresholdHeatmap.y_cluster,
                               xlab = input$thresholdHeatmap.xlab,ylab = input$thresholdHeatmap.ylab,legend.title = input$thresholdHeatmap.legend.title,
                               make.plotly = T
                   )
                 }
               }
             })
             
#              output$download.input.plot <- downloadHandler(
#                
#                filename = function() {
#                  paste0('Density plot.',input$save.input.type)
#                },
#                content = function(file,format=input$save.input.type) {
#                  if(format=='html')
#                    suppressWarnings(htmlwidgets::saveWidget(as.widget(
#                      rowDataPlot(y =  D.phenotype(),plot.type = input$plotInputType, make.plotly = T,
#                                 show.line = input$showDensityLine,bins = input$binsN,title = '')
#                    ), file=file,selfcontained=T))
#                  else ggsave(file,
#                              rowDataPlot(y =  D.phenotype(),plot.type = input$plotInputType, make.plotly = F,
#                                         show.line = input$showDensityLine,bins = input$binsN,title = ''),
#                              width = 16,height = 12,units = "cm")
#                })
             
           }
  )
}

