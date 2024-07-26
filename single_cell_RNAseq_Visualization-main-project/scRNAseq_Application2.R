#Unused packages:
#library(org.Hs.eg.db)
#library(org.Mm.eg.db)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(semantic.dashboard)
#library(RColorBrewer)



#Load packages:
library(shiny)
library(shinydashboard)
#library(shinydashboardPlus)
library(shinyjqui)
library(DT)
library(BiocManager)
library(pheatmap)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(Seurat)
library(DropletUtils)
library(dplyr)
library(scater)
library(cowplot)
library(biomaRt)
library(patchwork)
library(SCpubr)

Logo_name = "Bridge Informatics"
app_version = "0.25"
Application_name = "Single Cell RNAseq Data Analysis Tool"



###########################################################
ui <- dashboardPage(
  #This is setting the title of the application
  dashboardHeader(title = paste0(Logo_name),

                  #tags$img(src='BILogo.jpeg')), 
                  dropdownMenu(type="notifications",
                               taskItem("Project progress...", 10.00, color = "red")),
                  dropdownMenu(type = "tasks", badgeStatus = "success",
                               taskItem(value = 90, color = "green",
                                        "Documentation"
                               ),
                               taskItem(value = 17, color = "aqua",
                                        "Project X"
                               ),
                               taskItem(value = 75, color = "yellow",
                                        "Server deployment"
                               ),
                               taskItem(value = 80, color = "red",
                                        "Overall project")
                  ),
                  dropdownMenu(type = "notifications",
                               notificationItem(
                                 text = "5 new users today",
                                 icon("users")
                               ),
                               notificationItem(
                                 text = "12 items delivered",
                                 icon("truck"),
                                 status = "success"
                               ),
                               notificationItem(
                                 text = "Server load at 86%",
                                 icon = icon("exclamation-triangle"),
                                 status = "warning"
                               )
                  )
                  
  ), #End of dashboard header
  
  dashboardSidebar(width=NULL,
                   
             #Tags for centering text       
             tags$style(".skin-blue .sidebar
                .center {
                        text-align: center;
               }"),
            
             #Tags for borderbox         
             tags$style(".skin-blue .sidebar
                .borderbox {
                        border: 2px solid #666;
                        padding: 5px 5px 5px 5px;
                        margin: 2px;
               }"),
            
              #Tags for character color:
              ##Red
              tags$style(
                  ".first-p {
                       color: white;
                          }
                    #element {
                        color: white;
                            }"),
              ##Green
              tags$style(
                ".second-p {
                       color: Green;} #element {
                        color: Green;
                            }"),
              
             #Panel for data input
              conditionalPanel(condition = "input.program == 'Tab1'",
                     div(id ="samplesheetPage",
                         tags$br(),
                         tags$br(),  
                         tags$br(),
                         tags$br(),
                    tags$div('class'="borderbox",
                    tags$div('class'="center",
                    h5(HTML('<b>Load Data</b><br>'))
                            ),

                     sidebarMenu(
                     menuItem("load data (.txt)",icon=icon("th"),tabName=p(class="first-p","Load Data (Raw)"),
                              startExpanded = FALSE,
                              
                              fileInput(inputId  = "Text data",
                                        label    = "Choose .txt File",
                                        multiple = FALSE,
                                        accept   = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                              #textOutput('directory_name')
                              checkboxInput("header", "Header", TRUE),
                              #sets the sep parameter when uploading files into R
                              radioButtons("sep", "Separator",
                                           choices  = c(comma = ",",Semicolon = ";",Tab = "\t"),selected = ",")
                              ),
                     tags$br(),
                     tags$br(), 
                     menuItem("Select Folder", icon = icon("bell"), 
                              tags$div('class' = "center",
                              tags$br(),
                              shinyDirButton("dir", "Select Directory", "Upload"), 
                              h5(HTML('<br>This is your default directory. \
                                      Choose a directory that contains raw files generated from 10X Genomics. The output directory will be displayed below</b><br>'))
                                ),
                              verbatimTextOutput("dir", placeholder = TRUE)
                              ),
                     tags$br(),
                     tags$br()  
                           )
                          )
                         )
                        ),
             #Conditonal panel for gene annotation plots
             conditionalPanel(condition = "input.program == 'Tab3'",
                  div(id = 'ScatterPlot',
                    conditionalPanel(condition = "input.QC == 'Tab1'",      
                     div(id='ScatterPlot',
                       tags$div('class'= "center",
                                h4(HTML('<br><b>Settings</b><br>')),
                                ),
                       
                   tags$div('class'="borderbox",
                      menuItem("Axis Labels",icon=icon("th"),tabName=p(class="first-p","Load Data (Raw)"),
                       startExpanded = TRUE,
                       tags$br(),
                       textInput("main_title", label=h4("Main Title"),width="200",value="Empty Droplets",placeholder="Empty Droplets"),
                       textInput("yaxis_title", label=h4("Y axis"),width="200",value="Total UMI count",placeholder="Total UMI count"),
                       textInput("xaxis_title", label=h4("X axis"),width="200",value="-Log Probability",placeholder="-Log Probability")
                              ),
                      tags$br(),
                      tags$br(),
                      menuItem("Axis Font",icon=icon("th"),tabName=p(class="first-p","Load Data (Raw)"),
                           startExpanded = TRUE,
                           numericInput("Title_fontsize", label = "Main Title", value = 2),
                           numericInput("axis_fontsize", label = "Axis", value = 1),
                           numericInput("Lable_font", label = "Axis Label Font", value = 1)
                              ),
                             ),
                        tags$br(),
                        tags$br()
                          )
                         ), 
                  conditionalPanel(condition="input.QC== 'Tab2'",
                            div(id = 'Mitochondrial Percentage',
                              tags$br(),
                              tags$div('class'="borderbox",
                               menuItem("Axis Labels",icon=icon("th"),tabName=p(class="first-p","Load Data (Raw)"), startExpanded=TRUE,
                                textInput("main_mito", label=h5("Main Title"),width="200",value="Mitochondria",placeholder="Mitochondria"),
                                textInput("y_mito", label=h5("Y axis"),width="200",value="Y-axis",placeholder="Y-axis"),
                                textInput("x_mito", label=h5("X axis"),width="200",value="X-axis",placeholder="X-axis")
                               ),
                              tags$br(),
                             menuItem("Font Sizes",icon=icon("th"),tabName=p(class="first-p","Load Data (Raw)"), startExpanded=TRUE,
                                      numericInput("Axis_Font", label="Axis",value=10)
                                      
                                      
                             )
                               
                            )  
                         )
                        )
                       )
                      )
                   
                              
), #End Dashboard Sidebar



  dashboardBody(
    navbarPage(paste0(Application_name,", v. ", app_version),
    #fluidPage(
        tabsetPanel(id = "program",
           tabPanel(title = "Data Upload", id = "Tab1", value = "Tab1",
              HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
                 fluidRow(column(width=6,
                   tags$h3("Sample Information"),
                   tags$hr(),
                   p(h5(HTML('<b>The current selected directory is:</b>')), h5(textOutput(outputId="selected_dir",inline=TRUE))),
                   tags$br(),
                   tags$h5(HTML('Please note that once your files are read in successfully, the items below will be filled to display \
                                number of barcoded cells detected in your scRNAseq experiment.<br>')),
                   
                   p(h5(HTML('<b>Number of barcodes detected in your sample:</b>'))), 
                   h5(textOutput(outputId="barcode_count",inline=TRUE),
                  tags$hr(),
                  h5(HTML('<br><b>Select a reference genome:</b>')),
                  tags$h5("Currently, only mouse (mm10) and human (Hg19) genomes are configured for this program!"),
                  #select input for genome
                  tags$br(),
                  selectInput("ref_genome_organism", label = "Reference organism",
                              choices = list("Human" = 1, "Mouse" = 2),selected = 1),
                  
                  #Select the controls and treatment
                  uiOutput("control"),
                  uiOutput("treatment"),
                          )
                  ),column(width=5,
                     tags$h3("Main Settings"),
                     tags$hr(),
                
                     textInput("Mito_percent",label=h5("Filter by mitochondrial DNA (%)"),width="500",value="5",placeholder="5"),
                     textInput("Ribo_percent",label=h5("Filter by ribosomal reads (%)"),width="500",value="10",placeholder="10"),
                     textInput("Minimum_features",label=h5("Filter by minimum feature count"),width="500",value="200",placeholder="2500"),
                     textInput("Maximum_features",label=h5("Filter by maximum feature count"),width="500",value="2500",placeholder="2500"),
                     actionButton("Start_Analysis", 
                                  HTML('<b>Start My Analysis!</b>'),
                                  icon("Next Page"),
                                  style = "color: #fff;background-color: #27ae60; border-color: #fff;padding: 5px 14px 5px 14px;margin: 5px 5px 5px 5px; "),
                                ),
                               ),HTML('</p></div>'),
        
                              ),
     
                  #inputPanel(fileInput(inputId = "csv-upload", label = "File upload:", accept = ".csv"))
                    
           
            tabPanel(title = "Step 1: Pre-Processing", id = "Tab3", value = "Tab3",
               HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
               tabsetPanel(id = "QC",
                   tabPanel(title = "Scatter Plot", id = "Tab1", value = "Tab1",
                        tags$h3("Scatter Plot"),
                        tags$hr(),
                        fluidRow(width=40,
                           dropdown(label = "Save Scatter Plot",
                              selectInput("ScatterPlot_dpi", label = "Output dpi", width = "95",
                                          choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                          selected = 48),
                              downloadButton("saveSPpng","PNG"),
                              HTML('<br><br>'), #<br> creates breaks and moves item to next line
                              downloadButton("saveSPjpg","JPG"),
                              HTML('<br><br>'),
                              downloadButton("saveSPtiff","TIFF"),
                              size = "default",
                              icon = icon("download", class = ""), 
                              up = FALSE),
                               HTML('</p>'),
                             jqui_resizable(#jqui resizable canvass
                              plotOutput("Scatter_plot", height = "350", width = "350"),
                                   )
                                  )
                                ),#End scatter plot
                           tabPanel(title = "Percent Mitochondria", id = "Tab2", value = "Tab2",
                              #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                            fluidRow(width=40,
                               tags$h3("Percent Mitochondria"),
                               tags$hr(),
                               dropdown(label = "Save plot",
                                  selectInput("Heatmapo_dpi", label = "Output dpi", width = "95",
                                              choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                              selected = 48),
                                  downloadButton("saveHEATpng", "PNG"),
                                  HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                  downloadButton("saveHEATjpg", "JPG"),
                                  HTML('<br><br>'),
                                  downloadButton("saveHEATtiff", "TIFF"),
                                  size = "default",
                                  icon = icon("download", class = ""), 
                                  up   = FALSE),
                               HTML('</p>'),
                               jqui_resizable(#jqui resizable canvass
                                 plotOutput("HEATo_plot", height = "350", width = "350")
                               )#,HTML('</div>')
                              )
                            ),
                           tabPanel(title = "Percent Ribosomal Reads", id = "Tab3", value = "Tab3",
                              #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                              tags$h3("Percent Ribosomal Reads"),
                              tags$hr(),
                              fluidRow(width=40,
                                 dropdown(label = "Save Plot",
                                      selectInput("MAplot_dpi", label = "Output dpi", width = "95",
                                                  choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                  selected = 48),
                                      downloadButton("saveMApng", "PNG"),
                                      HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                      downloadButton("saveMAjpg", "JPG"),
                                      HTML('<br><br>'),
                                      downloadButton("saveMApdf", "PDF"),
                                      HTML('<br><br>'),
                                      downloadButton("saveMAtiff", "TIFF"),
                                      size = "default",
                                      icon = icon("download", class = ""), 
                                          up = FALSE),
                                   HTML('</p>'),
                                   jqui_resizable(
                                     #tagList(
                                     #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                     plotOutput("MA_plot", height = "300", width = "360"),
                                     #HTML('</p></div>')
                                       )
                                     )
                                   ),
                           tabPanel(title = "Features", id = "Tab4", value = "Tab4",
                              tags$h3("Features"),
                              tags$hr(),
                              fluidRow(width=40,
                                 dropdown(label = "Save Plot",
                                      selectInput("VOLCANO_plot_dpi", label = "Output dpi", width = "95",
                                                  choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                  selected = 48),
                                      downloadButton("saveVOLCANOpng", "PNG"),
                                      HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                      downloadButton("saveVOLCANOjpg", "JPG"),
                                      HTML('<br><br>'),
                                      downloadButton("saveVOLCANOpdf", "PDF"),
                                      HTML('<br><br>'),
                                      downloadButton("saveVOLCANOtiff", "TIFF"),
                                      size = "default",
                                      icon = icon("download", class = ""), 
                                      up = FALSE),
                                   HTML('</p>'),
                                   jqui_resizable(
                                     plotOutput("VOLCANO_plot", height = "300", width = "360")
                                       )
                                     )
                                    ),
                   tabPanel(title = "Summary Plots", id = "Tab5", value = "Tab5",
                    tags$h3("Summary Plots"),
                    tags$hr(),
                    fluidRow(width=40,
                       dropdown(label = "Save Plot",
                                selectInput("FourPlots_dpi", label = "Output dpi", width = "95",
                                            choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                            selected = 48),
                                downloadButton("saveFOURpng", "PNG"),
                                HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                downloadButton("saveFOURjpg", "JPG"),
                                HTML('<br><br>'),
                                downloadButton("saveFOURtiff", "TIFF"),
                                size = "default",
                                icon = icon("download", class = ""), 
                                      up = FALSE),
                             HTML('</p>'),
                             jqui_resizable(
                               plotOutput("four_plot", height = "300", width = "360")
                                       )
                                      )
                                    )
                                  ),HTML('</p></div>')
                                 ),
                  

           tabPanel(title = "Step 2: Processing", id = "Tab4", value = "Tab4",
                    HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
                    tabsetPanel(id = "QC",
                                tabPanel(title = "Feature Selection", id = "Tab1", value = "Tab1",
                                         tags$h3("Feature Selection"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("FeaturePlot_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveFEATUREPpng","PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveFEATUREjpg","JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveFEATUREtiff","TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(#jqui resizable canvass
                                                    plotOutput("FEATURE_plot", height = "350", width = "350"),
                                                  )
                                         )
                                ),#End scatter plot
                                tabPanel(title = "Feature Violin Plot", id = "Tab2", value = "Tab2",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         fluidRow(width=40,
                                                  tags$h3("Feature Violin Plot"),
                                                  tags$hr(),
                                                  dropdown(label = "Save plot",
                                                           selectInput("FeaturePlotViolin_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveFeatureViolinpng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveFeatureViolinjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveFeatureViolintiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up   = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(#jqui resizable canvass
                                                    plotOutput("FeatureViolin_plot", height = "350", width = "350")
                                                  )#,HTML('</div>')
                                         )
                                ),
                                tabPanel(title = "PCA Heatmaps", id = "Tab3", value = "Tab3",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         tags$h3("PCA Heatmaps"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("PCAHeatmap_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("savePCAHeatmappng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("savePCAHeatmapjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("savePCAHeatmappdf", "PDF"),
                                                           HTML('<br><br>'),
                                                           downloadButton("savePCAHeatmaptiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(
                                                    #tagList(
                                                    #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                                    plotOutput("PCAHeatmap_plot", height = "300", width = "360"),
                                                    #HTML('</p></div>')
                                                  )
                                         )
                                ),
                                tabPanel(title = "PCA Plot", id = "Tab4", value = "Tab4",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         tags$h3("PCA Plot (PC1 vs PC2)"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("VizDimPlot_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveVizDimPlotpng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveVizDimPlotjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveVizDimPlotpdf", "PDF"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveVizDimPlottiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(
                                                    #tagList(
                                                    #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                                    plotOutput("JackStraw_plot", height = "300", width = "360"),
                                                    #HTML('</p></div>')
                                                  )
                                         )
                                ),
                                tabPanel(title = "Elbow Plot", id = "Tab5", value = "Tab5",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         tags$h3("Elbow Plot"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("Elbowplot_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveELBOWpng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveELBOWjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveELBOWpdf", "PDF"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveELBOWtiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(
                                                    #tagList(
                                                    #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                                    plotOutput("ELBOW_plot", height = "300", width = "360"),
                                                    #HTML('</p></div>')
                                                  )
                                         )
                                ),
                                
                                tabPanel(title = "Resolution Parameter", id = "Tab6", value = "Tab6",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         tags$h3("Resolution Parameter"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("umapResolion_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveumapResolionpng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveumapResolionjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveumapResolionpdf", "PDF"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveumapResoliontiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(
                                                    #tagList(
                                                    #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                                    plotOutput("UMAP_Resolution_plot", height = "300", width = "360"),
                                                    #HTML('</p></div>')
                                                  )
                                         )
                                )
                    ),HTML('</p></div>')
                  ),
           
           tabPanel(title = "Step 3: Data Analysis", id = "Tab5", value = "Tab5",
                    HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
                    tabsetPanel(id = "Data_Analysis",
                                tabPanel(title = "UMAP Plot", id = "Tab1", value = "Tab1",
                                         tags$h3("UMAP Plot"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("UMAP_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveUMAPpng","PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveUMAPjpg","JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveUMAPtiff","TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(#jqui resizable canvass
                                                    plotOutput("UMAP_plot", height = "350", width = "350"),
                                                  )
                                         )
                                ),#End scatter plot
                                
                                tabPanel(title = "Feature Plots", id = "Tab2", value = "Tab2",
                                         #HTML('<div style="padding:5px 5px 5px 20px;"<p>'),
                                         tags$h3("Feature Violin Plots"),
                                         tags$hr(),
                                         fluidRow(width=40,
                                                  dropdown(label = "Save Plot",
                                                           selectInput("FeatureCluster_dpi", label = "Output dpi", width = "95",
                                                                       choices = list("24 dpi" = 24, "48 dpi" = 48, "96 dpi" = 96, "192 dpi" = 192),
                                                                       selected = 48),
                                                           downloadButton("saveFeatureClusterpng", "PNG"),
                                                           HTML('<br><br>'), #<br> creates breaks and moves item to next line
                                                           downloadButton("saveFeatureClusterjpg", "JPG"),
                                                           HTML('<br><br>'),
                                                           downloadButton("saveFeatureClustertiff", "TIFF"),
                                                           size = "default",
                                                           icon = icon("download", class = ""), 
                                                           up = FALSE),
                                                  HTML('</p>'),
                                                  jqui_resizable(
                                                    #tagList(
                                                    #HTML('<div style="border:1.5px solid grey;padding:5px 5px 0px 0px;width:5;background-color: #FFFFFF;float:left;display:block;"><p>'),
                                                    plotOutput("Featurecluster_plot", height = "300", width = "360"),
                                                    #HTML('</p></div>')
                                                  )
                                         )
                                )
     
                    ),HTML('</p></div>')
                  ),
           
           
           tabPanel(title = "Help Page", id = "Tab6", value = "Tab6",
                    HTML('<div style="padding:10px 10px 10px 10px;"><p>'),
                    fluidRow(
                      tags$h3("Frequently Asked Questions (FAQ)"),
                      tags$hr(),
                    ),HTML('</p></div>')
                  ),
       tags$script("
           $('body').mouseover(function() {
            list_tabs=[];
            $('#program li a').each(function(){
              list_tabs.push($(this).html())
                 });
              Shiny.onInputChange('List_of_tab', list_tabs);})"
                  )
               ),
    
          navbarMenu("Advanced Settings",
                
                     tabPanel("Normalization Settings"
                              ),
                     "-----",
                     tabPanel("Feature Selection"),
                     "-----",
                     tabPanel("Clustering Settings"),
                     "-----",
                     tabPanel("Biomarker Discovery Settings")
                )
    
    
), #End to navbarPage
  

#Link to CSS file
tags$head(tags$link(rel = "stylesheet",type = "text/css", href = "style.css")),



    
), #End DashboardBody
        
  
  skin = "blue"
  
) #End of UI component
    

###########################################################
  
shinyServer <- function(input, output, session) {
    # ---- DATA UPLOAD / input --------- #
    #reading in the user input and saving it as sample

    
    #Choosing and displaying directory:
    shinyDirChoose(
      input,
      'dir',
      roots = c(home = '~'),
      filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
       )
    global <- reactiveValues(datapath = getwd())
    dir <- reactive(input$dir)
    output$dir <- renderText(
      {
      path = global$datapath
      return(path)
      }
    )
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                   input$dir
                             },
                 handlerExpr = {
                   if (!"path" %in% names(dir())) 
                     return()
                   home <- normalizePath("~")
                   global$datapath <-
                     file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                               }
                             )
  
    
    #Display directory:
    output$selected_dir <- renderText({
      path_display <- global$datapath
      return(path_display)
      }
    )
    
    #Load in the file for DropletUtils:
    dropletutils <- reactive({
      #Add ProgressBar
      total_steps  = 3
      current_step = 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Initializing...."))
      path<-global$datapath
      
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Reading files into DropletUtils...."))
      
      sce <-  DropletUtils::read10xCounts(path, col.names=T)
      
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Finished loading files into DropletUtils...."))
      return(sce)
      }
    )
    
    #Annotate the genes:
    annotate <- reactive({
      sce<-dropletutils()
      ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
      
      attr.string = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name')
      attr.string = c(attr.string, 'start_position', 'end_position', 'strand')
      attr.string = c(attr.string, 'description', 'percentage_gene_gc_content')
      attr.string = c(attr.string, 'gene_biotype')
      rowData(sce)[1:2,]
      gene.annotation = getBM(attributes=attr.string, 
                              filters =  'ensembl_gene_id', 
                              values = rowData(sce)$ID, 
                              mart = ensembl)
      #gene.annotation[which(gene.annotation$ensembl_gene_id %in% names(t2)),]
      gene.annotation = distinct(gene.annotation, ensembl_gene_id, 
                                 .keep_all = TRUE)
      gene.missing = setdiff(rowData(sce)$ID, gene.annotation$ensembl_gene_id)
      w2kp = match(gene.annotation$ensembl_gene_id, rowData(sce)$ID)
      sce  = sce[w2kp,]
      rowData(sce)  = gene.annotation
      rownames(sce) = uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, 
                                           rowData(sce)$hgnc_symbol)
      return(sce)
    })
    
    #Apply emptydrop:
    emptydrop<-reactive({
      sce<-annotate()
      set.seed(100)
      e.out = emptyDrops(counts(sce))
      return(e.out)
    })
    
    #is cell:
    cell<- reactive({
      sce<-annotate()
      e.out<-emptydrop()
      is.cell <- (e.out$FDR <= 0.01) 
      filtered.counts<-counts(sce)[,which(is.cell)]
      return(filtered.counts)
    })
    
    #Create Seurat Object
    seurat<- reactive({
      filtered<-cell()
      object <- CreateSeuratObject(counts=filtered)
      return(object)
    })
    
    
    #Percent_ribo<-
    percent_ribo_mito<-reactive({
      object<-seurat()
      rb.genes <- rownames(object)[grep("^RP[SL]",rownames(object))]
      mt.genes <- rownames(object)[grep("^MT-",rownames(object))]
      C<-GetAssayData(object = object, slot = "counts")
      percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100 #calculates percent mitochondrial reads per barcode/ cell
      percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
      object2 <- AddMetaData(object, percent.ribo, col.name = "percent.ribo") #Add metadata
      object3 <- AddMetaData(object2, percent.mito, col.name = "percent.mito")
      return(object3)
    })
    
    #Apply filters:
    applyfilter<-reactive({
      object<-percent_ribo_mito()
      
      object2 <- subset(object, subset = nFeature_RNA > strtoi(input$Minimum_features) & nFeature_RNA < strtoi(input$Maximum_features) & percent.mito < strtoi(input$Mito_percent) & percent.ribo > strtoi(input$Mito_percent))
      return(object2)
    })
    
    #normalize dataset:
    scTransform <- reactive({
      object2<- applyfilter()
      object3<- SCTransform(object2)
      return(object3)
    })
    
    #Normalize data using log normalization:
    log_normalize<-reactive({
      object2<- applyfilter()
      object3 <- NormalizeData(object2, normalization.method = "LogNormalize", scale.factor = 10000)
      return(object3)
    })
    
    
    
    #feature selection:
    featureSelect <- reactive({
      object2<- scTransform()
      object2 <- FindVariableFeatures(object2, selection.method = "vst", nfeatures = 2000)
      return(object2)
    })
    
    #Scaling data:
    scaleData<- reactive({
      object2<-featureSelect()
      object2 <- ScaleData(object2) 
      return(object2)
    })
    
    ##compute PCA values:
    PCA_compute<-reactive({
      object2<-scaleData()
      object3 <- RunPCA(object2, features = VariableFeatures(object = object2))
      return(object3)
    })
    
  
    #Computing clusters:
    compute_neighbors<-reactive({
      object2<-PCA_compute()
      object3 <- FindNeighbors(object2, dims = 1:40)
      return(object3)
          })
    
    
    
    #Set resolution 
    compute_resolution<-reactive({
      object2<-compute_neighbors()
      object3 <- FindClusters(object2, resolution = 0.8) 
      object4 <- RunUMAP(object3, dims = 1:10)
      return(object4)
    })
    
    
    #Optimize resolution parameter
    compute_resolution_variable<-reactive({
      object2<-compute_neighbors()
      object3 <- FindClusters(object = object2,
                              resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
      return(object3)
    })
    
    set_resolution_ident<-reactive({
      object2<-compute_resolution_variable()
      Idents(object = object2) <- "RNA_snn_res.0.8"
      return(object2)
    })
    
    
    #Compute UMAP
    UMAP_per_resolution<-reactive({
      object2<-set_resolution_ident()
      object3 <- RunUMAP(object2, dims = 1:10)
      return(object3)
    })

    
    #Find cluster
    clusters<-reactive({
      object2<-UMAP_per_resolution()
      object2.markers <- FindAllMarkers(object2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #Finds all markers
      return(obect2)
    })
    
    
    ##################### OUTPUTS MAIN PAGE #####################

    #display barcodes:
    output$barcode_count <- reactive({
      withProgress(message = 'Reading files...', value=1, min=1,max=100,{
        object<- dropletutils()
        barcodes<-ncol(object)
        return(barcodes)
          }
        )
      }
    )
    
    
    #################### FUNCTIONS #####################
    
    VlnPlot_2 <- function(object, features.plot, xlab) {
      
      # Main function
      main_function <- function(object = object, features.plot = features.plot, xlab = xlab) {
        VlnPlot(object = object, features.plot = features.plot) + 
          labs(x = xlab)
      }
      
      # Apply main function on all features
      p <- lapply(X = features.plot, object = object, xlab = xlab,
                  FUN = main_function)
      
      # Arrange all plots using cowplot
      # Adapted from Seurat
      # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
      # ncol argument adapted from Josh O'Brien
      # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
      cowplot::plot_grid(plotlist = p, ncol = ceiling(sqrt(length(features.plot))))
    }
    
    ##################### QC PLOTS #####################
    
    output$Scatter_plot = renderPlot({
      withProgress(message = 'Rendering a scatter Plot...', value=1, min=1,max=100, {
        do_Scatter_plot()
              }
            )
          }
        )
    do_Scatter_plot <- reactive({
      e.out <- emptydrop()
      is.cell <- (e.out$FDR <= 0.01) 
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Converting to dba matrix."))
      #Scatter Plot
      plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
           main=input$main_title,
           cex.main=input$Title_fontsize,
           cex.lab=input$Lable_font,
           xlab=input$xaxis_title, 
           ylab=input$yaxis_title,
           cex.axis=input$axis_fontsize,
           cex=0.2, 
           xlim=c(0,2000), 
           ylim=c(0,2000)) #the black points need to be filtered out!
      legend("bottomright", legend=c("Empty Droplets", "Cell"), bty="n", 
             col=c("black", "red"), lty=1, cex=1.2)
          }
        )
    
    
    output$HEATo_plot = renderPlot({
      withProgress(message = 'Rendering a violin Plot...', value=1, min=1,max=100, {
        do_mitoplot()
         }
        )
      }
     )
    do_mitoplot <- reactive({
      object <- applyfilter()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Converting to dba matrix."))
      p <- SCpubr::do_ViolinPlot(object, 
              features = "percent.mito", 
              pt.size = 0, legend.position = "bottom") & theme(axis.title = element_text(size = input$Axis_Font))
      p + labs(y=input$y_mito,title=input$main_mito,x=input$x_mito)

      
      
   
    }
    )
    
    output$MA_plot = renderPlot({
      withProgress(message = 'Rendering a violin Plot...', value=1, min=1,max=100, {
        do_ribo()
       }
      )
     }
    )
    do_ribo <- reactive({
      object <- applyfilter()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Converting to dba matrix."))
      VlnPlot(object, features = "percent.ribo", pt.size = 0) + NoLegend()
     }
    )
    
    output$VOLCANO_plot = renderPlot({
      withProgress(message = 'Rendering a violin Plot...', value=1, min=1,max=100, {
        do_all()
       }
      )
     }
    )
    do_all <- reactive({
      object <- applyfilter()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Converting to dba matrix."))
      VlnPlot(object, features = "nFeature_RNA", pt.size = 0) + NoLegend()
     }
    )
    
    output$four_plot = renderPlot({
      withProgress(message = 'Rendering a scatter Plots...', value=1, min=1,max=100, {
        do_four_plot()
      }
      )
    }
    )
    do_four_plot<- reactive({
      object2 <- applyfilter()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      par(mfrow=c(1,2), bty="n")
      smoothScatter(log10(object2$nCount_RNA), log10(object2$nFeature_RNA), 
                    xlab="log10(Library sizes)", ylab="log10(# of expressed genes)", 
                    nrpoints=500, cex=0.5) #Graph of ncount RNA vs nFeature RNA 
      smoothScatter(object2$percent.ribo, object2$percent.mito,ylim=c(0,100), xlim=c(0,100),xlab="Ribosome prop. (%)", ylab="Mitochondrial prop. (%)",nrpoints=500, cex=0.5) #Graph of percent mito and percent ribo in your features
      abline(h=5,  lty=1)
      abline(v=10, lty=1)
    }
    )
    
    
    
    ##################### PROCESS PLOTS #####################
    
    #Feature Selection
    output$FEATURE_plot = renderPlot({
      withProgress(message = 'Rendering a feature seleciton Plots...', value=1, min=1,max=100, {
        do_FEATURE_plot()
      }
      )
    }
    )
    do_FEATURE_plot<- reactive({
      object2 <- featureSelect()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      top10 <- head(VariableFeatures(object2), 10) #Top 10 most variable genes
      
      plot1=VariableFeaturePlot(object2) 
      plot2<-LabelPoints(plot = plot1, points = top10, repel = TRUE) #label the top 10 HVG
      plot2
    }
    )
    
    #Feature violin plots
    output$FeatureViolin_plot = renderPlot({
      withProgress(message = 'Rendering a feature selection violin Plots...', value=1, min=1,max=100, {
        do_FeatureViolin_plot()
      }
      )
    }
    )
    do_FeatureViolin_plot<- reactive({
      object2 <- featureSelect()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      top5 <- head(VariableFeatures(object2), 5) #Top 10 most variable genes

      VlnPlot(object2, features = top5, pt.size = 0) + NoLegend()
    }
    )
    
    #PCA Heatmaps:
    output$PCAHeatmap_plot = renderPlot({
      withProgress(message = 'Rendering PCA Heatmaps Plots...', value=1, min=1,max=100, {
        do_PCAHeatmap_plot()
      }
      )
    }
    )
    do_PCAHeatmap_plot<- reactive({
      object2 <- PCA_compute()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      DimHeatmap(object2, dims = 1:5, cells = 1000, balanced = TRUE) 
    }
    )
    
    #Elbow plot
    output$ELBOW_plot = renderPlot({
      withProgress(message = 'Rendering PCA Heatmaps Plots...', value=1, min=1,max=100, {
        do_ELBOW_plot()
      }
      )
    }
    )
    do_ELBOW_plot<- reactive({
      object2 <- PCA_compute()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      ElbowPlot(object2, ndims=40)     
      }
    )
    
    
    #JackStraw Plot:
    output$JackStraw_plot = renderPlot({
      withProgress(message = 'Rendering PCA plot...', value=1, min=1,max=100, {
        do_JackStraw_plot()
      }
      )
    }
    )
    do_JackStraw_plot<- reactive({
      object2 <- PCA_compute()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      
      DimPlot(object2, reduction = "pca")
      }
    )
    
    #Resolution UMAPs
    output$UMAP_Resolution_plot = renderPlot({
      withProgress(message = 'Rendering resolution parameter specific UMAP Plots...', value=1, min=1,max=100, {
        do_UMAP_Resolution_plot()
      }
      )
    }
    )
    do_UMAP_Resolution_plot<- reactive({
      object2 <- UMAP_per_resolution()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      DimPlot(object2,
              reduction = "umap",
              label = TRUE,
              label.size = 6)
    }
    )
    
    ##################### ANALYSIS PLOTS #####################
    
    #UMAP
    output$UMAP_plot = renderPlot({
      withProgress(message = 'Rendering UMAP non-linear clustering and plotting the data...', value=1, min=1,max=100, {
        do_UMAP_plot()
      }
      )
    }
    )
    do_UMAP_plot<- reactive({
      object2 <- compute_resolution()
      #Add ProgressBar
      total_steps  = 2
      current_step = 1
      current_step = current_step + 1
      incProgress(current_step/total_steps*100,
                  detail = paste("Generating the plot."))
      DimPlot(object2,
              reduction = "umap",
              label = TRUE,
              label.size = 6)
    }
    )
    
    
    
    
    
    
    
    
    ##################### NAVIGATION BUTTONS #####################
    
    #Previous and next buttons in dashboard body:
    Previous_Button=tags$div(actionButton("Prev_Tab",HTML('
  <div class="col-sm-4"><i class="fa fa-angle-double-left fa-2x"></i></div>')))
    
    Next_Button=div(actionButton("Next_Tab",HTML('
  <div class="col-sm-4"><i class="fa fa-angle-double-right fa-2x"></i></div>')))
    
    output$Next_Previous=renderUI({
      div(column(1,offset=1,Previous_Button),column(1,offset=8,Next_Button))
    })
    observeEvent(input$Prev_Tab,
                 {
                   tab_list=input$List_of_tab
                   current_tab=which(tab_list==input$program)
                   updateTabsetPanel(session,"program",selected=tab_list[current_tab-1])
                 })
    observeEvent(input$Next_Tab,
                 {
                   tab_list=input$List_of_tab
                   current_tab=which(tab_list==input$program)
                   updateTabsetPanel(session,"program",selected=tab_list[current_tab+1])
                 })
    
    
    ##################### BUTTONS #####################
    observeEvent(
      input$Start_Analysis, {
        updateTabsetPanel(session = session, 
                          inputId = "program", 
                          selected = "Tab3")
      }
    )
    
} #End of server shiny
  

# Run the application 
shinyApp(ui = ui, server = shinyServer)
