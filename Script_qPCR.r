#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Script for analyzing qPCRs"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        fileInput("file1","Choose data file",accept = "text/csv"),
        fileInput("file2","Choose design file",accept = "text/csv"),
        textInput("normalized_group","Choose group to which normalize the RQ"),
        
        actionButton("GO","Analizar datos"),
        downloadButton("downloadData","Descargar datos graficos"),
        downloadButton("downloadReport","Descargar grafica"),
        downloadButton("downloadStats","Descargar analisis estadistico")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        
        tabsetPanel(
          tabPanel("Data", dataTableOutput("data")),
          tabPanel("Graphic Data", dataTableOutput("graphicdata")),
          tabPanel("Endogen Plot",plotOutput("endogenPlot",height = "800")),
          tabPanel("Outlier Plot", plotOutput("outlierPlot",height = "800")),
          tabPanel("Expression Plot", plotOutput("expressionPlot",height = "800")),
          tabPanel("Statistical Analysis", verbatimTextOutput("tests"))
        )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  expressionPlot <- geom_blank()
  endogenPlot <- geom_blank()
  outlierPlot <- geom_blank()
  observeEvent(input$GO, {
    
    library(dplyr)
    library(ggplot2)
    
    file1 <- input$file1
    datafile <- file1$datapath
    
    file2 <- input$file2
    designfile <- file2$datapath
    
    data = read.csv2(file= datafile,header=TRUE)
    design = read.csv2(file=designfile,header=TRUE)
    
    samples = unique(data$Sample)
    genes = unique(data$Detector)
    
    avgCtTable = aggregate(data$Ct, by = list(data$Sample, data$Detector), FUN = mean,na.rm=TRUE )
    colnames(avgCtTable) <- c("Sample","Detector","avgCt")
    
    sdCtTable = aggregate(data$Ct, by = list(data$Sample, data$Detector), FUN = sd,na.rm=TRUE )
    colnames(sdCtTable) <- c("Sample","Detector","sd")
    
    for (i in 1:length(data$Ct)) {
      sample = as.character(data[[i,"Sample"]])
      detector = as.character(data[[i,"Detector"]])
      Ct = data[i,"Ct"]
      avgCt = avgCtTable[which(avgCtTable$Sample== sample & avgCtTable$Detector==detector),"avgCt"]
      sd = sdCtTable[which(sdCtTable$Sample== sample & sdCtTable$Detector==detector),"sd"]
      
      if(!is.na(Ct) & !is.na(sd)) {
        if(Ct>(avgCt + sd) | Ct<(avgCt - sd)) {
          data[i,"Ct"] = NA
        }
      }
    }
    
    avgCtTable = aggregate(data$Ct, by = list(data$Sample, data$Detector), FUN = mean,na.rm=TRUE )
    colnames(avgCtTable) <- c("Sample","Detector","avgCt")
    
    gene_info = unique(select(data,Detector,Task))
    avgCtTable=merge(avgCtTable,gene_info)
    
    
    endoData = avgCtTable[which(avgCtTable$Task=="ENDO"),]
    endoData = endoData %>% select(Sample,avgCt)
    colnames(endoData) <- c("Sample","EndoCt")
    
    deltaCtTable = merge(avgCtTable,endoData)
    deltaCtTable = mutate(deltaCtTable,deltaCt = avgCt - EndoCt)
    deltaCtTable = mutate(deltaCtTable,RQ = 2^(-deltaCt))
    
    
    deltaCtTable = merge(deltaCtTable,design)
    
    output$data <- renderDataTable(deltaCtTable)
    
    std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))
    
    graphic_data =  filter(deltaCtTable,Group != "") ##Quitamos asi los blancos
    graphic_data =  filter(graphic_data,Task != "ENDO") ##Y ahora el endogeno
    gd1 = aggregate(graphic_data$RQ, by = list(graphic_data$Group,graphic_data$Detector), FUN =mean, na.rm = TRUE )
    colnames(gd1) <- c("Group","Detector","RQ")
    
    gd2 = aggregate(graphic_data$RQ, by = list(graphic_data$Group,graphic_data$Detector), FUN =std )
    colnames(gd2) <- c("Group","Detector","SEM")
    
    graphic_data = merge(gd1,gd2)
    
    normal_group = input$normalized_group
    WTData = graphic_data[which(graphic_data$Group==normal_group),]
    WTData = WTData %>% select(Detector,RQ)
    colnames(WTData) <- c("Detector","WT_RQ")
    
    graphic_data = merge(graphic_data,WTData)
    graphic_data = mutate(graphic_data,N_RQ = RQ/WT_RQ)
    graphic_data = mutate(graphic_data,N_SEM = SEM/WT_RQ)
    
    output$graphicdata <- renderDataTable(graphic_data)
    
    output$downloadData <- downloadHandler(
      
      filename = "graph_data.csv",
      
      content = function(file) {
        write.csv2(graphic_data,file,row.names = FALSE)
      }
    )
    
    avgCtTable =  filter(avgCtTable,Task == "ENDO") ##Y ahora el endogeno
    
    endogenPlot <<- ggplot(avgCtTable,aes(x=Detector,y=avgCt,label=Sample)) + 
      geom_point() +
      geom_label_repel() +
      labs(y = "avgCt",title = "Engenous gene expression",caption= Sys.Date()) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    deltaCtTable =  filter(deltaCtTable,Group != "") ##Quitamos asi los blancos/outliers
    deltaCtTable =  filter(deltaCtTable,Task != "ENDO") ##Y ahora el endogeno
    
    outlierPlot <<- ggplot(deltaCtTable,aes(x=Group,y=RQ,colour=Group,label=Sample)) + 
      geom_point(position = "dodge") +
      geom_text(check_overlap = TRUE) +
      labs(y = "RQ",title = "Gene expression for visualizing outliers",caption= Sys.Date()) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~Detector,scales = "free")
  
    dodge <- position_dodge(width=0.9)
    
    expressionPlot <<- ggplot(graphic_data,aes(x=Detector,y=N_RQ,fill=Group)) + 
      geom_col(position = dodge) + 
      geom_errorbar(aes(x= Detector,ymin=N_RQ,ymax=N_RQ + N_SEM),width = 0.25,position = dodge) + 
      labs(y = "Normalized RQ",title = "Relative gene expression",caption= Sys.Date()) +
      theme(plot.title = element_text(hjust = 0.5))
  
    output$endogenPlot <- renderPlot({
      endogenPlot
    }) 
    
    output$outlierPlot <- renderPlot({
      outlierPlot
    })
    
    output$expressionPlot <- renderPlot({
       expressionPlot
    })
    
    out <- character()
   
    output$tests <- renderPrint({
       
       endo = as.character(gene_info[gene_info$Task=="ENDO",1])
       genes = genes[genes != endo]
       
       out<<-capture.output(for (i in genes) {
         test <-  dplyr::filter(deltaCtTable,deltaCtTable$Group != "") ## Se quitan los blancos
         test <-  dplyr::filter(test,Task != "ENDO") ## Se quita el endogeno
         test <-  dplyr::filter(test,Detector == i)
         
         print(i)
         anova <- aov(data=test,deltaCt~Group)
         print(summary(anova))
         print(TukeyHSD(anova))
       })
       
      out
    })
    
    output$downloadReport <- downloadHandler(
      
      filename = "expressionPlot.png",
      
      content = function(file) {
        4
        ggsave(filename = file,plot = expressionPlot,device="png",width = 50,height = 25, units = "cm", dpi=300)
      }
    )
    
    output$downloadStats <- downloadHandler(
      
      filename = "stats.txt",
      
      content = function(file) {
        write.table(as.data.frame(out),file=file,sep="\n")
      }
    )

   
  
  })  
}

# Run the application 
shinyApp(ui = ui, server = server)

