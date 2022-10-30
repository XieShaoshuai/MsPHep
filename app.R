library(shiny)
library(shinydashboard)
library(markdown)
library(ggplot2)
library(cowplot)
library(IsoSpecR)
library(magick)
library(reactable)
library(pracma)
library(dplyr)
library(ggpmisc)
library(smoother)
library(DT)

source('structure_searching.R')
source('iso_prob.R')
#source('fragment.R')
#source('fragment_cal.R')
source('hepQuan.R')
source("get_summary.R")


options(shiny.maxRequestSize = 300*1024^2)

#database <- read.csv("mds/fragment database.csv",row.names=1, sep=";")

#header----
header <- dashboardHeader(title = ("MsPHep") )


#sidebar----
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName='home',icon = icon("home")),
    tags$hr(style="border-color: grey;"),
    
    menuItem("HepQual",tabName='MS',icon=icon("search")),
    
    menuItem("HepQuan",tabName='HepQuan',icon=icon("search")),
    
    #menuItem("MSMS",tabName='MSMS',icon=icon("cut")),
    tags$hr(style="border-color: grey;"),
    
    menuItem('Contact',tabName='help')
  )
  
)


#body----
body <-  dashboardBody(
  tabItems(
    #home-----
    tabItem(tabName = "home",
            fluidRow(
              column(width = 10,
                     includeMarkdown("mds/welcome.md"))
            )
              ),
    #????????????----
    tabItem(tabName = "MS",
            fluidRow(
              box(title=strong('User Guide'),status='warning',width=10,
                  solidHeader = FALSE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  closable = FALSE,
                  includeMarkdown("mds/structure find.md")
                 ),
            ),
            fluidRow(
              box(title=strong("HepQual"),status = "info",width=3,
                  h4('Step 1: Select Database'),
                  selectInput(
                    'lmwh', '', choices = c('Enoxaparin'="Eno",'Dalteparin'="Dal_Nadro",'Nadroparin'="Dal_Nadro",'Heparin/HS'='HS'),
                    selectize = FALSE
                  ),
                  h4("Step2: Adductive form"),
                  selectInput('method','',choice=c('NH3'='_hilic.txt',"No"='_sec.txt')),
                  h4("Step3: Input m/z and charge"),
                  splitLayout(numericInput("mz", "m/z:",min = 0, max = 5000, value =581.6388),
                              numericInput("charge", "charge:",min = 1, max = 15, value = 3)),
                  h4("Step4: Parameters"),
                  numericInput("ppm", "ppm:",
                              min = 0, max = 30,
                              value = 10, ),
                  
                  selectizeInput('isopeak', 'Is input m/z Monoisotopic peak?', choices = c('Yes' = 'yes', 'No' = 'no')
                  , multiple = FALSE),
                  actionButton("search", "Search!"),
                  
              ),
              box(title="Result:",width=7,
                  solidHeader = FALSE,status = "success",
                  #tableOutput("table1"),
                  reactableOutput("table2"),
              ),
            ),
    ),
    #??????????????????----
    tabItem(tabName ="HepQuan",
            fluidRow(
              box(title=strong('User Guide'),status='warning',width=12,
                  solidHeader = FALSE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  closable = FALSE,
                  includeMarkdown("mds/HepQuan.md")
              ),
            ),
            fluidRow(
              box(title=strong("HepQuan"),status = "info", width=3,
                  fileInput("scan", "Upload scan.csv",
                            multiple = FALSE,
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                  fileInput("iso", "Upload iso.csv",
                            multiple = FALSE,
                            accept=c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")),
                  selectInput(
                    'lmwh2', 'LMWH type', choices = c('Enoxaparin'="Eno",'Dalteparin'="Dal_Nadro",'Nadroparin'="Dal_Nadro",'Heparin/HS'='HS'),
                    selectize = FALSE
                  ),
                  #h4("Adductive"),
                  selectInput('add2','Adductive',choice=c('NH3'='_hilic.txt',"No adductive"='_sec.txt')),
                  numericInput("ppm2", "ppm:",min = 0, max = 30, value =15),
                  splitLayout(
                    textInput("dp_range", "dp range", value = "dp4-dp30"),
                    
                    numericInput("minscan", "Min Scan number",min = 0, max = 200, value = 15)),
              
                  
                  actionButton("search3", "Search!"),
                  
              ),
              box(title=strong("Result:"),width=9,
                  downloadButton("downloadData", "Download raw data"),
                  downloadButton("downloadData2", "Download merge data"),
                  DT::dataTableOutput("summary"),
                  
                  DT::dataTableOutput("details"))
            ),
    ),
    
    
    #????????????-----
    tabItem(tabName = "MSMS",h2("Widgets tab content"),)
    )
)


#ui----
ui <- dashboardPage(header,sidebar,body,
                    skin="purple")


#server----
server <- function(input, output,session) {
  #????????????----
  dataset <- reactive({
    file <- paste0("database/",input$lmwh,input$method)
    read.csv(file, sep=";")
  })
  res <- eventReactive(input$search,{
    data <-NULL
    if(input$isopeak=='yes'){
      data <- structure_searching(input$mz,input$charge,input$ppm, dataset())
    }else{
      data <- structure_searching2(input$mz,input$charge,input$ppm, dataset())
    }
     return(data)
  })

  output$table2 <- renderReactable(
    if(input$isopeak=='yes'){
      reactable(res()[,c("Theoretical.MW","Experimental.MW","Structure","Adductive","dp","PPM")],
                sortable = FALSE, pagination = FALSE, showPageInfo = FALSE,
                onClick = "expand",
                details = colDef(details = function(index) {
                  if (index >0) {
                    tabsetPanel(
                      tabPanel(plotOutput("plot",width = "400px", height = "250px")),
                    )
                  } 
                }, name = "Iso.Peak", width = 70))
    }else{
      reactable(res()[,c("Theoretical.MW","Experimental.MW","isotopic.peak","Structure","Adductive","dp","PPM")],rownames = FALSE,
                sortable = FALSE, pagination = FALSE, showPageInfo = FALSE,
                )
            
    }
    
    )
  output$plot <-renderPlot(iso_prob(res(),input$charge))
  
  #??????????????????-----
  dataset2 <- reactive({
    file <- paste0("database/",input$lmwh2,input$add2)
    d <- read.csv(file, sep=";")
    min <- as.numeric(substring(strsplit(input$dp_range,"-")[[1]][1],3,6))
    max <- as.numeric(substring(strsplit(input$dp_range,"-")[[1]][2],3,6))
    d <- d[which(d$dp>=min&d$dp<=max),]
    return(d)
  })
  quan_res <- eventReactive(input$search3,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:10) {
                     req(input$scan)
                     scan <-read.csv(input$scan$datapath,header=TRUE)
                     req(input$iso)
                     iso <- read.csv(input$iso$datapath,header=TRUE)
                     resx <- hepQuan(scan,iso,input$ppm2,dataset2(),input$minscan)
                     return(resx )
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
  }
  )
  
  #output$summary <- DT::renderDataTable(DT::datatable(quan_res(),rownames = FALSE))
  
  proxy_summary <- dataTableProxy("summary")
  proxy_details <- dataTableProxy("details")
  
  get_current_slice <- reactive({
    my_data <- req(quan_res())
    #my_data <-my_data[,c("Structure","Adductive","Exep_mz","charge","Inten","mw.count")]
    my_data %>%
      filter(Structure == get_summary(my_data) %>%
               slice(req(input$summary_rows_selected)) %>%
               pull(Structure))
  })
  

  
  output$summary <-DT::renderDataTable({
    datatable(
      #quan_res(),
      get_summary(quan_res()), 
      extensions = "Buttons",
      rownames   = FALSE,
      #filter     = "top",
      selection  = "single",
      editable   = FALSE,
    )
  })
  
  output$details <- DT::renderDataTable({
    req(input$summary_rows_selected)
    datatable(
      req(isolate(get_current_slice())),
      #extensions = "Buttons",
      rownames   = FALSE,
      #filter     = "top",
      selection  = "single",
      editable   = list(target = "cell", disable = list(columns = c(0:1, 3))),
    )
  })
  
  output$downloadData <- downloadHandler(
    filename="raw_result.csv",
    content=function(fname){
      write.table(quan_res(), fname,sep=";",row.names=FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename="merge_result.csv",
    content=function(fname){
      write.table(get_summary(quan_res()), fname,sep=";",row.names=FALSE)
    }
  )
  
}

shinyApp(ui, server)
