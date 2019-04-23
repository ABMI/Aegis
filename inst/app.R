
##install&require packages
# packages<-function(x){
#   x<-as.character(match.call()[[2]])
#   if (!require(x,character.only=TRUE)){
#     install.packages(pkgs=x,repos="http://cran.r-project.org")
#     require(x,character.only=TRUE)
#   }
#   else
#   {
#     require(x,character.only=TRUE)
#   }
# }
Sys.setlocale(category="LC_CTYPE", locale="C")
options(expressions=500000)

library(maps)
library(mapdata)
library(mapproj)
library(leaflet)
library(dplyr)
library(sqldf)
library(shiny)
library(bindrcpp)
library(pkgconfig)
library(shinyjs)
library(shinydashboard)
library(maptools)
library(SpatialEpi)
library(lubridate)
library(rgdal)
library(gpclib)
library(rgeos)
library(devtools)
devtools::install_github("ohdsi/aegis")
library(AEGIS)
library(googledrive)
library(raster)
library(DatabaseConnector)
library(SqlRender)
library(ggmap)
library(DCluster)


gpclibPermit()
##########################################################################
##########################################################################
##########################################################################
###################### function page #####################################

# if (!require(googledrive)){
#   devtools::install_github("tidyverse/googledrive")
#   require(googledrive)
# } else {
#   require(googledrive)
# }
#
# if (!require(devtools)){
#   install.packages("devtools")
# }
#
# if (!require(raster)){
#   devtools::install_github("cran/raster", ref="2.6-7")
#   require(raster)
# } else {
#   require(raster)
# }
#
# if (!require(DatabaseConnector)){
#   devtools::install_github("ohdsi/DatabaseConnector")
#   require(DatabaseConnector)
# } else {
#   require(DatabaseConnector)
# }
#
# if (!require(SqlRender)){
#   devtools::install_github("ohdsi/SqlRender")
#   require(SqlRender)
# } else {
#   require(SqlRender)
# }
#
# if (!require(AEGIS)){
#   devtools::install_github("ohdsi/aegis")
#   require(AEGIS)
# } else {
#   require(AEGIS)
# }
#
# if (!require(ggmap)){
#   devtools::install_github("dkahle/ggmap")
#   require(ggmap)
# } else {
#   require(ggmap)
# }
#
# if (!require(DCluster)){
#   devtools::install_github("cran/DCluster")
#   require(DCluster)
# } else {
#   require(DCluster)
# }




ggmap::register_google(key = 'AIzaSyAhtAT7_quTf0F1TTT2_rYnScYjO-hRA9Q', account_type = "standard", client = NA, signature = NA, day_limit = 2500, second_limit = 2500)

# temp <- tempfile(fileext = ".rds")
# dl <- drive_download(
#   as_id("1maOgnGqDypjIvH5K6EevY4mcSjhh0RhY"), path = temp, overwrite = TRUE)
# df <- readRDS(temp)
# colnames(df) <- tolower(colnames(df))
#
#
#
# temp <- tempfile(fileext = ".rds")
# dl <- drive_download(
#   as_id("1Cu1MhgI0XSBvFtfmUjXcRSTdS4rwLV5p"), path = temp, overwrite = TRUE)
# Country_list <- readRDS(temp)

# df <- "D:/git/AEGIS_demopage/data/df.rds"
# Country_list <- "D:/git/AEGIS_demopage/data/Country_list.rds"

df <- readRDS("df.rds")
Country_list <- readRDS("Country_list.rds")


Call.Cohortlist_demo<-function(){

  name <- c("All-person", "Vascular disorder")
  id <- c(1:2)

  cohortlist <- data.frame(cbind(name, id))
  cohortlist$cohortname <- paste(cohortlist$id, cohortlist$name, sep="; ")

  return(cohortlist)
}

###################### function page #####################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################




shinyApp(
  # Define UI for dataset viewer application
  ui <- dashboardPage(
    dashboardHeader(title = "AEGIS"),
    dashboardSidebar(sidebarMenu(menuItem("DB connection",tabName= "db" ),
                                 menuItem("Cohorts", tabName = "Cohorts" ),
                                 menuItem("Disease mapping", tabName = "Disease_mapping" ),
                                 menuItem("Clustering",tabName = "Clustering" )
    )
    ),
    dashboardBody(tabItems(
      tabItem(tabName = "db",
              fluidRow(
                titlePanel("Database Connection"),
                sidebarPanel(
                  textInput("ip","IP","0.0.0.0")
                  ,uiOutput("sqltype")
                  ,textInput("CDMschema","CDM Database schema","CDM_DB")
                  ,textInput("Resultschema","CDM Result schema","CDM_result_DB")
                  ,textInput("usr","USER","GUEST")
                  ,passwordInput("pw","PASSWORD","GUEST")
                  #input text to db information
                  ,actionButton("db_load","Load DB")
                  #,hr()
                  #,uiOutput("db_conn")
                  #,actionButton("cohort_load","Load cohort")
                  ,width=2
                ),
                mainPanel(
                  verbatimTextOutput("txt"),
                  tableOutput("view"),
                  h6("This version is only for demo, so it may not support all functions")
                )
              )
      ),
      tabItem(tabName = "Cohorts",
              fluidRow(
                titlePanel("Cohort selection"),
                sidebarPanel(
                  uiOutput("cohort_tcdi")
                  ,uiOutput("cohort_ocdi")
                  ,hr()
                  ,dateRangeInput(inputId = "dateRange", label = "Select Windows",  start = "2002-01-01", end = "2013-12-31")
                  ,hr()
                  ,radioButtons("GIS.Age","Age and gender adjust",choices = c("No" = "no", "Yes"="yes"))
                  #,numericInput("GIS.timeatrisk_startdt","Define the time-at-risk window start, relative to target cohort entry:", 0, min=0)
                  #,numericInput("GIS.timeatrisk_enddt","Define the time-at-risk window end:", 0, min=0)
                  #,selectInput("GIS.timeatrisk_enddt_panel","", choices =
                  #              c("from cohort start date" = "cohort_start_date","from cohort end date" = "cohort_start_date"),selected = "cohort_end_date")
                  #,textInput("fraction","fraction", 100000)
                  ,uiOutput("country_list")
                  ,actionButton("submit_table","submit")
                  ,width=2
                ),
                mainPanel
                (dataTableOutput('GIS.table'),
                  h6("This version is only for demo, so it may not support all functions"))
              )
      ),
      tabItem(tabName = "Disease_mapping",
              fluidRow(
                titlePanel("Disease mapping setting"),
                sidebarPanel(
                  radioButtons("GIS.level","Administrative level",choices = c("Level 2" = 1, "Level 3" = 2),selected = 2)
                  #radioButtons("GIS.level","Administrative level",choices = c("Level 1" = 0, "Level 2" = 1, "Level 3" = 2),selected = 1)
                  ,radioButtons("GIS.distribution","Select distribution options", choices = c("Count of the target cohort (n)" = "count","Propotion" = "proportion", "Standardized Incidence Ratio"="SIR"),selected = "count")
                  #,radioButtons("distinct","Select distinct options", c("Yes" = "distinct","No" = "" ),inline= TRUE)
                  ,textInput("plot.title","title"," ")
                  ,textInput("plot.legend","legend"," ")
                  ,actionButton("submit_plot","submit") #Draw plot button
                  ,width=2
                ),
                mainPanel(
                  #verbatimTextOutput("test")
                  #,
                  plotOutput("GIS.plot"),
                  h6("This version is only for demo, so it may not support all functions")
                )
              )
      ),

      tabItem(tabName = "Leaflet(beta)",
              fluidRow(
                titlePanel("Interactive disease map(beta)"),
                mainPanel(
                  leafletOutput("mappingLeaflet")
                )
              )
      ),


      tabItem(tabName ="Clustering",
              fluidRow(
                titlePanel("Disease clustering"),
                sidebarPanel(
                  radioButtons("Cluster.method","Cluster Method",choices = c("Kulldorff's method" = "kulldorff"))
                  ,textInput("Cluster.parameter","Kulldorff's method parameter", "0.05")
                  ,actionButton("submit_cluster","submit") #Draw plot button
                  ,width=2
                ),
                mainPanel(
                  dataTableOutput('Cluster.table')
                  ,plotOutput("Cluster.plot")
                  ,textOutput("Cluster.test"),
                  h6("This version is only for demo, so it may not support all functions")
                )
              )
      )
    )
    )
  ),


  server <- function(input, output,session)
  {

    output$sqltype <- renderUI({
      selectInput("sqltype", "Select DBMS",
                  choices = c(
                    "sql server" = "sql server",
                    "PostgreSQL" = "postresql",
                    "Amazon Redshift" = "redshift",
                    "Microsoft Parallel Data Warehouse" = "pdw",
                    "IBM Netezza" = "netezza",
                    "Google BigQuery" = "bigquery"
                              )
                  )
    })


    output$cohort_tcdi <- renderUI({
      cohort_list <- Call.Cohortlist_demo()
      if (length(cohort_list)>1) {
        selectInput("tcdi", "Select target cohort", choices = cohort_list[1,3])
      } else {
        selectInput("tcdi", "Select target cohort", choices = cohort_list[1,1])
      }
    })


    output$cohort_ocdi <- renderUI({
      cohort_list <-Call.Cohortlist_demo()
      #selectInput("ocdi", "Select outcome cohort", choices = cohort_list[,3])
      if (length(cohort_list)>1) {
          selectInput("ocdi", "Select target cohort", choices = cohort_list[2,3])
        } else {
          selectInput("ocdi", "Select target cohort", choices = cohort_list[2,1])
        }
    })


    output$country_list <- renderUI({
      country_list <- Country_list
      selectInput("country", "Select country", choices = country_list[,1])
    })

    cohort_listup <- eventReactive(input$db_load, {
      cohort_list <- cohort
    })


    render.table <- eventReactive(input$submit_table,{
      isolate({
        country_list <<- Country_list
        country <<- "United States"
        MAX.level <<- 2
        GADM <<- readRDS("GADM_2.8_United States_adm_total.rds")
        GADM.table <<- GADM[[3]]@data

        #Conditional input cohort number
        # if (length(cohort_list)>1) {
        #   tcdi <- cohort_list[which(cohort_list[,3] %in% input$tcdi == TRUE),1]
        #   ocdi <- cohort_list[which(cohort_list[,3] %in% input$ocdi == TRUE),1]
        # } else {
        #   tcdi <- input$tcdi
        #   ocdi <- input$ocdi
        # }
        tcdi <- 1
        ocdi <- 2


        colnames(df) <- tolower(colnames(df))
        df[, c("outcome_count", "target_count")][is.na(df[, c("outcome_count", "target_count")])] <- 0


        CDM.table <<- GIS.Indirect.AgeGenderadjust(df, 100000)


        table <- dplyr::left_join(CDM.table, GADM.table, by=c("gadm_id" = "ID_2"))
        switch(input$GIS.Age,
               "no"={
                 table <- table[, c("OBJECTID","ID_0", "ISO", "NAME_0", "ID_1", "NAME_1", "NAME_2",
                                    "target_count", "outcome_count", "crd_expected", "crd_prop", "crd_sir"
                 )]#"ID_2"

               },
               "yes"={
                 table <- table[, c("OBJECTID","ID_0", "ISO", "NAME_0", "ID_1", "NAME_1",  "NAME_2",
                                    "target_count", "outcome_count", "std_expected", "std_prop", "std_sir"
                 )]#"ID_2",

               }

        )
      })
      table
    })

    output$GIS.table <- renderDataTable(
      render.table()
    )


    draw.plot <- eventReactive(input$submit_plot,{

        countdf_level <<- GIS.calc0(input$GIS.level, input$GIS.distribution, input$GIS.Age)
        mapdf <<- GIS.calc3(input$GIS.level, input$fraction)
        plot <- GIS.plot2(mapdf,input$GIS.distribution, input$plot.legend, input$plot.title, input$GIS.Age)
      plot
    })


    output$GIS.plot <- renderPlot ({
      draw.plot()
    }, width = 1280, height = 1024, res = 100)


     output$mappingLeaflet <- renderLeaflet({
      leafletMapping()
     })

    #testing.cluster <- eventReactive(input$submit_cluster,{
    #  isolate({
    #    CDM.table$Observed <- CDM.table$outcome_count
    #    test.summ <- DCluster::achisq.stat(CDM.table, lambda=1)
    #  })
    #  test.summ[[1]]
    #})

    #output$Cluster.test <- renderText({
    #  testing.cluster()
    #})

    #finding.cluster <- eventReactive(input$submit_cluster,{
    #  isolate({
    #      table <- Cluster.find(input$Cluster.method, input$Cluster.parameter)
    #  })
    #  table
    #})

    #output$Cluster.table <- renderDataTable(
    #  finding.cluster()
    #)

    plotting.cluster <- eventReactive(input$submit_cluster,{
      isolate({
        plot <- Cluster.plot2(input$Cluster.method, input$Cluster.parameter, input$GIS.Age, input$country)
      })
      plot
    })



    output$Cluster.plot <- renderPlot ({
      plotting.cluster()
    }, width = 1024, height = 800, res = 100)


    ## End of server
  }, options = list(height = 1000)
)
