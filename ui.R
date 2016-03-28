# ui.R
require(shinythemes)
#require(networkD3)

makeSideBar <- function(){
  absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                draggable = FALSE, top = 60, left = 20, right = "auto", bottom = "auto",
                width = 300, height = "auto",
                #sidebarPanel(width = 3,
                
                selectInput("mutdata", "Mutation Data:",
                                   choices = c("GBM"="gbm", "AML"="aml"),
                                   selected = "GBM"
                ),
                sliderInput("CoMEtPval", "CoMEt pvalue cutoff:", 0, 1, 0.1, 0.1)
  )
  
}

makeHome <- function(){
  tabPanel('Home', value = 0,
           sidebarLayout(
             sidebarPanel(width=3),
             mainPanel(
               fluidRow(wellPanel(uiOutput('geneSetUI'), 
                                  style = "padding: 5px;")),
               fluidRow(column(6, wellPanel(uiOutput('buttonToGo'), 
                                            style = "padding: 5px;")), 
                        column(6, wellPanel(uiOutput('cometExactTest'), 
                                            style = "padding: 5px;"))),
               tags$hr(),
               tags$h4("Oncoprint:"),
               plotOutput("oncoprint")
             )
           )
  )
}

shinyUI(
  fluidPage(theme = shinytheme('cerulean'), 
            
            tags$head(
              tags$script(src ="https://code.highcharts.com/highcharts.js"),
              tags$script(src = 
                            "https://code.highcharts.com/highcharts-more.js"),
              tags$script(src = 
                            "https://code.highcharts.com/modules/exporting.js"),
              tags$script(src = 
                            "https://code.highcharts.com/modules/heatmap.js"),
              tags$link(rel="stylesheet", type="text/css",href="style.css"),
              tags$script(type="text/javascript", src = "js/busy.js")
            ),            
            navbarPage("Mutation and Expression data analysis", 
                       id="conditionedPanels",
                       makeHome()
                       )  ,
            div(class = "busy",                  
                img(src="loading2.gif", width="300", height="300")
            ),
            makeSideBar()
  )
)


