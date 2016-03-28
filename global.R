require(shiny)
require(dplyr)
require(tidyr)
require(cometExactTest)

select <- dplyr::select
slice <- dplyr::slice

source("src/data.R")
source("src/data_selection.R")
source("src/oncoprint.R")
#source("src/de_analysis.R")
source("src/utils.R")
#source("src/hn2.R")
#source("src/gene_set_analysis.R")

select <- dplyr::select
rename <- dplyr::rename

clean_NA <- function(tt){
  tt[complete.cases(tt) ,]
}

format_2 <- function(expr){
  as.numeric(format(round(expr,4),nsmall=4))
}

returnTextInput <- function(inputId, label, value = "") {
  tagList(
    singleton(tags$head(tags$script(src = "js/returnTextInputBinding.js"))),
    tags$label(label, `for` = inputId),
    tags$input(id = inputId, type = "text", value = value, class = "returnTextInput")
  )
}

# selectable data table from
# https://www.google.com/search?q=selDataTableOutput&oq=
# selDataTableOutput&aqs=chrome..69i57.181j0j7&sourceid=chrome&es_sm=91&
# ie=UTF-8
selDataTableOutput <- function (outputId) 
{
  tagList(singleton(
    tags$head(
      tags$link(
        rel = "stylesheet", 
        type = "text/css", 
        href = "shared/datatables/css/DT_bootstrap.css"),
      tags$style(type="text/css", ".rowsSelected td{background-color: 
                 rgba(112,164,255,0.2) !important}"),
      tags$style(type="text/css", ".selectable div table tbody tr{cursor: 
                 hand; cursor: pointer;}"),
      tags$style(type="text/css",".selectable div table tbody tr td{
               -webkit-touch-callout: none;
               -webkit-user-select: none;
               -khtml-user-select: none;
               -moz-user-select: none;
               -ms-user-select: none;
               user-select: none;}"),                          
      tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"), 
      tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
      tags$script(src = "js/DTbinding.js")
    )
  ), 
  div(id = outputId, class = "shiny-datatable-output selectable"))
}


# readPathway <-  function(){       
#   pathway<-readRDS("data/gsetdbs/kegg_austin.rds")
#   for (i in 1:length(pathway)){
#     if (class(pathway[[i]]) == "factor") {
#       pathway[[i]] <- levels(droplevels(pathway[[i]]))
#     }
    
#     pathway[[i]] <- pathway[[i]][!pathway[[i]]%in%""]
#   }
#   pathway
# }

########################################
# GBM
########################################
# read expression
# read mutations
gbmMUT <- readRDS("data/gbm/mut_data/GBM-EGFRA.r.matrix.rds")

amlMUT <- readRDS("data/aml/mut_data/aml.r.matrix.rds")

mutDataAll <- list(gbm=gbmMUT, aml=amlMUT)
#rawDataAll <- list(gbm=gbmRAW, aml=amlRAW)
########################################
# BRCA
########################################


#humanMart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# color scales
# category20
color_cat20 <- "'#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', 
'#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', 
'#8c564b', '#e377c2', '#bcbd22', '#17becf'"

color_cat20_plus1 <- "'#DDDDDD', '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', 
'#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', 
'#c5b0d5', '#8c564b', '#e377c2', '#bcbd22', '#17becf'"

color_blue_white_red <- "'#3232FF','#FFF','#FF3333'"

color_red <- "'#FFCBCB', '#FF3333'"

# RegTf_n <- 
#   read.table(file = "data/networks/RegTf/RegTf_viz_nodes.tsv", 
#              header=T) %>% 
#   clean_NA()
# RegTf_e <- 
#   read.table(file = "data/networks/RegTf/RegTf_viz_edges.tsv", 
#              header=T) %>% 
#   clean_NA()

#pathway_dbs <- readPathway() 
