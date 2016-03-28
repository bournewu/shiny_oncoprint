
#require(limma)
#require(DESeq2)
#require(limma)
#require(edgeR)
# require(genefilter)

shinyServer(
  function(input, output, session) {
    
    readExpr <-  reactive({
      # reformed the datasets based on users selection on samples    
      read_expr_data(input$mutdata, input$exprdata)      
    })
    
    readMut <-  reactive({
      read_mut_data(input$mutdata)
    })
    
    get_click_genes <- reactive({
      validate(
        need(length(input$exprTable) > 1, "Please click at least one gene."
        ))
      input$exprTable
    })
    
    getInputGenes <- reactive({
      validate(
        need(input$genesets != "", "Please type at least one mutated gene.")
      )
      input$genesets
    })
    

    oncoprint_dataprep <- eventReactive(input$goButtonOnco,{
      if (getInputGenes() != ""){
        readMut() %>% 
          filter(gene %in% getInputGenes())
      }
    })
    
   
    
    output$oncoprint <- renderPlot({              
      plot_oncoprint(oncoprint_dataprep(), input$showCoMEt, input$CoMEtPval)            
    })
    
   
    
    
    
    #source("src/generateData.R", local=T)
    source("src/renderUI.R", local=T)
  }
)


