
# UI
output$buttonToGo <- renderUI({
  if (getInputGenes() != ""){
    actionButton("goButtonOnco", "Show results")
  }
})
output$geneSetUI <- renderUI({
  
  mut_data <- readMut()
  g2freq <- bind_cols(
    (mut_data %>%
       select(-gene) %>% 
       mutate(freq=rowSums(.)) %>% select(freq)), 
    (mut_data %>% select(gene))) %>%
    mutate(key=paste0(gene,"[",freq,"]")) %>% 
    arrange(desc(freq)) %>% select(c(key,gene))
  
  g2freq_index <- setNames(as.character(g2freq$gene), g2freq$key )
  selectInput('genesets', 'Input genes:', g2freq_index, 
              multiple=TRUE, selectize=TRUE)
})

output$cometExactTest <- renderUI({
  if (getInputGenes() != ""){
    checkboxInput('showCoMEt', 'Show CoMEt pvalue:', value=F)
  }
})

output$logFCUI_t <- renderUI({
    sliderInput('logFCCut_t', 'Fold change cutoff:', min = 0, max=3, value = 1, 
                step=0.1)
  
})
output$aPvalUI_t <- renderUI({
  sliderInput('aPvalCut_t', 'Adj.P.Val cutoff:', min = 0, max=1, 
              value = 0.1, step=0.01)
})


output$slider1 <- renderUI({
  if (input$networkview == "hotnet2")
    sliderInput("edgeC", "Threshold of edge weight: ",min = 0,
                max = 0.1, 
                value = 0.022,
                step= 0.1/50 )
  else
    sliderInput('logFCCut', 'Fold change cutoff:', min = 0, max=3, value = .5, 
                step=0.1)
  
})
output$slider2 <- renderUI({
  if (input$networkview == "hotnet2")
    sliderInput("sizeC", "Minimum subnetwork size k:", min = 2, 
                max = 10, 
                value = 2 ,
                step= 1)
  else
    sliderInput('aPvalCut', 'Adj.P.Val cutoff:', min = 0, max=1, 
                value = 0.1, step=0.01)
})

