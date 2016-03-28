require(igraph)

## function for finding connected component
find_cc <- function(lst) {
  lst <- apply(lst, 1, as.list)#browser()
  edg <- do.call("rbind", lapply(lst, function(x) {
    if (length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
  }))
  
  g <- graph.data.frame(edg)
  options(stringsAsFactors = FALSE)
  data.frame(group=clusters(g)$membership, name=as.character(V(g)$name))
}

## function for finding strongly connected component
find_strong_cc <- function(in_lst) {
  
  lst <- apply(in_lst, 1, as.list)#browser()
  edg <- do.call("rbind", lapply(lst, function(x) {
    if (length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
  }))
  
  g <- graph.data.frame(edg, directed = TRUE)    
  options(stringsAsFactors = FALSE)
  data.frame(group=clusters(g, mode = "strong")$membership, 
             name=as.character(V(g)$name))
}




### Find overlapping between two hotnet2 results ###
find_overlapping_hn2 <- function(data1, data2, 
                                 edgecut1, edgecut2, 
                                 size1, size2){
  hn2_e <- data1$hn2_e
  deExpr <- data1$meta_data
  plot_n <- data1$nodes %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  filtered_hn2_e <- hn2_e %>% filter(similarity >= edgecut1)    
  if (nrow(filtered_hn2_e) <= 0) hn2_nodes_1 <- NULL
  else{    
    plot_n <- find_hn2_nodes(plot_n, filtered_hn2_e, size1)                    
    hn2_nodes_1 <- plot_n  %>% filter(group > 0) %>% select(name)
  }
  
  hn2_e <- data2$hn2_e
  deExpr <- data2$meta_data
  
  plot_n <- data2$nodes  %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  filtered_hn2_e <- hn2_e %>% filter(similarity >= edgecut2)    
  if (nrow(filtered_hn2_e) <= 0) hn2_nodes_2 <- NULL
  else{    
    plot_n <- find_hn2_nodes(plot_n, filtered_hn2_e, size2)                    
    hn2_nodes_2 <- plot_n  %>% filter(group > 0) %>% select(name)
  }      
  
  inner_join(hn2_nodes_1, hn2_nodes_2, by="name")
  
}

### Find overlapping between two fold change results ###
find_overlapping_fc <- function(data1, data2, noded1, noded2, nodep1, nodep2){
  
  deExpr <- data1$meta_data
  
  plot_n <- data1$nodes %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  selected_nodes_1 <- plot_n %>% 
    filter(abs(FC) >= noded1 & pvalue <= nodep1) %>%
    select(name)
  
  deExpr <- data2$meta_data
  
  plot_n <- data2$nodes  %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  selected_nodes_2 <- plot_n %>% 
    filter(abs(FC) >= noded2 & pvalue <= nodep2) %>%
    select(name)
  
  inner_join(selected_nodes_1, selected_nodes_2, by="name")
  
}

### Find overlapping between fold change and hotnet2 results ###
find_overlapping_fc_hn2 <- function(data1, data2, noded1, 
                                    edgecut2, nodep1, size2){
  
  deExpr <- data1$meta_data
  
  plot_n <- data1$nodes %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  selected_nodes <- plot_n %>% 
    filter(abs(FC) >= noded1 & pvalue <= nodep1) %>%
    select(name)
  
  hn2_e <- data2$hn2_e
  deExpr <- data2$meta_data
  
  plot_n <- data2$nodes  %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  filtered_hn2_e <- hn2_e %>% filter(similarity >= edgecut2)    
  if (nrow(filtered_hn2_e) <= 0) hn2_nodes_2 <- NULL
  else{    
    plot_n <- find_hn2_nodes(plot_n, filtered_hn2_e, size2)                    
    hn2_nodes <- plot_n  %>% filter(group > 0) %>% select(name)
  }      
  
  inner_join(selected_nodes, hn2_nodes, by="name")
  
}

network_filter <- function(source_nodes, steps, plot_n, plot_e){
  
  if (steps == "1"){
    # node selection - 1 step further
    nodes_1step <- find_1step_nodes(source_nodes, plot_n, plot_e)    
    plot_n <- nodes_1step$outn    
    plot_e <- nodes_1step$oute
  }
  else if (steps == "2"){      
    nodes_2step <- find_2step_nodes(source_nodes, plot_n, plot_e)
    plot_n <- nodes_2step$outn    
    plot_e <- nodes_2step$oute
  }
  else if (steps == "0"){
    nodes_2step <- find_noexpand_nodes(source_nodes, plot_n, plot_e)
    plot_n <- nodes_2step$outn    
    plot_e <- nodes_2step$oute
  }
  else {  # show all
    rownames(plot_n) <- seq_len(nrow(plot_n))
    # re-create index based on nodes order
    plot_e$target <- match(plot_e$target_g, plot_n$name)-1
    plot_e$source <- match(plot_e$source_g, plot_n$name)-1
    
  }
  return(list(n=plot_n, e=plot_e))
}


find_1step_nodes <- function(source_nodes, nodes, edges){
  
  out_edges <- edges %>% 
    filter( target_g %in% source_nodes$name | 
              source_g %in% source_nodes$name)
  
  out_nodes <- full_join( (out_edges %>% 
                             select(target_g) %>% 
                             distinct() %>% 
                             mutate(name=target_g) ), 
                          (out_edges %>%
                             select(source_g) %>% 
                             distinct() %>% 
                             mutate(name=source_g) ),
                          by="name") %>% select(name)
  
  nodes <- nodes[nodes$name %in% out_nodes$name,]   
  rownames(nodes) <- seq_len(nrow(nodes))  
  # re-create index based on nodes order
  out_edges$target <- match(out_edges$target_g, nodes$name)-1
  out_edges$source <- match(out_edges$source_g, nodes$name)-1  
  
  
  return(list(oute=out_edges, outn=nodes))
}

find_noexpand_nodes <- function(source_nodes, nodes, edges){
  
  out_edges <- edges %>% 
    filter( target_g %in% source_nodes$name & source_g %in%source_nodes$name)
  
  nodes <- nodes[nodes$name %in% source_nodes$name,]   
  rownames(nodes) <- seq_len(nrow(nodes))  
  # create index based on nodes order
  out_edges$target <- match(out_edges$target_g, nodes$name)-1
  out_edges$source <- match(out_edges$source_g, nodes$name)-1  
  
  return(list(oute=out_edges, outn=nodes))
}


find_2step_nodes <- function(source_nodes, nodes, edges){
  
  out_edges <- edges %>% 
    filter( target_g %in% source_nodes$name | source_g %in% source_nodes$name)
  
  out_nodes <- full_join(out_edges %>% 
                           select(target_g) %>% 
                           distinct() %>% 
                           mutate(name=target_g), 
                         out_edges%>% 
                           select(source_g) %>% 
                           distinct() %>% 
                           mutate(name=source_g),
                         by="name") %>% select(name)
  
  source_nodes <- full_join(out_nodes, source_nodes) %>% select(name)
  out_edges <- edges %>% 
    filter( target_g %in% source_nodes$name | source_g %in% source_nodes$name)
  out_nodes <- full_join(out_edges %>% 
                           select(target_g) %>% 
                           distinct() %>% 
                           mutate(name=target_g), 
                         out_edges%>%
                           select(source_g) %>% 
                           distinct() %>% 
                           mutate(name=source_g),
                         by="name") %>% select(name)
  
  
  nodes <- nodes[nodes$name %in% out_nodes$name,]   
  rownames(nodes) <- seq_len(nrow(nodes))  
  # create index based on nodes order
  out_edges$target <- match(out_edges$target_g, nodes$name)-1
  out_edges$source <- match(out_edges$source_g, nodes$name)-1  
  
  return(list(oute=out_edges, outn=nodes))
}