

# Create output of nodes and edges table for Network Viz library
select_network_data <- function(useData, edgecut, nodecut, nodecutp, 
                                steps, methodSelection, 
                                overlaygenes, inputG,sizeofsubnetwork){
  
  plot_n <- useData$nodes
  plot_e <- useData$edges
  hn2_e <- useData$hn2_e
  deExpr <- useData$meta_data  
  
  # check nodes in omics data 
  # merge omics with node data 
  # select edges within selected nodes
  
  plot_n <- plot_n %>% 
    rename(name=gene) %>% 
    left_join(deExpr, by="name") %>% 
    filter(!is.na(pvalue)) %>% 
    select(-index)
  
  # extract edges with both nodes in plot_n
  plot_e <- plot_e %>% 
    filter( target_g %in% plot_n$name & 
              source_g %in% plot_n$name )
  
  # Show HotNet2 results 
  if (methodSelection == 'hotnet2'){   
    selectHotNet2(plot_n, plot_e, hn2_e, edgecut, 
                  sizeofsubnetwork, steps, overlaygenes)  
  }   
  else if (methodSelection == 'fc'){  
    # Show Fold change results
    selectFCnodes(plot_n, plot_e, hn2_e, nodecut, nodecutp,
                  steps, overlaygenes)
  }
  else{
    # Show network given a set of genes (use overlaygenes as variable to pass
    # genes.    
    selectGivenNodes(plot_n, plot_e, steps, overlaygenes)
  }
}
  
selectFCnodes <- function(plot_n, plot_e, hn2_e, nodecut, nodecutp,
                          steps, overlaygenes){
  # two groups for FC filter mode. 
  # 1: nodes with FC and pvalue pass cutoffs.
  # 0: others.
  N <- plot_n %>% nrow()
  
  plot_n <- plot_n %>% 
    mutate( group=ifelse(
      abs(FC) >= nodecut & pvalue <= nodecutp, 1, 0))
  
  # retaining gene name and heat for viz
  plot_n <- plot_n %>% select(-pvalue) %>% rename(heat = FC)    
  selected_nodes <- plot_n %>% filter(group == 1)
  #cat(nrow(selected_nodes), "\n")
  
  # traverse from selected_nodes (group: 1)
  results <- network_filter(selected_nodes, steps, plot_n, plot_e)  
  plot_e <- results$e
  plot_n <- results$n
  
  if (steps == "0" & nrow(plot_e) > 0 ){ # separate connected components if # nodes in cc >= 2
    # finding strongly connected components in the network 
    # given edges with weight larger than edgecut      
    grouping <- plot_e %>% 
      select(c(source_g, target_g)) %>%         
      find_cc()
    
    
    grouping$group <- grouping$group+1
    plot_n$group <- grouping[match(plot_n$name, grouping$name),]$group
    # gene set analysis and give name.
    grpName <- gene_set_analysis(plot_n %>% filter(group > 1), N)   
    plot_n <- plot_n %>% left_join(grpName, by="group") 
    # set nodes not in any c.c. 
    plot_n[is.na(plot_n$group),]$group <- 1
    plot_n[is.na(plot_n$grpname),]$grpname <- ""
    plot_n[is.na(plot_n$grppval),]$grppval <- 1
    
  }
  
  # create colors, highlight genes, type of edges for visualization
  coldv <- min(selected_nodes$heat)    
  hotv <- max(selected_nodes$heat)
  
  # node color scale for viz
  if (coldv < 0) node_scale <- paste0("d3.scale.linear()
                                     .domain([", coldv, ",0," , hotv ,"])
                                     .range([", color_blue_white_red , "])
                                     .nice()")    
  else node_scale <- paste0("d3.scale.linear()
                           .domain([", coldv, "," , hotv ,"])
                           .range([" ,color_red ,"]).nice()")    
  
  # create domain range for hull color scale.
  domain_v <- paste(plot_n$group %>% unique() %>% sort() , collapse=",")

  if (length(unique(plot_n$group)) == 1 ) 
    hull_scale <- "d3.scale.ordinal().domain([1]).range(['#FFF'])"
  else#else if(steps == "0")  
    hull_scale <- paste0("d3.scale.ordinal().domain([", domain_v, 
                         "]).range([", color_cat20_plus1 ,"])")
  #else 
  #  hull_scale <- paste0("d3.scale.ordinal().domain([", domain_v, "]).
  #                         range([", color_cat20_plus1, "])")
  
  if (!("type" %in% colnames(plot_e)))  plot_e$type <- 'none'
  new_group_id <- max(plot_n$group) + 1
  plot_n <- plot_n %>% 
    mutate(newGroup = ifelse( !(name %in% selected_nodes$name ) &
                                (name %in% overlaygenes), 
                              new_group_id, group)) %>%
    mutate(group=newGroup) %>% select(-newGroup)
  
  # create new columns for highlighted genes
  plot_n <- plot_n %>% mutate(highlightedGene = 
                                ifelse(name %in% overlaygenes, 1, 0 ))
  # plot all edges with same stroke-width    
  if(nrow(plot_e) > 0) plot_e$stroke <- 1
  else  plot_e$stroke <- as.numeric()
  # fix node size of node in the viz
  plot_n$nodesize <- 0.01       
  if ("grpname" %in% colnames(plot_n))
    plot_n <- plot_n %>% mutate(grpname=ifelse(group==0, 'background', grpname))
  else
    plot_n <- plot_n %>% mutate(grpname=ifelse(group==0, 'background',
                                               '-'))
  
  return(list(nodes=plot_n, edges=plot_e, 
              nodeColorScale=node_scale, hullColorScale=hull_scale))
}


find_hn2_nodes <- function(plot_n, filtered_hn2_e, sizeofsubnetwork){
  # finding strongly connected components in the network 
  # given edges with weight larger than edgecut      
  grouping <- filtered_hn2_e %>% 
    select(c(source, target)) %>%         
    find_strong_cc()
  
  # assign group of hotnet2 to nodes
  plot_n$group <- grouping[match(plot_n$name, grouping$name),]$group              
  # NA, group 0 (background nodes)
  plot_n[is.na(plot_n)] <- 0               
  
  # filter subnetworks by size
  plot_n <- plot_n %>%
    group_by(group) %>%
    mutate(tmpG = ifelse(n()>=sizeofsubnetwork, group, 0) ) %>% 
    ungroup() %>%
    mutate(group=tmpG) %>% select(-tmpG)
  
  plot_n %>% select(-pvalue) %>% rename(heat = FC)                    
  
}


selectHotNet2 <- function(plot_n, plot_e, hn2_e, edgecut, 
                          sizeofsubnetwork, steps, overlaygenes){
  
  # create group of nodes based on delta
  
  filtered_hn2_e <- hn2_e %>% filter(similarity >= edgecut)   
  if (nrow(filtered_hn2_e) <= 0) NULL
  else{    
    plot_n <- find_hn2_nodes(plot_n, filtered_hn2_e, sizeofsubnetwork)                    
    hotnet2_nodes <- plot_n  %>% filter(group > 0)
    
    # gene set analysis and give name.
    grpName <- gene_set_analysis(hotnet2_nodes, plot_n%>%nrow())   
    plot_n <- plot_n %>% left_join(grpName, by="group") 
    
    # create key first (better way to do this?)    
    filtered_hn2_e <- filtered_hn2_e %>% 
      mutate(edge=ifelse(as.character(source)<as.character(target), 
                         paste(source,target,sep="_"), 
                         paste(target,source,sep="_"))) %>%
      group_by(edge) %>% summarise(stroke=max(similarity)) %>% ungroup()
            
    if (nrow(hotnet2_nodes) == 0) return(NULL)
    else{   
      # nodes edges filter      
      results <- network_filter(hotnet2_nodes, steps, plot_n, plot_e)  
      plot_e <- results$e
      plot_n <- results$n
    
      coldv <- min(hotnet2_nodes$heat)    
      hotv <- max(hotnet2_nodes$heat)  
      
      # assign colors
      if (coldv < 0) 
        node_scale <- paste("d3.scale.linear()
                            .domain([", coldv, ",0," , hotv ,"])
                            .range([",color_blue_white_red,"])
                            .nice()", sep="")    
      else node_scale <- paste("d3.scale.linear()
                               .domain([", coldv, "," , hotv ,"])
                               .range([", color_red, "])
                               .nice()", sep="")    
            
      domain_v <- paste(plot_n$group %>% unique() %>% sort() , collapse=",")
      if (length(unique(plot_n$group)) == 1 ) 
        hull_scale <- "d3.scale.ordinal().domain([1]).range(['#1f77b4'])"
      else if(steps == "0") # no background node, color start from blue  
        hull_scale <- paste0("d3.scale.ordinal()
                             .domain([", domain_v ,"])
                             .range([", color_cat20, "])")
      else 
        hull_scale <- paste0("d3.scale.ordinal()
                             .domain([", domain_v ,"])
                             .range([", color_cat20_plus1 , "])")   
      
      
      plot_n$nodesize <- 0.01
    
      # assign edge weight from hotnet2 edges: filtered_hn2_e    
      plot_e <- plot_e %>% 
        mutate(edge=ifelse(
          (paste(source_g, target_g, sep="_") %in% filtered_hn2_e$edge), 
          paste(source_g, target_g, sep="_"), 
          paste(target_g, source_g, sep="_"))) 
      plot_e <- plot_e %>% 
        left_join(filtered_hn2_e, by='edge') %>% 
        select(-edge)
      
      # edge with NA weight
      plot_e[is.na(plot_e)] <- 0.01               
      plot_e$stroke <- plot_e$stroke*100
      plot_n <- plot_n %>% 
        mutate(highlightedGene = ifelse(name %in% overlaygenes, 1,0 ))
      
      # assign new group for overlay genes not in foreground
      new_group_id <- max(plot_n$group) + 1
      plot_n <- plot_n %>% 
        mutate(newGroup = 
                 ifelse( !(name %in% hotnet2_nodes$name )&
                           (name %in% overlaygenes), 
                         new_group_id, group)) %>%
        mutate(group=newGroup) %>% select(-newGroup)
      
      if (!("type" %in% colnames(plot_e))) { plot_e$type <- 'none'}        
      return(list(nodes=plot_n, edges=plot_e, 
                  nodeColorScale=node_scale, hullColorScale=hull_scale))
    }
  }
  
  
  selectGivenNodes <- function(plot_n, plot_e, steps, groupGenes){
    
    plot_n <- plot_n %>% select(-pvalue) %>% rename(heat = FC)    
    plot_n$selectGroup <- 
      groupGenes[match(plot_n$name, groupGenes$gene),]$groupid
    
    plot_n[is.na(plot_n)] <- 0
    
    plot_n$group <- plot_n$selectGroup
    
    selected_nodes <- plot_n  %>% filter(selectGroup > 0)  
    # nodes edges filter
    results <- network_filter(selected_nodes, steps, plot_n, plot_e)  
    plot_e <- results$e
    plot_n <- results$n
    
    
    # create colors, highlight genes, type of edges for visualization
    coldv <- min(selected_nodes$heat)    
    hotv <- max(selected_nodes$heat)
    
    # node color scale for viz
    if (coldv < 0) node_scale <- paste0("d3.scale.linear()
                                        .domain([", coldv, ",0," , hotv ,"])
                                        .range([", color_blue_white_red , "])
                                        .nice()")    
    else node_scale <- paste0("d3.scale.linear()
                              .domain([", coldv, "," , hotv ,"])
                              .range([" ,color_red ,"]).nice()")    
    
    # create domain range for hull color scale.
    domain_v <- paste(plot_n$group %>% unique() %>% sort() , collapse=",")
    
    if (length(unique(plot_n$group)) == 1 ) 
      hull_scale <- "d3.scale.ordinal().domain([1]).range(['#FFF'])"
    else#else if(steps == "0")  
      hull_scale <- paste0("d3.scale.ordinal().domain([", domain_v, 
                           "]).range([", color_cat20_plus1 ,"])")
    
    if (!("type" %in% colnames(plot_e)))  plot_e$type <- 'none'
    
    # plot all edges with same stroke-width    
    if(nrow(plot_e) > 0) plot_e$stroke <- 5
    else  plot_e$stroke <- as.numeric()
    # fix node size of node in the viz
    plot_n$nodesize <- 0.01
    groupGenes <- groupGenes %>% 
      rename(name=gene) %>% 
      select(-group) %>% 
      select(-groupid) %>% rename(grpname=pname)
    plot_n <- plot_n %>% left_join(groupGenes, by="name")
    plot_n[is.na(plot_n)] <- 'background'
    
    return(list(nodes=plot_n, edges=plot_e, 
                nodeColorScale=node_scale, 
                hullColorScale=hull_scale))
  }
  
  
}