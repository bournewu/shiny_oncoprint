

read_expr_data <- function(mut_data, expr_data){
  # RAW[[mut_data]]
  rawDataAll[[mut_data]]
}

read_mut_data <- function(mut_data){
  mutDataAll[[mut_data]]
}

get_intersect_samples <- function(d1, d2){
  intersect(colnames(d1), colnames(d2))
}

expr_data <- function(data){

  # for multiple hits, chose the one with smallest pvalue
  data <- data %>% 
    rename(name = gene) %>% 
    rename(FC=logFC) %>% 
    rename(pvalue=adj.P.Val) %>%
    select(c(name, pvalue, FC)) 
  
    data_with_key <- data %>% 
    mutate(key=paste(name, pvalue, sep="_"))
  
  key <- data %>% 
    group_by(name) %>% 
    summarise(pvalue=min(pvalue) ) %>% 
    unite(key, name, pvalue)
  
  key %>% 
    left_join(data_with_key, by="key") %>%
    select(-key)
  
}