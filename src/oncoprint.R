# plot oncoprint
# non-interactive: https://gist.github.com/armish/564a65ab874a770e2c26

plot_oncoprint <- function(data, cometTag, pcut){
  # put gene column as rownames
  rownames(data) <- data[,1]
  # remove gene column
  data <- data[,-1]
  # make data matrix for oncoPrint
  M <- data.matrix(data)
  # create numbers (coverage, n, m) for oncoPrint
  if ( nrow(data) > 1 & cometTag ) cp <- get_CoMEt_pvalue(get_numbers(M), pcut)
  else cp <- pcut
  oncoPrint(get_numbers(M), cp, pcut)
  
}

get_mutated_samples <- function(data){
  # put gene column as rownames
  rownames(data) <- data[,1]
  # remove gene column
  data <- data[,-1]
  # make data matrix for oncoPrint
  M <- data.matrix(data)
  # return a list of samples with mutated and non-mutated status
  memoSort(M) %>% data.frame() %>% summarise_each(funs(max))
}

get_numbers <- function(M){
  alts <- memoSort(M)
  if (!is.vector(alts)){
    ngenes <- nrow(alts);
    nsamples <- ncol(alts);
    coverage <- sum(alts > 0)
  }
  else{
    coverage <- sum(alts);
    ngenes <- 1;
    nsamples <- length(alts);
  }
  return(list(M=M, alts=alts, ngenes=ngenes, nsamples=nsamples, 
              coverage=coverage))
}

# This function sort the matrix to better visualize Mutual Exclusivity
memoSort <- function(M) {
  
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  
  if (nrow(M) == 1){
    scores <- M[1,]
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[1, sampleOrder])
  }
  else{
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
    scores <- apply(M[geneOrder, ], 2, scoreCol);
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[geneOrder, sampleOrder]);
  }
  
}

get_CoMEt_pvalue <- function(nM, pval=0.1){
  M <- nM$M
  alts <- nM$alts
  ngenes <- nM$ngenes
  nsamples <- nM$nsamples
  coverage <- nM$coverage
  CoMEtInput <- rep(0, 2^ngenes)
  for (i in 1:nsamples) {
    index_i <- 0
    for (j in 1:ngenes){
      index_i <- index_i + alts[j,i] * 2^(j-1)
    }
    CoMEtInput[index_i+1] <- CoMEtInput[index_i+1] + 1
  }
  #cat("Comet input:", CoMEtInput, "\n")
  comet_exact_test(CoMEtInput, pval, mutmatplot = F)
}

# This is the plotting function
oncoPrint <- function(nM, cp, pcut) {
  
  # retrive all numbers and matrices
  M <- nM$M
  alts <- nM$alts
  ngenes <- nM$ngenes
  nsamples <- nM$nsamples
  coverage <- nM$coverage
  
  ### OncoPrint
  numOfOncos <- ngenes*nsamples;
  oncoCords <- matrix( rep(0, numOfOncos * 6), nrow=numOfOncos );
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", 
                           "ytop", "altered", "cooccurrence");
  
  # create rows for indicating co-occurrence or not
  if (!is.vector(alts)) {
    tmp_alts <- as.data.frame(alts)
    sumColumnsCooc <- tmp_alts %>% summarise_each(funs(sum))
  }
  else{
    sumColumnsCooc <- alts
  }
  coverage <- sum(sumColumnsCooc > 0)
  
  xpadding <- .01;
  ypadding <- .01;
  cnt <- 1;
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      xleft <- j-1 + xpadding;
      ybottom <- ((ngenes-i+1) -1) + ypadding;
      xright <- j - xpadding;
      ytop <- (ngenes-i+1) -ypadding;
      cooc <- sumColumnsCooc[[j]]
      
      if (!is.vector(alts)) altered <- alts[i, j]
      else altered <- alts[j]
      
      oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered, cooc);
      cnt <- cnt+1;
    }
  }
  
  if (!is.vector(alts)) agene <- rownames(alts)
  else agene <- rownames(M)
  
  tmp_M <- as.data.frame(M)
  alts_freq <- bind_cols(
      (tmp_M %>%
       mutate(freq=rowSums(.)) %>% select(freq)),
      (tmp_M %>% mutate(gene=agene) %>% select(gene) )
      ) %>% mutate(key=paste0(gene,"[",freq,"]")) %>% select(key) %>% unlist()
  
  colors <- rep("lightgray", cnt);
  colors[ which(oncoCords[, "altered"] == 1 & 
                  oncoCords[, "cooccurrence"] == 1 ) ] <- "steelblue1";
  colors[ which(oncoCords[, "altered"] == 1 & 
                  oncoCords[, "cooccurrence"] > 1 ) ] <- "darkorange2";
  plot(c(0, nsamples), c(0, ngenes), type="n", 
       main=ifelse(cp < pcut, 
              sprintf("Gene set altered in %.2f%%: %d of %d cases. 
                       CoMEt pvalue: %s.", coverage/nsamples*100, coverage, 
                      nsamples, format(cp, scientific = T, digits = 3)), 
              sprintf("Gene set altered in %.2f%%: %d of %d cases."
                      , coverage/nsamples*100, coverage, nsamples)), 
       xlab="Samples", ylab="", yaxt="n");
  rect(oncoCords[, "xleft"], 
       oncoCords[, "ybottom"],
       oncoCords[, "xright"], 
       oncoCords[, "ytop"], 
       col=colors, border="white");
  axis(2, at=(ngenes:1)-.5, labels=alts_freq, las=3, cex.axis = 0.5);
}

