pval.hist<-function(p)
  
{
  m<-length(p)                            # count the p-values 
  h.cdf<-h.pdf<-p.edf<-rep(NA,length(p))  # initialize h.pdf and p.edf vectors
  ord<-order(p)                           # order the p-values
  p.edf[ord]<-(1:m)/(m+1)                 # define the EDF estimator
  hh<-NULL                                # initialize the histogram heights 
  F0<-p0<-b.edf<-p.brks<-1                # initialize vectors with F0, p0, EDF break-points, and p-value break-points
  
  while(any(p<p0))  # Loop until all p-values are included in the histogram
  {
    slp<-min(1,(2*sum(p[p<p0])+1)/(1+p0*m))  # deterime the average abundance of p-values (under conditional uniform)
    slp.line<-(1-slp*(p0-p)/p0)*F0           # represent this line as a straight line in the EDF  
    ok<-(p<p0)&(p.edf>=slp.line)             # find the p-values that are "OK"
    if (any(ok))                             # Insert a new histogram break-point
    {
      p.cut<-min(p[ok])                   # Find the minimum p-value that is OK
      p.brks<-c(p.cut,p.brks)             # Update the p-value break-point vector
      hh<-c(slp,hh)                       # Update the histogram heights vector 
      p0<-p.cut                           # Update p0
      F0<-min(p.edf[p>=p0])               # Update F0
      b.edf<-c(F0,b.edf)                  # Update the break-points in the EDF histogram estimator
    }
    else # inserts a (p,EDF) point at (0,0)
    {
      p.brks<-c(0,min(p),max(p[p<p0]),p.brks) 
      b.edf<-c(0,min(p.edf[p==min(p)]),min(p.edf[p==max(p[p<p0])]),b.edf)
      p0<-0
    }
  }
  
  last.pt<-min((1:length(p.brks))[p.brks==1])
  p.brks<-p.brks[1:last.pt]
  b.edf<-b.edf[1:last.pt]
  
  p.diffs<-diff(p.brks)
  zero.diffs<-(1:length(p.diffs))[p.diffs==0]
  if (length(zero.diffs)>0)
  {
    p.brks<-p.brks[-(zero.diffs+1)]
    b.edf<-b.edf[-(zero.diffs+1)]   
  }
  
  hh<-exp(log(diff(b.edf))-log(diff(p.brks)))  
  if (any(is.na(hh)))
  {
    print(summary(p))
    print(p.brks)
    print(b.edf)
    print(hh)
    stop("Error: undefined histogram heights.")
    
  } 
  hd<-diff(hh)  # Difference in histogram heights
  p1<-1
  while(any(hd>0))  # Enforce monotonicity
  {
    k<-max((1:length(hd))[hd>0])                        # Find largest histogram bar index k such that bar k is shorter than bar k+1
    p.brks<-c(p.brks[1:k],p.brks[(k+2):length(p.brks)]) # combine histogram bars k and k+1
    b.edf<-c(b.edf[1:k],b.edf[(k+2):length(b.edf)])     # update EDF break-points with monotonocity constraint      
    hh<-exp(log(diff(b.edf))-log(diff(p.brks)))         # compute new histogram heights
    hd<-diff(hh)                                        # compare consecutive histogram heights
  }
  
  h.cdf<-approx(p.brks,b.edf,xout=p)$y                  # Get the histogram-based CDF estimates for each p-value
  p.hist<-exp(log(diff(b.edf))-log(diff(p.brks)))       # Get the height for each histogram bar
  pi0.hat<-min(p.hist)                                  # Get the pi0 estimate from the histogram
  h.ebp<-approx(p.brks,pi0.hat/c(p.hist,p.hist[length(p.hist)]),xout=p)$y # Get the conservative EBP interpolation
  h.fdr<-exp(log(min(pi0.hat))+log(p)-log(h.cdf))                         # Get the histogram-based FDR estimate
  h.ebp[p==0]=0
  h.fdr[p==0]=0
  
  return(list(p=p,              # input p-value vector
              edf=p.edf,        # the p-value edf
              h.cdf=h.cdf,      # the histogram-based cdf estimate
              h.fdr=h.fdr,      # the histogram-based FDR estimate
              h.ebp=h.ebp,      # the histogram-based EBP estimate
              p.brks=p.brks,    # the p-value break-points of the histogram
              p.hist=p.hist,    # the heights of each histogram bar
              edf.brks=b.edf,   # the break-points in the EDF of the histogram estimator
              pi0.hat=pi0.hat)) # the histogram-based estimate of the proportion of tests with a true null hypothesis
}

SortData <- function(counts, treatment, replic = NULL, sort.method, offset = NULL){
  #Trims and sorts matrix of counts into appropriate format to be used
  #in data.sim
  
  #Check inputs
  counts <- as.matrix(counts)
  n.col <- dim(counts)[2]
  
  if (!is.null(replic)){
    if(length(replic) != n.col) 
      stop("Error: Length of replic vector must equal number of 
           columns in counts matrix")
    if (any(tabulate(replic) > 2)) 
      stop("Error: Number of observations per replicate in counts 
           matrix must not be greater than 2.")
  }
  
  if(length(treatment) != n.col)
    stop("Error: Length of treatment vector must equal number of 
           columns in counts matrix")
  if (length(unique(treatment)) != 2) 
    stop("Error: Number of treatment groups in counts matrix 
         must be equal to two.")
  if (!sort.method %in% c("paired", "unpaired")) 
    stop("Error: sort.method must be set to either 'paired' or 'unpaired'.")
  if (sort.method == "paired" && is.null(replic)) 
    stop("Error: Must specify replic vector when sort.method equals 'paired'.")
  
  if(!is.null(offset)){
    if (!is.numeric(offset)) 
      stop("Error: offset must be a positive numeric vector with 
           length equal to the number of columns in the counts matrix.")
    if (is.numeric(offset)){
      if (any(offset <= 0) || length(offset) != n.col) 
        stop("Error: offset must be a positive numeric vector with 
             length equal to the number of columns in the counts matrix.")
    }
  }
  
  if (sort.method == "paired"){ 
    #sort replic by replic then by trt
    sorting <- order(replic, treatment, decreasing = FALSE)
    replic <- replic[sorting]
    
    #remove replicates in inputs that are non-paired and apply sorting
    #to all inputs
    sorting <- sorting[replic %in% unique(replic)[tabulate(replic) == 2]]
    
    counts <- counts[, sorting, drop = FALSE]
    replic <- factor(replic[sorting])
    treatment <- factor(treatment[sorting])
    offset <- offset[sorting]
  }
  
  if (sort.method == "unpaired"){
    #sort inputs by treatment
    sorting <- order(treatment)
    replic <- replic[sorting]
    treatment <- treatment[sorting]
    counts <- counts[, sorting, drop = FALSE]
    offset <- offset[sorting]
    
    #no need to remove any columns (replicates)
  }
  
  return(list(counts = counts, replic = replic, treatment = treatment, 
              offset = offset, sorting = sorting))
  }

SimData <- function(counts, treatment, replic = NULL, sort.method, 
                     k.ind, n.genes = NULL, n.diff = NULL, 
                     n1.genes = NULL, n1.diff = NULL,
                     offset = NULL, samp.independent = TRUE
                     genes.select = NULL, genes.diff = NULL, 
                     genes1.select = NULL, genes1.diff = NULL, 
                     probs = NULL, EBP = NULL, exact = FALSE, power = 1){
  #given matrix of gene expression counts with two treatment groups, 
  #simulates matrix of gene expression data with known list of DE
  #and EE genes. See help file for more info.
  
  #check inputs
  counts <- as.matrix(counts)
  if (any(!is.numeric(counts))) 
    stop("Error: counts matrix contains non-numeric values.")
  if (any(round(counts) != counts)) 
    stop("Error: counts matrix contains non-integer values.")

  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  
  if (is.null(n.genes) && is.null(genes.select)) 
    stop("Error: Both n.genes and genes.select are set to NULL.")
  
  if (!is.null(genes.select)){
    if (!is.numeric(genes.select) && !is.logical(genes.select)) 
      stop("Error: genes.select must be a NULL, numeric or logical vector.")
    if (is.logical(genes.select)){
      if (length(genes.select) != n.row) 
        stop("Error: When genes.select is a logical vector, 
             its length must equal the number of rows of the counts matrix.") 
      else {
        genes.select <- which(genes.select)
      }
    }
    if (!is.null(n.genes)){
      if (length(genes.select) != n.genes) 
        stop("Error: n.genes must equal number of genes selected in 
                 simulation through genes.select vector.")
    }
    if (is.null(n.genes)) n.genes <- length(genes.select) 
  }
  
  if (is.null(genes.select) && !is.null(genes.diff)) 
    stop("Error: genes.diff vector cannot be specified without 
         first specifying genes.select vector.")
  
  if (is.null(n.diff) && is.null(genes.diff)) 
    stop("Error: Both n.diff and genes.diff are set to NULL.")
  
  if (!is.null(genes.select) && !is.null(genes.diff)){
    if (!is.numeric(genes.diff) && !is.logical(genes.diff)) 
      stop("Error: genes.diff must be a NULL, numeric or logical vector.")
    if (is.numeric(genes.diff)){
      if (any(!genes.diff %in% genes.subset)) 
        stop("Error: genes.diff vector must be a subset of genes.select.")
    } 
    if (is.logical(genes.diff)){
      if (length(genes.diff) != length(genes.select)) 
        stop("Error: When genes.diff is a logical vector, length of genes.diff must equal 
             number of genes selected in genes.select vector.")
      else {
        genes.diff <- genes.select[genes.diff]
      }
    }
    if (!is.null(n.diff)){
      if (length(genes.diff) != n.diff) 
        stop("Error: n.diff must equal number of genes to be differentially expressed
              as selected through genes.diff vector.")
    }
    if (is.null(n.diff)) n.diff <- length(genes.diff)
  }
  
  if (!is.numeric(n.genes)) 
    stop("Error: Number of genes selected, n.genes, must be a positive integer 
          less than the total number of rows in the counts matrix.")
  else if (round(n.genes) != n.genes || n.genes > n.row || n.genes <= 0) 
    stop("Error: Number of genes selected, n.genes, must be a positive integer 
         less than the total number of rows in the counts matrix.")

  if (!is.numeric(n.diff)) 
    stop("Error: Number of DE genes, n.diff, must be a positive integer 
          less than n.genes.")
  else if (round(n.diff) != n.diff || n.diff > n.genes || n.diff <= 0) 
    stop("Error: Number of DE genes, n.diff, must be a positive integer 
          less than n.genes.")
  
  if (!is.numeric(k.ind)) 
    stop("Error: Number of replicates simulated per treatment group, k.ind, 
          must be a positive integer less than total the number of 
          columns in the counts matrix divided by two.")
  else if (round(k.ind) != k.ind || k.ind > n.col/2 || k.ind <= 0) 
    stop("Error: Number of replicates in each treatment group, k.ind, 
          must be a positive integerless than the total number of 
          columns in counts matrix divided by two.")

  if (!is.null(probs)){
    if (!is.numeric(probs))
      stop("Error: probs must be a numeric vector of length equal to the 
           number of rows in the counts matrix with each entry between 0 and 1.") 
    if (any(probs > 1) || any(probs < 0) || length(probs) != n.row) 
      stop("Error: probs must be a numeric vector of length equal to the 
           number of rows in the counts matrix with each entry between 0 and 1.") 
  }
  
  if (!is.null(EBP)){
    if (!is.numeric(EBP))
      stop("Error: EBP must be a numeric vector of length equal to the 
           number of rows of the counts matrix with each entry between 0 and 1.")
    if (any(EBP > 1) || any(EBP < 0) || length(EBP) != n.row) 
      stop("Error: EBP must be a numeric vector of length equal to the 
           number of rows of the counts matrix with each entry between 0 and 1.") 
  }
  
  if (!is.null(replic)){
    if(length(replic) != n.col) 
      stop("Error: Length of replic vector must equal number of 
           columns in counts matrix")
    if (any(tabulate(replic) > 2)) 
      stop("Error: Number of observations per replicate in counts 
           matrix must not be greater than 2.")
  }
  
  if(length(treatment) != n.col)
    stop("Error: Length of treatment vector must equal number of 
         columns in counts matrix")
  if (length(unique(treatment)) != 2) 
    stop("Error: Number of treatment groups in counts matrix 
         must be equal to two.")
  if (!sort.method %in% c("paired", "unpaired")) 
    stop("Error: sort.method must be set to either 'paired' or 'unpaired'.")
  if (sort.method == "paired" && is.null(replic)) 
    stop("Error: Must specify replic vector when sort.method equals 'paired'.")
  
  if(!is.null(offset)){
    if (!is.numeric(offset)) 
      stop("Error: offset must be a positive numeric vector with 
           length equal to the number of columns in the counts matrix.")
    if (is.numeric(offset)){
      if (any(offset <= 0) || length(offset) != n.col) 
        stop("Error: offset must be a positive numeric vector with 
             length equal to the number of columns in the counts matrix.")
      }
    }
  
  if (is.null(offset)) offset <- rep(1, n.col) 

  if (!is.null(genes1.select)){
    if (is.null(genes.select) || is.null(genes.diff))
      stop("Error: genes.diff vector cannot be specified without 
            first specifying genes.select and genes.diff vectors.")
  }
  
  if (is.null(genes1.select) && is.null(n1.genes)) 
    stop("Error: Both n1.diff and genes1.select are set to NULL.")

  if (!is.null(genes1.select)){
    if (!is.numeric(genes1.select) && !is.logical(genes1.select)) 
      stop("Error: genes1.select must be a NULL, numeric or logical vector.")
    if (is.numeric(genes1.select)){
      if (any(!genes1.select %in% genes.select)) 
        stop("Error: genes1.select must be a subset of genes.select.") 
    }
    if (is.logical(genes1.select)){
      if (length(genes1.select) != length(genes.select)) 
        stop("Error: If specifying genes1.select as a logical vector, the length of 
             genes1.select must equal the number of genes specified in the 
             genes.select vector.")
      else{
        genes1.select <- genes.select[genes1.select]
      }
    }
    if (!is.null(n1.genes)){
      if (length(genes1.select) != n1.genes) 
        stop("Error: n1.genes must equal number of genes selected in 
                 simulation through genes1.select vector.")
    }
    if (is.null(n1.genes)) n1.genes <- length(genes1.select)
  }
  
  if (!is.null(genes1.diff)){
    if (is.null(genes.select) || is.null(genes.diff) || is.null(genes1.select))
      stop("Error: genes1.diff vector cannot be specified without 
           first specifying genes.select, genes.diff and 
           genes1.select vectors.")
  }
  
  if (is.null(genes1.diff) && is.null(n1.diff)) 
    stop("Error: Both n1.diff and genes1.diff are set to NULL.")
  
  if (!is.null(genes1.diff)){
    if (!is.numeric(genes1.diff) && !is.logical(genes1.diff)) 
      stop("Error: genes1.diff must be a NULL, numeric or logical vector.")
    if (is.numeric(genes1.diff)){
      if (any(!genes1.diff %in% genes.diff) || any(!genes1.diff %in% genes1.select)) 
        stop("Error: genes1.diff vector must be a subset of both 
             genes.diff and genes1.select vectors.")
    }
    if (is.logical(genes1.diff)){
      if (length(genes1.diff) != genes.diff) 
        stop("Error: If specifying genes1.select as a logical vector, 
             the length of genes1.diff must equal the number of genes
             selected in genes.diff vector.")
      else{ 
        genes1.diff <- genes.diff[genes1.diff] 
      }
    }
    if (!is.null(n1.diff)){
      if (length(genes1.diff) != n1.diff) 
        stop("Error: n1.diff must equal number of genes selected in 
                 simulation through genes1.diff vector.")
    }
    if (is.null(n1.diff)) n1.diff <- length(genes1.diff)
  }
  
  if (is.null(n1.genes)) n1.genes <- n.genes
  if (is.null(n1.diff)) n1.diff <- n.diff  
  
  if (!is.numeric(n1.genes)) 
    stop("Error: n1.genes must be a positive integer less than n.genes.")
  else if (round(n1.genes) != n1.genes || n1.genes > n.genes || n1.genes < 0) 
    stop("Error: n1.genes must be a positive integer less than n.genes.")
  
  if (!is.numeric(n1.diff)) 
    stop("Error: n1.diff must be a positive integer less than the 
         minimum of n.diff and n1.genes.")
  else if (round(n1.diff) != n1.diff || n1.diff > n.diff || 
             n1.diff > n1.genes || n1.diff < 0) 
    stop("Error: n1.diff must be a positive integer less than the 
         minimum of n.diff and n1.genes.")

  if (is.numeric(n1.diff)){
    if (n.diff - n1.diff > n.genes - n1.genes) 
      stop("Error:n.diff - n1.diff, must be a positive integer less than the total
           number of genes available in the second group, n.genes - n1.genes.")
  }
  
  
  #sort necessary inputs
  sort.list <- SortData(counts, treatment, replic, sort.method, offset)
  counts <- sort.list[[1]]
  treatment <- sort.list[[3]]
  offset <- sort.list[[4]]
  
  #reset n.row and n.col variables since SortData function may have trimmed
  #counts matrix
  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  
  #calculate p-value of differential depression for each gene
  #if sort.method == paired, use signed rank test
  #if sort.method == unpaired, use ranksum test
  if (is.null(probs) && is.null(EBP) && is.null(genes.diff)){
    if (sort.method == "paired"){
      odds <- seq(from = 1, to = n.col, by = 2)
      diff1 <- t(counts[, odds, drop = FALSE] + 1)/offset[odds] 
      diff2 <- t(counts[, (odds + 1), drop = FALSE] + 1)/offset[(odds + 1)]
      diff <- log(t(diff1/diff2))
      probs <- apply(diff, 1 , function(x) wilcox.test(x, exact = exact)$p.value)        
    }
    if (sort.method == "unpaired"){
      seq.trt1 <- 1:table(treatment)[1]
      trt1 <- log(t(t(counts[, seq.trt1, drop = FALSE] + 1)/offset[seq.trt1]))
      trt2 <- log(t(t(counts[, -seq.trt1, drop = FALSE] + 1)/offset[-seq.trt1]))
      probs <- numeric(n.row)
      for (jj in 1:n.row) probs[jj] <-
        wilcox.test(trt1[jj, ], trt2[jj, ], exact = exact)$p.value
    }
  }  
  
  #calulate empirical bayes probability of differential expression for each gene
  if (is.null(EBP)){
    if (is.null(genes.select) | is.null(genes.diff)) EBP <- 1 - pval.hist(probs)$h.ebp
  } 
  
  SampGenes <- function(genes.select, genes.diff, genes1.select, genes1.diff, 
                         n.genes, n.diff, n1.genes, n1.diff, 
                         EBP, power, exact){
    #sample genes to be used in simulation and sample which genes are
    #to be differentially expressed
    
    genes <- 1:n.row
    if (is.null(genes.diff)){
      if(is.null(genes.select)){
        genes.diff <- 
          sample(genes, n.diff, prob = (EBP)^power, replace = F)
      } else {
        genes.diff <- 
          sample(genes.select, n.diff, prob = (EBP[genes.select])^power, replace = F)
      }
    } 
    if (is.null(genes1.diff)){
      if (n1.diff > 0){
        genes.diff1 <- genes.diff[1:n1.diff]
        genes.diff2 <- genes.diff[-(1:n1.diff)]
      } else{
        genes.diff1 <- NULL
        genes.diff2 <- genes.diff
      }
    } else{
      genes.diff1 <- genes1.diff
      genes.diff2 <- genes.diff[!genes.diff %in% genes1.diff]
    }
    
    if (is.null(genes.select)) genes.subset <- 
      c(sample(genes[-genes.diff], n.genes - n.diff, replace = F), genes.diff)
    else genes.subset <- genes.select
    if (is.null(genes1.select)){
      if (n1.genes > 0){
        if (n1.genes - n1.diff > 0){
          genes.subset1 <- 
            c(genes.subset[!genes.subset %in% genes.diff][1:(n1.genes - n1.diff)], genes.diff1)
          genes.subset2 <- 
            c(genes.subset[!genes.subset %in% genes.diff][-(1:(n1.genes - n1.diff))], genes.diff2)
        } else{
          genes.subset1 <- genes.diff1
          genes.subset2 <- genes.subset[!genes.subset %in% genes.diff1]
        }
      } else{
        genes.subset1 <- NULL
        genes.subset2 <- genes.subset
      }
    } else{
      genes.subset1 <- genes1.select
      genes.subset2 <- genes.select[!genes.select %in% genes1.select]
    }
    genes.subset1 <- sort(genes.subset1)
    genes.subset2 <- sort(genes.subset2)
    genes.subset <- c(genes.subset1,genes.subset2)
    genes.diff1 <- sort(genes.diff1)
    genes.diff2 <- sort(genes.diff2)
    genes.diff <- c(genes.diff1,genes.diff2)
    DE.genes1 <- genes.subset1 %in% genes.diff1
    DE.genes2 <- genes.subset2 %in% genes.diff2
    DE.genes <- c(DE.genes1, DE.genes2)
    
    return(list(genes.subset = genes.subset, genes.subset1 = genes.subset1, 
                genes.subset2 = genes.subset2, genes.diff = genes.diff, 
                genes.diff1 = genes.diff1, genes.diff2 = genes.diff2, 
                DE.genes = DE.genes, DE.genes1 = DE.genes1, DE.genes2 = DE.genes2))
  }
  
  #sample EE and DE genes
  samp.genes.list <- SampGenes(genes.select, genes.diff, genes1.select, genes1.diff, 
                                n.genes, n.diff, n1.genes, n1.diff,
                                EBP, power, exact)
  genes.subset <- samp.genes.list[[1]]
  genes.subset1 <- samp.genes.list[[2]]
  genes.subset2 <- samp.genes.list[[3]]
  genes.diff <- samp.genes.list[[4]]
  genes.diff1 <- samp.genes.list[[5]]
  genes.diff2 <- samp.genes.list[[6]]
  DE.genes <- samp.genes.list[[7]]
  DE.genes1 <- samp.genes.list[[8]]
  DE.genes2 <- samp.genes.list[[9]]
  
  
  SampCol <- function(n.col, k.ind, sort.method, treatment){
    #sample columns (replicates) to be used
    
    if (sort.method == "paired"){
      evens <- seq(2, n.col, by = 2)
      samp <- sample(evens, 2*k.ind, replace = F)
      samp <- c(samp, samp - 1)
    }
    if (sort.method == "unpaired"){  
      trt1 <- 1:table(treatment)[1]
      trt2 <- (table(treatment)[1] + 1):length(treatment)
      samp <- c(sample(trt1, 2*k.ind), sample(trt2, 2*k.ind))
    }
    samp <- as.numeric(samp)
    return(samp)
  }
  
  #swap in appropriate data to create matrix with known DE and EE genes
  #multiply by appropriate offsets and round
  if (samp.independent == FALSE){
    samp <- SampCol(n.col, k.ind, sort.method, treatment)
    offset.col <- offset[samp]
    
    #perform swaps for group 1
    data1 <- counts[genes.subset1, samp, drop = FALSE]
    data.table1 <- data1[, (2*k.ind + 1):(4*k.ind), drop = FALSE]
    data.table1[DE.genes1, 1:k.ind] <- 
      round(t(t(data1[DE.genes1, 1:k.ind])/offset.col[1:k.ind]*offset.col[(2*k.ind + 1):(3*k.ind)]))
    
    #perform swaps for group 2
    data2 <- counts[genes.subset2, samp, drop = FALSE]
    data.table2 <- data2[, 1:(2*k.ind), drop = FALSE]
    data.table2[DE.genes2, 1:k.ind] <- 
      round(t(t(data2[DE.genes2, (2*k.ind + 1):(3*k.ind)])/offset.col[(2*k.ind + 1):(3*k.ind)]*
                offset.col[1:k.ind]))
    
    #combine matrices from both groups
    data.table <- rbind(data.table1, data.table2)
  }
  
  #swap in appropriate data to create matrix with known DE and EE genes
  #multiply by appropriate offsets and round
  #under independent method, sample columns separately and independently
  #for each gene
  if (samp.independent == TRUE){
    samp <- t(replicate(n.genes, SampCol(n.col, k.ind, sort.method, treatment)))
    offset.col <- apply(samp, 2, function(x) offset[x])
    
    #perform swaps for group 1
    if (n1.genes > 0){
      samp1 <- samp[1:n1.genes, , drop = FALSE]
      offset.col1 <- offset.col[DE.genes1, , drop = FALSE]
      data.temp1 <- counts[genes.subset1, , drop = FALSE]
      data1 <- matrix(mapply(function(x, y) data.temp1[x, y], x = 1:n1.genes, 
                             y = samp1, SIMPLIFY = TRUE), ncol = 4*k.ind)
      data.table1 <- data1[, (2*k.ind + 1):(4*k.ind), drop = FALSE]
      if (n1.diff > 0){
        if (n1.diff == 1){
          data.table1[1:k.ind] <- 
            round(data1[1:k.ind]/offset.col1[1:k.ind]*offset.col1[(2*k.ind + 1):(3*k.ind)])
        } else{
            for(jj in 1:n1.diff){
              data.table1[DE.genes1, 1:k.ind][jj, ] <- 
                round(t(t(data1[DE.genes1, 1:k.ind][jj, ])/offset.col1[jj, 1:k.ind]*
                          offset.col1[jj, (2*k.ind + 1):(3*k.ind)]))
          }
        }
      }
    }
    if (n1.genes == 0) data.table1 <- NULL
    
    #perform swaps for group 2
    if (n1.genes < n.genes){
      if (n1.genes > 0) samp2 <- samp[-(1:n1.genes), , drop = FALSE]
      if (n1.genes == 0) samp2 <- samp
      offset.col2 <- offset.col[DE.genes2, , drop = FALSE]
      data.temp2 <- counts[genes.subset2, , drop = FALSE]
      data2 <- 
        matrix(mapply(function(x, y) data.temp2[x, y], x = 1:(n.genes - n1.genes), 
                      y = samp2, SIMPLIFY = TRUE), ncol = 4*k.ind)
      data.table2 <- data2[, 1:(2*k.ind), drop = FALSE]
      if (n.diff - n1.diff > 0){
        if (n.diff - n1.diff == 1){
          data.table2[1:k.ind] <- 
            round(data2[(2*k.ind + 1):(3*k.ind)]/offset.col2[(2*k.ind + 1):(3*k.ind)]*
                    offset.col2[1:k.ind])
        } else{
            for(jj in 1:(n.diff - n1.diff)){
              data.table2[DE.genes2, 1:k.ind][jj, ] <- 
                round(t(t(data2[DE.genes2, (2*k.ind + 1):(3*k.ind)][jj, ])/offset.col2[jj, (2*k.ind + 1):(3*k.ind)]*
                          offset.col2[jj, 1:k.ind]))
          }
        }
      }
    }
    if (n1.genes == n.genes) data.table2 <- NULL
    
    #combine matrices from both groups
    data.table <- rbind(data.table1, data.table2)
  }
  return(list(counts = data.table, genes.subset = genes.subset, 
              genes.subset1 = genes.subset1, genes.subset2 = genes.subset2, 
              DE.genes = genes.diff, DE.genes1 = genes.diff1, 
              DE.genes2 = genes.diff2))
  }
