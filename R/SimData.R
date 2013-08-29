SimData <-
function(counts, treatment, replic = NULL, sort.method, 
                     k.ind, n.genes = NULL, n.diff = NULL, 
                     n1.genes = NULL, n1.diff = NULL,
                     offset = NULL, samp.independent = FALSE,
                     genes.select = NULL, genes.diff = NULL, 
                     genes1.select = NULL, 
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
  
  genes1.diff <- NULL
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
      if (any(!genes.diff %in% genes.select)) 
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
  
  if (is.null(offset)) offset <- calcNormFactors(counts, method = "TMM") 

  if (!is.null(genes1.select)){
    if (is.null(genes.select) || is.null(genes.diff))
      stop("Error: genes.diff vector cannot be specified without 
            first specifying genes.select and genes.diff vectors.")
  }

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
    
    if (!is.null(n1.diff)){
      if (sum(genes.diff %in% genes1.select) != n1.diff) 
        stop("Error: n1.diff must equal number of genes selected in 
                 simulation through genes1.select vector.")
    }
    if (is.null(n1.diff)) n1.diff <- sum(genes.diff %in% genes1.select)
    genes1.diff <- genes.diff[genes.diff %in% genes1.select]
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
  if (is.null(probs) && is.null(EBP) && is.null(genes.diff)){
    probs <- CalcPvalWilcox(counts, treatment, sort.method, sorted = TRUE, offset)
  }  
  
  #calulate empirical bayes probability of differential expression for each gene
  if (is.null(EBP)){
    if (is.null(genes.select) | is.null(genes.diff)) EBP <- 1 - pval.hist(probs)$h.ebp
  } 
  
  SampGenes <- function(genes.select, genes.diff, genes1.select, genes1.diff = NULL,
                         n.genes, n.diff, n1.genes, n1.diff, 
                         EBP, power, exact){
    #sample genes to be used in simulation and sample which genes are
    #to be differentially expressed
    
    genes <- 1:n.row
    if (is.null(genes.diff)){
      if(is.null(genes.select)){
        genes.diff <- 
          sample(genes, n.diff, prob = (EBP)^power, replace = FALSE)
      } else {
        genes.diff <- 
          sample(genes.select, n.diff, prob = (EBP[genes.select])^power, replace = FALSE)
      }
    } 
    if (is.null(genes1.diff)){
      if (n1.diff > 0){
        genes.diff1 <- sample(genes.diff, n1.diff, replace = FALSE)
        genes.diff2 <- genes.diff[!genes.diff %in% genes.diff1]
      } else{
        genes.diff1 <- NULL
        genes.diff2 <- genes.diff
      }
    } else{
      genes.diff1 <- genes1.diff
      genes.diff2 <- genes.diff[!genes.diff %in% genes1.diff]
    }
    
    if (is.null(genes.select)) genes.subset <- 
      c(sample(genes[-genes.diff], n.genes - n.diff, replace = FALSE), genes.diff)
    else genes.subset <- genes.select
    if (is.null(genes1.select)){
      if (n1.genes > 0){
        if (n1.genes - n1.diff > 0){
          genes.subset1 <- c(
            sample(genes.subset[!genes.subset %in% genes.diff], n1.genes - n1.diff, replace = FALSE), genes.diff1)
          genes.subset2 <- genes.subset[!genes.subset %in% genes.subset1]
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
    genes.subset <- sort(c(genes.subset1,genes.subset2))
    genes.diff1 <- sort(genes.diff1)
    genes.diff2 <- sort(genes.diff2)
    genes.diff <- sort(c(genes.diff1,genes.diff2))
    DE.genes1 <- genes.subset1 %in% genes.diff1
    DE.genes2 <- genes.subset2 %in% genes.diff2
    DE.genes <- genes.subset %in% genes.diff
    
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
    
    #order by genes.subset
    data.table <- data.table[order(c(genes.subset1, genes.subset2)), ]
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
      offset.col1 <- offset.col[1:n1.genes, , drop = FALSE]
      data.temp1 <- counts[genes.subset1, , drop = FALSE]
      data1 <- matrix(mapply(function(x, y) data.temp1[x, y], x = 1:n1.genes, 
                             y = samp1, SIMPLIFY = TRUE), ncol = 4*k.ind)
      normalized1 <- data1/offset.col1
      data.table1 <- normalized1[, (2*k.ind + 1):(4*k.ind), drop = FALSE]
      if (n1.diff > 0){
        data.table1[DE.genes1, 1:k.ind] <- normalized1[DE.genes1, 1:k.ind, drop = FALSE]
      }
    }
    if (n1.genes == 0) data.table1 <- NULL
    
    #perform swaps for group 2
    if (n1.genes < n.genes){
      if (n1.genes > 0) {
        samp2 <- samp[-(1:n1.genes), , drop = FALSE]
        offset.col2 <- offset.col[-(1:n1.genes), , drop = FALSE]
      }
      if (n1.genes == 0) {
        samp2 <- samp
        offset.col2 <- offset.col
      }
      data.temp2 <- counts[genes.subset2, , drop = FALSE]
      data2 <- 
        matrix(mapply(function(x, y) data.temp2[x, y], x = 1:(n.genes - n1.genes), 
                      y = samp2, SIMPLIFY = TRUE), ncol = 4*k.ind)
      normalized2 <- data2/offset.col2
      data.table2 <- normalized2[, 1:(2*k.ind), drop = FALSE]
      if (n.diff - n1.diff > 0){
        data.table2[DE.genes2,1:k.ind] <- normalized2[DE.genes2, (2*k.ind+1):(3*k.ind), drop = FALSE]
      }
    }
    if (n1.genes == n.genes) data.table2 <- NULL
    
    #combine matrices from both groups
    data.table <- rbind(data.table1, data.table2)
    #unnormalize dataset
    data.table = round(t(t(data.table)*offset[sample(1:n.col, 2*k.ind, replace = FALSE)]))
    #order by genes.subset
    data.table <- data.table[order(c(genes.subset1, genes.subset2)), ]
  }
  
  return(list(counts = data.table, genes.subset = genes.subset, 
              genes.subset1 = genes.subset1, genes.subset2 = genes.subset2, 
              DE.genes = genes.diff, DE.genes1 = genes.diff1, 
              DE.genes2 = genes.diff2))
}
