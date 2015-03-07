\name{SimData}
\alias{SimData}
\title{SimData}
\description{
Given a matrix of RNA-seq data with large samples sizes in at least two treatment groups, SimData simulates a new matrix of RNA-seq data with a known list of differentially expressed (DE) and equivalently expressed (EE) genes from the original data.
}
\usage{
SimData(counts, treatment, replic = NULL, sort.method, k.ind, 
  n.genes = NULL, n.diff = NULL, norm.factors = NULL, 
  samp.independent = FALSE, genes.select = NULL, genes.diff = NULL, switch.trt = FALSE,
  probs = NULL, weights = NULL, exact = FALSE, power = 1)
}
\arguments{
  \item{counts}{
A matrix of counts where each row specifies a gene and each column specifies a replicate.
}
  \item{treatment}{
A vector specifying the treatment group for each column of the counts matrix. Only two treatment groups of either paired or unpaired data are allowed.
}
  \item{replic}{
A vector specifying the replicate for each column of the counts matrix when there is paired data; optional if data is unpaired.
}
  \item{sort.method}{
One of either "paired" or "unpaired", depending on the structure of the counts matrix.
}
  \item{k.ind}{
The number of experimental units to be simulated for each treatment group.
}
  \item{n.genes}{
The number of genes to be subsetted from the counts matrix in the simulation. Must be less than the total number of rows in the counts matrix and greater than n.diff. Optional if genes.select vector is specified.
}
  \item{n.diff}{
The number of genes simulated to be differentially expressed. Must be less than n.genes. Optional if genes.diff vector is specified.
}
  \item{norm.factors}{
A positive numeric vector of multiplicative normalization factors for each column of the counts matrix.  Will default to CalcNormFactors function from the package edgeR with method = "TMM".
}
  \item{samp.independent}{
A logical vector specifying whether the columns sampled for each gene should sampled separately and independently for each gene. Defaults to FALSE.
}
  \item{genes.select}{
A vector specifying genes to be subsetted from the counts matrix in the simulation. Can be either a logical vector or a numeric vector indicating the rows of the counts matrix to be used. Optional if n.genes is specified.
}
  \item{genes.diff}{
A vector specifying genes to be differentially expressed in the simulation. Genes selected to be differentially expressed must be a subset of the genes selected in the genes.select vector. Can be either logical vector with length equal to genes selected or a numeric vector indicating the rows of the counts matrix to be used. Optional if n.diff is specified.
}
  \item{switch.trt}{
Logical specifying which treatment group should be sampled from for EE genes. Default is first treatment group.  
  }
  \item{probs}{
Optional vector specifying the p-value of differential expression for each gene to be used in the estimate of empirical Bayes probability for each gene. If not provided and weights are not specified, SimData will perform either a signed rank test (paired case) or a rank sum test (unpaired) case for each gene in the counts matrix.
}
  \item{weights}{
Optional vector specifying weights to be used for sampling which genes are to be differentially expressed in the simulation. If null, weights will be calculated using the fdrtool function from the package 'fdrtool' to calculate one minus local fdr. If desired, the sampling of differentially expressed genes can be done without respect to any weights by providing a vector of ones.
}
  \item{exact}{
Specifies whether an exact signed rank test (paired) or exact ranksum test (unpaired) should be used.
}
  \item{power}{
Transorms the weights for each gene by raising each weights to some power. Must be greater than 0. Default is set to 1.

}
}
\value{
List containing:
\item{counts}{matrix of simulated RNA-seq data with known list of DE and EE genes.}
\item{treament}{a numeric vector specifying the treatment group structure in the new simulated counts matrix.}
\item{genes.subset}{all genes included in the simulated matrix. The i-th entry corresponds to the i-th row of the simulated counts matrix.}
\item{DE.genes}{vector of genes used for differential expression in the simulated counts matrix.}
\item{DE.ind}{logical vector indicating which genes are differentially expressed in the simulated counts matrix.}
\item{col}{vector of lenght 3*k.ind specifying which columns were used in the simulation. The first k.ind columns were used to simulate the first treatment group. The last 2*k.ind columns were used to simulate the second treatment group. If switch.trt equals TRUE, then the first 2*k.ind columns were used to simulate the first treatment group and the last k.ind columns were used to simulate the second treatment group.}
}

\details{
SimData simulates an RNA-seq matrix of counts from a large source RNA-sequence dataset using a resampling based approach. If you use the Kidney Cancer (KIRC) dataset from the The Cancer Genome Atlas as the source RNA-seq counts matrix in the SimData function, please cite The Cancer Genome Atlas Research Network (2013) in any publications.
}

\author{
Samuel Benidt <sbenidt@iastate.edu>
}

\source{
Benidt, S., and Nettleton, D. (2015). SimSeq: a nonparametric approach to simulation of RNA-sequence datasets. Bioinformatics. doi: 10.1093/bioinformatics/btv124. 

The Cancer Genome Atlas Research Network (2013). Comprehensive molecular characterization of clear cell renal cell carcinoma. Nature, 499(7456), 43-49.

https://tcga-data.nci.nih.gov/tcga/
}

\examples{
data(kidney)
counts <- kidney$counts # Matrix of read counts from KIRC dataset
replic <- kidney$replic # Replic vector indicating paired columns
treatment <- kidney$treatment # Treatment vector indicating Non-Tumor or Tumor columns

nf <- apply(counts, 2, quantile, 0.75)

require(fdrtool)

\dontrun{

  ### Example 1: Simulate Matrix with 1000 DE genes and 4000 EE genes
  data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                      sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                      norm.factors = nf)
  
  ### Example 2: Calculate weights vector beforehand to save run time in
  ### repeated simulations
  sort.list <- SortData(counts = counts, treatment = treatment, replic = replic,
                        sort.method = "paired", norm.factors = nf)
  counts <- sort.list$counts
  replic <- sort.list$replic
  treatment <- sort.list$treatment
  nf <- sort.list$norm.factors
  
  probs <- CalcPvalWilcox(counts, treatment, sort.method = "paired", 
                          sorted = TRUE, norm.factors = nf, exact = FALSE)
  weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 
  
  data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                      sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                      weights = weights, norm.factors = nf)
  
  ### Example 3: Specify which genes you want to use in the simulation
  
  # Randomly sample genes or feed in the exact genes you wish to use
  genes.diff <- sample(1:nrow(counts), size = 1000, prob = weights)
  genes <- c(sample(1:nrow(counts)[-genes.diff], 4000), genes.diff)
  
  data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                      sort.method = "paired", k.ind = 5, genes.select = genes,
                      genes.diff = genes.diff, weights = weights, norm.factors = nf)
  
  ### Example 4: Simulate matrix with DE genes having log base 2 fold change greater than 1
  
  # add one to counts matrix to avoid infinities when taking logs
  tumor.mean <- rowMeans(log2((counts[, treatment == "Tumor"] + 1) \%*\% 
    diag(1/nf[treatment == "Tumor"])))
  nontumor.mean <- rowMeans(log2((counts[, treatment == "Non-Tumor"] + 1) \%*\% 
    diag(1/nf[treatment == "Non-Tumor"])))
  
  lfc <- tumor.mean - nontumor.mean
  weights.zero <- abs(lfc) < 1
  weights[weights.zero] <- 0

  data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                      sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                      weights = weights, norm.factors = nf)
  
  ### Example 5: Simulate three treatment groups:
  ### 3 Different types of Differential Expression Allowed
  ### First Group Diff, Second and Third group Equal
  ### Second Group Diff, First and Third group Equal
  ### Third Group Diff, First and Second group Equal
  
  k <- 5 # Sample Size in Each treatment group
  
  ### Sample DE genes beforehand
  N <- nrow(counts)
  genes.de <- sample(1:N, size = 1000, prob = weights) # Sample all DE genes
  DE1 <- genes.de[1:333] # Sample DE genes with first trt diff
  DE2 <- genes.de[334:666] # Sample DE genes with sec trt diff
  DE3 <- genes.de[667:1000] # Sample DE genes with third trt diff
  EE <- sample( (1:N)[-genes.de], size = 4000) #Sample EE genes
  
  genes.tot <- c(EE, genes.de)
  genes.de1 <- union(DE2, EE) #Assign DE genes for first sim
  genes.de2 <- union(DE2, DE3) #Assign DE genes for second sim
  
  data.sim1 <- SimData(counts = counts, replic = replic, treatment = treatment, 
                       sort.method = "paired", k.ind = k, genes.select = genes.tot,
                       genes.diff = genes.de1, weights = weights, norm.factors = nf)
  
  #remove pairs of columns used in first simulation
  cols.rm <- c(data.sim1$col[1:(2*k)], data.sim1$col[1:(2*k)] + 1)
  counts.new <- counts[, -cols.rm]
  nf.new <- nf[-cols.rm]
  replic.new <- replic[-cols.rm]
  treatment.new <- treatment[-cols.rm]
  
  ### Set switch.trt = TRUE for second sim
  data.sim2 <- SimData(counts = counts.new, replic = replic.new, treatment = treatment.new, 
                       sort.method = "paired", k.ind = k, genes.select = genes.tot,
                       genes.diff = genes.de2, weights = weights, norm.factors = nf.new,
                       switch.trt = TRUE)
  
  ### Remove first k.ind entries from first sim and combine two count matrices
  counts.sim <- cbind(data.sim1$counts[, -(1:k)],  data.sim2$counts)
  
  ### treatment group levels for simulated matrix
  trt.grp <- rep(NA, 5000)
  trt.grp[is.element(data.sim1$genes.subset, DE1)] <- "DE_First_Trt"
  trt.grp[is.element(data.sim1$genes.subset, DE2)] <- "DE_Second_Trt"
  trt.grp[is.element(data.sim1$genes.subset, DE3)] <- "DE_Third_Trt"
  trt.grp[is.element(data.sim1$genes.subset, EE)] <- "EE"
}

\dontshow{
### Shorter running code for CRAN checks
### subsetted counts matrix down to 1500 genes

data(kidney)
attach(kidney)

nf <- apply(counts, 2, quantile, 0.75)

### subset counts matrix 1500 genes to make CRAN check faster
samp.rows <- sample(1:nrow(counts), 1500)
counts <- counts[samp.rows, ]

require(fdrtool)
data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    norm.factors = nf)

##Save run time in repeated simulations
sort.list <- SortData(counts = counts, treatment = treatment, replic = replic,
                      sort.method = "paired", norm.factors = nf)
counts <- sort.list$counts
replic <- sort.list$replic
treatment <- sort.list$treatment
nf <- sort.list$norm.factors

probs <- CalcPvalWilcox(counts, treatment, sort.method = "paired", 
                        sorted = TRUE, norm.factors = nf, exact = FALSE)
weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 

data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    weights = weights, norm.factors = nf)

#specify exactly which genes to use in simulation
genes.select <- sample(1:nrow(counts), 500)
genes.diff <- sample(genes.select, 100)

data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, genes.select = genes.select,
                    genes.diff = genes.diff, weights = weights, norm.factors = nf)

### Simulate three treatment groups:
### 3 Different types of Differential Expression Allowed
### First Group Diff, Second and Third group Equal
### Second Group Diff, First and Third group Equal
### Third Group Diff, First and Second group Equal

k <- 5 # Sample Size in Each treatment group

### Sample DE genes beforehand
N <- nrow(counts)
genes.de <- sample(1:N, size = 100, prob = weights) # Sample all DE genes
DE1 <- genes.de[1:33] # Sample DE genes with first trt diff
DE2 <- genes.de[34:66] # Sample DE genes with sec trt diff
DE3 <- genes.de[67:100] # Sample DE genes with third trt diff
EE <- sample( (1:N)[-genes.de], size = 400) #Sample EE genes

genes.tot <- c(EE, genes.de)
genes.de1 <- union(DE2, EE) #Assign DE genes for first sim
genes.de2 <- union(DE2, DE3) #Assign DE genes for second sim

data.sim1 <- SimData(counts = counts, replic = replic, treatment = treatment, 
                     sort.method = "paired", k.ind = k, genes.select = genes.tot,
                     genes.diff = genes.de1, weights = weights, norm.factors = nf)

#remove pairs of columns used in first simulation
cols.rm <- c(data.sim1$col[1:(2*k)], data.sim1$col[1:(2*k)] + 1)
counts.new <- counts[, -cols.rm]
nf.new <- nf[-cols.rm]
replic.new <- replic[-cols.rm]
treatment.new <- treatment[-cols.rm]

### Set switch.trt = TRUE for second sim
data.sim2 <- SimData(counts = counts.new, replic = replic.new, treatment = treatment.new, 
                     sort.method = "paired", k.ind = k, genes.select = genes.tot,
                     genes.diff = genes.de2, weights = weights, norm.factors = nf.new,
                     switch.trt = TRUE)

### Remove first k entries from first sim and combine two count matrices
counts.sim <- cbind(data.sim1$counts[, -(1:k)],  data.sim2$counts)

### treatment group levels for simulated matrix
trt.grp <- rep(NA, 500)
trt.grp[is.element(data.sim1$genes.subset, DE1)] <- "DE_First_Trt"
trt.grp[is.element(data.sim1$genes.subset, DE2)] <- "DE_Second_Trt"
trt.grp[is.element(data.sim1$genes.subset, DE3)] <- "DE_Third_Trt"
trt.grp[is.element(data.sim1$genes.subset, EE)] <- "EE"
}
}