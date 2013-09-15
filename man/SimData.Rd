\name{SimData}
\alias{SimData}
\title{SimData}
\description{
Given a matrix of RNA-Seq data, SimData simulates a new matrix of RNA-Seq data with a known list of differentially expressed (DE) and equivalently expressed (EE) genes from the original data.
}
\usage{
SimData(counts, treatment, replic = NULL, sort.method, k.ind, n.genes = NULL, 
	n.diff = NULL, norm.factors = NULL, samp.independent = FALSE, genes.select = NULL, 
	genes.diff = NULL, probs = NULL, weights = NULL, exact = FALSE, power = 1)
}
%- maybe also 'usage' for other objects documented here.
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
Character vector specifying one of "unpaired" or "paired", depending on the structure of the counts matrix.
}
  \item{k.ind}{
The number of replicates to be simulated for each treatment group.
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
A logical vector specifying whether the columns sampled for each gene should sampled separately and independently. Defaults to FALSE.
}
  \item{genes.select}{
A vector specifying genes to be subsetted from the counts matrix in the simulation. Can be either logical or numeric. Optional if n.genes is specified.
}
  \item{genes.diff}{
A vector specifying genes to be differentially expressed in the simulation. Genes selected to be differentially expressed must be a subset of the genes selected in the genes.select vector. Can be either logical or numeric. Optional if n.diff is specified.
}
  \item{probs}{
Optional vector specifying the p-value of differential expression for each gene to be used in the estimate of EBP for each gene. If not provided and weights are not specified, SimData will perform either a signed rank test (paired case) or a rank sum test (unpaired) case for each gene in the counts matrix.
}
  \item{weights}{
Optional vector specifying weights to be used for sampling which genes are to be differentially expressed in the simulation. If null, weights will be calculated using the fdrtool function from the package 'fdrtool' to calculate one minus local fdr. If desired, the sampling of differentially expressed genes can be done without respect to weights by providing a vector of ones.
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
\item{counts}{matrix of simulated gene expression data with known list of DE and EE genes}
\item{treament}{a numeric vector specifying the treatment group structure in the new simulated counts matrix}
\item{genes.subset}{all genes included in the simulated matrix}
\item{DE.genes}{vector of differentially expressed genes}
\item{samp.col}{if samp.indepent equals FALSE, gives columns sampled from the counts matrix}

}

\details{
SimData simulates RNA-Sequence from a given RNA-Sequence dataset by the following method:
\enumerate{
\item Sample n.diff genes to be differentially expressed according to a vector of weights supplied.
\item Sample n.genes - n.diff genes to be equivalentlly expressed according to equal weights on the remaining genes.
\item Sample one column of data from each of the two treatment groups and denote these two columns as N and T. Let o_N and o_T denote their normalization factors respectively.
\item Subset these two columns to the number of genes chosen to simulated.
\item Create the column Y_1 in the following manner: for each simulated gene i in the subsetted columns where i = 1,...,n.genes, if gene i is chosen to be differentially expressed, set Y_{1,i} = N_i * o_T/o_N. Else set Y_{1,i} = T_i.
\item Sample one additional column Y_2 from the second treatment group. Then we have simulated a pair of columns Y_1 and Y_2 with a known list of DE and EE genes between them. 
\item Repeat steps 3-6, k.ind times with different columns sampled each time.
}


}

\author{
Samuel Benidt <sbenidt@iastate.edu>
}

\examples{
data(kidney)
attach(kidney)

library(edgeR)
nf <- calcNormFactors(kidney.counts, method = "TMM")

library(fdrtool)
data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    norm.factors = nf)

##Save run time in repeated simulations
sort.list <- SortData(counts = kidney.counts, treatment = treatment, replic = replic,
                      sort.method = "paired", norm.factors = nf)
kidney.counts <- sort.list$counts
replic <- sort.list$replic
treatment <- sort.list$treatment
nf <- sort.list$norm.factors

probs <- CalcPvalWilcox(kidney.counts, treatment, sort.method = "paired", 
                        sorted = TRUE, norm.factors = nf, exact = FALSE)
weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 500, n.diff = 100,
                    weights = weights, norm.factors = nf)

#specify exactly which genes to use in simulation
genes.select <- sample(1:nrow(kidney.counts), 500)
genes.diff <- sample(genes.select, 100)

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, genes.select = genes.select,
                    genes.diff = genes.diff, weights = weights, norm.factors = nf)



}