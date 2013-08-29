\name{SimData}
\alias{SimData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SimData}
\description{
Given a matrix of RNA-Seq data, data.sim simulates a new RNA-Seq dataset with a known list of differentially expressed (DE) and equivalently expressed (EE) genes by sampling two groups of k.ind replicates and performing appropriate swaps between the two treatment group for each differentially expressed gene. One with a combination of the two treatment groups swapped in to ensure a set list of DE genes.
}
\usage{
SimData(counts, treatment, replic = NULL, sort.method, k.ind, n.genes = NULL, 
	n.diff = NULL, n1.genes = NULL, n1.diff = NULL, offset = NULL, 
	samp.independent = FALSE, genes.select = NULL, genes.diff = NULL, 
	genes1.select = NULL, probs = NULL, EBP = NULL, exact = FALSE, power = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
A matrix of counts where each row specifies a gene and each column specifies a replicate; counts matrix must contain data from an experiment with two treatment groups of either paired or unpaired data.
}
  \item{treatment}{
A vector specifying the treatment group for each column of the counts matrix; should match number of columns of counts matrix.
}
  \item{replic}{
A vector specifying the replicate for each column of the counts matrix when there is paired data; optional if data is unpaired.
}
  \item{sort.method}{
Character vector specifying one of "unpaired" or "paired", depending on the structure of the data.
}
  \item{k.ind}{
The number of replicates in the simulated matrix for each treatment group.
}
  \item{n.genes}{
The number of genes subsetted from the counts matrix. Must be less than the total number of rows in the counts matrix but greater than n.diff. Optional if genes.select vector is specified.
}
  \item{n.diff}{
The number of genes simulated to be differentially expressed. Must be less than n.genes. Optional if genes.diff vector is specified.
}
  \item{n1.genes}{
The number of genes simulated to be equivalently expressed by sampling only from the first treatment group plus the number of genes simulated to be differentially expressed by sampling the first column from the first treatment group and the second column from the second treatment group. Optional.
}
  \item{n1.diff}{
The number of genes simulated to be differentially expressed by sampling the first column from the first treatment group and the second column from the second treatment group. The other n.diff - n1.diff differentially expressed genes are simulated by sampling the first column from the second treatment group and the second column from the first treatment group. Optional. 
}
  \item{offset}{
A positive numeric vector with length equal to the number of columns of counts matrix. Will default to CalcNormFactors function from the package edgeR with method = "TMM".
}
  \item{samp.independent}{
A logical vector specifying whether the columns sampled for each gene should sampled separately and independently. Defaults to FALSE.
}
  \item{genes.select}{
A vector specifying genes subsetted from the counts matrix in the simulation. Can be either logical or numeric. Optional if n.genes is specified.
}
  \item{genes.diff}{
A vector specifying genes to be differentially expressed in the simulation. Genes selected to be differentially expressed must be a subset of the genes selected in thhe genes.select vector. Can be either logical or numeric. Optional if n.diff is specified.
}
  \item{genes1.select}{
A vector specifying genes to be equivalently expressed through sampling only through the first treatment group or to be differentially expressed from sampling first from the first treatment group and sampling second from the second treatment group.
}
  \item{probs}{
Optional vector specifying the p-value of differential expression for each gene to be used in the estimate of EBP for each gene. If null, data.sim will perform a either a signed rank test (paired case) or a rank sum test (unpaired) case for each gene in the counts matrix.
}
  \item{EBP}{
Optional vector specifying the Empirical Bayes Probability (EBP) of differential expression for each gene in the counts matrix.
}
  \item{exact}{
Specifies whether an exact signed rank test (paired) or exact ranksum test (unpaired) should be used.
}
  \item{power}{
Transorms EBP for each gene by raising each EBP estimate to some power. Must be greater than 0. Default is set to 1.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
List containing:
\item{counts}{Matrix of simulated gene expression data with known list of DE and EE genes}
\item{genes.subset}{all genes included in the simulated matrix}
\item{genes.subset1}{DE and EE genes simulated from group 1}
\item{genes.subset2}{DE and EE genes dimulated from group 2}
\item{DE.genes}{Logical vector specicying which genes in the simulated matrix are Differentially Expressed}
\item{DE.genes1}{DE genes simulated from group 1}
\item{DE.genes2}{DE genes simulated from group 2}

}
\references{
Pounds SB, Gao CL, and Zhang H (2012). Empirical Bayesian Selection of Hypothesis Testing Procedures for Analysis of Sequence Count Expression Data. Statistical Applications in Genetics and Molecular Biology. 11-5, 7.
}
\author{
Samuel Benidt <sbenidt@iastate.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(kidney)
kidney.counts <- kidney[[1]]
replic <- kidney[[2]]
treatment <- kidney[[3]]
offset <- calcNormFactors(kidney.counts, method = "TMM")

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
			sort.method = "paired", k.ind = 5, n.genes = 100, n.diff = 50)

##OR
sort.list <- SortData(counts = kidney.counts, treatment = treatment, replic = replic,
                      sort.method = "paired", offset = offset)
kidney.counts <- sort.list[[1]]
replic <- sort.list[[2]]
treatment <- sort.list[[3]]
offset <- sort.list[[4]]

probs <- CalcPvalWilcox(kidney.counts, treatment, sort.method = "paired", 
                        sorted = TRUE, offset = offset, exact = FALSE)
EBP <- pval.hist(probs)$h.ebp

data.sim <- SimData(counts = kidney.counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 100, n.diff = 50,
                    EBP = EBP, offset = offset)



}