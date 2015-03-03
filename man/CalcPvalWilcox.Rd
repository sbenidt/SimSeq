\name{CalcPvalWilcox}
\alias{CalcPvalWilcox}
\title{
Calculate P-values of Differential Expression
}
\description{
A function called within SimData. Calculates the P-value of differential expression for each gene under either Wilcoxin signed rank test (paired data) or rank sum test (unpaired data). The vector of differences used in the test is based off of the log of the counts for each gene plus 1 divided by their respective multiplicative normalization factors provided by the user.
}
\usage{
CalcPvalWilcox(counts, treatment, replic = NULL, sort.method, 
               sorted = FALSE, norm.factors, exact = FALSE)
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
Character vector specifying one of "unpaired" or "paired", depending on the structure of the data.
}
  \item{sorted}{
logical vector specifying whether data has already gone through SortData function. Defaults to FALSE.
}
  \item{norm.factors}{
A positive numeric vector of multiplicative normalization factors for each column of the counts matrix.
}
  \item{exact}{
Specifies whether an exact signed rank test (paired) or exact ranksum test (unpaired) should be used.
}
}
\value{
probs: a vector of p-values of differential expression for each gene.
}
\author{
Samuel Benidt <sbenidt@iastate.edu>
}
