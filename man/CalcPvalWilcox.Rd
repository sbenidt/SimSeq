\name{CalcPvalWilcox}
\alias{CalcPvalWilcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate P-values of Differential Expression
}
\description{
A function called within SimData. Calculates the P-value of differential expression for each gene under either Wilcoxin signed rank test (paired data) or rank sum test (unpaired data). Will base test of the log of the counts for each gene plus 1 divided by their respective column offsets.
}
\usage{
CalcPvalWilcox(counts, treatment, replic = NULL, sort.method, sorted = FALSE, offset, exact = FALSE)
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
  \item{sorted}{
logical vector specifying whether data has already gone through SortData function. Defaults to FALSE.
}
  \item{offset}{
A positive numeric vector with length equal to the number of columns of counts matrix.
}
  \item{exact}{
Specifies whether an exact signed rank test (paired) or exact ranksum test (unpaired) should be used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
probs: A vector of p-values of differential expression for each gene.
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
