\name{pval.hist}
\alias{pval.hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pval.hist
}
\description{
A function called in SimData used to compute the Empirical Bayes Probability of differential expression for each gene given a vector of p-values.
}
\usage{
pval.hist(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
a vector of p-values of differential expression for each gene
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
\item{p}{input p-value vector}
\item{edf}{the p-value edf}
\item{h.cdf}{the histogram-based cdf estimate}
\item{h.fdr}{the histogram-based FDR estimate}
\item{h.ebp}{the histogram-based EBP estimate}
\item{p.brks}{the p-value break-points of the histogram}
\item{p.hist}{the heights of each histogram bar}
\item{p.brks}{the p-value break-points of the histogram}
\item{edf.brks}{the break-points in the EDF of the histogram estimator}
\item{pi0.hat}{the histogram-based estimate of the proportion of tests with a true null hypothesis}
}
\references{
Pounds SB, Gao CL, and Zhang H (2012). Empirical Bayesian Selection of Hypothesis Testing Procedures for Analysis of Sequence Count Expression Data. Statistical Applications in Genetics and Molecular Biology. 11-5, 7.
}
\author{
Pounds SB, Gao CL, and Zhang H
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
