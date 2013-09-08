\name{kidney}
\alias{kidney}
\docType{data}
\title{
Kidney Renal Clear Cell Carcinoma [KIRC] RNA Seq data
}
\description{
A list of RNA Seq data downloaded from The Cancer Genome Atlas. Provides a matrix of paired data with rows corresponding to genes and columns corresponding to replicates; replic vector specifies replicates and treatment vector specifies non-tumor and tumor group samples respectively within replicate.
}
\usage{data(kidney)}
\format{
List containing:
\itemize{
\item kidney.counts: Matrix of RNA Seq data for 100 sampled genes and 69 paired columns from individuals with Kidney Renal Clear Cell Carcinoma
\item replic: vector detailing which column in kidney.counts matrix belongs to each individual
\item treatment: vector detailing whether each column in kidney.counts matrix is a non-tumor or tumor sample
}
}
\source{
https://tcga-data.nci.nih.gov/tcga/
}
\examples{
data(kidney)
}
\keyword{datasets}
