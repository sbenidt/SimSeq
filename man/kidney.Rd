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
  The format is:
List of 3
 $ kidney.counts: num [1:100, 1:138] 0 6 6 140 2042 ...
 $ replic       : Factor w/ 69 levels "3358","3387",..: 1 1 2 2 3 3 4 4 5 5 ...
 $ treatment    : Factor w/ 2 levels "0","1": 1 2 1 2 1 2 1 2 1 2 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
https://tcga-data.nci.nih.gov/tcga/
}
\examples{
data(kidney)
}
\keyword{datasets}
