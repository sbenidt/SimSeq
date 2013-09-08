\name{SortData}
\alias{SortData}
\title{
SortData
}
\description{
A function called in SimData used to trim and sort the matrix of counts provided
}
\usage{
SortData(counts, treatment, replic = NULL, sort.method, offset = NULL)
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
  \item{offset}{
A positive numeric vector with length equal to the number of columns of counts matrix.
}
}
\value{
List containing:
\item{counts}{sorted and trimmed matrix of counts}
\item{replic}{sorted and trimmed replic vector}
\item{treatment}{sorted and trimmed treatment vector}
\item{offset}{sorted and trimmed offset vector}
\item{sorting}{sorting vector used to sort and trim}
}
\author{
Samuel Benidt <sbenidt@iastate.edu>
}