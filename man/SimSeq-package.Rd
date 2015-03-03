\name{SimSeq-package}
\alias{SimSeq-package}
\alias{SimSeq}
\docType{package}
\title{Nonparametric Simulation of RNA-Seq Data}
\description{RNA sequencing analysis methods are often derived by relying on hypothetical parametric models for read counts that are not likely to be precisely satisfied in practice. Methods are often tested by analyzing data that have been simulated according to the assumed model. This testing strategy can result in an overly optimistic view of the performance of an RNA-seq analysis method. We develop a data-based simulation algorithm for RNA-seq data. The vector of read counts simulated for a given experimental unit has a joint distribution that closely matches the distribution of a source RNA-seq dataset provided by the user. Users control the proportion of genes simulated to be differentially expressed (DE) and can provide a vector of weights to control the distribution of effect sizes. The algorithm requires a matrix of RNA-seq read counts with large sample sizes in at least two treatment groups. Many datasets are available that fit this standard.}
\details{
\tabular{ll}{
Package: \tab SimSeq\cr
Type: \tab Package\cr
Version: \tab 1.2.0\cr
Date: \tab 2014-12-09\cr
License: \tab GPL-2\cr
}
SimSeq performs data based simulation of RNA-Seq data creating a dataset with a known list of DE and EE genes. The core function that implements of the methodology of SimSeq is the SimData function.
}
\author{
Samuel Benidt

Maintainer: Samuel Benidt <sbenidt@iastate.edu>
}