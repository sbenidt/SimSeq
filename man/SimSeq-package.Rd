\name{SimSeq-package}
\alias{SimSeq-package}
\alias{SimSeq}
\docType{package}
\title{Data Based Simulation of RNA-Seq Data}
\description{SimSeq is a package that creates a simulated RNA-Seq dataset with a known set of differentially expressed (DE) genes and equivalently expressed (EE) genes. This simulation assumes no underlying statistical model, but rather is based on a provided matrix of real gene expression data with large sample sizes in each treatment group. Such a simulation is useful to evaluate statistical methodology for RNA-Seq in that false discovery proportions (FDP's) are known exactly due to knowledge of which genes are simulated to be DE and which are simulated to be EE.}
\details{
\tabular{ll}{
Package: \tab SimSeq\cr
Type: \tab Package\cr
Version: \tab 0.4.0\cr
Date: \tab 2013-09-07\cr
License: \tab GPL-2\cr
}
SimSeq performs data based simulation of RNA-Seq data creating a dataset with a known list of DE and EE genes. The core function that implements of the methodology of SimSeq is the SimData function.
}
\author{
Samuel Benidt

Maintainer: Samuel Benidt <sbenidt@iastate.edu>
}