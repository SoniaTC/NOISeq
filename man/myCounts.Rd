\name{myCounts}
\alias{myCounts-class}
\alias{myCounts}
\docType{class}
\title{Class myCounts}

\description{
This is the main class which contains the information needed to do the different analyses.
}

\section{Extends}{
Class \code{eSet} (package 'Biobase').
}

\section{Quick View}{
This object will contain the expression data and further information needed to do the exploratory analysis or the normalization such as the length, GC content, biotypes, chromosomes and positions for each feature.

Internally, the data is stored as follows:

As \code{myCounts} derives from \code{eSet}, we have used the slot \code{assayData} to store all the expression data, \code{phenoData} to store the factors with the conditions, \code{featureData} which will contain the variables \code{Length}, \code{GCcontent}, \code{Biotype}, \code{Chromosome}, {Start Position}, \code{End Position} for each feature. It has been used the slot \code{experimentData} derived from \code{MIAME-class} which will contain the type of replicates (biological replicates, technical replicates or no replicates at all).
}
\seealso{

If you need further information to know the methods that can be used, see \code{eSet}, \code{AnnotatedDataFrame-class}, \code{MIAME-class}. 
}

\author{Sonia Tarazona}

\keyword{classes}
