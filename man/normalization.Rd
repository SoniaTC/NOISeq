\name{Normalization}
\alias{uqua}
\alias{rpkm}
\alias{tmm}

\title{
Normalization methods
}
\description{
Normalization procedures such as RPKM (Mortazavi et al., 2008),
Upper Quartile (Bullard et al., 2010) and TMM (Trimmed Mean of M) (Robinson and Oshlack, 2010). 
These normalization functions are used within the \code{noiseq} or \code{noiseqbio} functions but may be also used by themselves to normalize a dataset.
}
\usage{
uqua(datos, long = 1000, lc = 1, k = 0)
rpkm(datos, long = 1000, lc = 1, k = 0)
tmm(datos, long = 1000, lc = 1, k = 0, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
}

\arguments{
  \item{datos}{
Matrix containing the read counts for each sample.
}
  \item{long}{
Numeric vector containing the length of the features. If long == 1000, no length correction is applied.
}
  \item{lc}{
Correction factor for length normalization. This correction is done by dividing the counts 
vector by length^lc. By default, lc = 1.
}
  \item{k}{
Counts equal to 0 are changed to k in order to avoid indeterminations when applying logarithms, for instance. By default, k = 0.5
}
\item{refColumn}{
Column to use as reference (only needed for \code{tmm} function).
}
\item{logratioTrim}{
Amount of trim to use on log-ratios ("M" values) (only needed for \code{tmm} function).
}
\item{sumTrim}{
Amount of trim to use on the combined absolute levels ("A" values) (only needed for \code{tmm} function).
}
\item{doWeighting}{
Logical, whether to compute (asymptotic binomial precision) weights (only needed for \code{tmm} function).
}
\item{Acutoff}{
Cutoff on "A" values to use before trimming (only needed for \code{tmm} function).
}
}
\details{
\code{tmm} normalization method was taken from \emph{edgeR} package (Robinson et al., 2010). 

Although \code{Upper Quartile} and \code{TMM} methods themselves do not correct for the length of the features, these functions in
\code{NOISeq} allow users to combine the normalization procedures with an additional length correction whenever 
the length information is available.
}


\examples{

## Simulate some count data and the features length
datasim = matrix(sample(0:100, 2000, replace = TRUE), ncol = 4)
lengthsim = sample(100:1000, 500)

## RPKM normalization
myrpkm = rpkm(datasim, long = lengthsim, lc = 1, k = 0)

## Upper Quartile normalization, dividing normalized data by the square root of the features length and replacing counts=0 by k=1
myuqua = uqua(datasim, long = lengthsim, lc = 0.5, k = 1)

## TMM normalization with no length correction
mytmm = tmm(datasim, long = 1000, lc = 1, k = 0)

}


\references{

Bullard J.H., Purdom E., Hansen K.D. and Dudoit S. (2010) Evaluation of statistical methods for 
normalization and differential expression in mRNA-seq experiments. \emph{BMC Bioinformatics} 11(1):94+.

Mortazavi A., Williams B.A., McCue K., Schaeer L. and Wold B. (2008) Mapping and quantifying 
mammalian transcriptomes by RNA-seq. \emph{Nature Methods} 5(7):621-628.

Robinson M.D. and Oshlack A. (2010) A scaling normalization method for differential expression
analysis of RNA-seq data. \emph{Genome Biology} 11(3):R25+.

Robinson M.D., McCarthy D.J. and Smyth G.K. (2010) edgeR: a Bioconductor package for differential 
expression analysis of digital gene expression data. \emph{Bioinformatics} 26(1):139-140.

}
\author{
Sonia Tarazona
}

