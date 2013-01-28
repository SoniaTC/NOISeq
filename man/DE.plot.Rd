\name{Differential expression plots}
\alias{DE.plot}

\title{
Plot to compare expression values for two conditions or to represent differential expression statistics
}
\description{
Function to represent three possible plots: one to compare the expression values in each condition, one also to compare the expression values in each condition but according to the chromosome and position within the chromosome, and the other one with the (M,D) statistics. In all three plots, the differentially expressed features may be highlighted.
}
\usage{
DE.plot(output, q = NULL, graphic = c("MD","expr","chrom","distr"), pch = 20, cex = 0.5, col = 1, pch.sel = 1, cex.sel = 0.6, col.sel = 2,log.scale = TRUE, chromosomes = NULL, join = FALSE,...)
}

\arguments{

  \item{output}{
Object of class \code{\link{Output}}.
}
  \item{q}{
Probability of differential expression threshold to determine differentially expressed features.
}
  \item{graphic}{
String indicating which kind of plot is to be generated. If "expr", the feature expression values are depicted. If "MD", the values for the (M,D) statistics when comparing both conditions are used. If "chrom", the feature expression values are depicted across their positions in the chromosomes (if chromosome information has been provided).
}
  \item{pch, cex, col,...}{
Graphical parameters as in any other R plot. See \code{\link{par}}. They do not apply for graphic="chrom".
}
  \item{pch.sel, cex.sel, col.sel}{
\code{pch}, \code{cex} and \code{col}, respectively, to represent differentially expressed features. They do not apply for graphic="chrom".
}
  \item{log.scale}{
If TRUE, log2(data+K) values are depicted instead of the expression data in the \code{\link{Output}} object. K is an appropriate constant to avoid negative values. It does not apply for graphic="MD".
}
  \item{chromosomes}{
Character vector indicating the chromosomes to be plotted. If NULL, all chromosomes are plotted. It only applies for graphic="chrom" and graphic="distr".  For graphic="chrom", the chromosomes are plotted in the given order. In some cases (e.g. chromosome names are character strings), it is very convenient to specify the order although all chromosomes are being plotted. For graphic="distr", the chromosomes are plotted according to the number of features they contain (from the highest number to the lowest).
}
  \item{join}{
If FALSE, each chromosome is depicted in a separate line. If TRUE, all the chromosomes are depicted in the same line, consecutively (useful for prokaryote organisms). It only applies for graphic="chrom".
}
}

\author{
Sonia Tarazona
}



\seealso{
       \code{\link{readData}}, \code{\link{noiseq}}, \code{\link{degenes}}.
     }

\examples{
## We load the object generated after running noiseq on Marioni's data
data(noiseq)

## Third, plot the expression values for all genes and highlighting the differentially expressed genes
DE.plot(mynoiseq, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseq, q = 0.8, graphic = "MD", ylim = c(0,50))
DE.plot(mynoiseq, chromosomes = c(1,2), log.scale = TRUE,join = FALSE, q = 0.8, graphic = "chrom")
DE.plot(mynoiseq, chromosomes = NULL, q = 0.8, graphic = "distr")
}


\references{

Bullard J.H., Purdom E., Hansen K.D. and Dudoit S. (2010) Evaluation of statistical methods for 
normalization and differential expression in mRNA-seq experiments. \emph{BMC Bioinformatics} 11(1):94+.

Mortazavi A., Williams B.A., McCue K., Schaeer L. and Wold B. (2008) Mapping and quantifying 
mammalian transcriptomes by RNA-seq. \emph{Nature Methods} 5(7):621-628.

Robinson M.D. and Oshlack A. (2010) A scaling normalization method for differential expression
analysis of RNA-seq data. \emph{Genome Biology} 11(3):R25+.

Marioni, J.C. and Mason, C.E. and Mane, S.M. and Stephens, M. and Gilad, Y. (2008) RNA-seq: an assessment of
technical reproducibility and comparison with gene expression arrays. \emph{Genome Research}, \bold{18}:
1509--1517.

}
