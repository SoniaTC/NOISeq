\name{Differential expression plots}
\alias{DE.plot}

\title{
Plotting differential expression results
}
\description{
Function to generate plots showing different aspects of differential expression results. 
Expression plot is to compare the expression values in each condition for all features. Differentially expressed features can be highlighted.
Manhattan plot is to compare the expression values in each condition across all the chromosome positions. Differentially expressed features can also be highlighted.
MD plot shows the values for (M,D) statistics. Differentially expressed features can also be highlighted.
Distribution plot displays the percentage of differentially expressed features per chromosome and biotype (if this information is provided by the user).
}
\usage{
DE.plot(output, q = NULL, graphic = c("MD","expr","chrom","distr"), pch = 20, cex = 0.5, col = 1, pch.sel = 1, cex.sel = 0.6, col.sel = 2,log.scale = TRUE, 
chromosomes = NULL, join = FALSE,...)
}

\arguments{

  \item{output}{
Object of class \code{\link{Output}}.
}
  \item{q}{
Probability of differential expression threshold to determine differentially expressed features.
}
  \item{graphic}{
String indicating which kind of plot is to be generated. If "expr", the feature expression values are depicted. 
If "MD", the values for the (M,D) statistics when comparing both conditions are used. 
If "chrom", the feature expression values are depicted across their positions in the chromosomes (if chromosome information has been provided).
If "distr", two plots showing the percentage of differentially expressed features per both chromosome and biotype are generated 
(only if this information is available).
}
  \item{pch, cex, col,...}{
Graphical parameters as in any other R plot. See \code{\link{par}}. They do not apply for graphic="chrom".
}
  \item{pch.sel, cex.sel, col.sel}{
\code{pch}, \code{cex} and \code{col}, respectively, to represent differentially expressed features. They do not apply for graphic="chrom".
}
  \item{log.scale}{
If TRUE, log2(data+K) values are depicted instead of the expression data in the \code{\link{Output}} object. K is an appropriate constant to avoid negative values. 
It does not apply for graphic="MD" and graphic="distr".
}
  \item{chromosomes}{
Character vector indicating the chromosomes to be plotted. If NULL, all chromosomes are plotted. It only applies for graphic="chrom" and graphic="distr". For graphic="chrom", the chromosomes are plotted in the given order. In some cases (e.g. chromosome names are character strings), it is very convenient to specify the order although all chromosomes are being plotted. For graphic="distr", the chromosomes are plotted according to the number of features they contain (from the highest number to the lowest).
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
DE.plot(mynoiseq, q = 0.8, graphic = "MD")
DE.plot(mynoiseq, chromosomes = c(1,2), log.scale = TRUE,join = FALSE, q = 0.8, graphic = "chrom")
DE.plot(mynoiseq, chromosomes = NULL, q = 0.8, graphic = "distr")
}
