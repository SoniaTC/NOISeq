\name{noiseqbio}
\alias{noiseqbio}

\title{
Differential expression method for biological replicates
}
\description{
\code{noiseqbio} computes differential expression between two experimental conditions from read count data (e.g. RNA-seq).
}
\usage{
noiseqbio(input, k = 0.5, norm = c("rpkm","uqua","tmm","n"), nclust = 15, plot = FALSE,
                      factor=NULL, conditions = NULL, lc = 0, r = 50, adj = 1.5,
                      a0per = 0.9, random.seed = 12345, filter = 1, depth = NULL,
                      cv.cutoff = 500, cpm = 1)
}

\arguments{
  \item{input}{
Object of eSet class coming from \code{\link{readData}} function or other R packages such as DESeq.
}
  \item{k}{
Counts equal to 0 are replaced by k. By default, k = 0.5.
}
  \item{norm}{
Normalization method. It can be one of "rpkm" (default), "uqua" 
(upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
}
  \item{factor}{
A string indicating the name of factor whose levels are the conditions to be compared.
}
  \item{conditions}{
A vector containing the two conditions to be compared by the differential expression algorithm (needed when the \code{factor} contains more than 2 different conditions).
}
  \item{lc}{
Length correction is done by dividing expression by length^lc. By default, lc = 0.
}
  \item{r}{
Number of permutations to generate noise distribution by resampling.
}
  \item{adj}{
Smoothing parameter for the Kernel Density Estimation of noise distribution. Higher values produce smoother curves.
}
  \item{nclust}{
Number of clusters for the K-means algorithm. Used when the number of replicates per condition is less than 5.
}
\item{plot}{
  If TRUE, a plot is generated showing the mixture distribution (f) and the noise distribution (f0) of theta values.
}
\item{a0per}{
  M and D values are corrected for the biological variability by being divided by S + a0, where S is the standard error of the corresponding statistic 
and a0 is determined by the value of a0per parameter. If a0per is NULL, a0 = 0. If a0per is a value between 0 and 1, 
a0 is the a0per percentile of S values for all features. If a0per = "B", a0 takes the highest value given by 100*max(S).
}
\item{random.seed}{
  Random seed. In order to get the same results in different runs of the method (otherwise the resampling procedure would produce different resulst), 
the random seed is set to this parameter value. 
}
\item{filter}{
  Method to filter out low count features before computing differential expression analysis. If filter=0, no filtering is performed. 
If 1, CPM method is applied. 
If 2, Wilcoxon test method (not recommended when the number of replicates per condition is less than 5),
If 3, proportion test method. Type \code{?filtered.data} for more details.
}
\item{depth}{
  Sequencing depth of each sample to be used by filtering method. It must be data provided when the data is already normalized 
and filtering method 3 is to be applied.
}
\item{cv.cutoff}{
  Cutoff for the coefficient of variation per condition to be used in filtering method 1.
}
\item{cpm}{
  Cutoff for the counts per million value to be used in filtering methods 1 and 3.
}
}
\value{
The function returns an object of class \code{\link{Output}}
}

\author{
Sonia Tarazona
}

\seealso{
       \code{\link{readData}}.
     }

\examples{
## Load the input object from Marioni's data as returned by readData()
data(myCounts)

## Computing differential expression probability by NOISeqBIO using factor "Tissue" (data will be RPKM-normalized)
mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "rpkm", factor="Tissue", lc = 1, r = 50, adj = 1.5, plot = FALSE,
                        a0per = 0.9, random.seed = 12345, filter = 1, cv.cutoff = 500, cpm = 1)
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
