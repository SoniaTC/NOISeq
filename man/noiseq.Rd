\name{noiseq}
\alias{noiseq}

\title{
Differential expression method for technical replicates or no replicates at all
}
\description{
\code{noiseq} computes differential expression between two experimental conditions from read count data (e.g. RNA-seq).
}
\usage{
noiseq(input, k = 0.5, norm = c("rpkm","uqua","tmm","n"), replicates = c("technical","biological","no"), 
factor=NULL, conditions=NULL, pnr = 0.2, nss = 5, v = 0.02, lc = 1)
}

\arguments{
  \item{input}{
Object of eSet class coming from \code{\link{readData}} function or other R packages such as DESeq.
}
  \item{factor}{
A string indicating the name of factor whose levels are the conditions to be compared.
}
  \item{conditions}{
A vector containing the two conditions to be compared by the differential expression algorithm (needed when the \code{factor} contains more than 2 different conditions).
}
\item{replicates}{
  In this argument, the type of replicates to be used is defined. Technical, biological or none. By default, technical replicates option is chosen.
}
  \item{k}{
Counts equal to 0 are replaced by k. By default, k = 0.5.
}
  \item{norm}{
Normalization method. It can be one of "rpkm" (default), "uqua" 
(upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).
}
  \item{lc}{
Length correction is done by dividing expression by length^lc.
By default, lc = 1.
}
  \item{pnr}{
Percentage of the total reads used to simulated each sample when no replicates are available. By default, pnr = 0.2.
}
  \item{nss}{
Number of samples to simulate for each condition (nss>= 2). By default, nss = 5.
}
  \item{v}{
Variability in the simulated sample total reads. 
By default, v = 0.02. Sample total reads is computed as a
random value from a uniform distribution in the interval
[(pnr-v)*sum(counts), (pnr+v)*sum(counts)]
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

## Computing differential expression probability on RPKM-normalized data by NOISeq-real using factor "Tissue"
mynoiseq = noiseq(mydata, k = 0.5, norm = "rpkm", replicates = "technical", factor="Tissue", pnr = 0.2, nss = 5, v = 0.02, lc = 1)


## Computing differential expression probability on Upper Quartile normalized data by NOISeq-real using factor "TissueRun" and comparing samples in Run 1 (levels "Kidney_1" and "Liver_1")
mynoiseq.uqua = noiseq(mydata, k = 0.5, norm = "uqua", replicates = "technical", factor="TissueRun", conditions = c("Kidney_1","Liver_1"), 
pnr = 0.2, nss = 5, v = 0.02, lc = 1)

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
