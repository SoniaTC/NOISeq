\name{ARSyNseq}
\alias{ARSyNseq}
\alias{arsynseq}

\title{ASCA Removal of Systematic Noise on Seq data}

\description{
 \code{ARSyNseq} filters the noise associated to identified or not identified batch effects considering the experimental design 
and applying Principal Component Analysis (PCA) to the ANOVA parameters and residuals.
}
\usage{
ARSyNseq(data, factor = NULL, batch = FALSE, norm = "rpkm", logtransf = FALSE, Variability = 0.75, beta = 2)
}

\arguments{
  \item{data}{A Biobase's eSet object created with the \code{readData} function.}
  
  \item{factor}{Name of the factor (as it was given to the \emph{readData} function) to be used in the ARSyN model (e.g. the factor containing the batch information). When it is NULL, all the factors are considered.}
  
  \item{batch}{TRUE to indicate that the \emph{factor} argument indicates the batch information. In this case, the \emph{factor} argument must be used to specify the names of the onlu factor containing the information of the batch.}
  
  \item{norm}{Type of normalization to be used. One of ``rpkm'' (default), ``uqua'', ``tmm'' or ``n'' (if data are already normalized). If length was provided through the \emph{readData} function, it will be considered for the normalization (except for ``n''). Please note that if a normalization method if used, the arguments \emph{lc} and  \emph{k} are set to 1 and 0 respectively.}
  
  \item{logtransf}{If FALSE, a log-transformation will be applied on the data before computing ARSyN model 
  to improve the results of PCA on count data.}
  
  \item{Variability}{Parameter for Principal Componentents (PCs) selection of the ANOVA models effects. This is the desired proportion 
    of variability explained for the PC of the main effects (time and experimental group). Variability=0.75 by default.}
  
  \item{beta}{Parameter for PCs selection of the residual model. Components selected will be those that explain more than beta times 
		the average component variability computed as the total data variability divided by the rank of the matrix associated to the factor.
		Default beta=2. }
}

\details{
  
  When batch is identified with one of the factors described in the argument \code{factor} of the \code{data} object, \code{ARSyNseq}
	estimates this effect and removes it by estimating the main PCs of the ANOVA effects associated. Selected PCs will be those that explain 
  more than the variability proportion specified in \code{Variability}.
  
	When batch is not identified, the model estimates the effects associated to each factor of interest and analyses if there exists 
  systematic noise in the residuals. If there is batch effect, it will be identified with the main PCs of these residuals.
  Selected PCs will be those that explain more than \code{beta} times the average component variability.
   
}

\value{
The Biobase's eSet object created with the \code{readData} function that was given as input but replacing the expression data with the filtered expression data matrix.
}
\references{

Nueda, M.J.; Ferrer, A. and Conesa, A. (2012) ARSyN: a method for the identification and removal of
systematic noise in multifactorial time-course microarray experiments. \emph{Biostatistics} 13(3), 553-566.

}

\author{Maria Jose Nueda, \email{mj.nueda@ua.es} }

\examples{

# Generating an artificial batch effect from Marioni's data
data(Marioni)
set.seed(123)
mycounts2 = mycounts
mycounts2[,1:4] = mycounts2[,1:4] + runif(nrow(mycounts2)*4, 3, 5)
myfactors = data.frame(myfactors, "batch" = c(rep(1,4), rep(2,6)))
mydata2 = readData(mycounts2, factors = myfactors)

# Exploring batch effect with PCA
myPCA = dat(mydata2, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")

# Removing batch effect when the batch is identified for each sample and exploring results with PCA
mydata2corr1 = ARSyNseq(mydata2, factor = "batch", batch = TRUE, norm = "rpkm",  logtransf = FALSE)
myPCA = dat(mydata2corr1, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")

# If we consider that exist a batch but it is not identified (we do not know the batch information):
mydata2corr2 = ARSyNseq(mydata2, factor = "Tissue", batch = FALSE, norm = "rpkm",  logtransf = FALSE)
myPCA = dat(mydata2corr2, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Tissue")
explo.plot(myPCA, factor = "batch")
}

\keyword{ASCA, ANOVA, PCA, batch }


