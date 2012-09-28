\name{readData}
\alias{readData}
\alias{addData}

\title{
Create an object of eSet class
}
\description{
Create an object of eSet class to be used by NOISeq functions by taking R matrix or data.frame objects.
}
\usage{
readData(data, factors, length = NULL, biotype = NULL, chromosome = NULL)
addData(data, length = NULL, biotype = NULL, chromosome = NULL, factors = NULL)
}

\arguments{
  \item{data}{
The matrix or data.frame containing the counts (or expression data) for each feature and condition. Rows are features and columns are conditions.
}
  \item{length}{
Optional argument. Vector containing the length of each feature. The names of the vector must be the feature names or ids with the same type of identifier used in \code{data}.
}
  \item{factors}{
A data.frame containing the different conditions of each sample included in the \code{data} object.
}
  \item{biotype}{
Optional argument. Vector containing the biological classification of each feature. The names of the vector must be the feature names or ids with the same type of identifier used in \code{data}.
}
  \item{chromosome}{
Optional argument. A matrix or data.frame containing the chromosome, start position and end position of each feature. The row names must be the feature names or ids with the same type of identifier used in \code{data}.
}
}

\value{
It returns an object of eSet class \code{\link{myCounts}} with all the information defined and ready to be used.
}

\author{
Sonia Tarazona
}

\examples{
# Load an object containing the information explained above
data(Marioni)

# Create the object with the data        
mydata <- readData(data=mycounts, biotype=mybiotypes, chromosome=mychroms, factors=myfactors)

# Add length annotation to the existing data object
mydata <- addData(mydata, length=mylength)
}

\references{ 
    Marioni, J.C. and Mason, C.E. and Mane, S.M. and Stephens, M. and Gilad, Y. (2008) RNA-seq: an assessment of
technical reproducibility and comparison with gene expression arrays. \emph{Genome Research}, \bold{18}:
1509--1517.
}
