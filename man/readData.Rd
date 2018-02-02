\name{readData}
\alias{readData}
\alias{addData}

\title{
Creating an object of eSet class
}
\description{
This function is to create an object of eSet class to be used by NOISeq functions from matrix or data.frame R objects.
}
\usage{
readData(data, factors, length = NULL, biotype = NULL, chromosome = NULL, gc = NULL)
addData(data, length = NULL, biotype = NULL, chromosome = NULL, factors = NULL, gc = NULL)
}

\arguments{
  \item{data}{
Matrix or data.frame containing the counts (or expression data) for each feature and sample. Features must be in rows and samples must be in columns.
}

  \item{factors}{
A data.frame containing the experimental condition or group for each sample (columns in the \code{data} object).
}
  \item{biotype}{
    Optional argument.Vector, matrix or data.frame containing the biological group (biotype)
  for each feature. In case of giving a vector, the names of the
vector must be the feature names or ids  with the same type of
identifier used in \code{data}. If a matrix or a data.frame is provided, 
and it has two columns, it is expected that the feature names or ids are
in the first column and the biotypes of the features in the
second. If it only has one column containing the biotypes, the rownames of the
object must be the feature names or ids.
}
  \item{chromosome}{
Optional argument. A matrix or data.frame containing the chromosome, start position and end position of each feature. 
The rownames must be the feature names or ids with the same type of identifier used in \code{data}.
}

\item{gc}{
Optional argument.Vector, matrix or data.frame containing the GC content
  of each feature. In case of giving a vector, the names of the
vector must be the feature names or ids  with the same type of
identifier used in \code{data}. If a matrix or a data.frame is provided, 
and it has two columns, it is expected that the feature names or ids are
in the first column and the GC content of the features in the
second. If it only has one column containing the GC content, the rownames of the
object must be the feature names or ids.
}

  \item{length}{
    Optional argument.Vector, matrix or data.frame containing the length
  of each feature. In case of giving a vector, the names of the
vector must be the feature names or ids  with the same type of
identifier used in \code{data}. If a matrix or a data.frame is provided, 
and it has two columns, it is expected that the feature names or ids are
in the first column and the length of the features in the
second. If it only has one column containing the length, the rownames of the
object must be the feature names or ids.
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
