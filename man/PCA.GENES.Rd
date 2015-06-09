\name{PCA.GENES}
\alias{PCA.GENES}

\title{
Principal Component Analysis
}
\description{
Computes a Principal Component Analysis on any data matrix.
}
\usage{
PCA.GENES(X)
}

\arguments{
  \item{X}{
		Matrix or data.frame with variables (e.g. genes) in columns and observations (e.g. samples) in rows.
	}
}

\examples{

	## Simulate data matrix with 500 variables and 10 observations
	datasim = matrix(sample(0:100, 5000, replace = TRUE), nrow = 10)

	## PCA
	myPCA = PCA.GENES(datasim)

	## Extracting the variance explained by each principal component
	myPCA$var.exp

}


\author{
Maria Jose Nueda
}

