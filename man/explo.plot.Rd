\name{Exploratory_Plots}
\alias{explo.plot}

\docType{methods}

\title{
Exploratory plots for expression data.
}

\description{
Standard generic function. Different types of plots showing the biological classification for detected features, the expression distribution across samples or biological groups, the detection of technical bias such as length, GCcontent or RNA composition, the dependence of expression on sequencing depth, etc.
}


\usage{
explo.plot(object, ...)
}

\arguments{
  \item{object}{
Object generated with \code{dat()} function.
}

  \item{...}{
    Any argument from \code{par}.
  }
}

\value{
The explo.plot() function takes the object generated by dat() and draws the corresponding plot.
}


\author{
Sonia Tarazona
}


\seealso{
       \code{\link{Biodetection}},\code{\link{CD}},\code{\link{CountsBio}},\code{\link{GCbias}},\code{\link{lengthbias}},\code{\link{Saturation}}, \code{\link{PCA}}, \code{\link{readData}}, \code{\link{addData}}, \code{\link{dat}}.
     }

\examples{
## Load the input object with the expression data and the annotations
data(myCounts)

## Generating data for the plot "biodetection" and samples in columns 3 and 4 of expression data
mydata2plot = dat(mydata, type = "biodetection", k = 0)

## Generating the corresponding plot
explo.plot(mydata2plot)

}
