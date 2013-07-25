\name{CD}
\alias{CD-class}
\alias{CD}
\alias{show,CD-method}
\alias{dat2save,CD-method}
\alias{explo.plot,CD-method}
\docType{class}
\title{CD class}

\description{
CD class generated from dat() function with type="cd". This object contains the cumulative distributions of reads for the samples.
}

\usage{
\S4method{explo.plot}{CD}(object, samples = NULL, ...)
\S4method{dat2save}{CD}(object)
}

\arguments{
  \item{object}{
    Object generated from \code{dat()} function.
  }
  \item{samples}{
    Samples or conditions to be plotted. If NULL, the two first samples are plotted because the plot for this object only admit a maximum of two samples.
  }
  \item{...}{
    Any argument from \code{par}.
  }
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:

  \code{dat}: List containing the following elements:

  \code{data2plot}: Data frame where the first column is a number between 0 and 100 and the following columns are the corresponding percentiles for each sample.
  \code{DiagnosticTest}: FDR adjusted pvalues for Kolmogorov-Smirnov tests comparing each pair of samples.
  
}

\section{Methods}{
This class has an specific \code{show} method in order to print Kolmogorov-Smirnov tests results 
and a \code{dat2save} method to save the relevant information in the object in a user-friendly way. It also has an
\code{explo.plot} method to plot the data contained in the object.
}

\author{Sonia Tarazona}

\keyword{classes}
