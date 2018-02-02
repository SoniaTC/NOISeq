\name{CD}
\alias{CD-class}
\alias{CD}
\alias{show,CD-method}
\alias{dat2save,CD-method}
\alias{explo.plot,CD-method}
\docType{class}
\title{CD class}

\description{
CD class generated from dat() function with type="cd". This object contains the distributions of log-fold changes (M values) between each of the samples and a reference sample as well as confidence intervals for the median of these distributions that are used to detect a potential RNA composition bias in the data.
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
    Samples or conditions to be plotted. If NULL, the twelve first samples are plotted because the plot for this object only admit a maximum of twelve samples.
  }
  \item{...}{
    Any argument from \code{par}.
  }
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:

  \code{dat}: List containing the following elements:

  \code{data2plot}: Data frame where each column contains the M values obtained as the log2-ratio of each sample against the reference sample.
  \code{refColumn}: Column number in input data that is taken as the reference sample.
  \code{DiagnosticTest}: Data frame that contains the lower and upper limits of the confidence intervals for the median of M values per each sample. The last column indicates if the diagnostic test for that sample has been passed or failed (so normalization has to be applied).
  
}

\section{Methods}{
This class has an specific \code{show} method in order to show the confidence intervals for the M median 
and a \code{dat2save} method to save the relevant information in the object in a user-friendly way. It also has an
\code{explo.plot} method to plot the data contained in the object.
}

\author{Sonia Tarazona}

\keyword{classes}
