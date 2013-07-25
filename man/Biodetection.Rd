\name{Biodetection}
\alias{Biodetection-class}
\alias{Biodetection}
\alias{show,Biodetection-method}
\alias{explo.plot,Biodetection-method}
\alias{dat2save,Biodetection-method}
\docType{class}
\title{Biodetection class}

\description{
Biodetection class generated from dat() function with type="biodetection". This object contains the percentage of each biological 
class (e.g. biotype) in the genome (i.e. in the whole set of features provided), the corresponding percentage detected by the sample 
and the percentage of the biotype within the sample.
}

\usage{
\S4method{explo.plot}{Biodetection}(object, samples = c(1, 2), ...)
\S4method{dat2save}{Biodetection}(object)
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
An object of this class contains an element (dat) which is a list with the following components:

  \code{genome}: Vector containing the percentage of features per biotype in the genome.
  
  \code{biotables}: List with as many elements as samples or conditions. Each element of the list contains the percentage of features in the genome per biotype detected in that sample or condition features per biotype and the percentage of detected features in the sample or condition per biotype.  
}

\section{Methods}{
This class has an specific \code{show} method in order to work and print
a summary of the elements which are contained and a \code{dat2save} method
to save the relevant information in an object cleanly. It also has an
\code{explo.plot} method to plot the data contained in the object.
}

\author{Sonia Tarazona}

\keyword{classes}
