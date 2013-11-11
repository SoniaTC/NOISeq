\name{CountsBio}
\alias{CountsBio-class}
\alias{CountsBio}
\alias{show,CountsBio-method}
\alias{explo.plot,CountsBio-method}
\alias{dat2save,CountsBio-method}
\docType{class}
\title{CountsBio class}

\description{
CountsBio class generated from dat() function with
type="countsbio". This object contains the count distribution for each biological group 
and also the percentage of features with counts per million higher than 0, 1, 2, 5 or 10,
per each sample independently and in at least one of the samples (total).
}

\usage{
\S4method{explo.plot}{CountsBio}(object, samples = c(1,2), toplot = "global", plottype = c("barplot", "boxplot"),...)
\S4method{dat2save}{CountsBio}(object)
}

\arguments{
  \item{object}{
    Object generated with \code{dat()} function.
  }
  \item{toplot}{
    This parameter indicates which biological group is to be plotted.
    It may be a number or a text with the name of the biological group. 
    If toplot=1 (or "global"), a global plot with all the biological groups will be generated.

  }
  \item{samples}{
    Samples or conditions to be plotted. If NULL, the two first samples are plotted because the plot for this object only admit a maximum of two samples.
  }
  \item{plottype}{
    Type of plot to be generated for "countsbio" data. If "barplot", the plot indicates the percentage of features with counts per millior higher than 0, 1, 2, 5 or 10 counts or less. Above each bar, the sequencing depth (million reads) is shown. 
If "boxplot", a boxplot is drawn per sample or condition showing the count distribution for features with more than 0 counts. 
Both types of plot can be obtained for all features ("global") or for a specified biotype (when biotypes are available).
  }
  \item{...}{
    Any argument from \code{par}.
  }
}

\section{Slots/List Components}{
Objects of this class contain a list (\code{dat}) with the following components:
  
  \code{result}: Matrix containing the expression data for all the detected features and all samples or conditions.

  \code{bionum}: Vector containing the number of detected features per biological group (global indicates the total).
  
  \code{biotypes}: Vector containing the biological group (biotype) for each detected feature.

  \code{summary}: List with as many elements as number of biotypes and an additional element with the global information (for all features).
Each element is a data frame containing for each sample or condition the number of features with 0 counts, 1 count or less, 
2 counts or less, 5 counts or less and 10 counts or less, more than 10 counts, the total number of features and the sequencing depth.  
}

\section{Methods}{
This class has an specific \code{show} method in order to work and print
a summary of the elements which are contained and a \code{dat2save} method
to save the relevant information in an object cleanly. It also has an
\code{explo.plot} method to plot the data contained in the object.
}

\author{Sonia Tarazona}

\keyword{classes}
