\name{Data_Exploration}
\alias{dat}

\title{
Exploration of expression data.
}
\description{
Take the expression data and the feature annotations to generate the results that will be used for the exploratory plots (\code{explo.plot}) or saved by the user to perform other analyses.
}
\usage{
dat(input, type = c("biodetection","cd","countsbio","DLbio","saturation"), 
                    selection=c(1,2), k = 0, ndepth = 5, newdetections = TRUE)
}

\arguments{
  \item{input}{
Object of eSet class with expression data and optional annotation.
}
  \item{type}{
Type of plot for which the data are to be generated. It can be one of: "biodetection","cd","countsbio","DLbio","saturation".
}newdetections = TRUE
  \item{selection}{
Vector containing the number or names of the columns (samples) from the expression data to be plotted. 
Depending on the chosen \code{type}, a different maximum number of samples is allowed. 
For "biodetection", the number of samples may be 1 or 2; for "cd", 2; for "countsbio", "DLbio" and "saturation", 
the data may be generated for all the available samples.  
}
  \item{k}{
A feature is considered to be detected if the corresponding number of read counts is > k. By default, k = 0. This parameter is used by all types except "cd".
}
  \item{ndepth}{
Number of different sequencing depths to be simulated and plotted. By default, ndepth = 5. This parameter is only used by types "countsbio","DLbio" and "saturation".
}
  \item{newdetections}{
If TRUE, a second Y-axis is drawn for new detections per million of new sequenced reads. This parameter is only used by type "saturation". By default, newdetections = TRUE. If samples to be plotted are more than two, this option is disabled. 
}
}


\value{
\code{dat()} function returns an S4 object to be used by \code{explo.plot()} or to be converted into a more friendly formatted object by the \code{dat2save()} function.
}


\author{
Sonia Tarazona
}


\seealso{
       \code{\link{readData}},\code{\link{addData}},\code{dat2save},\code{explo.plot}
     }

\examples{
## Load the input object with the expression data and the annotations
data(myCounts)

## Generating data for the plot "biodetection" and samples in columns 3 and 4 of expression data
mydata2plot = dat(mydata, type = "biodetection", k = 0, selection = c(3,4))

## Generating the corresponding plot
explo.plot(mydata2plot)

}
