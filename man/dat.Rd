\name{Data_Exploration}
\alias{dat}

\title{
Exploration of expression data.
}
\description{
Take the expression data and the feature annotations to generate the results that will be used for the exploratory plots (\code{explo.plot}) or saved by the user to perform other analyses.
}
\usage{
dat(input, type = c("biodetection","cd","countsbio","GCbias","lengthbias","saturation"), 
                    k = 0, ndepth = 6, factor = NULL, norm = FALSE, refColumn = 1)
}

\arguments{
  \item{input}{
Object of eSet class with expression data and optional annotation.
}
  \item{type}{
Type of plot for which the data are to be generated. It can be one of: "biodetection","cd","countsbio","GCbias","lengthbias","saturation".
}
  \item{k}{
A feature is considered to be detected if the corresponding number of read counts is > k. By default, k = 0. This parameter is used by types "biodetection" and "saturation".
}
  \item{ndepth}{
Number of different sequencing depths to be simulated and plotted apart from the real depth. By default, ndepth = 6. This parameter is only used by type "saturation".
}
  \item{factor}{
If factor = NULL (default), the calculations are done for each sample independently. When the factor is specified, the calculations are done for each experimental condition. 
Samples within the same condition are summed up ("biodetection") or averaged and normalized by sequencing depth ("countsbio", "GCbias" and "lengthbias").
}
  \item{norm}{
To indicate if provided data are already normalized (TRUE) or they are raw data (FALSE), which is the default. This parameter is used by types "cd", "lengthbias", "GCbias" and "countsbio".
}
  \item{refColumn}{
Column number in input data that is taken as the reference sample to compute M values. This parameter is only used by type "cd".
}
}


\value{
\code{dat()} function returns an S4 object to be used by \code{explo.plot()} or to be converted into a more friendly formatted object by the \code{dat2save()} function.
}


\author{
Sonia Tarazona
}


\seealso{
       \code{\link{Biodetection}},\code{\link{CD}},\code{\link{CountsBio}},\code{\link{GCbias}},\code{\link{lengthbias}},\code{\link{Saturation}},\code{\link{readData}},\code{\link{addData}},\code{dat2save},\code{explo.plot}
     }

\examples{
## Load the input object with the expression data and the annotations
data(myCounts)

## Generating data for the plot "biodetection" and samples in columns 3 and 4 of expression data
mydata2plot = dat(mydata, type = "biodetection", k = 0)

## Generating the corresponding plot
explo.plot(mydata2plot, samples = c(3,4))

}
