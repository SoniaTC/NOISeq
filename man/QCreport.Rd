\name{QCreport}
\alias{QCreport}

\title{
Quality Control report for expression data
}
\description{
Generate a report with the exploratory plots for count data that can be generated from the biological information provided. This report is designed to compare two samples or two experimental conditions. 
}
\usage{
QCreport(input, file = NULL, samples = NULL, factor = NULL)
}

\arguments{
  \item{input}{
Object of eSet class coming from \code{\link{readData}} function or other R packages such as DESeq.
}
  \item{file}{
String indicating the name of the PDF file that will contain the report. It should be in this format: "filename.pdf". The default name is like this: "QCreport_2013Sep26_15:58:16.pdf".
}
  \item{samples}{
Vector with the two samples to be compared in the report when "factor" is NULL. If "factor" is not NULL and has more than two levels, samples has to indicate the two conditions to be compared. It can be numeric or character (when names of samples or conditions are provided).
}
  \item{factor}{
If NULL, individual samples indicated in "samples" are compared. Otherwise, it should be a string indicating the factor containing the experimental conditions to be compared in the report.
}
}
\value{
A pdf file.
}

\author{
Sonia Tarazona
}



\examples{
## Load the input object from Marioni's data as returned by readData()
data(myCounts)

## Generate the report
QCreport(mydata, samples = NULL, factor = "Tissue")

}

\references{
Marioni, J.C. and Mason, C.E. and Mane, S.M. and Stephens, M. and Gilad, Y. (2008) RNA-seq: an assessment of
technical reproducibility and comparison with gene expression arrays. \emph{Genome Research}, \bold{18}:
1509--1517.
}
