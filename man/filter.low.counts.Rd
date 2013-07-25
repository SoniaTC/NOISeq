\name{FilterLowCounts}
\alias{filtered.data}

\title{
Methods to filter out low count features
}
\description{
Function to filter out the low count features according to three different methods.
}
\usage{
filtered.data(dataset, factor, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1)
}

\arguments{
  \item{dataset}{
Matrix or data.frame containing the expression values for each sample (columns) and feature (rows).
}
  \item{factor}{
Vector or factor indicating which condition each sample (column) in dataset belongs to.
}
  \item{norm}{
Logical value indicating whether the data are already normalized (TRUE) or not (FALSE).
}
  \item{depth}{
Sequencing depth of samples (column totals before normalizing the data). When norm=FALSE, depth does not need to be provided. 
}
\item{method}{
Method must be one of 1,2 or 3. 
Method 1 (CPM) removes those features that have an average expression per condition less than cpm value  
and a coefficient of variation per condition higher than cv.cutoff in all the conditions.
Method 2 (Wilcoxon) performs a Wilcoxon test per condition and feature where in the null hypothesis the median expression is 0 and
in the alternative the median is higher than 0. Those features with p-value greater than 0.05 in all the conditions are removed.
Method 3 (Proportion test) performs a proportion test on the counts per condition and feature (or pseudo-counts if data were normalized) 
where null hypothesis is that the feature relative expression (count proportion) is equal to cpm/10^6 and higher than cpm/10^6 for the alternative.
Those features with p-value greater than 0.05 in all the conditions are removed.
}
\item{cv.cutoff}{
Cutoff for the coefficient of variation per condition to be used in method 1.
}
\item{cpm}{
Cutoff for the counts per million value to be used in methods 1 and 3.
}
}
\details{
}


\examples{

## Simulate some count data
datasim = matrix(sample(0:100, 2000, replace = TRUE), ncol = 4)

## Filtering low counts (method 1)
myfilt1 = filtered.data(datasim, factor = c("cond1", "cond1", "cond2", "cond2"), norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1)

## Filtering low counts (method 2)
myfilt2 = filtered.data(datasim, factor = c("cond1", "cond1", "cond2", "cond2"), norm = FALSE, method = 2)

## Filtering low counts (method 3)
myfilt3 = filtered.data(datasim, factor = c("cond1", "cond1", "cond2", "cond2"), norm = FALSE, method = 3, cpm = 1)

}



\author{
Sonia Tarazona
}

