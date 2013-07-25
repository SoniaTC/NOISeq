\name{Output}
\alias{Output-class}
\alias{Output}
\alias{show,Output-method}
\docType{class}
\title{Output class of NOISeq}

\description{
Output object containing the results from differential expression analysis by \code{noiseq} or \code{noiseqbio}.
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{comparison}: String indicating the two experimental conditions being compared and the sense of the comparison.

  \code{factor}: String indicating the factor chosen to compute the differential expression.
  
  \code{k}: Value to replace zeroes in orden to avoid inderminations when computing logarithms. 
  
  \code{lc}: Correction factor for length normalization. Counts are divided by length^lc.
  
  \code{method}: Normalization method chosen. It can be one of "rpkm" (default), "uqua" (Upper Quartile), "tmm" 
		 (Trimmed Mean of M) or "n" (no normalization).

  \code{replicates}: Type of replicates: "technical" for technical replicates and "biological" for biological ones.

  \code{results}: R data frame containing the differential expression results, where each row corresponds to a feature. The columns are: Expression values for each condition to be used by \code{noiseq} 
or \code{noiseqbio} (the columns names are the levels of the factor); differential expression statistics (columns "M" and "D" for \code{noiseq} or "theta" for \code{noiseqbio}); probability of differential expression ("prob"); "ranking", which is a summary statistic of "M" and "D" values equal to -sign(M)*sqrt(M^2 + D^2), than can be used for instance in gene set enrichment analysis (only when \code{noiseq} is used); "length" and "GC" of each feature (if provided); chromosome where the feature is ("Chrom"), if provided; start and end position of the feature within the chromosome ("GeneStart", "GeneEnd"), if provided.  

  \code{nss}: Number of samples to be simulated for each condition (only when there are not replicates available). 

  \code{pnr}: Percentage of the total sequencing depth to be used in each simulated replicate (only when there are not replicates available). If, for instance, pnr = 0.2 , each simulated replicate will have 20\% of the total reads of the only available replicate in that condition.

  \code{v}: Variability of the size of each simulated replicate (only used by NOISeq-sim).

}

\section{Methods}{
This class has an specific \code{show} method in order to work and print a summary of the elements which are contained.
}

\author{Sonia Tarazona}

\keyword{classes}
