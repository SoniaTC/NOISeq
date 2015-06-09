ARSyNseq <- function(data, factor = NULL, batch = FALSE, norm = "rpkm", logtransf = FALSE, Variability = 0.75, beta = 2)
{

	#   data:     A Biobase's eSet object created with the readData function
	#   factor:   Column name choosen from factors argument of data. 
	#             When it is NULL, all the factors (1,2 or 3) specified in data are considered.
	#   batch:    TRUE when the factor is an identified batch effect. This option can be run only with 1 factor.
	#   norm:     Normalization method. It can be one of "rpkm" (default), "uqua"
	#             upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).

	#------- Parameters for PCs selection: ----------------------
	#Variability: Parameter for PCs selection of the ANOVA models effects.
	#    beta:    Parameter for PCs selection of the residual model. 

	#------- Only used for 2 or 3 factors:---------------------
	#    Join:    Logical to indicate whether interaction Factor1xFactor2 must be analysed jointly with the second factor. 
	#Interaction: Logical to indicate whether interaction/s between factors should be analyzed.

	#------- Arguments for normalization: ----------------------
	#    k:       Counts equal to 0 are replaced by k. By default, no replacement, k=0.  
	#    lc:      Length correction is done by dividing expression by length^lc. By default, no correction, lc = 1.

	#----------------------------------
	# --- Compute Inputs for ARSyN
	#----------------------------------


	Join = TRUE
	Interaction = TRUE


	 dat <- as.matrix(assayData(data)$exprs)
	 long <- featureData(data)@data[rownames(dat), "Length"]
	 if (is.null(long)) long = 1000

	 if (norm == "rpkm") {     
	   dat <- rpkm(dat, long = long, k = 0, lc = 1)
	  }

	 if (norm == "uqua") {
	   dat <- uqua(dat, long = long, lc = 1, k = 0)
	    }
	   
	 if (norm == "tmm") {
	    dat <- tmm(dat, long = long, lc = 1, k = 0)      
	    }

	 if (!logtransf)   dat <- log(dat + 1)
	 
	X <- t(dat) 		#conditions x genes


	#----------------------------------

	 if(is.null(factor))
	{
	 Covariates <- t(pData(data))
	 Num.factors <- nrow(Covariates)
	 labels.factors <- rownames(Covariates)
	 Design <- list(NULL,NULL,NULL)
	 for (i in 1:Num.factors)
	 {
		x <- as.character(Covariates[i,])
		Design[[i]] <- make.ASCA.design(x)
	 } 
	}else{
	 Covariates <- pData(data)[,factor]
	 Num.factors <- 1
	 labels.factors <- factor
	 Design <- list(NULL,NULL,NULL)
		x <- as.character(Covariates)
		Design[[1]] <- make.ASCA.design(x)
	 }

	####################################
	### --- Execute ASCAmodel 
	####################################

	 my.asca <- ARSyNmodel(Factors=Num.factors,X=X,Designa=Design[[1]],Designb=Design[[2]],Designc=Design[[3]],Join=Join,Interaction=Interaction,Variability=Variability,beta=beta)

	#################################### 
	### --- Writing filtered matrix 
	####################################

	 X.filtered <- X
	 M<-length(my.asca)-1

	if(!batch)
	{
	# for (i in 1:(M-1))
	#  {
	#	X.filtered <- X.filtered-my.asca[[i]]$E 
	#  }
		X.filtered <- X.filtered-my.asca[[M]]$TP
	}

	if(batch)
	{
		X.filtered <- X.filtered-my.asca[[1]]$TP 
	}

	data.filtered <- t(X.filtered)

	 if (!logtransf)  data.filtered <- exp(data.filtered)+1
	                     
	exprs(data) = data.filtered

	return(data)
}

