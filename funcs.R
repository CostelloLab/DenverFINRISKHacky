## From: https://gist.github.com/variani/d6a42ac64f05f8ed17e6c9812df5492b
#
#' Inverse normal transformation (INT)
#'
#' https://www.biostars.org/p/80597/
#' See the supplement of Yang et al. Nature 2012. 
#'
#' @example
#' x1 <- 5:1
#' inormal(x1)
#'
#' x2 <- c(NA, 5:1, 5:1, NA)
#' inormal(x2)
inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}


## Collapse microbiome data to a desired level
#'
#' Take for example: "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_gallolyticus"
#' ... could be collapsed to levels k, p, c, o, f, g binding all rows underneath together
#' Presumably by summing for aligned read counts.
#'
collapseMicrobiome <- function(
	x, # Input matrix
	stratFUN = colnames, # The extractor function for strata, by default column names
	sep = ";", # Separator for the various strata
	strata = c("k", "p", "c", "o", "f", "g"), # Chosen level to collapse to; 's' == 'species' is the lowest level identified, thus omitted. Can also be integer, with k==1, p==2, c==3, ...
	collapseFUN = sum, # Function to apply over rows to collapse, applied across columns within chosen strata
	...
){
	stratas <- c("k", "p", "c", "o", "f", "g")
	if(is.integer(strata) | is.numeric(strata)){
		strata <- stratas[strata]
	}else if(!any(strata %in% stratas)){
		stop("Invalid strata; please indicate letter 'k', 'p', 'c', 'o', 'f', or 'g' or an integer corresponding to these in increasing order")
	}
	
	# Collapse matrix
	mat <- do.call("rbind", lapply(stratFUN(x), FUN=\(q){ strsplit(q, ";")[[1]] }))
	colnames(mat) <- c("k", "p", "c", "o", "f", "g", "s")
	
	# Collapsed names per columns
	nams <- apply(mat[,1:which(colnames(mat) == strata)], MARGIN=1, FUN=\(x){ paste0(x, collapse=";") })
	# Collapse matrix based on these indices
	cmat <- do.call("cbind", by(t(x), INDICES=nams, FUN=\(q){ apply(q, MARGIN=2, FUN=collapseFUN) }))
	
	#cmat
}
# Examples: 
# tmp_g <- collapseMicrobiome(t(train_raw), strata = "g")
# tmp_f <- collapseMicrobiome(t(train_raw), strata = "f")
# tmp_o <- collapseMicrobiome(t(train_raw), strata = "o")
# tmp_c <- collapseMicrobiome(t(train_raw), strata = "c")
# tmp_p <- collapseMicrobiome(t(train_raw), strata = "c")
# >  dim(tmp_g)
# [1] 3615 2019
# >  dim(tmp_f)
# [1] 3615  631
# >  dim(tmp_o)
# [1] 3615  342
# >  dim(tmp_c)
# [1] 3615  190
# >  dim(tmp_p)
#[1] 3615  190
# Sanity checking read counts:
#
#> train_raw[grep(grep("strep|Strep", colnames(tmp_g), value=TRUE)[21], rownames(train_raw)),1:5]
#                                                                                                                   Simulated_328 Simulated_1644 Simulated_1710 Simulated_1732 Simulated_1727
#k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__                                 0              0              0              0              0
#k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus_garvieae             5              2              0              0              3
#k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus_lactis              42              9              0              9             12
#k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus;s__Lactococcus_piscium              0              6              6              0              0
#> tmp_g[1:5,"k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Lactococcus"]
# Simulated_328 Simulated_1644 Simulated_1710 Simulated_1732 Simulated_1727 
#            47             17              6              9             15
# -> Matching


#' Relatively straight-forward filter based on threshold values and proportion of samples that need to qualify
filterMicrobiome <- function(
	x, # Data matrix x
	MARGIN = 2, # MARGIN to filter on (passed to apply etc; 1 = rows, 2 = columns)
	prop = 0.10, # Proportion of samples required to have at least <threshold> value for variable to be included
	threshold = 5, # Threshold value that samples ought to >= above
	...
){
	# 'which' rows/cols to keep
	w <- which(apply(x, MARGIN=MARGIN, FUN=\(q){
		(sum(q>=threshold) / length(q)) >= prop
	}))
	# Return filtered rows or columns
	switch(MARGIN, 
		"1" = x[w,,drop=FALSE],
		"2" = x[,w,drop=FALSE]
	)
}


