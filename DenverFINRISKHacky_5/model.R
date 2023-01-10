# Model to be submitted to DREAM FINRISK Challenge
# DenverFINRISKHacky

###
#
# Load libraries
#
###
library(survival)
library(glmnet)
library(microbiome)
library(phyloseq)
library(mia)
library(TreeSummarizedExperiment)
library(ecodist)
library(vegan)

###
#
# Define main model function
#
###
model <- function(
	# Clinical data matrix for both training and testing cohort
	train_clin,
	test_clin,
	# Microbiome (raw) read counts for both training and testing cohort
	train_micr,
	test_micr,
	# phyloseq-objects for the train and test cohorts
	train_phyloseq,
	test_phyloseq,
	# Which model submission version to use 
	v,
	# Additional parameters
	seed = 1, # Random seed for reproducibility
	...
){
	# Start constructing the output df which will be output as csv
	output_temp <- data.frame(SampleID = rownames(test_clin))
	output_final <- data.frame(SampleID = rownames(test_clin), Score = 0)

	# Discard negative or zero event times
	discard <- which(train_clin$Event_time <= 0)
	train_clin <- train_clin[-discard,]
	otu_table(train_phyloseq) <- otu_table(train_phyloseq)[,-discard]
	train_micr <- train_micr[,-discard]

	# Construct an ensemble training df
	ensemble_temp <- data.frame(SampleID = rownames(train_clin), Event = train_clin[,"Event"], Event_time = train_clin[,"Event_time"])

	# Naive median imputation to ages to remove NAs
	train_clin <- as.data.frame(imputationNaive(train_clin))
	test_clin <- as.data.frame(imputationNaive(test_clin))
	
	# Squared & shifted age
	train_clin[,"AgeShift2"] <- sqshift(train_clin[,"Age"])
	train_clin[,"AgeNlogNshift"] <- nlognshift(train_clin[,"Age"])

	test_clin[,"AgeShift2"] <- sqshift(test_clin[,"Age"])
	test_clin[,"AgeNlogNshift"] <- nlognshift(test_clin[,"Age"])

	# Create shifted & z-scored clinical variables, with all pairwise interactions incorporated
	train_clin2a <- apply(train_clin[,-which(colnames(train_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=shift)
	train_clin2b <- apply(train_clin[,-which(colnames(train_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=zscale)
	colnames(train_clin2a) <- paste0("s_", colnames(train_clin2a))
	colnames(train_clin2b) <- paste0("z_", colnames(train_clin2b))
	train_clin2 <- cbind(train_clin2a, train_clin2b, interact.all(train_clin2a), interact.all(train_clin2b))

	test_clin2a <- apply(test_clin[,-which(colnames(test_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=shift)
	test_clin2b <- apply(test_clin[,-which(colnames(test_clin) %in% c("Event", "Event_time"))], MARGIN=2, FUN=zscale)
        colnames(test_clin2a) <- paste0("s_", colnames(test_clin2a))
        colnames(test_clin2b) <- paste0("z_", colnames(test_clin2b))	
	test_clin2 <- cbind(test_clin2a, test_clin2b, interact.all(test_clin2a), interact.all(test_clin2b))

	# Omit combinations that produce NaN
	train_clin2 <- train_clin2[,-which(apply(train_clin2, MARGIN=2, FUN=\(x){ any(!is.finite(x)) }))]
	# Keep same columns
	test_clin2 <- test_clin2[,colnames(train_clin2)]

	# Construct response Surv and remove it from the pheno data
	train_y <- survival::Surv(time = train_clin$Event_time, event = train_clin$Event)	

	# Microbiome diversity metrics

	# Alpha
	print("Calculating microbiome::alpha for train cohort...")	
	train_alpha <- train_phyloseq |>
		microbiome::alpha() |>
		imputationNaive()

	print("Calculating microbiome::alpha for test cohort...")	
	test_alpha <- test_phyloseq |>
		microbiome::alpha() |>
		imputationNaive()

	# Introduce diversity/entropy metric interactions
	train_alpha <- cbind(train_alpha, interact.all(train_alpha))
	test_alpha <- cbind(test_alpha, interact.all(test_alpha))

	# Beta, examine PCoA vectors up to 10th
	print("Calculating beta diversity for train cohort using mia...")
	print("mia...")
	train_tse <- train_phyloseq |>
		mia::makeTreeSummarizedExperimentFromPhyloseq() |>
		mia::agglomerateByRank(x = _, rank = "Genus") |>
		mia::transformSamples(x = _, method = "relabundance") |>
		mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr")

	print("vegan & ecodist...")
	train_beta <- t(assays(train_tse)$relabundance) |>
		vegan::vegdist(x = _, method = "bray") |>
		ecodist::pco(x = _) |>
		(\(x) { x <- x$vectors[,1:10]; colnames(x) <- paste0("PCoA", 1:10); x })()

	print("Calculating beta diversity for test cohort using mia...")
	print("mia...")
        test_tse <- test_phyloseq |>
                mia::makeTreeSummarizedExperimentFromPhyloseq() |>
                mia::agglomerateByRank(x = _, rank = "Genus") |>
                mia::transformSamples(x = _, method = "relabundance") |>
                mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr")

	print("vegan & ecodist...")
        test_beta <- t(assays(test_tse)$relabundance) |>
                vegan::vegdist(x = _, method = "bray") |>
                ecodist::pco(x = _) |>
                (\(x) { x <- x$vectors[,1:10]; colnames(x) <- paste0("PCoA", 1:10); x })()

        # Collapse raw reads on suitable phyla level, filter, inverse normal rank transformation, and pick correct columns
        # genus
        print("Collapsing genus...")
        train_g <- collapseMicrobiome(x = t(train_micr), strata = "g") |>
                filterMicrobiome() |>
                (\(x) { apply(x, MARGIN=2, FUN=inormal) })()
        test_g <- collapseMicrobiome(x = t(test_micr), strata = "g") |>
                (\(x) { x[,colnames(train_g)] })() |>
                (\(x) { apply(x, MARGIN=2, FUN=inormal) })()
        # family
        print("Collapsing family...")
        train_f <- collapseMicrobiome(x = t(train_micr), strata = "f") |>
                filterMicrobiome() |>
                (\(x) { apply(x, MARGIN=2, FUN=inormal) })()
        test_f <- collapseMicrobiome(x = t(test_micr), strata = "f") |>
                (\(x) { x[,colnames(train_f)] })() |>
                (\(x) { apply(x, MARGIN=2, FUN=inormal) })()

	# Test mia/microbiome packages' suggested approaches
	print("Using previously generated TreeSummarizedExperiment-objects to extract variables on multiple levels...")
	# Little piping function
	pip <- function(phylo, level){
		phylo |>
                mia::makeTreeSummarizedExperimentFromPhyloseq() |>
                mia::agglomerateByRank(x = _, rank = level) |>
                mia::transformSamples(x = _, method = "relabundance") |>
                mia::transformSamples(x = _, abund_values = "relabundance", pseudocount = 1, method = "clr", name = "clr_transformation") |>
                (\(x) { assay(x, "clr_transformation") })() 
	}
	# Training relative abundances
	print("Training data relative abundances...")
	train_relabus <- t(rbind(
		#pip(train_phyloseq, level = "Genus"),
		pip(train_phyloseq, level = "Family"),
		pip(train_phyloseq, level = "Order"),
		pip(train_phyloseq, level = "Class"),
		pip(train_phyloseq, level = "Phylum")
	))
	# Test relative abundances
	print("Test data relative abundances...")
	test_relabus <- t(rbind(
		#pip(test_phyloseq, level = "Genus"),
		pip(test_phyloseq, level = "Family"),
		pip(test_phyloseq, level = "Order"),
		pip(test_phyloseq, level = "Class"),
		pip(test_phyloseq, level = "Phylum")
	))

	## Modular modelling of risk

	# A risk module capturing age effect, adjusted by sex (females tend to have more favorable outcome)
	# Using various shapes; linear, squared, and n*log(n) after shifting, with interactions to Sex
	module_agesex <- function(
		train, # Training clinical data
		test, # Testing clinical data
		debug = FALSE, # If debug output is needed
		...
	){
		cph <- coxph(Surv(event = Event, time = Event_time) ~ Sex*Age + Sex*AgeShift2 + Sex*AgeNlogNshift, data = train)
		pred <- predict(cph, newdata = test)
		pred
	}
	# Collapse microbiome to suitable phyla level, inverse normal rank transformation, and utilize glmnet 10-fold CV
	# ... or other generic use of LASSO such as odd combinations of metadata
	module_glmnet <- function(
		trainx,
                trainy,
		testx,
		debug = FALSE, # If debug output is needed		
		...
	){
		# Omit samples with NAs
		trainx <- trainx[which(!is.na(trainy)),]
		trainy <- trainy[which(!is.na(trainy))]

		# Fit, CV, predict
		fit <- glmnet(x = as.matrix(trainx), y = trainy, family = "cox")
		cv <- cv.glmnet(x = as.matrix(trainx), y = trainy, family = "cox")

		print("lambda.1se")			
		print(predict(fit, s = cv$lambda.1se, type = "nonzero"))
		print("lambda.min")
		print(predict(fit, s = cv$lambda.min, type = "nonzero"))

		if(v >= 5){
			pred <- predict(fit, newx = as.matrix(testx), s = cv$lambda.min, type = "response")
		}else{
			pred <- predict(fit, newx = as.matrix(testx), s = cv$lambda.1se, type = "response")
		}

		pred[,1]
	}

	# Compute different modules
	set.seed(seed)
	
	# Part Ia: Training data, individual modules
	print("Pt Ia")	
	print("module_agesex")
	ensemble_temp[,"module_agesex"] <- module_agesex(train = train_clin, test = train_clin)
	print("module_metamix")
	ensemble_temp[,"module_metamix"] <- module_glmnet(trainx = train_clin2, trainy = train_y, test = train_clin2)
	print("module_genus_glmnet")
	ensemble_temp[,"module_genus_glmnet"] <- module_glmnet(trainx = train_g, trainy = train_y, test = train_g)
	print("module_family_glmnet")
	ensemble_temp[,"module_family_glmnet"] <- module_glmnet(trainx = train_f, trainy = train_y, test = train_f)
	print("module_alpha_glmnet")
	ensemble_temp[,"module_alpha_glmnet"] <- module_glmnet(trainx = train_alpha, trainy = train_y, test = train_alpha)
	print("module_beta_glmnet")
	ensemble_temp[,"module_beta_glmnet"] <- module_glmnet(trainx = train_beta, trainy = train_y, test = train_beta)
	print("module_relabus_glmnet")
	ensemble_temp[,"module_relabus_glmnet"] <- module_glmnet(trainx = train_relabus, trainy = train_y, test = train_relabus)

	print("Ensemble head")
	print(head(ensemble_temp))

	# Part Ib: Test data, individual modules
	print("Pt Ib")
	print("module_agesex")
	output_temp[,"module_agesex"] <- module_agesex(train = train_clin, test = test_clin)
	print("module_metamix")
	output_temp[,"module_metamix"] <- module_glmnet(trainx = train_clin2, trainy = train_y, test = test_clin2)
	print("module_genus_glmnet")
	output_temp[,"module_genus_glmnet"] <- module_glmnet(trainx = train_g, trainy = train_y, test = test_g)
	print("module_family_glmnet")
	output_temp[,"module_family_glmnet"] <- module_glmnet(trainx = train_f, trainy = train_y, test = test_f)
	print("module_alpha_glmnet")
	output_temp[,"module_alpha_glmnet"] <- module_glmnet(trainx = train_alpha, trainy = train_y, test = test_alpha)
	print("module_beta_glmnet")
	output_temp[,"module_beta_glmnet"] <- module_glmnet(trainx = train_beta, trainy = train_y, test = test_beta)
	print("module_relabus_glmnet")
	output_temp[,"module_relabus_glmnet"] <- module_glmnet(trainx = train_relabus, trainy = train_y, test = test_relabus)
		
	print("Temp output head")
	print(head(output_temp))

	# Part II: Find coefficients that maximize ensemble modules' linear sum for coxph in training data
	print("Pt II")	

	# Submissions 1-3 were not penalized Cox ensembles
	if(v < 4){
		print("Regularized derived features combined into Cox PH")
		ensemble_cox <- coxph(Surv(event = Event, time = Event_time) ~ module_agesex + module_metamix + module_genus_glmnet + module_family_glmnet + module_alpha_glmnet + module_beta_glmnet, data = ensemble_temp)
		print("Identified ensemble coefficients (coxph)")
		print(summary(ensemble_cox))

		# Part III: Construct the predicted test data score as a combination of weighted ensemble components
		# Combine modules to final output
		print("Pt III")
		output_final[,"Score"] <- predict(ensemble_cox, newdata = output_temp)	
	# Submissions 4+ testing penalized Cox ensembles (nested basically); sub 4 is more conservative with lambda.1se, sub 5 is lambda.min
	}else{
		print("Pt III")
		print("Nested regularization")
		w1 <- grep("module", colnames(ensemble_temp), value = TRUE)
		w2 <- grep("module", colnames(output_temp), value = TRUE)
		output_final[,"Score"] <- module_glmnet(trainx = ensemble_temp[,w1], trainy = train_y, test = output_temp[,w2])
	}		
		
	# Scale within [0,1] as instructed
	# From instructions: "The predictions have to be between 0 and 1, with larger numbers 
	# being associated with higher likelihood of having HF"
	output_final[,"Score"] <- outscale(output_final[,"Score"])

	# Return scores df
	output_final		
}

###
#
# Helper functions
#
###

# Scale risk scores between [0,1]
outscale <- \(x){ (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) }
# Just left shift variables to start at zero
shift <- \(x){ (x-min(x, na.rm=TRUE)) }
# Scale variables with z transformation
zscale <- \(x){ (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE) }
# Left shifted and squared
sqshift <- \(x){ (x-min(x, na.rm=TRUE))^2 }
# Left shifted and n*log(n+1) multiplied
nlognshift <- \(x){ x <- x-min(x, na.rm=TRUE); x*log(x+1) }

# Modified helper reader function modified from microbiome::read_csv2phyloseq
# package 'microbiome' in bioconductor, also using 'phyloseq'
#
# Taxonomy table provided in 'taxtable.csv' missing row names, so needed to hack the function a bit
#
csv2phylo <- function (otu.file = NULL, taxonomy.file = NULL, metadata.file = NULL, sep = ",") 
{
    s.meta <- read.csv(metadata.file, row.names = 1, check.names = FALSE, sep = sep)
    s.sampledata <- phyloseq::sample_data(s.meta)
    s.otu <- read.csv(otu.file, row.names = 1, check.names = FALSE, sep = sep)
    if (any(rownames(s.otu) %in% rownames(s.meta))) {
        s.otu <- t(s.otu)
    }
    s.otu.table <- otu_table(s.otu, taxa_are_rows = TRUE)
    
    s.tax_table <- phyloseq::tax_table(as.matrix(read.csv(taxonomy.file, sep = sep)))
    rownames(s.tax_table) <- rownames(s.otu)

    if (!all(rownames(s.tax_table) == s.tax_table[, ncol(s.tax_table)])) {
        s.tax_table <- cbind(s.tax_table, OTU = rownames(s.tax_table))
        s.tax_table <- phyloseq::tax_table(s.tax_table)
    }
    pseq <- phyloseq::merge_phyloseq(s.otu.table, s.tax_table, s.sampledata)
    return(pseq)
}

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

## Take all pairwise interactions and create derived variables
#
# From old code: https://github.com/Syksy/ePCR/blob/master/R/int.R
interact.all <- function(input){
	output <- do.call("cbind", lapply(1:ncol(input), FUN=function(z){ 
		do.call("cbind", lapply(z:ncol(input), FUN=function(x){
			tmp <- data.frame(input[,z] * input[,x])
			colnames(tmp)[1] <- paste(colnames(input)[z], "x", colnames(input)[x], sep="") 
			tmp
		}))
	}))
	output
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
	nams <- apply(mat[,1:which(colnames(mat) == strata),drop=FALSE], MARGIN=1, FUN=\(x){ paste0(x, collapse=";") })
	# Collapse matrix based on these indices
	cmat <- do.call("cbind", by(t(x), INDICES=nams, FUN=\(q){ apply(q, MARGIN=2, FUN=collapseFUN) }))

	# Return collapsed matrix
	cmat
}

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

#' Impute median values to missing slots based on columns
imputationNaive <- function(
	x, # Data matrix 
	FUN = \(q){ q[is.na(q)] <- median(q, na.rm=TRUE); q },
	MARGIN = 2,
	...
){
	dimn <- dimnames(x)
	x <- apply(x, MARGIN=MARGIN, FUN=FUN)
	dimnames(x) <- dimn
	x
}

###
#
# Run script inside the environment
#
###

# Submission number
subv <- 5

args=(commandArgs(TRUE))
PARAM <- list()
# Path of input folder
PARAM$folder.R <- paste0(args[1]) 

# Once model actually runs
# Create team folder
dir.create(file.path(PARAM$folder.R, paste0("DenverFINRISKHacky_", subv)))
# Create team output folder
dir.create(file.path(PARAM$folder.R, paste0("DenverFINRISKHacky_", subv),"output"))
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
PARAM$folder.result <- paste0(PARAM$folder.data, paste0("DenverFINRISKHacky_", subv), "/output/")

# Pheno data (both meta as well as response)
test_p <- read.csv(file = paste0(PARAM$folder.data, "test/pheno_test.csv"), row.names=1)
train_p <- read.csv(file = paste0(PARAM$folder.data, "train/pheno_training.csv"), row.names=1)

# Read count raw data
test_r <- read.csv(file = paste0(PARAM$folder.data, "test/readcounts_test.csv"), row.names=1)
train_r <- read.csv(file = paste0(PARAM$folder.data, "train/readcounts_training.csv"), row.names=1)

# Read taxa along with raw reads and metadata into a phyloseq object
test_phylo <- csv2phylo(
	otu.file=paste0(PARAM$folder.data, "test/readcounts_test.csv"), 
	taxonomy.file=paste0(PARAM$folder.data, "test/taxtable.csv"),
	metadata.file=paste0(PARAM$folder.data, "test/pheno_test.csv")
)
train_phylo <- csv2phylo(
	otu.file=paste0(PARAM$folder.data, "train/readcounts_training.csv"), 
	taxonomy.file=paste0(PARAM$folder.data, "train/taxtable.csv"),
	metadata.file=paste0(PARAM$folder.data, "train/pheno_training.csv")
)

# Run model and obtain scores result
res <- model(
	train_clin = train_p, # Training clinical metadata
	test_clin = test_p, # Test clinical metadata
	train_micr = train_r, # Training microbiome raw read counts
	test_micr = test_r, # Testing microbiome raw read counts
	train_phyloseq = train_phylo, # Training phyloseq object
	test_phyloseq = test_phylo, # Train phyloseq object
	v = subv # Submission version
)

# Write the resulting scores.csv
write.csv(res, file=paste0(PARAM$folder.result, "scores.csv"), quote=FALSE, row.names=FALSE)

# Just in case for debugging
print("args provided:")
print(args)

