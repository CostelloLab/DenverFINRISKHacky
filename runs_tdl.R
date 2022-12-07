# Base working directory
setwd("C:\\Users\\Daniel\\DenverFINRISKHacky\\")

# Code specific to TDL
source("funcs.R")

# Read in training data
train <- list()
train[["pheno"]] <- read.table("data\\DreamHF\\train\\pheno_training.csv", sep=",", header=TRUE, row.names=1)
train[["readcounts"]] <- read.table("data\\DreamHF\\train\\readcounts_training.csv", sep=",", header=TRUE, row.names=1)

test <- list()
test[["pheno"]] <- read.table("data\\DreamHF\\test\\pheno_test.csv", sep=",", header=TRUE, row.names=1)
test[["readcounts"]] <- read.table("data\\DreamHF\\test\\readcounts_test.csv", sep=",", header=TRUE, row.names=1)

library(survival)

train_raw <- train[["readcounts"]]
# Expand train_raw by collapsing to different strata
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="k")) )
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="p")) )
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="c")) )
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="o")) )
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="f")) )
train_raw <- rbind(train_raw, t(collapseMicrobiome(t(train[["readcounts"]]), strata="g")) )

train_x <- apply(train[["readcounts"]], MARGIN=1, FUN=inormal)
train_p <- train[["pheno"]]
train_p <- train_p[,-which(colnames(train_p) %in% c("Event", "Event_time"))]

test_raw <- test[["readcounts"]]
# Expand test_raw by collapsing to different strata
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="k")) )
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="p")) )
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="c")) )
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="o")) )
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="f")) )
test_raw <- rbind(test_raw, t(collapseMicrobiome(t(test[["readcounts"]]), strata="g")) )

test_x <- apply(test[["readcounts"]], MARGIN=1, FUN=inormal)
test_p <- test[["pheno"]]
test_p <- test_p[,-which(colnames(train_p) %in% c("Event", "Event_time"))]

train_y <- survival::Surv(event = train[["pheno"]]$Event, time = train[["pheno"]]$Event_time)
test_y <- survival::Surv(event = test[["pheno"]]$Event, time = test[["pheno"]]$Event_time)


# Save only essential objects
save(
	# Train
	train_raw, # Raw counts
	train_x, # inormalized counts
	train_p, # phenotype
	train_y, # Survival response
	
	# Test
	test_raw, # Raw counts
	test_x, # inormalized counts
	test_p, # phenotype
	test_y, # Survival response

	# Filename
	file="train_test.RData"
)

#train[["taxtable"]] <- read.table("data\\DreamHF\\train\\taxtable.csv", sep=",", header=TRUE)

# Inspect distributions

# Based on taxa
apply(train[["readcounts"]], MARGIN=1, FUN=\(x) quantile(x, probs=seq(0,1,by=0.1)))
quantile(apply(train[["readcounts"]], MARGIN=1, FUN=\(x) median(x)), probs=seq(0,1,by=.1))

# Based on patients
apply(train[["readcounts"]], MARGIN=2, FUN=\(x) quantile(x, probs=seq(0,1,by=0.1)))
quantile(apply(train[["readcounts"]], MARGIN=2, FUN=\(x) max(x)), probs=seq(0,1,by=.1))

#> quantile(apply(train[["readcounts"]], MARGIN=2, FUN=\(x) max(x)), probs=seq(0,1,by=.1))
#       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
#   5340.0   36996.0   56735.6   75377.6   98197.0  124672.0  157342.6  203113.8  267962.0  407615.2 3302165.0

# Odd maximum

library(ComplexHeatmap)

as.matrix(train[["pheno"]])

library(survival)

train[["x"]] <- train[["pheno"]]
#train[["y"]] <- survival::Surv(

par(mfrow=c(2,2))
boxplot(Age ~ Sex, data = train[["pheno"]], range=0)
boxplot(BodyMassIndex ~ Sex, data = train[["pheno"]], range=0)
boxplot(Smoking ~ Sex, data = train[["pheno"]], range=0)
boxplot(PrevalentHFAIL ~ Sex, data = train[["pheno"]], range=0)

ComplexHeatmap::Heatmap()



# Distribution of follow-up times

plot(y=1:nrow(train_p), x=train_p[,"Event_time"], col=1+train_p[,"Event"], pch=15+train_p[,"Event"], xlab="Time (relative to 2002)", ylab="Patient index")
abline(v=0, col="grey", lwd=2)
legend("bottomleft", col=1:2, pch=15:16, legend = c("Censored", "HF"))



plot(y=1:nrow(train_p), x=train_p[,"Event_time"], col=1+train_p[,"Event"], pch=15+train_p[,"Event"], xlab="Time (relative to 2002)", ylab="Patient index", xlim=c(0,16))
abline(v=0, col="grey", lwd=2)
legend("bottomleft", col=1:2, pch=15:16, legend = c("Censored", "HF"))


colnames(train_x)[grep("Streptococcus|streptococcus", colnames(train_x))]


# The S. bovis group includes S. equinus, S. gallolyticus, S. infantarius, and other closely related species; they are the nonenterococcal group D streptococci. 

streppos <- 
	c(
		"k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__",
		"k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_gallolyticus",
		"k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus_infantis"
	)

train_raw[streppos,]


grep("k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus", colnames(train_x), value=TRUE)

streppo_x <- apply(train_x[,grep("k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus", colnames(train_x), value=TRUE)], MARGIN=1, FUN=\(x){ sum(x, na.rm=TRUE) })


par(mfrow=c(1,2))

plot(y=train_x[,streppos[2]], x=train_p[,"Event_time"], col=1+train_p[,"Event"], pch=15+train_p[,"Event"], xlab="Time (relative to 2002)", ylab="Gallolyticus")
#legend("bottomleft", col=1:2, pch=15:16, legend = c("Censored", "HF"))

plot(y=train_x[,streppos[3]], x=train_p[,"Event_time"], col=1+train_p[,"Event"], pch=15+train_p[,"Event"], xlab="Time (relative to 2002)", ylab="Infantis")
#legend("bottomleft", col=1:2, pch=15:16, legend = c("Censored", "HF"))


coxph(train_y ~ train_x[,streppos[2]])

coxph(train_y ~ streppo_x)


train_p <- apply(train_p, 


train_p <- apply(train_p, MARGIN=2, FUN=\(x){ x[is.na(x)] <- median(x, na.rm=TRUE); scale(x) })

pca_p <- prcomp(train_p)
plot(pca_p$x[,1:2], pch=16, col=1+train_y[,2])

plot(x=train_p[,"Age"], y=train_p[,"BodyMassIndex"], pch=16, col=1+train_y[,2], xlab="Age (scaled)", ylab="BMI (scaled)")

