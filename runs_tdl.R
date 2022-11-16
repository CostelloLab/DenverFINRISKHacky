# Code specific to TDL
source("funcs.R")

# Base working directory
setwd("C:\\Users\\Daniel\\DenverFINRISKHacky\\")

# Read in training data
train <- list()
train[["pheno"]] <- read.table("data\\DreamHF\\train\\pheno_training.csv", sep=",", header=TRUE, row.names=1)
train[["readcounts"]] <- read.table("data\\DreamHF\\train\\readcounts_training.csv", sep=",", header=TRUE, row.names=1)

test <- list()
test[["pheno"]] <- read.table("data\\DreamHF\\test\\pheno_test.csv", sep=",", header=TRUE, row.names=1)
test[["readcounts"]] <- read.table("data\\DreamHF\\test\\readcounts_test.csv", sep=",", header=TRUE, row.names=1)

library(survival)

train_x <- apply(train[["readcounts"]], MARGIN=1, FUN=inormal)
test_x <- apply(test[["readcounts"]], MARGIN=1, FUN=inormal)

train_y <- survival::Surv(event = train[["pheno"]]$Event, time = train[["pheno"]]$Event_time)
test_y <- survival::Surv(event = test[["pheno"]]$Event, time = test[["pheno"]]$Event_time)

save.image("train_test.RData")

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
train[["y"]] <- survival::Surv(

par(mfrow=c(2,2))
boxplot(Age ~ Sex, data = train[["pheno"]], range=0)
boxplot(BodyMassIndex ~ Sex, data = train[["pheno"]], range=0)
boxplot(Smoking ~ Sex, data = train[["pheno"]], range=0)
boxplot(PrevalentHFAIL ~ Sex, data = train[["pheno"]], range=0)

ComplexHeatmap::Heatmap()




