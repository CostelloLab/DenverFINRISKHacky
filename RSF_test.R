library(randomForestSRC)
setwd("~/Downloads/DenverFINRISKHacky")

#############################################################################
# Import clinical data from the files they provide directly
test<-rio::import("./data/DreamHF/train/pheno_training.csv")

test[sapply(test, is.character)] <- lapply(test[sapply(test, is.character)], 
                                       as.factor)

test$Event_time=as.numeric(test$Event_time)
test$Event=as.numeric(test$Event)
test2<-test[!test$Event_time<0,]
test2=test2[!is.na(test2$Event_time),]
test.obj <- rfsrc(Surv(Event_time,Event) ~ ., data = test2,importance="permute")
plot.variable(test.obj)
plot(test.obj)

# Tests
test3<-test2[,-2]
test.obj3 <- rfsrc(Surv(Event_time,Event) ~ ., data = test3,importance="permute")
test_age_bmi<- test2[,c(2,3,9,10)]
test_age_bmi.obj <- rfsrc(Surv(Event_time,Event) ~ ., data = test_age_bmi,importance="permute")

# Combine all microbiome data

# Load Processed data by Daniel
load("~/Downloads/DenverFINRISKHacky/train_test.RData")

train_x<-as.data.frame(train_x)
train_y<-as.matrix(train_y)
train_p<-cbind(train_p,train_y)
train_p[sapply(train_p, is.character)] <- lapply(train_p[sapply(train_p, is.character)], 
                                                 as.factor)
train_p<-train_p[!train_p$time<0,]
train_p=train_p[!is.na(train_p$time),]

combined<-merge(train_x,train_p,by=0)

combined[sapply(combined, is.character)] <- lapply(combined[sapply(combined, is.character)], 
                                           as.factor)

test.obj <- rfsrc(Surv(time,status) ~ ., data = combined,importance="permute")
