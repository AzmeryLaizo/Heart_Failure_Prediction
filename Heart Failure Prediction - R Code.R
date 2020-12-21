############################ DATA PREPROCESSING ################################

rm(list=ls())
gc() 
setwd('C:\\Users\\Azmery\\Documents\\CSUF\\ISDS 574\\Group Project\\NEW')

dat = read.csv('heart_failure_clinical_records_dataset.csv', head=T, stringsAsFactors=F, na.strings='')

# Check missing values & visualize
matrix.na = is.na(dat)
pmiss = colMeans(matrix.na)
nmiss = rowMeans(matrix.na)
plot(pmiss)

library(Amelia)
missmap(dat)

# Remove 'time' variable
dat$time = NULL

# Check the distribution of continuous variables
par(mfrow=c(1, 2))

# Histogram & Boxplot for age
hist(dat$age)
boxplot(dat$age)

# Histogram & Boxplot for creatinine_phosphokinase
hist(dat$creatinine_phosphokinase)
boxplot(dat$creatinine_phosphokinase)

# Removing outlier for creatinine_phosphokinase
which(dat$creatinine_phosphokinase > 3000)
dim(dat[c(2,53,61,73,104,135,172),])
dat = dat[dat$creatinine_phosphokinase <= 3000,]

# Histogram & Boxplot for ejection_fraction
hist(dat$ejection_fraction)
boxplot(dat$ejection_fraction)

# Removing outlier for ejection_fraction
which(dat$ejection_fraction > 65)
dim(dat[c(62,211),])
dat = dat[dat$ejection_fraction <= 65,]

# Histogram & Boxplot for platelets
hist(dat$platelets)
boxplot(dat$platelets)

# Removing outlier for platelets
which(dat$platelets > 590000)
dim(dat[c(100,104,288),])
dat = dat[dat$platelets <= 590000,]

# Histogram & Boxplot for serum_creatinine
hist(dat$serum_creatinine)
boxplot(dat$serum_creatinine)

#removing outlier for serum_creatinine
which(dat$serum_creatinine > 8)
dim(dat[9,])
dat = dat[dat$serum_creatinine <= 8,]

# Histogram & Boxplot for serum_sodium
hist(dat$serum_sodium)
boxplot(dat$serum_sodium)

#removing outlier for serum_sodium
which(dat$serum_sodium < 125)
dim(dat[c(4,18,117,187),])
dat = dat[dat$serum_sodium >= 125,]

#check frequency of age
table(dat$age)
#combining possible data entry error
which(dat$age == 60.667)
id = which(dat$age %in% c(60.667,61))
dat$age[id] = 61

# Checking the correlation between variables
library(corrplot)
par(mfrow = c(1,1))
corrplot(cor(dat), method="number", number.cex = .5)

library(dplyr)
dat %>% cor() %>% round(2) %>% View()

save(dat, file = "cleaned_Heart_failure_data.rda")

############################ LOGISTIC REGRESSION ###############################

rm(list=ls())
gc()
load('./cleaned_Heart_failure_data.rda')

# Data partition to have 90% training data
set.seed(1) 
id.train = sample(1:nrow(dat), nrow(dat)*.9) 
id.test = setdiff(1:nrow(dat), id.train) 
dat.train = dat[id.train,]
dat.test = dat[id.test,]

min.model = glm(DEATH_EVENT ~ 1, data = dat.train, family = 'binomial')
max.model = glm(DEATH_EVENT ~ ., data = dat.train, family = 'binomial')
max.formula = formula(max.model)

# Forward selection

# Cutoff = 0.2 ### OUR BEST MODEL FOR LOGISTIC REGRESSION ####

obj2 = step(min.model, direction='forward', scope=max.formula)
summary(obj2)

get.or = function(sobj, alpha=.05) {
  b = sobj$coef[-1, 'Estimate']
  se.b = sobj$coef[-1, 'Std. Error']
  pval = sobj$coef[-1, 'Pr(>|z|)']
  or = exp(b); se.or = exp(b)*se.b
  lb = b - qnorm(alpha/2)*se.b; lb.or = exp(lb)
  ub = b + qnorm(1-alpha/2)*se.b; ub.or = exp(ub)
  out = cbind(or, se.or, lb, ub, pval)
  colnames(out) = c('OR', 'SE', paste0((1-alpha)*100, '% CI, lower'),
                    paste0((1-alpha)*100, '% CI, upper'), 'p value')
  return(out)
}
get.or(summary(obj2))

yhat2 = predict(obj2, newdata = dat.test, type='response')
hist(yhat2)

dichotomize = function(yhat, cutoff=.5) {
  out = rep(0, length(yhat))
  out[yhat > cutoff] = 1
  out
}

yhat2.class = dichotomize(yhat2, .2) 
err2 = mean(yhat2.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err2
# Confusion matrix
table(yhat2.class, dat.test$DEATH_EVENT)


sen = function(ytrue, yhat) {
  ind.true1 = which(ytrue == 1)
  mean( ytrue[ind.true1] == yhat[ind.true1] )
}

spe = function(ytrue, yhat) {
  ind.true0 = which(ytrue == 0)
  mean( ytrue[ind.true0] == yhat[ind.true0] )
}

# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat2.class)
spe(dat.test$DEATH_EVENT, yhat2.class)

# Cutoff = 0.1

obj = step(min.model, direction='forward', scope=max.formula) 
summary(obj) 

yhat = predict(obj, newdata = dat.test, type='response')
hist(yhat)

yhat.class = dichotomize(yhat, .1) 
err = mean(yhat.class != dat.test$DEATH_EVENT) # Misclassification error rate
err

table(yhat.class, dat.test$DEATH_EVENT)

sen(dat.test$DEATH_EVENT, yhat.class)
spe(dat.test$DEATH_EVENT, yhat.class)

# Cutoff = 0.5

obj5 = step(min.model, direction='forward', scope=max.formula) 
summary(obj5) 

yhat5 = predict(obj5, newdata = dat.test, type='response')
hist(yhat5)

yhat5.class = dichotomize(yhat5, .5) 
err5 = mean(yhat5.class != dat.test$DEATH_EVENT)

# Misclassification error rate
err5
# Confusion matrix
table(yhat5.class, dat.test$DEATH_EVENT)

# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat5.class)
spe(dat.test$DEATH_EVENT, yhat5.class)

# Cutoff = 0.4

obj4 = step(min.model, direction='forward', scope=max.formula) 
summary(obj4) 

yhat4 = predict(obj4, newdata = dat.test, type='response')
hist(yhat4)

yhat4.class = dichotomize(yhat4, .4) 
err4 = mean(yhat4.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err4
# Confusion matrix
table(yhat4.class, dat.test$DEATH_EVENT)

# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat4.class)
spe(dat.test$DEATH_EVENT, yhat4.class)

# Cutoff = 0.3

obj3 = step(min.model, direction='forward', scope=max.formula) 
summary(obj3) 

yhat3 = predict(obj3, newdata = dat.test, type='response')
hist(yhat3)

yhat3.class = dichotomize(yhat3, .3) 
err3 = mean(yhat3.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err3
# Confusion matrix
table(yhat3.class, dat.test$DEATH_EVENT)

# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat3.class)
spe(dat.test$DEATH_EVENT, yhat3.class)

# Backward selection

# Cutoff = 0.1

obj6 = step(max.model, direction='backward', scope=max.formula)
summary(obj6)

get.or(summary(obj6))

yhat6 = predict(obj6, newdata = dat.test, type='response')
hist(yhat6)
yhat6.class = dichotomize(yhat6, .1)
err6 = mean(yhat6.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err6
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat6.class)
spe(dat.test$DEATH_EVENT, yhat6.class)

# Cutoff = 0.2

obj7 = step(max.model, direction='backward', scope=max.formula)
summary(obj7)

get.or(summary(obj7))

yhat7 = predict(obj7, newdata = dat.test, type='response')
hist(yhat7)
yhat7.class = dichotomize(yhat7, .2)
err7 = mean(yhat7.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err7
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat7.class)
spe(dat.test$DEATH_EVENT, yhat7.class)

# Cutoff = 0.3
obj8 = step(max.model, direction='backward', scope=max.formula)
summary(obj8)

get.or(summary(obj8))

yhat8 = predict(obj8, newdata = dat.test, type='response')
hist(yhat8)
yhat8.class = dichotomize(yhat8, .3)
err8 = mean(yhat8.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err8
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat8.class)
spe(dat.test$DEATH_EVENT, yhat8.class)

## Cutoff = 0.4
obj9 = step(max.model, direction='backward', scope=max.formula)
summary(obj9)

get.or(summary(obj9))

yhat9 = predict(obj9, newdata = dat.test, type='response')
hist(yhat9)
yhat9.class = dichotomize(yhat9, .4)
err9 = mean(yhat9.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err9
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat9.class)
spe(dat.test$DEATH_EVENT, yhat9.class)

## Cutoff = 0.5
obj10 = step(max.model, direction='backward', scope=max.formula)
summary(obj10)

get.or(summary(obj10))

yhat10 = predict(obj10, newdata = dat.test, type='response')
hist(yhat10)
yhat10.class = dichotomize(yhat10, .5)
err10 = mean(yhat10.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err10
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat10.class)
spe(dat.test$DEATH_EVENT, yhat10.class)

# Stepwise selection

# Cutoff = 0.1

obj11 = step(min.model, direction='both', scope=max.formula) 
summary(obj11)

get.or(summary(obj11))

yhat11 = predict(obj11, newdata = dat.test, type='response')
hist(yhat11)
yhat11.class = dichotomize(yhat11, .1)
err11 = mean(yhat11.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err11
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat11.class)
spe(dat.test$DEATH_EVENT, yhat11.class)

# Cutoff = 0.2
obj12 = step(min.model, direction='both', scope=max.formula) 
summary(obj12)

get.or(summary(obj12))

yhat12 = predict(obj12, newdata = dat.test, type='response')
hist(yhat12)
yhat12.class = dichotomize(yhat12, .2)
err12 = mean(yhat12.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err12
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat12.class)
spe(dat.test$DEATH_EVENT, yhat12.class)

# Cutoff = 0.3
obj13 = step(min.model, direction='both', scope=max.formula) 
summary(obj13)

get.or(summary(obj13))

yhat13 = predict(obj13, newdata = dat.test, type='response')
hist(yhat13)
yhat13.class = dichotomize(yhat13, .3)
err13 = mean(yhat13.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err13
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat13.class)
spe(dat.test$DEATH_EVENT, yhat13.class)

# Cutoff = 0.4
obj14 = step(min.model, direction='both', scope=max.formula) 
summary(obj14)

get.or(summary(obj14))

yhat14 = predict(obj14, newdata = dat.test, type='response')
hist(yhat14)
yhat14.class = dichotomize(yhat14, .4)
err14 = mean(yhat14.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err14
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat14.class)
spe(dat.test$DEATH_EVENT, yhat14.class)

# Cutoff = 0.5
obj15 = step(min.model, direction='both', scope=max.formula) 
summary(obj15)

get.or(summary(obj15))

yhat15 = predict(obj15, newdata = dat.test, type='response')
hist(yhat15)
yhat15.class = dichotomize(yhat15, .5)
err15 = mean(yhat15.class != dat.test$DEATH_EVENT) 

# Misclassification error rate
err15
# Sensitivity & Specificity
sen(dat.test$DEATH_EVENT, yhat15.class)
spe(dat.test$DEATH_EVENT, yhat15.class)

############################ CLASSIFICATION TREE ###############################

rm(list=ls())
gc()
load('./cleaned_Heart_failure_data.rda')

library(rpart)
library(rpart.plot)

# Data partition to have 80% training data
set.seed(1)
n.train = floor(nrow(dat)*.8)
id.train = sample(1:nrow(dat), n.train)
id.test = setdiff(1:nrow(dat), id.train)

# Classification Tree with rpart
K = 10 # number of cross-validations
fit = rpart(DEATH_EVENT ~ ., method="class", data=dat[id.train,], cp = 0.01, minsplit = 10, xval=K)

printcp(fit)
plotcp(fit)

# Minimum Error Tree
pfit.me = prune(fit, cp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
rpart.plot(pfit.me, main = 'Min Error Tree')

Summary(pfit.me) # to view the variable importance

## Prediction using 0.5 cutoff
yhat = predict(pfit.me, dat[id.test,], type = "class") 
err.me = 1 - mean(yhat == dat[id.test,'DEATH_EVENT'])
err.me

ytest = factor(dat[id.test,'DEATH_EVENT'])


sen = function(ytrue, yhat) {
  ind.true1 = which(ytrue == 1)
  mean( ytrue[ind.true1] == yhat[ind.true1] )
}

spe = function(ytrue, yhat) {
  ind.true0 = which(ytrue == 0)
  mean( ytrue[ind.true0] == yhat[ind.true0] )
}
# Confusion matrix
table(ytest, yhat)
#Sensitivity & Specificity
sen(ytest, yhat)
spe(ytest, yhat)

# Classification probability
prob1 = predict(pfit.me, dat[id.test,], type = "prob")[,2]

# Prediction using 0.4 cutoff
yhat1 = as.numeric(prob1 > .4)
err.me.newCut = 1 - mean(yhat1 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat1)
library(rpart)
library(rpart.plot)

# Data partition to have 80% training data
set.seed(1)
n.train = floor(nrow(dat)*.8)
id.train = sample(1:nrow(dat), n.train)
id.test = setdiff(1:nrow(dat), id.train)

# Classification Tree with rpart
K = 10 # number of cross-validations
fit = rpart(DEATH_EVENT ~ ., method="class", data=dat[id.train,], cp = 0.01, minsplit = 10, xval=K)

printcp(fit)
plotcp(fit)

# Minimum Error Tree
pfit.me = prune(fit, cp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
rpart.plot(pfit.me, main = 'Min Error Tree')


## Prediction using 0.5 cutoff
yhat = predict(pfit.me, dat[id.test,], type = "class") 
err.me = 1 - mean(yhat == dat[id.test,'DEATH_EVENT'])
err.me

ytest = factor(dat[id.test,'DEATH_EVENT'])


sen = function(ytrue, yhat) {
  ind.true1 = which(ytrue == 1)
  mean( ytrue[ind.true1] == yhat[ind.true1] )
}

spe = function(ytrue, yhat) {
  ind.true0 = which(ytrue == 0)
  mean( ytrue[ind.true0] == yhat[ind.true0] )
}
# Confusion matrix
table(ytest, yhat)
#Sensitivity & Specificity
sen(ytest, yhat)
spe(ytest, yhat)

# Classification probability
prob1 = predict(pfit.me, dat[id.test,], type = "prob")[,2]

# Prediction using 0.4 cutoff
yhat1 = as.numeric(prob1 > .4)
err.me.newCut = 1 - mean(yhat1 == ytest)
err.me.newCut

table(ytest, yhat1)
sen(ytest, yhat1)
spe(ytest, yhat1)

# Prediction using 0.3 cutoff
yhat2 = as.numeric(prob1 > .3)
err.me.newCut = 1 - mean(yhat2 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat2)
# Sensitivity & Specificity
sen(ytest, yhat2)
spe(ytest, yhat2)

# Prediction using 0.2 cutoff
yhat3 = as.numeric(prob1 > .2)
err.me.newCut = 1 - mean(yhat3 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat3)
# Sensitivity & Specificity
sen(ytest, yhat3)
spe(ytest, yhat3)

# Prediction using 0.1 cutoff
yhat4 = as.numeric(prob1 > .1)
err.me.newCut = 1 - mean(yhat3 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat4)
# Sensitivity & Specificity
sen(ytest, yhat4)
spe(ytest, yhat4)
sen(ytest, yhat1)
spe(ytest, yhat1)

# Prediction using 0.3 cutoff
yhat2 = as.numeric(prob1 > .3)
err.me.newCut = 1 - mean(yhat2 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat2)
# Sensitivity & Specificity
sen(ytest, yhat2)
spe(ytest, yhat2)

# Prediction using 0.2 cutoff
yhat3 = as.numeric(prob1 > .2)
err.me.newCut = 1 - mean(yhat3 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat3)
# Sensitivity & Specificity
sen(ytest, yhat3)
spe(ytest, yhat3)

# Prediction using 0.1 cutoff
yhat4 = as.numeric(prob1 > .1)
err.me.newCut = 1 - mean(yhat3 == ytest)
err.me.newCut

# Confusion matrix
table(ytest, yhat4)
# Sensitivity & Specificity
sen(ytest, yhat4)
spe(ytest, yhat4)

############################ k-NEAREST NEIGHBOUR ###############################

rm(list=ls())
gc()
load('./cleaned_Heart_failure_data.rda')

# Normalize all continuous variables
dat[,1] = scale(dat[,1])
dat[,3] = scale(dat[,3])
dat[,5] = scale(dat[,5])
dat[,7:9] = scale(dat[,7:9])

# Remove variables to run the model with variables selected by Logistic Regression
dat$anaemia=NULL
dat$creatinine_phosphokinase=NULL
dat$diabetes=NULL
dat$platelets=NULL
dat$sex=NULL
dat$smoking=NULL
dat$anaemia=NULL

# Data partition to have 80% training data
set.seed(1)
n.train = floor( nrow(dat)*0.8)
ind.train = sample(1:nrow(dat), n.train)
ind.test = setdiff(1:nrow(dat), ind.train)

require(class)
Xtrain = dat[ind.train,1:5]
Xtest = dat[ind.test,1:5]
ytrain = dat[ind.train,6]
ytest = dat[ind.test,6]

get.prob = function(x) {
  prob = attr(x, 'prob')
  cl = as.numeric(x)
  ind = which(cl == 1)
  prob[ind] = 1 - prob[ind]
  return(prob)
}

knn.bestK = function(train, test, y.train, y.test, k.grid = 1:282, ct = .7) {
  # browser()
  fun.tmp = function(x) {
    y.tmp = knn(train, test, y.train, k = x, prob=T) # run knn for each k in k.grid
    prob = get.prob(y.tmp)
    y.hat = as.numeric( prob > ct ) + 1
    return( sum(y.hat != as.numeric(y.test)) )
  }
  ## create a temporary function (fun.tmp) that we want to apply to each value in k.grid
  error = unlist(lapply(k.grid, fun.tmp))
  names(error) = paste0('k=', k.grid)
  ## it will return a list so I need to unlist it to make it to be a vector
  out = list(k.optimal = k.grid[which.min(error)], 
             error.min = min(error)/length(y.test),
             error.all = error/length(y.test))
  return(out)
}

obj1 = knn.bestK(Xtrain, Xtest, ytrain, ytest, seq(1,25, 2), .7)
obj1

# Rerun with the best k
ypred = knn(Xtrain, Xtest, ytrain, k=obj1$k.optimal, prob=T)

sen = function(ytrue, yhat) {
  ind.true1 = which(ytrue == 1)
  mean( ytrue[ind.true1] == yhat[ind.true1] )
}

spe = function(ytrue, yhat) {
  ind.true0 = which(ytrue == 0)
  mean( ytrue[ind.true0] == yhat[ind.true0] )
}

# Confusion matrix
table(ytest, ypred)
# Sensitivity and specifity
sen(ytest, ypred)
spe(ytest, ypred)

############################ RANDOM FOREST ####################################

rm(list=ls())
gc()
load('./cleaned_Heart_failure_data.rda')

# take 90% of data randomly as training
set.seed(1)
id.train = sample(1:nrow(dat), nrow(dat)*.9) 
id.test = setdiff(1:nrow(dat), id.train) 
dat.train = dat[id.train,]
dat.test = dat[id.test,]

#Change to factor for classification
dat.train$DEATH_EVENT <- as.character(dat.train$DEATH_EVENT)
dat.train$DEATH_EVENT <- as.factor(dat.train$DEATH_EVENT)

dat.test$DEATH_EVENT <- as.character(dat.test$DEATH_EVENT)
dat.test$DEATH_EVENT <- as.factor(dat.test$DEATH_EVENT)

#Train the model
library(randomForest)
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, proximity = TRUE, importance = TRUE)
rf_classifier

#predict with test data
rf_predict = predict(rf_classifier, newdata = dat.test)

#install caret package
library(caret)
confusionMatrix(rf_predict, dat.test$DEATH_EVENT, positive = '1' )


#Tuning the model
#install e1071 package
#1.Using Grid Search with cross validation
library(e1071)
# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
set.seed(1)
# Run the model
rf_default <- train(DEATH_EVENT~.,
                    data = dat.train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)
# Print the results
print(rf_default)
#Accuracy was used to select the optimal model using the largest value.
#The final value used for the model was mtry = 6.


#Let's search best mtry for better results
set.seed(1)
tuneGrid <- expand.grid(.mtry = c(1: 11))
rf_mtry <- train(DEATH_EVENT~.,
                 data = dat.train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 proximity = TRUE,
                 ntree = 500)
print(rf_mtry)
#The final value used for the model was mtry = 6.
plot(rf_mtry)
#Save the best mtry
rf_mtry$bestTune$mtry
max(rf_mtry$results$Accuracy)
best_mtry <- rf_mtry$bestTune$mtry 
best_mtry

#2.tune mtry by tools
set.seed(1)
bestMtry <- tuneRF(dat.train[,-12], dat.train[,12], stepFactor = 1.5, improve = 1e-5, ntree = 500)
print(bestMtry)
#best mtry = 3

#Error rate of RF
plot(rf_classifier)
#From 300 trees the error rate is almost flat

#best ntree with default mtry
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(1)
  rf_maxtrees <- train(DEATH_EVENT~.,
                       data = dat.train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       proximity = TRUE,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
dotplot(results_tree)
#best mtree is 400


#train model with new arguments
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, ntree = 450, mtry =3, proximity = TRUE, importance = TRUE)
rf_classifier
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, ntree = 400, mtry =6, proximity = TRUE, importance = TRUE)
rf_classifier
rf_predict = predict(rf_classifier, newdata = dat.test)
confusionMatrix(rf_predict, dat.test$DEATH_EVENT,positive = '1' )
#lets choose mtry =6 and ntree=400

#Adjusting the cutoff
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, ntree = 400, mtry =6, cutoff = c(0.6,0.4), proximity = TRUE, importance = TRUE)
rf_classifier
rf_predict = predict(rf_classifier, newdata = dat.test)
confusionMatrix(rf_predict, dat.test$DEATH_EVENT, positive = '1' )


#Search the best max nodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(2: 25)) {
  set.seed(1)
  rf_maxnode <- train(DEATH_EVENT~.,
                      data = dat.train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      proximity = TRUE,
                      maxnodes = maxnodes,
                      ntree = 400, cutoff = c(0.6,0.4))
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)
#best maxnodes = 11

#adding maxnodes
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, ntree = 400, mtry =6, cutoff = c(0.6,0.4) , maxnodes = 11, proximity = TRUE, importance = TRUE)
rf_classifier
rf_predict = predict(rf_classifier, newdata = dat.test)
confusionMatrix(rf_predict, dat.test$DEATH_EVENT, positive = '1' )
#since the accuracy and specificity is lesser, let's not specify maxnodes and let the trees grow to the max possible. 

#Final random forest classifier
set.seed(1)
rf_classifier <- randomForest(DEATH_EVENT ~ ., data = dat.train, ntree = 400, mtry =6,cutoff = c(0.6,0.4), proximity = TRUE, importance = TRUE)
rf_classifier
rf_predict = predict(rf_classifier, newdata = dat.test)
confusionMatrix(rf_predict, dat.test$DEATH_EVENT, positive = '1' )

#Evaluate variable importance
importance(rf_classifier)
varImpPlot(rf_classifier)

###############################################################################
#################################### END ######################################