setwd("C:/Users/Yash/Desktop/Sepsis Project/datasets/july 26")
test_20 <- read.csv("test_20_scaled.csv",header = T)
test_30 <- read.csv("test_30_scaled.csv",header = T)
train <- read.csv("train_scaled.csv",header = T)
validation <- read.csv("validation_scaled.csv",header = T)


train$X <- NULL

si <- train$sepsis_identifier #si is numeric
levels(si)

train$sepsis_identifier <- as.factor(train$sepsis_identifier)
levels(train$sepsis_identifier)

#Random Forest Feature Selection
require(randomForest)
fit=randomForest(sepsis_identifier~., data=train)

library(caret)
varImp(fit)
which(fit$importance>225)
varImpPlot(fit,type=2)

importanceOrder=order(-fit$importance)
names=rownames(fit$importance)[importanceOrder][1:30]
names

#####MODEL - BMA -#####
library(BMS)
library(ROCR)

#func. to calculate AUC 
calc_auc <- function(actual_y, pred_proba) {
  pred <- ROCR::prediction(pred_proba, actual_y)
  # AUC
  auc.perf <- performance(pred, "auc")
  return(auc.perf@y.values[[1]])
}

#func. to calculate optimal F scores
calc_optimal_threshold <- function(actual_y, pred_proba) {
  pred <- ROCR::prediction(pred_proba, actual_y)
  # Recall-Precision values        	
  RP.perf <- performance(pred, "prec", "rec")
  P <- RP.perf@y.values[[1]] #Precision
  R <- RP.perf@x.values[[1]] #Recall
  threshold <- RP.perf@alpha.values[[1]]
  df <- as.data.frame(cbind(P, R, threshold))
  df <- df[!is.na(df$P) & !is.na(df$R),]
  df <- mutate(df, f2_score = 5 * P * R / (4 * P + R))
  df <- mutate(df, f1_score = 2 * P * R / (P + R))
  df <- mutate(df, f1.5_score = 3.25 * P * R / (2.25 * P + R))
  optimal_f2_threshold <- df[which.max(df$f2_score),'threshold']
  optimal_f1_threshold <- df[which.max(df$f1_score),'threshold']
  optimal_f1.5_threshold <- df[which.max(df$f1.5_score),'threshold']
  return(c(optimal_f2_threshold, optimal_f1_threshold, optimal_f1.5_threshold))
}

# function to calculate cutoff
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    #returns the cutoff based on the youden index
    return(p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
#=======================================================
#
#data <- train
summary(train)
head(train)
names(train)
data_test <- train[c("GCS", "Respiration", "SpO2", "Pulse", "Anion_Gap", "Temp_Source", "Chloride","Alkaline_Phosphatase", "Platelet_Count", "Glucose",
                "Neutrophils_Percent", "Calcium", "Hemoglobin", "ALT_GPT", "WBCC","sepsis_identifier")]
#data <- NULL
data <- as.data.frame(data_test)
data <- train
summary(data)
#data$si <- "0"
#data$si <- as.factor(data$sepsis_identifier=="1")
#levels(data$si)[levels(data$si)=="TRUE"] <- "1"
#levels(data$si)[levels(data$si)=="FALSE"] <- "0"
#data$sepsis_identifier <- data$si
#data$si <- NULL
#-----data$sepsis_identifier is a factor--- check!

#set.seed(123)
data <- data[sample(nrow(data)),]
folds <- cut(seq(1,nrow(data)), breaks=10, labels=FALSE)
#10-fold cross-validation:
all_auc <- vector()
all_optimal_threshold_2 <- vector()
all_optimal_threshold_1 <- vector()
all_optimal_threshold_1.5 <- vector()
all_optimal_threshold_youden <- vector()
all_sens <- vector()
all_spec <- vector()
k <- vector()
i=1
for (i in 1:10) {
  test_indx <- which(folds==i, arr.ind=TRUE)
  test_data <- data[test_indx, ]
  train_data <- data[-test_indx, ]
  summary(train_temp)
  train_temp <- train_data
  train_temp$sepsis_identifier<- 0
  train_temp$sepsis_identifier<- as.numeric(train_data$sepsis_identifier==1)
  bma.model <-bms(train_temp$sepsis_identifier~.,data = train_data)
  
  #label<- train_data$sepsis_identifier
  #train.t<- data.frame(train_data[,-dim(train_data)])
  #bma.model<-bic.glm(train.t, label, strict=FALSE,
  #                  OR=10, glm.family="binomial", factor.type=FALSE)
  
  
  #  summary(bma.model)
  # test <- as.matrix(valset[,-dim(valset)])
  #  bma_pred <- predict( bma.glm, newdata = test,type = "response")
  
  
  #summary(pred_proba)
  
  #summary(xdat)
  #temp <- xdat$sepsis_identifier
  #xdat$sepsis_identifier <- as.numeric(xdat$sepsis_identifier==1)
  
  xdat <- test_data[,1:(ncol(test_data)-1)]
  pred_proba <- predict( bma.model, newdata = xdat,type = "response")
  
  actual_y <- test_data[,ncol(test_data)]
  all_auc[i] <- calc_auc(actual_y, pred_proba)
  best_thresholds <- calc_optimal_threshold(actual_y, pred_proba)
  
  all_optimal_threshold_2[i] <- best_thresholds[1]
  all_optimal_threshold_1[i] <- best_thresholds[2]
  all_optimal_threshold_1.5[i] <- best_thresholds[3]
  CUTOFF <- best_thresholds[2]
  
  #all_sens[i] <- calc_sens(actual_y, pred_proba,CUTOFF)[1]
  #all_spec[i] <- calc_sens(actual_y, pred_proba,CUTOFF)[2]
  
  t<- 0
  t <- as.numeric(test_data$sepsis_identifier==1)
  
  summary(t)
  
  
  sepsis.pr <- ROCR::prediction(pred_proba, test_data$sepsis_identifier)
  sepsis.roc_perf <- performance(sepsis.pr, measure ="tpr", x.measure = "fpr")
  all_optimal_threshold_youden[i] <- opt.cut(sepsis.roc_perf, sepsis.pr)
}

cat("Train data: \n")
cat(paste0('Mean CV sensitivity: ', as.character(mean(all_sens)), '\n'))
cat(paste0('Mean CV specificity: ', as.character(mean(all_spec)), '\n'))
cat(paste0('Mean CV auc: ', as.character(mean(all_auc)), '\n'))
cat(paste0('Best threshold based on F2 score: ', as.character(mean(all_optimal_threshold_2)), '\n'))
cat(paste0('Best threshold based on F1 score: ', as.character(mean(all_optimal_threshold_1)), '\n'))
cat(paste0('Best threshold based on F1.5 score: ', as.character(mean(all_optimal_threshold_1.5)), '\n'))
cat(paste0('Best threshold based on youden index: ', as.character(mean(all_optimal_threshold_youden)), '\n'))
cat(paste0('\n'))

#top - 5
#based on 
#F2 : 
#F1 :
#F1.5 :
#youden-index: 


#top - 10
#based on 
#F2 : 
#F1 :
#F1.5 :
#youden-index: 

#top -15
#based on 
#F2 : 
#F1 :
#F1.5 :
#youden-index: 

#all
#based on 
#F2 : 0.304294165512741
#F1 : 0.432113985939036
#F1.5 : 0.346455701243973
#youden-index: 0.352985536620346

