#load libraries
library(caret)
library(dplyr)
library(tibble)

#set working directory to parent folder of the current script file
myDir <- unlist(strsplit(dirname(this.path::this.path()), '/'))
setwd(paste0(myDir[-length(myDir)], collapse = '/'))
getwd()

#targets file
targets_ML <- read.csv(file="data/CCLE_targets_with_NMF_result.csv") #targets file

#metagenes files
genes_k5_all <- read.csv(file="data/NMF_metagenes_k5_all.csv", row.names=1)

#format data
dat <- merge(targets_ML, t(genes_k5_all), by.x=2, by.y=0)
dat <- dat[,-c(2,3,4,6)]

#create test/training partitions
set.seed(123)
folds <- createMultiFolds(dat$NMF_cluster, k = 4, times = 10)

#set training parameters
train_control <- trainControl(method = "repeatedcv", number = 4, repeats = 10, index = folds, savePredictions = TRUE)

#train model (random forrest)
model_rf <- train(dat[,3:210], dat[,2], trControl = train_control, preProcess = c("center","scale"), method = 'rf')
#extract per class sensitivity and specificity
confusion_rf <- model_rf$pred
metrics_rf <- confusionMatrix(confusion_rf$pred, confusion_rf$obs)
metrics_rf_per_class <- as.data.frame(metrics_rf[["byClass"]])
metrics_rf_per_class$model <- "RF"
metrics_rf_per_class <- metrics_rf_per_class %>%
  rownames_to_column(var = "Subtype") %>%
  mutate(Subtype = gsub("^Class: ", "", Subtype)) %>%
  dplyr::select(Subtype, Sensitivity, Specificity, model)

#train model (KNN)
model_knn <- train(dat[,3:210], dat[,2], trControl = train_control, preProcess = c("center","scale"), method = 'knn')
#extract per class sensitivity and specificity
confusion_knn <- model_knn$pred
metrics_knn <- confusionMatrix(confusion_knn$pred, confusion_knn$obs)
metrics_knn_per_class <- as.data.frame(metrics_rf[["byClass"]])
metrics_knn_per_class$model <- "KNN"
metrics_knn_per_class <- metrics_knn_per_class %>%
  rownames_to_column(var = "Subtype") %>%
  mutate(Subtype = gsub("^Class: ", "", Subtype)) %>%
  dplyr::select(Subtype, Sensitivity, Specificity, model)

#train model (KNN)
model_svm <- train(dat[,3:210], dat[,2], trControl = train_control, preProcess = c("center","scale"), method = 'svm')
#extract per class sensitivity and specificity
confusion_knn <- model_knn$pred
metrics_knn <- confusionMatrix(confusion_knn$pred, confusion_knn$obs)
metrics_knn_per_class <- as.data.frame(metrics_knn[["byClass"]])
metrics_knn_per_class$model <- "KNN"
metrics_knn_per_class <- metrics_knn_per_class %>%
  rownames_to_column(var = "Subtype") %>%
  mutate(Subtype = gsub("^Class: ", "", Subtype)) %>%
  dplyr::select(Subtype, Sensitivity, Specificity, model)

#train model (KNN)
model_svm <- train(dat[,3:210], dat[,2], trControl = train_control, preProcess = c("center","scale"), method = "svmLinear")
#extract per class sensitivity and specificity
confusion_svm <- model_svm$pred
metrics_svm <- confusionMatrix(confusion_svm$pred, confusion_svm$obs)
metrics_svm_per_class <- as.data.frame(metrics_svm[["byClass"]])
metrics_svm_per_class$model <- "SVM"
metrics_svm_per_class <- metrics_svm_per_class %>%
  rownames_to_column(var = "Subtype") %>%
  mutate(Subtype = gsub("^Class: ", "", Subtype)) %>%
  dplyr::select(Subtype, Sensitivity, Specificity, model)

#bind all metrics together for plotting
metrics_per_class_all <- rbind(metrics_knn_per_class,metrics_rf_per_class,metrics_svm_per_class)

#sensitivity and specificity per class and model plot
ggplot(metrics_per_class_all, aes(x=Sensitivity, y=Specificity, colour=Subtype)) +
  geom_point(size=4, shape=1) +
  facet_grid(model ~ .) +
  theme_bw()

#export RF model
saveRDS(model_rf, file = "data/RF_model.rds")
