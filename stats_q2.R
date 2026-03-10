# === PACKAGES USED ========
library(readxl)
library(glmnet)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)

# === 1 | DATA CLEANING AND EDA ========
COV <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
View(COV)
str(COV)

#dataset contains measurements of 28 analytes (cytokines, chemokines, and acute phase proteins) from plasma samples of 76 COVID-19-positive patients. Response variable: disease progression/severity (binary, mild vs. severe)

table(COV$Severirty)

## DO SOME EDA AND CHECK FOR MISSING VALUES (there should be none)

# Convert severity to a factor variable
#Severity <- ifelse(COV$Severirty == "Severe", 1, 0)
#table(Severity)

COV$Severity <- ifelse(COV$Severirty == "Severe", 1, 0)
table(COV$Severity) # 0 = mild, 1 = severe

# Remove Patient Number and Severity
COV <- COV[, !(names(COV) %in% c("Severirty", "Patient Number"))]
str(COV)

# Convert SEX to numeric
COV$SEX <- ifelse(COV$SEX == "F", 1, 0)
table(COV$SEX)  

dim(COV)       # 76 x 31
colnames(COV)  


# Feature Engineering ?????

# === 2.1 | MODEL MATRIX AND TRAINING/TEST SETS ========
# Create model matrix for glmnet including all effects (even quadratic?)
#f0 <- lm(Severity ~ ., data = COV)
#f0 <- lm(Severity ~ .^2, data = COV)

X <- model.matrix(Severity ~., COV)[, -1]
#X <- model.matrix(Severity ~.^2, COV)[, -1]
dim(X)
colnames(X)
# X <- cbind(X,X[, 1:29]^2)

Y <- COV$Severity

# predictor standard deviations
plot(as.ts(apply(X, 2, sd)), xlab = "Covariate", ylab = "SD",
     main = "Standard Deviations of Predictors")
# very variable so thats why standardization is needed


# Create train (75%) & test (25%) sets
set.seed(1717)

test_index <- sample.int(nrow(X), round(nrow(X) * 0.25), replace = FALSE)
test_index

X_train <- X[-test_index, ]
Y_train <- Y[-test_index]
X_test <- X[test_index, ]
Y_test <- Y[test_index]

## Confirm class balance in train and test
table(Y_train)
table(Y_test)

# === 2.2 | GRID SEARCH FOR ALPHA ========
# Function to perform grid search for optimal alpha value with different n-folds 
alpha_search <- function (nfolds){
  alpha_grid <- seq(0, 1, .01)
  results <- data.frame()
  set.seed(1717)
  for (a in alpha_grid) {
  cv_model <- cv.glmnet(X_train, 
                        Y_train, 
                        alpha = a, 
                        nfolds = nfolds,
                        family = "binomial",
                        #foldid check if needed
                        type.measure = "deviance",
                        keep = T)  # Accuracy measure is deviance because AUC could not be used; number of observations is too small per fold
                        #keep = TRUE) # check if keep should be T or F
  
  # Store results
  results <- rbind(results, data.frame(
    alpha = a,
    lambda_min = cv_model$lambda.min,
    lambda_1se = cv_model$lambda.1se,
    cvm_min = min(cv_model$cvm),
    cvm_1se = cv_model$cvm[which(cv_model$lambda == cv_model$lambda.1se)]
  ))
  }
  return(results)
}
## Add part to function that collects and prints best alpha values automatically

alpha_10 <- alpha_search(10)
best_min_10 <- alpha_10$alpha[which.min(alpha_10$cvm_min)]
best_1se_10 <- alpha_10$alpha[which.min(alpha_10$cvm_1se)]

alpha_20 <- alpha_search(20)
best_min_20 <- alpha_20$alpha[which.min(alpha_20$cvm_min)]
best_1se_20 <- alpha_20$alpha[which.min(alpha_20$cvm_1se)]

# Visualize alpha grid search
ggplot(alpha_10, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (MSE)",
       color = "Lambda Type") +
  theme_minimal()

ggplot(alpha_20, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (MSE)",
       color = "Lambda Type") +
  theme_minimal()

## Decide between alpha values

## Create functions to streamline comparison between 10 and 20 folds
# Elastic net model with 10 folds
cv_10 <- cv.glmnet(X_train, Y_train, nfolds = 10, family = "binomial", alpha = best_min_10, type.measure = "deviance")
plot(cv_10)

prd_train <- predict(cv_10, newx = X_train, type = "response", s = cv_10$lambda.min)[, 1]
prd_test <- predict(cv_10, newx = X_test, type = "response", s = cv_10$lambda.min)[, 1]

# Testing prediction accuracy (Change plots to be made with ggplot and overlay curves)
par(mfrow=c(1,2))
AUC_train <- roc(Y_train, prd_train)
plot(AUC_train)
AUC_test <- roc(Y_test, prd_test)
plot(AUC_test)

# Check for performance of test set
max(cv_10$cvm)
which.max(cv_10$cvm)
cv_10$cvlo[which.max(cv_10$cvm)];cv_10$cvup[which.max(cv_10$cvm)]

# Determining cutoffs for sensitivity, specificity, and classification; make function to test both; return dataframe with all
cutoffs_snsp <- function(AUC) {
  snsp <- cbind(AUC$sensitivities, AUC$specificities)
  indx <- which.max(apply(snsp, 1, min))
  cutoff <- AUC$threshold[indx]
  return(cutoff)
}

cutoffs_snsp(AUC_train)
cutoffs_snsp(AUC_test)
# Make confusion matrix?

# Training set
snsp_train <- cbind(AUC_train$sensitivities, AUC_train$specificities)
snsp_train

indx_train <- which.max(apply(snsp_train, 1, min))
indx_train
snsp_train[indx_train, ]

train_cutoff <- AUC_train$thresholds[indx_train]
train_cutoff

# Test set
snsp_test <- cbind(AUC_test$sensitivities, AUC_test$specificities)
snsp_test

indx_test <- which.max(apply(snsp_test, 1, min))
indx_test
snsp_test[indx_test, ]

test_cutoff <- AUC_test$thresholds[indx_test]
test_cutoff

# Make ggplot versions
par(mfrow=c(1,2))
plot(AUC_train)
abline(h=snsp_train[indx_train,1],v=snsp_train[indx_train,2], col='blue', lty=2)
plot(AUC_test)
abline(h=snsp_test[indx_test,1],v=snsp_test[indx_test,2], col='blue', lty=2)
par(mfrow=c(1,1))

## A function to compute Sensitivity and Specificity (probably remake function or extract SN and SP in a different way)
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

sn.sp(table(Y_train, yhat = as.numeric(prd_train > train_cutoff))) ## TRAIN SET
sn.sp(table(Y_test, yhat = as.numeric(prd_test > train_cutoff))) ## TEST SET
sn.sp(table(Y_test, yhat = as.numeric(prd_test > test_cutoff)))

# Variable selection
coef(cv_10, s = cv_10$lambda.min)
coef.min <- coef(cv_10,s = cv_10$lambda.min)[-1,1]
coef.min[coef.min!=0] 
names(coef.min[coef.min!=0])

# Provide ranking of top predictors

# Elastic net model with 20 folds
## Make function that does all the same steps as 10-fold
cv_20 <- cv.glmnet(X_train, Y_train, nfolds = 20, family = "binomial", alpha = best_min_20, type.measure = "deviance")
plot(cv_20)

coef(cv_20, s = cv_20$lambda.min)
coef.min <- coef(cv_20,s = cv_20$lambda.min)[-1,1]
coef.min[coef.min!=0]

## FOR MS. MAHNOOR NIZAMANI: Make some plots and do some analysis for how AGE correlates with severity and also which cytokines are most effective and most significant