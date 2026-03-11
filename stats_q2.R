# === PACKAGES USED ========
library(readxl)
library(glmnet)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)

# === 1 | DATA CLEANING AND EDA ========
COV <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
str(COV)

# Checking for missing values
sum(is.na(COV))
# no missing values

table(COV$Severirty)

# Convert severity to a factor variable
COV$Severity <- ifelse(COV$Severirty == "Severe", 1, 0)
table(COV$Severity) # 0 = mild, 1 = severe

# Remove Patient Number and Severity
COV <- COV[, !(names(COV) %in% c("Severirty", "Patient Number"))]
str(COV)
# SEX is still chr

# Convert SEX to numeric
COV$SEX <- ifelse(COV$SEX == "F", 1, 0)
table(COV$SEX)  

dim(COV)       # 76 x 31
colnames(COV)  

# EDA Age distribution by severity


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
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1) +
  geom_vline(xintercept = best_min_10, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = best_1se_10, linetype = "dashed", color = "red") +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 2) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 2) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (MSE)",
       color = "Lambda Type") +
  theme_minimal()

ggplot(alpha_20, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 2) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 2) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "CV Error (MSE)",
       color = "Lambda Type") +
  theme_minimal()

## Decide between alpha values

## Create functions to streamline comparison between 10 and 20 folds

# this one gets sensitivity and specificity from confusion matrix (made later down)
sn.sp <- function(mat) {
  sn <- mat[2, 2] / sum(mat[2, ])
  sp <- mat[1, 1] / sum(mat[1, ])
  acc <- (mat[1,1] + mat[2,2]) / sum(mat)

  return(unlist(list(sensitivity = sn, specificity = sp, accuracy = acc)))
}

# getting cutoffs function (same as yours but I added return sensitivity and specificity too)
get_cutoff <- function(auc_obj) {
  snsp <- cbind(auc_obj$sensitivities, auc_obj$specificities)
  indx <- which.max(apply(snsp, 1, min))
  
  return(list(cutoff = auc_obj$thresholds[indx],
              sensitivity = snsp[indx, 1],
              specificity = snsp[indx, 2]))}



# Elastic net model with 10 folds
set.seed(1717) #added set seed here again idk if its nessesary but the graph was coming out different for me every time lol
cv_10 <- cv.glmnet(X_train, Y_train, nfolds = 10, family = "binomial", alpha = best_min_10, type.measure = "deviance")
plot(cv_10)
title(main = paste0("10-Fold CV (alpha = ", best_min_10, ")"), line = 2.3)


#prediction (i changed the variable names so prd_train_10 to differenciate from prd_train_20, )
prd_train_10 <- predict(cv_10, newx = X_train, type = "response", s = cv_10$lambda.min)[, 1]
prd_test_10 <- predict(cv_10, newx = X_test, type = "response", s = cv_10$lambda.min)[, 1]

# Testing prediction accuracy (Change plots to be made with ggplot and overlay curves)
AUC_train_10 <- roc(Y_train, prd_train_10)
AUC_test_10 <- roc(Y_test, prd_test_10)
# removed the plots for now, will add them with optimal cuttoffs after


auc(AUC_train_10)
auc(AUC_test_10)

# Check for performance of test set
#I think this part would be min(cv_10$cvm) no? maybe im wrong but we used min earlier multiple times since deviance measures error right, so lower deviance would be better?

#max(cv_10$cvm)
#which.max(cv_10$cvm)
#cv_10$cvlo[which.max(cv_10$cvm)];cv_10$cvup[which.max(cv_10$cvm)]

#heres a version w min instead

#best cv deviance (lamda min)
min(cv_10$cvm)
cv_10$cvlo[which.min(cv_10$cvm)]
cv_10$cvup[which.min(cv_10$cvm)]

# Determining cutoffs for sensitivity, specificity, and classification; make function to test both; return dataframe with all
#cutoffs_snsp <- function(AUC) {
#  snsp <- cbind(AUC$sensitivities, AUC$specificities)
#  indx <- which.max(apply(snsp, 1, min))
#  cutoff <- AUC$threshold[indx]
#  return(cutoff)
#}

#cutoffs_snsp(AUC_train_10)
#cutoffs_snsp(AUC_test_10)

#using function that gives sensitivity & specificity
get_cutoff(AUC_test_10)
get_cutoff(AUC_train_10)

cut_10 <- get_cutoff(AUC_test_10)

# Make confusion matrix? Oui for test set
conf_10 <- table(y = Y_test, yhat = as.numeric(prd_test_10 > cut_10$cutoff))
print(conf_10)
print(sn.sp(conf_10))


# This is all done in the function from earlier, I commented it out incase you wanted to keep it over the function.

# Training set
snsp_train_10 <- cbind(AUC_train_10$sensitivities, AUC_train_10$specificities)
#snsp_train

indx_train_10 <- which.max(apply(snsp_train_10, 1, min))
#indx_train
#snsp_train[indx_train, ]

#train_cutoff <- AUC_train$thresholds[indx_train]
#train_cutoff

# Test set
snsp_test_10 <- cbind(AUC_test_10$sensitivities, AUC_test_10$specificities)
#snsp_test_10

indx_test_10 <- which.max(apply(snsp_test_10, 1, min))
#indx_test
#snsp_test[indx_test, ]

#test_cutoff <- AUC_test$thresholds[indx_test]
#test_cutoff

# ROC curve with optimal cutoff

# Make ggplot versions
par(mfrow=c(1,2))
plot(AUC_test_10, main = "ROC Curve - Test Set (10-Fold) with Optimal Cutoff")
abline(h = snsp_test_10[indx_test_10, 1],
       v = snsp_test_10[indx_test_10, 2], col = "blue", lty = 2)

plot(AUC_train_10, main = "ROC Curve - Train Set (10-Fold) with Optimal Cutoff")
abline(h = snsp_train_10[indx_train_10, 1],
       v = snsp_train_10[indx_train_10, 2], col = "blue", lty = 2)

par(mfrow=c(1,1))

# Variable selection (lambda.min)
coef(cv_10, s = cv_10$lambda.min)
coef_10_min <- coef(cv_10,s = cv_10$lambda.min)[-1,1]
coef_10_min[coef_10_min!=0] 
names(coef_10_min[coef_10_min!=0])


# Variable selection (lambda.1se)
coef(cv_10, s = cv_10$lambda.1se)
coef_10_1se <- coef(cv_10,s = cv_10$lambda.1se)[-1,1]
#all reduced to 0 

selected_10_1se  <- sort(abs(coef_10_1se[coef_10_1se != 0]), decreasing = TRUE)

selected_10  <- sort(abs(coef_10_min[coef_10_min != 0]), decreasing = TRUE)

#these are the top ranked predictors from highest coef to lowest
names(selected_10)


# ranking of top predictors
par(mar = c(9, 4, 4, 2))
barplot(selected_10, las = 2, col = "navy",
        main = "Top Predictors - 10-Fold Elastic Net (lambda.min)",
        ylab = "Coefficient", cex.names = 0.75)
par(mar = c(5, 4, 4, 2))


# Elastic net model with 20 folds
## Make function that does all the same steps as 10-fold
set.seed(1717)
cv_20 <- cv.glmnet(X_train, Y_train, nfolds = 20, family = "binomial", alpha = best_min_20, type.measure = "deviance")
plot(cv_20)
title(main = paste0("20-Fold CV (alpha = ", best_min_20, ")"), line = 2.5)

# Predictions
prd_train_20 <- predict(cv_20, newx = X_train, type = "response", s = cv_20$lambda.min)[, 1]
prd_test_20  <- predict(cv_20, newx = X_test,  type = "response", s = cv_20$lambda.min)[, 1]


#ROC & AUC
AUC_train_20 <- roc(Y_train, prd_train_20)
AUC_test_20  <- roc(Y_test,  prd_test_20)
auc(AUC_train_20)
auc(AUC_test_20)
min(cv_20$cvm)

# optimal cutoff
cut_20 <- get_cutoff(AUC_test_20)
cut_20

# Confusion matrix (test set)
conf_20 <- table(Y = Y_test, Yhat = as.numeric(prd_test_20 > cut_20$cutoff))
sn.sp(conf_20)

#ROC Curve
snsp_train_20 <- cbind(AUC_train_20$sensitivities, AUC_train_20$specificities)
indx_train_20 <- which.max(apply(snsp_train_20, 1, min))

snsp_test_20 <- cbind(AUC_test_20$sensitivities, AUC_test_20$specificities)
indx_test_20 <- which.max(apply(snsp_test_20, 1, min))


par(mfrow=c(1,2))
plot(AUC_test_20, main = "ROC Curve - Test Set (20-Fold) with Optimal Cutoff")
abline(h = snsp_test_20[indx_test_20, 1],
       v = snsp_test_20[indx_test_20, 2], col = "blue", lty = 2)

plot(AUC_train_20, main = "ROC Curve - Train Set (20-Fold) with Optimal Cutoff")
abline(h = snsp_train_10[indx_train_20, 1],
       v = snsp_train_10[indx_train_20, 2], col = "blue", lty = 2)

par(mfrow=c(1,1))


# variable selection
coef(cv_20, s = cv_20$lambda.min)
coef_20_min <- coef(cv_20,s = cv_20$lambda.min)[-1,1]
coef_20_min[coef_20_min!=0] 
names(coef_20_min[coef_20_min!=0])

selected_20  <- sort(abs(coef_20_min[coef_20_min != 0]), decreasing = TRUE)

#for lamda.1se
coef_20_1se <- coef(cv_20,s = cv_20$lambda.1se)[-1,1]
coef_20_1se[coef_20_1se!=0] 
names(coef_20_1se[coef_20_1se!=0])
selected_20_1se <- sort(abs(coef_20_1se[coef_20_1se != 0]), decreasing = TRUE)
# retained 3 variables, PTX3, IL-6, and TNF-α

par(mar = c(9, 4, 4, 2))
barplot(selected_20, las = 2, col = "maroon",
        main = "Top Predictors - 20-Fold Elastic Net (lambda.min)",
        ylab = "Coefficient", cex.names = 0.75)
par(mar = c(5, 4, 4, 2))


## FOR MS. MAHNOOR NIZAMANI: Make some plots and do some analysis for how AGE correlates with severity and also which cytokines are most effective and most significant


# plot for camparison, 10 fold vs 20 fold

plot(AUC_test_10, col = "navy", main = "Test Set ROC: 10-Fold vs 20-Fold")
plot(AUC_test_20, col = "maroon", add = TRUE)
legend("bottomright",
       legend = c(paste0("10-Fold AUC = ", round(auc(AUC_test_10), 3)),
                  paste0("20-Fold AUC = ", round(auc(AUC_test_20), 3))),
       col = c("navy", "maroon"), lwd = 2, bty = "n")

#omg the same literally identical I hope i didnt fuck something up

# summary comparison table
# ill do that after


#analysis for how AGE correlates with severity

# Wilcoxon test 
wilcox.test(AGE ~ Severity, data = COV)

# i used this instead of a T test bc age is a wee bit skewed more older people in severeity = 0

#can see using this
#boxplot(AGE ~ Severity, data = COV)

#p-value = 0.02017

#correlation of age w predicted probability in 10 fold model
prd_all_10 <- predict(cv_10, newx = X, type = "response", s = cv_10$lambda.min)[, 1]

cor(COV$AGE, prd_all_10)
#-0.3660866


#scatterplot of age vs probability

ggplot(data.frame(Age = COV$AGE, Prob = prd_all_10, Severity = factor(COV$Severity, labels = c("Mild","Severe"))), aes(x = Age, y = Prob, color = Severity)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", color = "black") +
  scale_color_manual(values = c("Mild" = "navy", "Severe" = "maroon")) +
  labs(title = "Age vs Probability of Severe COVID-19",
       x = "Age (years)", y = "Probability of Severity") +
  theme_bw()


#which cytokines are most effective and most significant

wilcox_results <- data.frame()

for ( i in names(selected_10)) {if ( i %in% colnames(X)) {
  test <- wilcox.test(X[, i] ~ COV$Severity)
  
  wilcox_results <- rbind(wilcox_results, data.frame(
    Predictor = i,
    W_statistic = test$statistic,
    p_value = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No")))}}

print(wilcox_results[order(wilcox_results$p_value), ])

