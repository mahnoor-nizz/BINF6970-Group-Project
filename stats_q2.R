# === PACKAGES USED ========
library(readxl)
library(glmnet)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)
library(caret)

### TO-DO: print or cat important statistics after each analysis
### TO-DO: need analysis or comparison between 10 and 20-fold (e.g. their predictive performance metrics + which cytokines they selected) --> put into combined visualization

# === 1.1| DATA CLEANING ========
COV <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")
str(COV) # dataset contains measurements of 28 analytes from 76 COVID-19-positive patients; response variable: severity (binary, mild vs. severe)

# Checking for missing values
sum(is.na(COV)) # no missing values

# Convert severity to a factor variable
# COV$Severity <- ifelse(COV$Severirty == "Severe", 1, 0)
COV$Severity <- as.factor(COV$Severirty)
table(COV$Severity) # 0 = mild, 1 = severe

# Remove Patient Number and Severirty
COV <- COV[, !(names(COV) %in% c("Severirty", "Patient Number"))]
str(COV) # SEX is still character

# Convert SEX to numeric
COV$SEX <- ifelse(COV$SEX == "F", 1, 0)
table(COV$SEX)  

dim(COV)       # 76 x 31
colnames(COV)  

# === 1.2| EXPLORATORY DATA ANALYSIS ========


# === 2.1 | MODEL MATRIX AND TRAINING/TEST SETS ========
# Create model matrix for glmnet including all effects 
X <- model.matrix(Severity ~., COV)[, -1]
dim(X) # 30 predictor variables
colnames(X)

Y <- COV$Severity # severity set as response variable

# Checking predictor standard deviations
plot(as.ts(apply(X, 2, sd)), xlab = "Covariate", ylab = "SD",
     main = "Standard Deviations of Predictors")
# Very variable so standardization is needed; glmnet automatically does standardization


# Create train (75%) & test (25%) sets
set.seed(1717)
test_index <- sample.int(nrow(X), round(nrow(X) * 0.25), replace = FALSE)
X_train <- X[-test_index, ]
Y_train <- Y[-test_index]
X_test <- X[test_index, ]
Y_test <- Y[test_index]

# Change to data partition maybe
set.seed(1717)
train_index <- createDataPartition(Y, p = 0.75, list = FALSE)
X_train2 <- X[train_index, ]
Y_train2 <- Y[train_index]
X_test2 <- X[-train_index, ]
Y_test2 <- Y[-train_index]

## Confirm class balance in train and test
table(Y_train)
table(Y_test)

table(Y_train2)
table(Y_test2)

# === 2.2 | GRID SEARCH FOR ALPHA ========
# Function to perform grid search for optimal alpha value with different n-folds 
alpha_search <- function (X, Y, nfolds){
  alpha_grid <- seq(0, 1, .01)
  results <- data.frame()
  set.seed(1717)
  for (a in alpha_grid) {
  cv_model <- cv.glmnet(x = X, 
                        y = Y, 
                        alpha = a, 
                        nfolds = nfolds,
                        family = "binomial",
                        keep = T,
                        standardize = T,
                        foldid = sample(rep(1:nfolds, length.out = length(Y))),
                        type.measure = "deviance")  # Accuracy measure is deviance because AUC could not be used; number of observations is too small per fold

  # TO-DO: change to keep = T, maybe add foldid
  
  # Store results
  results <- rbind(results, data.frame(
    alpha = a,
    lambda_min = cv_model$lambda.min,
    lambda_1se = cv_model$lambda.1se,
    cvm_min = min(cv_model$cvm),
    cvm_1se = cv_model$cvm[which(cv_model$lambda == cv_model$lambda.1se)]
  ))
  }
  # Extract optimal alpha and lambda
  best_min <- results$alpha[which.min(results$cvm_min)]
  best_1se <- results$alpha[which.min(results$cvm_1se)]
  
  return(list(results = results,
              best_min = best_min,
              best_1se = best_1se))
}

# Perform grid search and obtain optimal alpha for each k-fold CV
#alpha_10 <- alpha_search(X_train, Y_train, 10)
#best_min_10 <- alpha_10$best_min
#best_1se_10 <- alpha_10$best_1se

alpha_10 <- alpha_search(X_train2, Y_train2, 10)
best_min_10 <- alpha_10$best_min
best_1se_10 <- alpha_10$best_1se

#alpha_20 <- alpha_search(20)
#best_min_20 <- alpha_20$best_min
#best_1se_20 <- alpha_20$best_1se

alpha_20 <- alpha_search(X_train2, Y_train2, 20)
best_min_20 <- alpha_20$best_min
best_1se_20 <- alpha_20$best_1se


# Visualize alpha grid search to decide between alpha values
ggplot(alpha_10$results[-1, ], aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), linewidth = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), linewidth = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0 = Ridge, 1 = Lasso)",
       y = "CV Error (Deviance)",
       color = "Lambda Type") +
  theme_minimal()

ggplot(alpha_20$results, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), linewidth = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), linewidth = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha",
       x = "Alpha (0 = Ridge, 1 = Lasso)",
       y = "CV Error (Deviance)",
       color = "Lambda Type") +
  theme_minimal()
# alpha_min will be used as it has lowest CV error and will thus maximize predictive accuracy

# Helper functions for next steps
# Function to determine cutoffs for sensitivity, specificity, and classification
get_cutoff <- function(AUC) {
  snsp <- cbind(AUC$sensitivities, AUC$specificities)
  indx <- which.max(apply(snsp, 1, min))
  return(data.frame(cutoff = AUC$threshold[indx],
                    sensitivity = snsp[indx, 1],
                    specificity = snsp[indx, 2]))
} # get_metrics already gets sn and sp, so might be redundant for this function
# I agree I had it in there before I changed the get metrics function

# Function to obtain metrics from confusion matrix
get_metrics <- function(mat) {
  sn <- mat[2, 2] / sum(mat[2, ])
  sp <- mat[1, 1] / sum(mat[1, ])
  acc <- (mat[1,1] + mat[2,2]) / sum(mat)
  return(unlist(list(sensitivity = sn, specificity = sp, accuracy = acc)))
}

# === 3.1 | ELASTIC NET MODEL (10-FOLD) ========
set.seed(1717)
cv_10 <- cv.glmnet(X_train2, Y_train2, nfolds = 10, family = "binomial", alpha = best_min_10, type.measure = "deviance")
plot(cv_10)
title(main = paste0("10-Fold CV (alpha = ", best_min_10, ")"), line = 2.3)

prd_train_10 <- predict(cv_10, newx = X_train2, type = "response", s = cv_10$lambda.min)[, 1]
prd_test_10 <- predict(cv_10, newx = X_test2, type = "response", s = cv_10$lambda.min)[, 1]

# Test prediction accuracy with AUC-ROC curves(Change plots to be made with ggplot and overlay curves)
AUC_train_10 <- roc(Y_train2, prd_train_10)
AUC_test_10  <- roc(Y_test2,  prd_test_10)
auc(AUC_train_10)
auc(AUC_test_10)

# Determine cutoffs for model with 10-folds
cut_train_10 <- get_cutoff(AUC_train_10)
cut_test_10 <- get_cutoff(AUC_test_10)

cut_train_10
cut_test_10

# Make confusion matrix for model
conf_10 <- table(y = Y_test2, yhat = as.numeric(prd_test_10 > cut_test_10$cutoff))
conf_10
get_metrics(conf_10)

# Plot ROC curve with optimal cutoff; TO-DO: make ggplot versions (change color and linewidth and stuff)
plot(AUC_train_10, main  = "ROC Curve - Training Set (10-Fold) with Optimal Cutoff")
abline(h = cut_train_10$sensitivity, v = cut_train_10$specificity, col = "blue", lty = 2)

plot(AUC_test_10, main = "ROC Curve - Test Set (10-Fold) with Optimal Cutoff")
abline(h = cut_test_10$sensitivity, v = cut_test_10$specificity, col = "blue", lty = 2)

ggroc(AUC_train_10, linewidth = 1.5) +
  labs(title = "ROC Curve - Training Set (10-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
               linetype = "dashed",
               color = "grey") +
  geom_hline(yintercept = cut_train_10$specificity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_train_10$sensitivity,
             linetype = "dashed",
             color = "blue") +
  theme_minimal() 
# add labels, title, either put cutoffs on graph or legend for sn and sp, maybe even put area under curve; could also put both curves in one graph

# Obtain most significant coefficients as predicted by model
coef(cv_10, s = cv_10$lambda.min)
coef_10_min <- coef(cv_10,s = cv_10$lambda.min)[-1,1]
coef_10_min[coef_10_min!=0] 

selected_10  <- sort(abs(coef_10_min[coef_10_min != 0]), decreasing = TRUE)
selected_10 
# With 10-fold CV and using min for alpha and lambda, elastic net regression selected 12 relevant variables
# Notes: AGE has the most influence; all are positive so increasing these by 1 unit increases severity by the corresponding coefficient value 

# Ranking of top predictors TO-DO: make ggplot version
par(mar = c(9, 4, 4, 2))
barplot(selected_10, las = 2, col = "navy",
        main = "Top Predictors - 10-Fold Elastic Net (lambda.min)",
        ylab = "Coefficient", cex.names = 0.75)
par(mar = c(5, 4, 4, 2))

# === 3.2 | ELASTIC NET MODEL (20-FOLD) ========
set.seed(1717)
cv_20 <- cv.glmnet(X_train2, Y_train2, nfolds = 20, family = "binomial", alpha = best_min_20, type.measure = "deviance")
plot(cv_20)
title(main = paste0("20-Fold CV (alpha = ", best_min_20, ")"), line = 2.5)

prd_train_20 <- predict(cv_20, newx = X_train2, type = "response", s = cv_20$lambda.min)[, 1]
prd_test_20  <- predict(cv_20, newx = X_test2,  type = "response", s = cv_20$lambda.min)[, 1]

# Test prediction accuracy with AUC-ROC curves
AUC_train_20 <- roc(Y_train2, prd_train_20)
AUC_test_20  <- roc(Y_test2,  prd_test_20)
auc(AUC_train_20)
auc(AUC_test_20)

# Determine cutoffs for model with 10-folds
cut_train_20 <- get_cutoff(AUC_train_20)
cut_test_20 <- get_cutoff(AUC_test_20)
cut_train_20
cut_test_20

# Make confusion matrix for model
conf_20 <- table(y = Y_test2, yhat = as.numeric(prd_test_20 > cut_test_20$cutoff))
conf_20
get_metrics(conf_20)

# Plot ROC curve with optimal cutoff
ggroc(AUC_train_20, linewidth = 1.5) +
  labs(title = "ROC Curve - Training Set (20-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
              linetype = "dashed",
              color = "grey") +
  geom_hline(yintercept = cut_train_20$specificity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_train_20$sensitivity,
             linetype = "dashed",
             color = "blue") +
  theme_minimal() 

# Obtain most significant coefficients as predicted by model
coef(cv_20, s = cv_20$lambda.min)
coef_20_min <- coef(cv_20,s = cv_20$lambda.min)[-1,1]
coef_20_min[coef_20_min!=0] 

selected_20  <- sort(abs(coef_20_min[coef_20_min != 0]), decreasing = TRUE)
selected_20 
# With 20-fold CV and using min for alpha and lambda, elastic net regression selected 13 relevant variables
# Notes: AGE has the most influence; all are positive so increasing these by 1 unit increases severity by the corresponding coefficient value

# Ranking of top predictors TO-DO: make ggplot version
par(mar = c(9, 4, 4, 2))
barplot(selected_20, las = 2, col = "maroon",
        main = "Top Predictors - 20-Fold Elastic Net (lambda.min)",
        ylab = "Coefficient", cex.names = 0.75)
par(mar = c(5, 4, 4, 2))

# Comparison plot for 10-fold vs 20-fold
plot(AUC_test_10, col = "navy", main = "Test Set ROC: 10-Fold vs 20-Fold")
plot(AUC_test_20, col = "maroon", add = TRUE)
legend("bottomright",
       legend = c(paste0("10-Fold AUC = ", round(auc(AUC_test_10), 3)),
                  paste0("20-Fold AUC = ", round(auc(AUC_test_20), 3))),
       col = c("navy", "maroon"), lwd = 2, bty = "n")

# === 4 | STATISTICAL ANALYSIS OF PREDICTORS ======== #TO-DO: check if this is correct
# Wilcoxon test for age 
wilcox.test(AGE ~ Severity, data = COV) #p-value = 0.02017
boxplot(AGE ~ Severity, data = COV)

# Model with age as only predictor
model_age <- glm(Severity ~ AGE, data = COV, family = "binomial")
summary(model_age)
### 4.3% decrease in odds per 1 unit increase in age??? 
# yeah that kinda makes sense based on the graph below, I think its because most ages in the mild range were over 50 while severe had ones in their 20s. the correlation thing under also gives a negative value 

# Correlation of age with predicted probability in 10 and 20 fold model
prd_all_10 <- predict(cv_10, newx = X, type = "response", s = cv_10$lambda.min)[, 1]
cor(COV$AGE, prd_all_10)

prd_all_20 <- predict(cv_20, newx = X, type = "response", s = cv_20$lambda.min)[, 1]
cor(COV$AGE, prd_all_20)

# Scatterplot of age vs probability
ggplot(data.frame(Age = COV$AGE, Prob = prd_all_10, Severity = factor(COV$Severity, labels = c("Mild","Severe"))), aes(x = Age, y = Prob, color = Severity)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", color = "black") +
  scale_color_manual(values = c("Mild" = "navy", "Severe" = "maroon")) +
  labs(title = "Age vs Probability of Severe COVID-19",
       x = "Age (years)", y = "Probability of Severity") +
  theme_bw()
### I don't get this graph, please explain later
# its like the one we did in another class I made it based on that but idk if it is actually meaningful lol, it just shows that when age increases your probability of severeity decreases, 
#which is supported by both the correlation thing and your model with age as the only predictor


### I don't think we really need to do this part since she never taught us but I'll leave this here for now
#I'm also good with getting rid of it i dont think its necessary
# Statistical significance of cytokines
wilcox_results <- data.frame()

for (i in names(selected_10)) {if (i %in% colnames(X)) {
  test <- wilcox.test(X[, i] ~ COV$Severity)
  
  wilcox_results <- rbind(wilcox_results, data.frame(
    Predictor = i,
    W_statistic = test$statistic,
    p_value = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No")))}}

print(wilcox_results[order(wilcox_results$p_value), ])

# Odds ratio of cytokines (doing with 10 first before 20 to validate)
odds_10 <- cbind(X[, names(selected_10)], Severity = COV$Severity) 
odds_10 <- as.data.frame(odds_10)
model_10 <- glm(Severity ~., data = odds_10, family = "binomial")
summary(model_10)

or_table <- exp(cbind(OR = coef(model_10), confint(model_10)))
print(round(or_table, 3))


