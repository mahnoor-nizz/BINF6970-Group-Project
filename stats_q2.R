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
cat("Number of missing values:", sum(is.na(COV))) # no missing values

# Convert severity to a factor variable
COV$Severity <- as.factor(COV$Severirty)
cat("Severity class distribution:")
table(COV$Severity) # 32 mild, 44 severe

# Remove Patient Number and Severirty
COV <- COV[, !(names(COV) %in% c("Severirty", "Patient Number"))]
str(COV) # SEX is still character

# Convert SEX to numeric
COV$SEX <- ifelse(COV$SEX == "F", 1, 0)
table(COV$SEX)  

cat("Dataset dimensions:", dim(COV))
cat("Column names:")
print(colnames(COV))

# === 1.2| EXPLORATORY DATA ANALYSIS ========
# Age distribution by severity boxplot 

ggplot(data.frame(Age = COV$AGE, Severity = COV$Severity), aes(x = Severity, y = Age, fill = Severity)) +
  geom_boxplot() +
  labs(title = "Age Distribution by COVID-19 Severity", x = "Severity Group", y = "Age (years)") +
  theme_bw() +
  theme(legend.position = "none")

 #Cytokine distribution by severity, single combined plot to be added

# === 1.3 Assumption checks ========

#Class Balance
sev_tab <- table(COV$Severity)
print(sev_tab)
cat("Class ratio (mild:severe):", round(sev_tab[1] / sev_tab[2], 2), "\n")
#Classes reasonably balanced, so no resampling required?

#Multicollinearity check maybe a correlation heatmap of the predictors to be added

# === 2.1 | MODEL MATRIX AND TRAINING/TEST SETS ========
# Create model matrix for glmnet including all effects 
X <- model.matrix(Severity ~., COV)[, -1]
cat("Predictor matrix dimensions:", dim(X))
cat("Predictor names:"); print(colnames(X))


Y <- COV$Severity # severity set as response variable

# Checking predictor standard deviations
plot(as.ts(apply(X, 2, sd)), xlab = "Covariate", ylab = "SD",
     main = "Standard Deviations of Predictors")
# Very variable so standardization is needed; glmnet automatically does standardization


# Create train (75%) & test (25%) sets
set.seed(1717)
train_index <- createDataPartition(Y, p = 0.75, list = FALSE)
X_train <- X[train_index, ]
Y_train <- Y[train_index]
X_test <- X[-train_index, ]
Y_test <- Y[-train_index]

## Confirm class balance in train and test
cat("Training set class distribution:"); print(table(Y_train)) # 57 total, 42% in mild, 58% in severe
cat("Test set class distribution:"); print(table(Y_test)) # 19 total, 42% in mild, 58% in severe

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
                        type.measure = "deviance")  # Accuracy measure is deviance because AUC could not be used; number of observations is too small per fold
  
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
alpha_10 <- alpha_search(X_train, Y_train, 10)
best_min_10 <- alpha_10$best_min
best_1se_10 <- alpha_10$best_1se
cat("Best 10 fold alpha (lambda.min):", best_min_10)
cat("Best 10 foldalpha (lambda.1se):", best_1se_10)

alpha_20 <- alpha_search(X_train, Y_train, 20)
best_min_20 <- alpha_20$best_min
best_1se_20 <- alpha_20$best_1se
cat("Best 20 fold alpha (lambda.min):", best_min_20)
cat("Best 20 fold alpha (lambda.1se):", best_1se_20)

# Visualize alpha grid search to decide between alpha values
ggplot(alpha_10$results, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), linewidth = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), linewidth = 1.2) +
  geom_vline(xintercept = best_min_10, linetype = "dashed") +
  labs(title = "CV Error vs Alpha (10-Fold CV)",
       x = "Alpha (0 = Ridge, 1 = Lasso)",
       y = "CV Error (Deviance)",
       color = "Lambda Type") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

ggplot(alpha_20$results, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), linewidth = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), linewidth = 1.2) +
  geom_vline(xintercept = best_min_20, linetype = "dashed") +
  labs(title = "CV Error vs Alpha (20-Fold CV)",
       x = "Alpha (0 = Ridge, 1 = Lasso)",
       y = "CV Error (Deviance)",
       color = "Lambda Type") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 
# alpha_min will be used as it has the lowest CV error and will thus maximize predictive accuracy

# Helper functions for next steps
# Function to determine cutoffs for sensitivity, specificity, and classification
get_cutoff <- function(AUC) {
  snsp <- cbind(AUC$sensitivities, AUC$specificities)
  indx <- which.max(apply(snsp, 1, min))
  return(data.frame(sensitivity = snsp[indx, 1],
                    specificity = snsp[indx, 2],
                    cutoff = AUC$threshold[indx]))
} 

# Function to obtain metrics from confusion matrix; positive is set to second factor
get_metrics <- function(mat, pos = 2) {
  neg <- if (pos == 2) 1 else 2
  sn <- mat[pos, pos] / sum(mat[pos, ])
  sp <- mat[neg, neg] / sum(mat[neg, ])
  acc <- (mat[neg, neg] + mat[pos, pos]) / sum(mat)
  return(unlist(list(sensitivity = sn, specificity = sp, accuracy = acc)))
}

# === 3.1 | ELASTIC NET MODEL (10-FOLD) ========
set.seed(1717)
cv_10 <- cv.glmnet(X_train, Y_train, nfolds = 10, family = "binomial", alpha = best_min_10, type.measure = "deviance")
plot(cv_10)
title(main = paste0("10-Fold CV (alpha = ", best_min_10, ")"), line = 2.3)

prd_train_10 <- predict(cv_10, newx = X_train, type = "response", s = cv_10$lambda.min)[, 1]
prd_test_10 <- predict(cv_10, newx = X_test, type = "response", s = cv_10$lambda.min)[, 1]

# Test prediction accuracy with AUC-ROC curves(Change plots to be made with ggplot and overlay curves)
AUC_train_10 <- roc(Y_train, prd_train_10)
AUC_test_10  <- roc(Y_test,  prd_test_10)
cat("Training AUC:", auc(AUC_train_10))
cat("Test AUC:    ", auc(AUC_test_10))


# Determine cutoffs for model with 10-folds
cut_train_10 <- get_cutoff(AUC_train_10)
cut_test_10 <- get_cutoff(AUC_test_10)
cat("Training set optimal cutoff:"); print(cut_train_10)
cat("Test set optimal cutoff:"); print(cut_test_10)

# Make confusion matrix for model
conf_10 <- table(y = Y_test, yhat = as.numeric(prd_test_10 > cut_test_10$cutoff))
metrics_10 <- get_metrics(conf_10)
cat("Confusion Matrix (10-Fold, Test Set):") 
print(conf_10)
cat("Performance Metrics (10-Fold, Test Set):") 
print(metrics_10)

# Plot ROC curve with optimal cutoff
# Training set
ggroc(AUC_train_10, linewidth = 1.1) +
  labs(title = "ROC Curve - Training Set (10-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
               color = "grey") +
  geom_hline(yintercept = cut_train_10$sensitivity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_train_10$specificity,
             linetype = "dashed",
             color = "blue") +
  annotate("text", 
           x = 0.84,
           y = 0.25,
           label = "paste(Specificity, \" = 0.917\")", parse = T) +
  annotate("text", 
           x = 0.25,
           y = 0.93,
           label = "paste(Sensitivity, \" = 0.970\")", parse = T) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

# Test set
ggroc(AUC_test_10, linewidth = 1.1) +
  labs(title = "ROC Curve - Test Set (10-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
              color = "grey") +
  geom_hline(yintercept = cut_test_10$sensitivity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_test_10$specificity,
             linetype = "dashed",
             color = "blue") +
  annotate("text", 
           x = 0.8,
           y = 0.3,
           label = "paste(Specificity, \" = 0.875\")", parse = T) +
  annotate("text", 
           x = 0.125,
           y = 0.78,
           label = "paste(Sensitivity, \" = 0.818\")", parse = T) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

# Obtain most significant coefficients as predicted by model
coef(cv_10, s = cv_10$lambda.min)
coef_10_min <- coef(cv_10,s = cv_10$lambda.min)[-1,1]

# With 10-fold CV and using min for alpha and lambda, elastic net regression selected 9 relevant variables; 7 are negative and 2 are positive

select_10  <- sort(abs(coef_10_min[coef_10_min != 0]), decreasing = TRUE)
cat("Selected predictors (10-Fold, lambda.min):") 
print(round(coef_10_min[coef_10_min != 0], 4))
cat("Number of non-zero predictors:", length(select_10))
cat("Most influential predictor:", names(select_10)[1])

selected_10 <- data.frame(Predictor = names(select_10),
                          Value = select_10)
selected_10$Predictor <- factor(selected_10$Predictor, levels = selected_10$Predictor)
rownames(selected_10) <- NULL
# Notes: TNF-a has the most influence

# Ranking of top predictors 
ggplot(selected_10, aes(x = Predictor, y = Value)) +
  geom_col(fill = "navy") +
  geom_text(aes(label = signif(Value, 3)), vjust = -0.5) +
  labs(title = "Top Selected Predictors - 10-Fold Elastic Net (lambda.min)",
       x = "Predictor",
       y = "Coefficient score") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1)) 


# what are your thoughts on this kinda graph instead
df_coef_10 <- data.frame(Predictor = names(select_10), Coefficient = select_10, Direction = ifelse(coef_10_min[names(select_10)] > 0, "Positive", "Negative"))

ggplot(df_coef_10, aes(x = reorder(Predictor, Coefficient), y = Coefficient, fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "navy", "Negative" = "maroon")) +
  labs(title = "Top Predictors — 10-Fold Elastic Net (lambda.min)", x = "Predictor", y = "Coefficient(absolute value)", caption = "Red = associated with higher severity; Blue = associated with lower severity") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#it highlights what was pos and negative >:) 

# === 3.2 | ELASTIC NET MODEL (20-FOLD) ========
set.seed(1717)
cv_20 <- cv.glmnet(X_train, Y_train, nfolds = 20, family = "binomial", alpha = best_min_20, type.measure = "deviance")
plot(cv_20)
title(main = paste0("20-Fold CV (alpha = ", best_min_20, ")"), line = 2.5)

prd_train_20 <- predict(cv_20, newx = X_train, type = "response", s = cv_20$lambda.min)[, 1]
prd_test_20  <- predict(cv_20, newx = X_test,  type = "response", s = cv_20$lambda.min)[, 1]

# Test prediction accuracy with AUC-ROC curves
AUC_train_20 <- roc(Y_train, prd_train_20)
AUC_test_20  <- roc(Y_test,  prd_test_20)
cat("Training AUC:", auc(AUC_train_20))
cat("Test AUC:    ", auc(AUC_test_20))

# Determine cutoffs for model with 20-folds
cut_train_20 <- get_cutoff(AUC_train_20)
cut_test_20 <- get_cutoff(AUC_test_20)
cat("Training set optimal cutoff:"); print(cut_train_20)
cat("Test set optimal cutoff:"); print(cut_test_20)


# Make confusion matrix for model
conf_20 <- table(y = Y_test, yhat = as.numeric(prd_test_20 > cut_test_20$cutoff))
conf_20
metrics_10 <- get_metrics(conf_20)
cat("Confusion Matrix (20-Fold):\n"); print(conf_20)
cat("Performance Metrics (20-Fold):\n"); print(get_metrics(conf_20))

# Plot ROC curve with optimal cutoff
# Training set
ggroc(AUC_train_20, linewidth = 1.1) +
  labs(title = "ROC Curve - Training Set (20-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
              color = "grey") +
  geom_hline(yintercept = cut_train_20$sensitivity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_train_20$specificity,
             linetype = "dashed",
             color = "blue") +
  annotate("text", 
           x = 0.84,
           y = 0.25,
           label = "paste(Specificity, \" = 0.917\")", parse = T) +
  annotate("text", 
           x = 0.25,
           y = 0.84,
           label = "paste(Sensitivity, \" = 0.879\")", parse = T) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

# Test set
ggroc(AUC_test_20, linewidth = 1.1) +
  labs(title = "ROC Curve - Test Set (20-Fold) with Optimal Cutoff",
       x = "False Positive Rate (Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  geom_abline(slope = 1, intercept = 1, 
              color = "grey") +
  geom_hline(yintercept = cut_test_20$sensitivity,
             linetype = "dashed",
             color = "blue") + 
  geom_vline(xintercept = cut_test_20$specificity,
             linetype = "dashed",
             color = "blue") +
  annotate("text", 
           x = 0.8,
           y = 0.3,
           label = "paste(Specificity, \" = 0.875\")", parse = T) +
  annotate("text", 
           x = 0.125,
           y = 0.78,
           label = "paste(Sensitivity, \" = 0.818\")", parse = T) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) 

# Obtain most significant coefficients as predicted by model
coef(cv_20, s = cv_20$lambda.min)
coef_20_min <- coef(cv_20,s = cv_20$lambda.min)[-1,1]
coef_20_min[coef_20_min!=0] 


# With 20-fold CV and using min for alpha and lambda, elastic net regression selected 9 relevant variables; 8 are negative and 1 is positive

select_20  <- sort(abs(coef_20_min[coef_20_min != 0]), decreasing = TRUE)

cat("Selected predictors (20-Fold):")
print(round(coef_20_min[coef_20_min != 0], 4))
cat("Number of non-zero predictors:", length(select_20))
cat("Most influential predictor:", names(select_20)[1])

selected_20 <- data.frame(Predictor = names(select_20),
                          Value = select_20)
selected_20$Predictor <- factor(selected_20$Predictor, levels = selected_20$Predictor)
rownames(selected_20) <- NULL
# Notes: TNF-a has the most influence

# Ranking of top predictors 
ggplot(selected_20, aes(x = Predictor, y = Value)) +
  geom_col(fill = "maroon") +
  geom_text(aes(label = signif(Value, 3)), vjust = -0.5) +
  labs(title = "Top Selected Predictors - 20-Fold Elastic Net (lambda.min)",
       x = "Predictor",
       y = "Coefficient score") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

# Comparison plot for 10-fold vs 20-fold
ggroc(list("10-Fold" = AUC_test_10, "20-Fold" = AUC_test_20), linewidth = 1) +
  scale_color_manual(values = c("10-Fold" = "navy", "20-Fold" = "maroon")) +
  geom_abline(slope = 1, intercept = 1, color = "grey") +
  annotate("text", x = 0.25, y = 0.15, label = paste0("10-Fold AUC = ", round(auc(AUC_test_10), 3)), color = "navy", fontface = "bold") +
  annotate("text", x = 0.25, y = 0.07, label = paste0("20-Fold AUC = ", round(auc(AUC_test_20), 3)), color = "maroon", fontface = "bold") +
  labs(title = "Test Set ROC Curves: 10-Fold vs 20-Fold", x = "Specificity", y = "Sensitivity", color = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# For MN: please compare the predictors in the 10-fold and 20-fold models; idk how but show which ones wre exclusive in one and shared in both
# would a table work here 

comparison_df <- data.frame(
  Metric = c("Alpha (lambda.min)", "Train AUC", "Test AUC", "Test Sensitivity", "Test Specificity", "Test Accuracy", "Predictors Selected"),
  `10-Fold` = c(best_min_10, round(auc(AUC_train_10), 4), round(auc(AUC_test_10), 4), round(metrics_10["sensitivity"], 4), round(metrics_10["specificity"], 4), round(metrics_10["accuracy"], 4), length(selected_10)),
  `20-Fold` = c(best_min_20, round(auc(AUC_train_20), 4), round(auc(AUC_test_20), 4), round(metrics_20["sensitivity"], 4), round(metrics_20["specificity"], 4), round(metrics_20["accuracy"], 4), length(select_20)))
print(comparison_df, row.names = FALSE)

#
#Shared vs exclusive predictors i feel like we have enough graphs ;-;
vars_10 <- names(select_10)
vars_20 <- names(select_20)
shared <- intersect(vars_10, vars_20)
only_10 <- setdiff(vars_10, vars_20)
only_20 <- setdiff(vars_20, vars_10)

cat("Shared predictors (in both models):", paste(shared,  collapse = ", "))
cat("Exclusive to 10-fold model:", paste(only_10, collapse = ", "))
cat("Exclusive to 20-fold model:", paste(only_20, collapse = ", "))


#combined coefficient data frame for if i can figure out how i wanna plot this
all_preds <- union(vars_10, vars_20)
df_compare <- data.frame(Predictor = all_preds, Coef_10 = coef_10_min[all_preds], Coef_20 = coef_20_min[all_preds], 
                         Group = ifelse(all_preds %in% shared, "Shared", ifelse(all_preds %in% only_10, "10-Fold Only", "20-Fold Only")))


# === 4 | STATISTICAL ANALYSIS OF PREDICTORS ======== #TO-DO: check if this is correct
# Wilcoxon test for age 
wilcox.test(AGE ~ Severity, data = COV) #p-value = 0.02017
wtest_age <- wilcox.test(AGE ~ Severity, data = COV)
cat("Wilcoxon test for Age ~ Severity: p-value =", round(wtest_age$p.value, 4))

# Model with age as only predictor
model_age <- glm(Severity ~ AGE, data = COV, family = "binomial")
summary(model_age)
### 4.3% decrease in odds per 1 unit increase in age??? 
# yeah that kinda makes sense based on the graph below, I think its because most ages in the mild range were over 50 while severe had ones in their 20s. the correlation thing under also gives a negative value 

# Correlation of age with predicted probability in 10 and 20 fold model
prd_all_10 <- predict(cv_10, newx = X, type = "response", s = cv_10$lambda.min)[, 1]
prd_all_20 <- predict(cv_20, newx = X, type = "response", s = cv_20$lambda.min)[, 1]

cat("Correlation of Age with predicted severity probability:")
cat("10-Fold model:", round(cor(COV$AGE, prd_all_10), 4))
cat("20-Fold model:", round(cor(COV$AGE, prd_all_20), 4))

# Scatterplot of age vs probability
ggplot(data.frame(Age = COV$AGE, Prob = prd_all_10, Severity = factor(COV$Severity, labels = c("Mild","Severe"))), aes(x = Age, y = Prob, color = Severity)) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", color = "black") +
  scale_color_manual(values = c("Mild" = "navy", "Severe" = "maroon")) +
  labs(title = "Age vs Probability of Severe COVID-19",
       x = "Age (years)", y = "Probability of Severity") +
  theme_bw()

### I don't think we really need to do this part since she never taught us but I'll leave this here for now
#I'm also good with getting rid of it i dont think its necessary
# Statistical significance of cytokines
wilcox_results <- data.frame()

for (i in names(select_10)) {if (i %in% colnames(X)) {
  test <- wilcox.test(X[, i] ~ COV$Severity)
  
  wilcox_results <- rbind(wilcox_results, data.frame(
    Predictor = i,
    W_statistic = test$statistic,
    p_value = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "Yes", "No")))}}

print(wilcox_results[order(wilcox_results$p_value), ])



