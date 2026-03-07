library(readxl) 
library(glmnet)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)

COV <- read_excel("Immunologic profiles of patients with COVID-19.xlsx")

View(COV)
str(COV)

#dataset contains measurements of 28 analytes (cytokines, chemokines,and acute phase proteins) from plasma samples of 76 COVID-19-positive patients. Response variable: disease progression/severity (binary, mild vs. severe)

table(COV$Severirty)

Severity <- ifelse(COV$Severirty == "Severe", 1, 0)
table(Severity)


# Remove Patient Number and Severity
X_raw <- COV[, !(names(COV) %in% c("Severirty", "Patient Number"))]

str(X_raw)

# Convert SEX to numeric
X_raw$SEX <- ifelse(X_raw$SEX == "F", 1, 0)
table(X_raw$SEX)  

#model matrix 
X <- as.matrix(X_raw)

dim(X)       # 76 x 30
colnames(X)  



# Feature Engineering



# model matrix including all effects
f0 <- lm(Severity ~ ., data = data.frame(Severity = Severity, X_raw))

X <- model.matrix(f0)[, -1]

dim(X)
colnames(X)

# predictor standard deviations
plot(as.ts(apply(X, 2, sd)), xlab = "Covariate", ylab = "SD",
     main = "Standard Deviations of Predictors")
# very variable so thats why standerdization is needed




#Train & Test sets


set.seed(1717)


n <- nrow(X)
n_test <- round(n * 0.25)

test.index <- sample.int(n, n_test, replace = FALSE)

cat("Test set size:", n_test)

## Confirm class balance in train and test
table(y[-test.index])
table(y[test.index])


# Elastic Net 10-Fold Cross-Validation

