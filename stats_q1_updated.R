# Stats 6970 - Problem 1
# PCA Analysis of T-Cell Gene Expression Data
# Reference: Holmes et al. (2005)

# Background:
# This analysis investigates how T cells differentiate after thymic selection.
# Two competing models exist:
#   - Linear model: NAI -> EFFE -> MEM
#   - Parallel model: NAI -> EFFE and MEM (simultaneously)

# Dataset: 156 differential expressed genes (Holmes et al., 2005)
# Row name format: [Status]_[CellType]_[ID]
#   Status:   HEA = Healthy, MEL = Melanoma
#   CellType: NAI = Naïve, EFFE = Effector, MEM = Memory
# Each column = a gene's expression measurement (a numeric value indicating how highly/lowly that gene is expressed in that sample)

# ---------------------- Import Libraries ----------------------------
library(ggplot2)
library(ggfortify)

# Load the gene expression file and explore the data
load("geneexpression2.rda")
head(dat)
dim(dat) # Rows = samples, Columns = 156 differentially expressed genes

# Extract cell type counts
cell_types <- ifelse(grepl("NAI", rownames(dat)), "Naive",
  ifelse(grepl("EFFE", rownames(dat)), "Effector", "Memory")
)

table(cell_types)

#--------------------- Part 0: Pre-processing ---------------------
##### Distribution of raw data #####
# Histogram of all expression values across all genes
hist(as.matrix(dat),
  breaks = 50,
  main = "Distribution of Raw Gene Expression Values",
  xlab = "Expression Value",
  ylab = "Frequency"
)

##### Outliers Detection #####
# Standardize the raw gene expression data
dat_scaled <- scale(dat)

# Compute a summary outlier score for each sample
# Using the maximum absolute z-score across all genes
sample_outlier_score <- apply(abs(dat_scaled), 1, max)

# Identify samples with extreme expression values
# Threshold: |z| > 3
outlier_samples <- sample_outlier_score > 3

# Create results table
outlier_df <- data.frame(
  Sample = rownames(dat),
  Max_Zscore = sample_outlier_score,
  Outlier = outlier_samples
)

# Print outliers
print("Samples flagged as potential outliers (|z| > 3):")
print(outlier_df[outlier_df$Outlier == TRUE, ])

# Plot outliers
ggplot(outlier_df, aes(x = reorder(Sample, Max_Zscore), y = Max_Zscore, fill = Outlier)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Outlier Detection of Raw Data",
    x = "Sample",
    y = "Max |Z-score|"
  )


##### Multicollinearity Check #####
# With 156 genes and only 30 samples, traditional VIF is not appropriate
# (more variables than observations). Instead we assess multicollinearity
# via: (1) correlation heatmap, (2) eigenvalue analysis, (3) pairwise correlation distribution.

# --- Method 1: Correlation Matrix Heatmap ---
# Compute pairwise Pearson correlations between all 156 genes
# producing a 156x156 matrix where each cell is a value between -1 and +1
cor_matrix <- cor(dat, use = "pairwise.complete.obs")

# Visualize as heatmap — strong off-diagonal blocks indicate multicollinearity
heatmap(cor_matrix,
        symm = TRUE,
        col  = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Gene-Gene Correlation Heatmap",
        labRow = NA, labCol = NA)

legend("topright",
       legend = c("High Positive (+1)", "No Correlation (0)", "High Negative (-1)"),
       fill   = c("red", "white", "blue"),
       border = "grey50",
       bty    = "n",
       cex    = 0.8)

# --- Method 2: Eigenvalue Analysis of Correlation Matrix ---
# Decompose the correlation matrix into eigenvalues.
# Each eigenvalue represents the amount of variance in a particular direction
# of the gene expression space. Near-zero eigenvalues signal multicollinearity.
eigenvalues <- eigen(cor_matrix)$values

# Note: because we have more variables (156 genes) than observations (30 samples),
# the correlation matrix is singular (rank 29 at most). This means up to 127
# eigenvalues will be essentially zero or slightly negative due to floating point
# rounding — these are numerical noise, not real directions in the data.
# We therefore only use positive eigenvalues for the condition number calculation.
eigenvalues_pos <- eigenvalues[eigenvalues > 0]

n_small <- sum(eigenvalues_pos < 0.01)
print(paste0("Near-zero eigenvalues (< 0.01): ", n_small, " out of ", length(eigenvalues_pos)))

# Condition number = ratio of largest to smallest positive eigenvalue
# A value > 1000 is a strong indicator of multicollinearity
condition_number <- max(eigenvalues_pos) / min(eigenvalues_pos)
print(paste0("Condition number: ", round(condition_number, 2)))

# Note: a high condition number is expected in gene expression data due to
# co-regulated gene networks. PCA handles this by producing orthogonal components,
# eliminating multicollinearity by construction.


# --- Method 3: Distribution of pairwise correlations ---
# Extract upper triangle of correlation matrix to get all unique gene pairs
# For 156 genes this gives 156×155/2 = 12,090 unique pairs.
upper_tri <- cor_matrix[upper.tri(cor_matrix)]

# Plots all 12,090 correlation values as a histogram.
hist(upper_tri,
  breaks = 50,
  col    = "steelblue",
  border = "white",
  main   = "Distribution of Pairwise Gene Correlations",
  xlab   = "Pearson Correlation Coefficient",
  ylab   = "Frequency"
)
abline(v = c(-0.8, 0.8), col = "red", lty = 2, lwd = 1.5)
mtext("Red dashed lines = ±0.8 threshold", 
      side = 3,    # top of plot
      line = 0.3,  # just below the main title
      cex  = 0.8,  # smaller than main title
      col  = "grey40") # Draws vertical red dashed lines at ±0.8 as a conventional threshold for highly correlated. Anything beyond these lines is considered strongly collinear.

# Report proportion of highly correlated pairs
prop_high <- mean(abs(upper_tri) > 0.8) # Calculates what fraction of all gene pairs exceed the ±0.8 threshold
# Gives you a single interpretable number summarising the extent of multicollinearity.
print(paste0("Proportion of gene pairs with |r| > 0.8: ", round(prop_high * 100, 1), "%"))


# --------------------- Part 1: PCA ---------------------
# Scale = TRUE is appropriate for gene expression data because genes are measured on different scales; standardizing ensures no single gene dominates variance to prevent unequal weights impacting PCA.
# This matters because Different genes have genuinely different expression ranges.
# Ensures PCA reflects patterns of co-variation across genes, not just which genes happen to have large absolute values.
# It ensures all 156 genes contribute fairly to the principal components.

##### Sample PCA #####
pca_res <- prcomp(dat, scale. = TRUE)
summary(pca_res)

# Calculate variance explained by each PC
# Eigenvalues = squared standard deviations of each PC
eigenvalues_samples <- pca_res$sdev^2

# Proportion of total variance explained (as percentages)
var_explained <- eigenvalues_samples / sum(eigenvalues_samples) * 100

# Print variance explained for key PCs
print(paste0("PC1 variance explained: ", round(var_explained[1], 2), "%"))
print(paste0("PC2 variance explained: ", round(var_explained[2], 2), "%"))
print(paste0("PC3 variance explained: ", round(var_explained[3], 2), "%"))
print(paste0("Cumulative variance (PC1+PC2): ", round(sum(var_explained[1:2]), 2), "%"))


##### Gene PCA #####
# Transpose: genes become rows, samples become columns
pca_genes <- prcomp(t(dat), scale. = TRUE)

# Variance explained
eigenvalues_genes <- pca_genes$sdev^2
var_genes_explained <- eigenvalues_genes / sum(eigenvalues_genes) * 100

# Build plot data frame
gene_df <- as.data.frame(pca_genes$x)
gene_df$gene <- rownames(gene_df)

# Print key PCs
print(paste0("Gene PC1 variance explained: ", round(var_genes_explained[1], 2), "%"))
print(paste0("Gene PC2 variance explained: ", round(var_genes_explained[2], 2), "%"))
print(paste0("Gene PC3 variance explained: ", round(var_genes_explained[3], 2), "%"))
print(paste0("Cumulative variance (PC1+PC2): ", round(sum(var_genes_explained[1:2]), 2), "%"))

##### Scree Plot for Samples #####
# Build a data frame with 2 columns:
# PC — the PC number
# Eigenvalue — the raw eigenvalue for each PC (sdev^2)
# Capped at 20 PCs so the plot isn't overcrowded

scree_samples_df <- data.frame(
  PC           = 1:min(20, length(eigenvalues_samples)),
  Eigenvalue   = eigenvalues_samples[1:min(20, length(eigenvalues_samples))]
)

# Each point/line height = the eigenvalue for that PC
# Look for the elbow — PCs before it are the most important
ggplot(scree_samples_df, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 1:20) +
  theme_minimal(base_size = 12) +
  labs(
    title   = "Scree Plot for Samples",
    x       = "PC Number",
    y       = "Eigenvalue"
  )


##### Scree Plot for gene #####
# Scree plot data frame — capped at 30 PCs
# (gene PCA can only have min(156, 30)-1 = 29 meaningful PCs
# since we only have 30 samples as variables)
scree_genes_df <- data.frame(
  PC         = 1:min(30, length(eigenvalues_genes)),
  Eigenvalue = eigenvalues_genes[1:min(30, length(eigenvalues_genes))]
)

ggplot(scree_genes_df, aes(x = PC, y = Eigenvalue)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 3) +
  scale_x_continuous(breaks = 1:30) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Scree Plot for Genes",
    x        = "PC Number",
    y        = "Eigenvalue"
  )

#--------------- Part 2: Visualize and Investigate PCs --------------

# Build a data frame with PC scores and extracted metadata
# pca_res$x contains the PC scores — the coordinates of each sample in the new PC space. Converting to a data frame and saving row names (e.g. HEA26_EFFE_1) allows the extraction of metadata from them.
pca_df <- as.data.frame(pca_res$x)
pca_df$sample_names <- rownames(pca_df)

# Extract Status from row names (HEA = Healthy, MEL = Melanoma)
# grepl() searches each row name for the pattern "HEA" — if found, labels it "Healthy", otherwise "Melanoma".
pca_df$Status <- ifelse(grepl("HEA", pca_df$sample_names), "Healthy", "Melanoma")

# Extract Cell Type from row names
pca_df$CellType <- "Unknown"
pca_df$CellType[grepl("NAI", pca_df$sample_names)] <- "Naïve"
pca_df$CellType[grepl("EFFE", pca_df$sample_names)] <- "Effector"
pca_df$CellType[grepl("MEM", pca_df$sample_names)] <- "Memory"

# Order cell types biologically: Naïve first as the precursor cell type
pca_df$CellType <- factor(pca_df$CellType, levels = c("Naïve", "Effector", "Memory"))


# --- PCA Scatter Plot of Samples: PC1 vs PC2 ---
# Each point = one sample, colour = cell type, shape = disease status
ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType, shape = Status)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    name = "Cell Type",
    values = c("Naïve" = "#2ca02c", "Effector" = "#d62728", "Memory" = "#1f77b4")
  ) +
  scale_shape_manual(
    name = "Subject Status",
    values = c("Healthy" = 15, "Melanoma" = 17) # square = healthy, triangle = melanoma
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "Sample PCA: T-Cell Gene Expression (PC1 vs PC2)",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)")
  )


# --- PCA Scatter Plot of Genes: PC1 vs PC2 ---
# Each point = one gene
gene_df <- as.data.frame(pca_genes$x)

ggplot(gene_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, alpha = 0.6, color = "steelblue") +
  theme_minimal(base_size = 12) +
  labs(
    title   = "Gene PCA: Expression Patterns Across Samples (PC1 vs PC2)",
    x       = paste0("PC1 (", round(var_genes_explained[1], 1), "%)"),
    y       = paste0("PC2 (", round(var_genes_explained[2], 1), "%)")
  )

