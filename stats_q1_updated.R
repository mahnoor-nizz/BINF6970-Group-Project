# Stats 6970 - Problem 1
# PCA Analysis of T-Cell Gene Expression Data
# Reference: Holmes et al. (2005)
#
# Background:
# This analysis investigates how T cells differentiate after thymic selection.
# Two competing models exist:
#   - Linear model: NAI -> EFFE -> MEM
#   - Parallel model: NAI -> EFFE and MEM simultaneously
#
# Dataset: 156 most differentially expressed genes (Holmes et al., 2005)
# Row name format: [Status]_[CellType]_[ID]
#   Status:   HEA = Healthy, MEL = Melanoma
#   CellType: NAI = Naïve, EFFE = Effector, MEM = Memory
# Each column = a gene's expression measurement (a numeric value indicating how highly/lowly that gene is expressed in that sample)

##### Import Libraries ####
library(ggplot2)
library(ggfortify)

# Load the gene expression file and explore the data
load("geneexpression2.rda")
head(dat)
dim(dat)  # Rows = samples, Columns = 156 differentially expressed genes


##### Question 1: PCA #####

# Scale = TRUE is appropriate for gene expression data because genes are measured on different scales; standardizing ensures no single gene dominates variance to prevent unequal weights impacting PCA.
# - This matters because Different genes have genuinely different expression ranges.
# - Ensures PCA reflects patterns of co-variation across genes, not just which genes happen to have large absolute values.
# - It ensures all 156 genes contribute fairly to the principal components.
pca_res <- prcomp(dat, scale. = TRUE)
summary(pca_res)

## Calculate variance explained by each PC
# Eigenvalues = squared standard deviations of each PC
variance <- pca_res$sdev^2

# Proportion of total variance explained (as percentages)
var_explained <- variance / sum(variance) * 100

# Print variance explained for key PCs
print(paste0("PC1 variance explained: ", round(var_explained[1], 2), "%"))
print(paste0("PC2 variance explained: ", round(var_explained[2], 2), "%"))
print(paste0("Cumulative variance (PC1+PC2): ", round(sum(var_explained[1:2]), 2), "%"))

# --- Scree Plot: to justify how many PCs to retain ---
# This creates a table with 3 columns:
# PC — the PC number
# VarianceExplained — the % variance each individual PC explains
# Cumulative — the running total of variance explained, using cumsum(). So row 1 = PC1's %, row 2 = PC1+PC2's %.
# Note: capped at 20 PCs so the plot isn't overcrowded.
scree_df <- data.frame(
  PC = 1:min(20, length(var_explained)),
  VarianceExplained = var_explained[1:min(20, length(var_explained))],
  Cumulative = cumsum(var_explained)[1:min(20, length(var_explained))]
)

# Each bar height = the variance explained by that individual PC. 
# stat = "identity" tells ggplot to use the actual values in the data rather than counting anything.
# This overlays a second layer using the Cumulative column on the same y-axis as the bars. The line climbs from left to right, showing how quickly the PCs collectively explain the total variance — if it flattens early, the first few PCs are sufficient.
ggplot(scree_df, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_line(aes(y = Cumulative), color = "darkred", linewidth = 0.8) +
  geom_point(aes(y = Cumulative), color = "darkred", size = 2) +
  scale_x_continuous(breaks = 1:20) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Scree Plot: Variance Explained by Each Principal Component",
    x = "Principal Component",
    y = "Variance Explained (%)",
    caption = "Red line = cumulative variance explained"
  )


##### Question 2: Visualize and Investigate PCs #####

# Build a data frame with PC scores and extracted metadata
#pca_res$x contains the PC scores — the coordinates of each sample in the new PC space. Converting to a data frame and saving row names (e.g. HEA26_EFFE_1) allows the extraction of metadata from them.
pca_df <- as.data.frame(pca_res$x)
pca_df$sample_names <- rownames(pca_df)

# Extract Status from row names (HEA = Healthy, MEL/other = Melanoma)
#grepl() searches each row name for the pattern "HEA" — if found, labels it "Healthy", otherwise "Melanoma".
pca_df$Status <- ifelse(grepl("HEA", pca_df$sample_names), "Healthy", "Melanoma")

# Extract Cell Type from row names
pca_df$CellType <- "Unknown"
pca_df$CellType[grepl("NAI",  pca_df$sample_names)] <- "Naïve"
pca_df$CellType[grepl("EFFE", pca_df$sample_names)] <- "Effector"
pca_df$CellType[grepl("MEM",  pca_df$sample_names)] <- "Memory"

# Forces the legend and any downstream grouping to follow biological order rather than alphabetical — Naïve comes first since it's the precursor cell type.
pca_df$CellType <- factor(pca_df$CellType, levels = c("Naïve", "Effector", "Memory"))

# --- Primary PCA Scatter Plot: PC1 vs PC2 ---
#x/y axes → PC1 and PC2 scores (the two directions of most variance)
#color → cell type
#shape → disease status
ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType, shape = Status)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    name = "Cell Type",
    values = c("Naïve" = "#2ca02c", "Effector" = "#d62728", "Memory" = "#1f77b4")
  ) +
  scale_shape_manual(
    name = "Subject Status",
    values = c("Healthy" = 15, "Melanoma" = 17)  # square = healthy, triangle = melanoma
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "PCA of T-Cell Gene Expression (PC1 vs PC2)",
    subtitle = "Color = Cell Type | Shape = Subject Status",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance explained)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance explained)"),
    caption = "Data: 156 differentially expressed genes from Holmes et al. (2005)"
  )

# --- Secondary Plot: PC1 vs PC3 (to check additional structure) ---
ggplot(pca_df, aes(x = PC1, y = PC3, color = CellType, shape = Status)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    name = "Cell Type",
    values = c("Naïve" = "#2ca02c", "Effector" = "#d62728", "Memory" = "#1f77b4")
  ) +
  scale_shape_manual(
    name = "Subject Status",
    values = c("Healthy" = 15, "Melanoma" = 17)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCA of T-Cell Gene Expression (PC1 vs PC3)",
    subtitle = "Color = Cell Type | Shape = Subject Status",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC3 (", round(var_explained[3], 1), "%)"),
    caption = "Data: 156 differentially expressed genes from Holmes et al. (2005)"
  )


##### Multicollinearity Check #####
# With 156 gene expression variables, traditional VIF is not appropriate
# (more variables than observations). Instead, we assess multicollinearity
# via: (1) correlation heatmap, (2) eigenvalue analysis of the correlation matrix.

# --- Method 1: Correlation Matrix Heatmap (subset for readability) ---
# Compute pairwise correlations between all 156 genes. This computes the Pearson correlation between every possible pair of genes — producing a 156×156 matrix where each cell contains a value between -1 and +1.
cor_matrix <- cor(dat, use = "pairwise.complete.obs")

# Visualize as heatmap — strong off-diagonal blocks indicate multicollinearity
# Using base R heatmap to avoid extra dependencies
#symm = TRUE — tells R the matrix is symmetric (correlation of A→B equals B→A), so it's displayed correctly.
# colorRampPalette(c("blue", "white", "red"))(100) — creates a 100-colour gradient where:
#  Red = strong positive correlation (+1)
# White = no correlation (0)
#  Blue = strong negative correlation (-1)
heatmap(cor_matrix,
        symm = TRUE,
        col  = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Gene-Gene Correlation Heatmap\n(Red = high positive, Blue = high negative correlation)",
        labRow = NA, labCol = NA)  # suppress labels — 156 genes too dense to label

# --- Method 2: Proportion of near-zero eigenvalues ---
# Eigenvalues of the correlation matrix: near-zero values signal multicollinearity
eigenvalues <- eigen(cor_matrix)$values

# Count eigenvalues below threshold (< 0.01 indicates near-singular directions)
# When the 156×156 correlation matrix is decomposed, you get 156 eigenvalues. Each one represents how much variance exists in a particular direction of the gene expression space.
n_small <- sum(eigenvalues < 0.01)
print(paste0("Eigenvalues < 0.01: ", n_small, " out of ", length(eigenvalues)))

# Condition number = ratio of largest to smallest eigenvalue
# > 1000 is a strong indicator of multicollinearity
condition_number <- max(eigenvalues) / min(eigenvalues)
print(paste0("Condition number: ", round(condition_number, 2)))
# Note: High condition number is EXPECTED in gene expression data due to co-regulated gene networks. PCA handles this by working in PC space, which is orthogonal (zero multicollinearity by construction).

# --- Method 3: Distribution of pairwise correlations ---
# Extract upper triangle of correlation matrix (avoid duplicates/diagonal)
#  For 156 genes this gives 156×155/2 = 12,090 unique pairs.
upper_tri <- cor_matrix[upper.tri(cor_matrix)]

#Plots all 12,090 correlation values as a histogram.
hist(upper_tri,
     breaks = 50,
     col    = "steelblue",
     border = "white",
     main   = "Distribution of Pairwise Gene Correlations",
     xlab   = "Pearson Correlation Coefficient",
     ylab   = "Frequency")
abline(v = c(-0.8, 0.8), col = "red", lty = 2, lwd = 1.5)
legend("topright", legend = "±0.8 threshold", col = "red", lty = 2) #Draws vertical red dashed lines at ±0.8 as a conventional threshold for "highly correlated." Anything beyond these lines is considered strongly collinear.

# Report proportion of highly correlated pairs
prop_high <- mean(abs(upper_tri) > 0.8) #Calculates what fraction of all gene pairs exceed the ±0.8 threshold — gives you a single interpretable number summarising the extent of multicollinearity.
print(paste0("Proportion of gene pairs with |r| > 0.8: ", round(prop_high * 100, 1), "%"))


##### Outlier Detection #####
# We assess outliers in two complementary ways:
# (1) In PCA space — samples with extreme PC scores
# (2) Mahalanobis distance — multivariate distance from group centroid

# --- Method 1: PC Score Outliers (boxplots on PC1 and PC2) ---
par(mfrow = c(1, 2)) #Splits the plot window into 1 row and 2 columns so both boxplots appear side by side.

boxplot(pca_df$PC1,
        main  = "PC1 Scores",
        ylab  = "Score",
        col   = "lightblue",
        outline = TRUE) #makes the boxplot show individual points beyond the whiskers as circles — these are the potential outliers.
stripchart(pca_df$PC1, vertical = TRUE, method = "jitter",
           pch = 16, col = "darkblue", add = TRUE)

boxplot(pca_df$PC2,
        main  = "PC2 Scores",
        ylab  = "Score",
        col   = "lightcoral",
        outline = TRUE)
stripchart(pca_df$PC2, vertical = TRUE, method = "jitter",
           pch = 16, col = "darkred", add = TRUE) #overlays all individual data points with slight horizontal jitter so you can see every sample, not just the flagged ones.

par(mfrow = c(1, 1))

# Flag samples whose PC1 or PC2 score exceeds 3 SD from the mean
# Computes the mean and standard deviation for each PC to define the outlier threshold.
pc1_mean <- mean(pca_df$PC1); pc1_sd <- sd(pca_df$PC1)
pc2_mean <- mean(pca_df$PC2); pc2_sd <- sd(pca_df$PC2)

#is this sample more than 3 standard deviations away from the centre? The 3 SD rule comes from the normal distribution where 99.7% of data falls within 3 SD — anything beyond is unusual.
outliers_pc <- pca_df[
  abs(pca_df$PC1 - pc1_mean) > 3 * pc1_sd |
  abs(pca_df$PC2 - pc2_mean) > 3 * pc2_sd, 
  c("sample_names", "CellType", "Status", "PC1", "PC2")
]
print("Samples flagged as outliers (>3 SD on PC1 or PC2):")
print(outliers_pc)

# --- Method 2: Mahalanobis Distance (per cell type group) ---
# Mahalanobis distance accounts for correlation structure within each group, stretching and squishing distances based on how the data is actually spread.
mah_dist <- numeric(nrow(pca_df))

# Use first 3 PCs (captures majority of variance, avoids noise PCs)
pc_cols <- paste0("PC", 1:3)

#Loops through each cell type separately — this is important to find outliers within each group, not samples that are just far from the overall centre because they belong to a different cell type.
for (ct in unique(pca_df$CellType)) {
  idx   <- which(pca_df$CellType == ct)
  group <- as.matrix(pca_df[idx, pc_cols])
  
  # Need at least p+1 observations to invert covariance matrix
  if (nrow(group) > ncol(group)) {
    cov_mat      <- cov(group) #the covariance matrix describing how the 3 PCs relate to each other within the group
    center       <- colMeans(group) #the centroid (centre point) of the group
    mah_dist[idx] <- mahalanobis(group, center = center, cov = cov_mat)
  } # computes how far each sample is from that centroid, adjusted for the covariance structure
}

pca_df$MahalanobisD <- mah_dist

# Chi-squared critical value: df = 3 PCs, alpha = 0.001 (conservative). Any sample exceeding this threshold is flagged 
chi_crit <- qchisq(0.999, df = 3)
pca_df$Outlier <- pca_df$MahalanobisD > chi_crit

# Plot Mahalanobis distances with threshold
ggplot(pca_df, aes(x = seq_len(nrow(pca_df)), y = MahalanobisD,
                   color = CellType, shape = Outlier)) +
  geom_point(size = 3) +
  geom_hline(yintercept = chi_crit, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("Naïve" = "#2ca02c", "Effector" = "#d62728", "Memory" = "#1f77b4")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8), labels = c("Normal", "Outlier")) + #normal samples get a filled circle (16) and outliers get a star/asterisk shape (8)
  theme_minimal(base_size = 12) +
  labs(
    title    = "Mahalanobis Distance per Sample (within Cell Type groups)",
    subtitle = paste0("Red dashed line = χ²(df=5, p=0.001) threshold = ", round(chi_crit, 2)),
    x        = "Sample Index",
    y        = "Mahalanobis Distance",
    shape    = "Outlier Status"
  )

# Print flagged outliers
outliers_mah <- pca_df[pca_df$Outlier == TRUE,
                        c("sample_names", "CellType", "Status", "PC1", "PC2", "MahalanobisD")]
print("Samples flagged as multivariate outliers (Mahalanobis, p < 0.001):")
print(outliers_mah)


##### Question 3: Interpretation #####
#
# The PCA plot reveals clear clustering by cell type along PC1 and PC2, with the three T-cell populations (Naïve, Effector, Memory) forming distinct, well-separated groups.
#
# PC1 (63.8%) — Differentiation Axis
# PC1 explains the overwhelming majority of variance and clearly separates cell types. Naïve cells sit on the right (positive PC1), while Effector cells cluster on the far left (negative PC1), with Memory cells sitting between them but closer to Effector. This axis represents the primary transcriptional shift that occurs when a Naïve T cell encounters an antigen and becomes activated.

# PC2 (9.6%) — Effector vs Memory Axis
# PC2 separates Effector and Memory cells from each other vertically, despite both being on the same side of PC1. Memory cells sit higher (positive PC2) and Effector cells lower (negative PC2). This suggests that while both populations share a common activation signature relative to Naïve cells, they maintain distinct transcriptional programs.

# PC3 (4%) — Minor Structure
# PC1 vs PC3 reveals no additional meaningful separation beyond what PC2 already showed. The cell type clusters remain intact but PC3 adds little biological insight, confirming that the first two PCs capture the dominant structure in the data.
#
#Differentiation Model
# Critically, both Effector and Memory cells are separated from Naïve on PC1 simultaneously — Memory cells do not sit intermediate between Naïve and Effector as the linear model (NAI → EFFE → MEM) would predict. Instead, both populations diverge from Naïve in parallel, supporting the parallel differentiation model proposed by Holmes et al. (2005).

# Healthy vs Melanoma
# Within each cell type cluster, Healthy (squares) and Melanoma (triangles) samples overlap considerably, indicating that disease status has minimal influence on the overall transcriptional landscape compared to cell type identity.
# 
# Multicollinearity
# The correlation heatmap shows clear block structure, confirming that genes do not act independently but rather in co-regulated clusters, which is expected given the coordinated nature of T-cell gene networks. The pairwise correlation histogram reveals a bimodal distribution, suggesting two distinct gene co-expression modules — one group of genes that are strongly positively co-expressed with each other (right peak, ~0.6–0.8), and another group that behaves more independently (left peak, near 0). This bimodal structure is consistent with the two dominant axes of variation captured by PCA: the strongly correlated gene module likely drives the Naïve vs. activated cell separation along PC1, while the more independent gene module contributes to the Effector vs. Memory separation along PC2. Importantly, this multicollinearity does not affect the validity of the PCA results — the principal components are orthogonal by construction, fully eliminating redundancy in the reduced space.
# 
# Outliers NEED TO FINISH