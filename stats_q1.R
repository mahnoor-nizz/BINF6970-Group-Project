# Stats 6970
# Question 1

# Overview
# EFFE = effector
# MEM = memory
# Do Native T cells differentiate into EFFE first before becoming MEM?
# Or do they differentiate into EFFE and MEM simultaneously?

##### Import Libraries ####
library(ggplot2)
library(ggfortify)

# Load the gene expression file and explore the data
load("geneexpression2.rda")
head(dat)
dim(dat)



##### PCA #####
# Perform PCA (scaling is recommended for gene expression)
pca_res = prcomp(dat, scale. = TRUE)
summary(pca_res)


## Calculate variance explained
# Calculate eigen values
variance = pca_res$sdev^2

# Calculate proportion of variance explained
var_explained <- variance / sum(variance) * 100

# Print variance
print(paste0("PC1 variance: ", round(var_explained[1], 2), "%"))
print(paste0("PC2 variance: ", round(var_explained[2], 2), "%"))

# Create a dataframe for plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$sample_names <- rownames(pca_df)

# Extract metadata using regex
pca_df$Status <- ifelse(grepl("HEA", pca_df$sample_names), "Healthy", "Melanoma")

# Extract Cell Type
pca_df$CellType <- "Unknown"
pca_df$CellType[grepl("NAI", pca_df$sample_names)] <- "Naive"
pca_df$CellType[grepl("EFFE", pca_df$sample_names)] <- "Effector"
pca_df$CellType[grepl("MEM", pca_df$sample_names)] <- "Memory"



##### Visualization #####

# Plotting
ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType, shape = Status)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c("Healthy" = 15, "Melanoma" = 17)) +
  theme_minimal() +
  labs(title = "PCA of T-cell Gene Expression",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  scale_color_manual(values = c("Effector" = "red", "Memory" = "blue", "Naive" = "green"))

