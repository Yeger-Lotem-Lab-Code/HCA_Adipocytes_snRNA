
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)

input_path <- ".."
output_path <- "..."

obj_seurat <- readRDS(paste0(output_path,'adipo_vissc_diet.rds', sep=""))
gene_data <- read.csv('github/HCA_Adipocytes_snRNA/Figures/FIG.S11/mito_genes.csv', 
                      na.strings = c("", "NA"))  # Treat empty values as NA

# Remove rows where all columns are NA (i.e., completely empty rows)
gene_data <- gene_data[complete.cases(gene_data), ]
# Assuming the column name for genes is 'Gene'
# Exclude genes starting with "MT-"
genes_of_interest <- gene_data %>%
  filter(!startsWith(Symbol, "MT-")) %>%
  pull(Symbol)

# Check the first few entries of the filtered gene list
head(genes_of_interest)
genes_in_obj <- intersect(genes_of_interest, rownames(obj_seurat))

obj_seurat[["percent.nuc_mito"]] <- PercentageFeatureSet(obj_seurat, features = genes_in_obj)
obj_seurat@meta.data$adipo_ann<- factor(x = obj_seurat@meta.data$adipo_ann, levels = c("SA1", "SA2", "SA3", "SA4","SA5", "SA6", "SA7", "VA1", "VA2", "VA3", "VA4", "VA5", "VA6", "VA7", "VA8"))
seurat_metadata <- obj_seurat@meta.data

##First part check VAT adipocytes compared to SAT adipocytes:
# Calculate the median percent.nuc_mito for each cluster
median_percent_nuc_mito <- seurat_metadata %>%
  group_by(adipo_ann) %>%  # Replace seurat_clusters with your cluster column name
  summarize(median_percent = median(percent.nuc_mito, na.rm = TRUE))


# Split data into two groups: clusters starting with 'S' and 'V'
group_S <- subset(median_percent_nuc_mito, grepl("^S", adipo_ann))$median_percent
group_V <- subset(median_percent_nuc_mito, grepl("^V", adipo_ann))$median_percent

# Check if data is normally distributed using Shapiro-Wilk test
shapiro.test(group_S)
shapiro.test(group_V)

##Since the p-value (0.9211) is greater than the common significance level of 0.05, we fail to reject the null hypothesis. This suggests that group_S is normally distributed.
# Perform T-test (parametric) since both groups are normally distributed
t_test_result <- t.test(group_S, group_V)
print(t_test_result)

# If the data is not normally distributed, use the Wilcoxon rank-sum test (non-parametric)
# wilcox_test_result <- wilcox.test(group_S, group_V)
# print(wilcox_test_result)


##Second part compare VA6 adipocytes compared to all other VAT adipocytes:
# Filter for VA6 and other VA clusters
va6_data <- seurat_metadata %>%
  filter(adipo_ann == "VA6") %>%
  pull(percent.nuc_mito)

other_va_data <- seurat_metadata %>%
  filter(grepl("VA", adipo_ann) & adipo_ann != "VA6") %>%
  pull(percent.nuc_mito)

# Step 1: Check normality of VA6 and other VA clusters using Shapiro-Wilk test
# shapiro_va6 <- shapiro.test(va6_data)
# shapiro_other_va <- shapiro.test(other_va_data)

# # Step 2: Perform appropriate test based on normality results
# if (shapiro_va6$p.value > 0.05 & shapiro_other_va$p.value > 0.05) {
#   # If both groups are normally distributed, perform t-test
#   t_test_result <- t.test(va6_data, other_va_data)
#   test_type <- "t-test"
# } else {
#   # If either group is not normally distributed, perform Mann-Whitney U test
#   t_test_result <- wilcox.test(va6_data, other_va_data)
#   test_type <- "Mann-Whitney U test"
# }

# Q-Q plot for VA6 cluster
qqnorm(va6_data, main = "Q-Q Plot for VA6 Cluster")
qqline(va6_data)

# Q-Q plot for other VA clusters
qqnorm(other_va_data, main = "Q-Q Plot for Other VA Clusters")
qqline(other_va_data)

# Perform Mann-Whitney U test
mw_test_result <- wilcox.test(va6_data, other_va_data)

# Print the Mann-Whitney U test result
print(mw_test_result)

# Calculate medians for VA6 and other VA clusters
va6_median <- median(va6_data, na.rm = TRUE)
other_va_median <- median(other_va_data, na.rm = TRUE)

# Print medians for a comparison view
cat("VA6 Median percent.nuc_mito:", va6_median, "\n")
cat("Other VA Clusters Median percent.nuc_mito:", other_va_median, "\n")

# Determine and print which group has a higher median
if (va6_median > other_va_median) {
  cat("The median percent.nuc_mito of the VA6 cluster is higher than that of the other VA clusters.\n")
} else {
  cat("The median percent.nuc_mito of the VA6 cluster is not higher than that of the other VA clusters.\n")
}


## Third part create a ridge plot for all the nuclear  encoded genes:

p1 <-  ggplot(seurat_metadata, aes(x = percent.nuc_mito , y = adipo_ann,  fill=adipo_ann)) +  
  scale_y_discrete(limits = rev) +
  geom_density_ridges(scale = 0.9, alpha=0.3) +
  theme_ridges(center_axis_labels = TRUE) +
  scale_x_continuous(limits = c(0, 10))  # Set x-axis range from 0 to 5
#scale_fill_manual(values = colors10) #+
# facet_wrap(~orig.ident)
jpeg(paste0(input_path, "try_adipo_all_mt_nuc_ridge_levels.jpeg"), width = 3000, height = 3000, bg = "transparent",res=300)
print(p1)
dev.off()