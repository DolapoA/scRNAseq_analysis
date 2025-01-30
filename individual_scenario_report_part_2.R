# Install required packages
#BiocManager::install("limma")
#install.packages('ggrepel')
#install.packages("WebGestaltR")

library(dplyr)
library(ggplot2)
library(tidyr)
library(limma)
library(ggrepel)
library(WebGestaltR)


# Read in data
exp_data <- read.csv('/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/X.csv', header=FALSE)
metadata <-read.csv('/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/obs.csv')
var_data<- read.csv('/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/var.csv')

head(exp_data)
head(metadata)
head(var_data)

# Make the row names of exp_data be CellIDs, and the column names be gene names
# CellIDs are being used for rownames to serve as a way of referring to specific genes
rownames(exp_data)<-metadata$CellID
colnames(exp_data)<-var_data$Gene

dim(metadata)
dim(var_data)
dim(exp_data)
head(exp_data)

# Calculate the total number of cells per patient]# by grouping by 
# patient and then counting the cells associated with it

cells_per_patient <- metadata %>%
  group_by(Patient) %>%
  tally(name="Total_Cells")

# Use ggplot to plot this as a bar graph
ggplot(cells_per_patient, aes(x = reorder(Patient, -Total_Cells), y = Total_Cells)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Total Number of Cells per Patient", x = "Patient ID", y = "Total Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Calculate and plot total cell counts per disease stage
cells_per_stage <- metadata %>%
  group_by(Disease_stage) %>%
  tally(name="Total_Cells")
head(cells_per_stage)

ggplot(cells_per_stage, aes(x = Disease_stage, y = Total_Cells, fill = Disease_stage)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Number of Cells per Disease Stage", x = "Disease Stage", y = "Total Number of Cells") +
  theme_minimal()

# Calculate and plot counts of cell types per disease stage
cell_counts <- metadata %>%
  group_by(Celltype, Disease_stage) %>%
  tally(name="Total_Cells")
head(cell_counts)

# Plotting 3 columns of info
ggplot(cell_counts, aes(x = Disease_stage, y = Total_Cells, fill = Celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Type Distribution by Disease Stage", x = "Disease Stage", y = "Cell Count") +
  theme_minimal()

# Perform a statistical test to show if the differences in cell type proportions
# between disease stage are significant

filtered_data <- cell_counts[cell_counts$Disease_stage %in% c("IIIc", "Benign"),]

ttest_result <- t.test(Total_Cells ~ Disease_stage, data = filtered_data)
ttest_result


################################
# Differential gene expression #
################################


# Subset the metadata and expression data to the two secretory epithelial cell types

#secretorymeta<- subset(metadata, metadata$Celltype %in%c('Secretory Epithelial-1','Secretory Epithelial-2'))
secretorymeta<- subset(metadata, metadata$Celltype %in%c('High Grade Carcinoma','STIC lesion'))
secretoryexp<- subset(exp_data, rownames(exp_data) %in% secretorymeta$CellID)
unique(metadata$Celltype)
head(secretorymeta$CellID)
head(exp_data)

# Use dim to get row and column number
dim(secretorymeta)
dim(secretoryexp)

# Sanity check on IDs
secretoryexp <- t(secretoryexp)
identical(colnames(secretoryexp), secretorymeta$CellID)

# Perform differential gene expression test #

# This creates a design matrix in preparation for differential expression
# The relationship to be modelled is "the effect of Celltype on gene expression"
# The tilde separates the response (gene expression) from the predictor (Celltype)
# The end result is a table with "CelltypeSecretory Epithelial-2" or "treatment" with 0s and 1s
# The 0s mean the sample (row) doesn't belong to that cell type and the opposite
# is true for the 1s
design <- model.matrix(~ Celltype, data = secretorymeta)

# The lmFit (limma) function tries to find a mathematical relationship between the 
# predictor and the response, in other words, the gene expression and the experimental design (celltype)
fit <- lmFit(secretoryexp, design)

# A beneficial step for dealing with small sample sizes or high variability
fit <- eBayes(fit, trend = TRUE)
head(design)

# Print top 20 differentially expressed genes
top20genes <- topTable(fit, coef = 2, number = 20)
print(top20genes)

# Plot a volcano plot of top 100 differentially expressed genes
top_genes_df <- topTable(fit, coef = 2, number = 100)
top_genes_df$significance <- ifelse(top_genes_df$P.Value < 0.05 & abs(top_genes_df$logFC) > 1,
                                    "Significant", "Not Significant")
head(top_genes_df[4:7])

# Add a column for labeling genes with -log10(P.Value) > 150
top_genes_df$label <- ifelse(-log10(top_genes_df$P.Value) > 125, rownames(top_genes_df), "")


pdf(file = "~/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/r_plots/volcano_plot_HGC_vs_STIC.pdf",
    height = 12,
    width = 12)
ggplot(top_genes_df, aes(x = logFC, y = -log10(P.Value), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 10) +  # Add labels
  scale_color_manual(values = c("Not Significant" = "lightgrey", "Significant" = "dodgerblue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression - High Grade Carcinoma vs STIC Lesion",
       x = "Log2 Fold Change", y = "-Log10(p-value)") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top")

dev.off()

################################
# Gene Set Enrichment Analysis #
################################

# Import the list of ranked genes from cell types we generated using Scanpy yesterday

rankedgenes <- read.csv("/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/celltypemarkergenes-2.csv")
head(rankedgenes)

# We will look at High Grade Carcinoma, STIC lesion and Immune cell clusters
HGC <- rankedgenes[, colnames(rankedgenes) %in% c('High.Grade.Carcinoma_n', 'High.Grade.Carcinoma_s'), drop = FALSE]
STIC <- rankedgenes[, colnames(rankedgenes) %in% c('STIC.lesion_n', 'STIC.lesion_s'), drop = FALSE]
I2 <- rankedgenes[, colnames(rankedgenes) %in% c('Immune_n', 'Immune_s'), drop = FALSE]


head(HGC)
head(STIC)
head(I2)


# Save as rank files for input into Webgestalt
write.table(HGC, file = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/HGCrankedgenes.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(STIC, file = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/STICrankedgenes.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(I2, file = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/I2rankedgenes.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run Webgestalt GSEA for each rank file
HGCResult <- WebGestaltR(
  enrichMethod = "GSEA", organism = "hsapiens",
  enrichDatabase = "pathway_KEGG", interestGeneFile = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/HGCrankedgenes.rnk",
  interestGeneType = "genesymbol", sigMethod = "top", topThr = 10, minNum = 5)

STICResult <- WebGestaltR(
  enrichMethod = "GSEA", organism = "hsapiens",
  enrichDatabase = "pathway_KEGG", interestGeneFile = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/STICrankedgenes.rnk",
  interestGeneType = "genesymbol", sigMethod = "top", topThr = 10, minNum = 5)

I2Result <- WebGestaltR(
  enrichMethod = "GSEA", organism = "hsapiens",
  enrichDatabase = "pathway_KEGG", interestGeneFile = "/Users/dajayi/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/raw_data/I2rankedgenes.rnk",
  interestGeneType = "genesymbol", sigMethod = "top", topThr = 10, minNum = 5)

# To plot a dotplot of all the results combined
HGCResult$Celltype = 'High Grade Carcinoma'
STICResult$Celltype = 'STIC Lesion'
I2Result$Celltype = 'Immune'
gsea<- HGCResult %>% bind_rows(STICResult) %>% bind_rows(I2Result)
head(gsea)

# Plot the dot plot
pdf(file = "~/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/r_plots/GSEA_on_HGC_STIC_Immune_cells.pdf",
    height = 12,
    width = 12)
ggplot(gsea, aes(x = Celltype, y = description)) +
  geom_point(aes(size= FDR, color = pValue)) +
  scale_color_gradient(low = "darkgoldenrod", high = "dodgerblue", name = "pValue") +
  scale_size_continuous(range = c(0.2, 15)) +
  theme_minimal() +
  labs(
    title = "Gene Set Enrichment Analysis for High Grade Carcinoma, STIC lesions and Immune cells",
    x = "Cell Types",
    y = "Pathway Description"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.8))

dev.off()

#####################
# Linear Regression #
#####################


# Subset the data to High Grade Carcinoma and STIC Lesion again
secretorymeta<- subset(metadata, metadata$Celltype %in%c('High Grade Carcinoma','STIC lesion'))
secretoryexp<- subset(exp_data, rownames(exp_data) %in% secretorymeta$CellID)

# Check the highest expressed genes
head(HGC)
head(STIC)
head(I2)

# Create a new column in the metadata with HMGA1, HSPB8 and SRGN expression for each cell
secretorymeta$HMGA1<-exp_data[row.names(secretoryexp),"HMGA1"]
secretorymeta$HSPB8<-exp_data[row.names(secretoryexp),"HSPB8"]
secretorymeta$SRGN<-exp_data[row.names(secretoryexp),"SRGN"]

# Subset the secretory metadata to contain only cells with HMGA1, HSPB8 and SRGN expression equal to or greater than 0
filtered_meta <- secretorymeta %>%
  filter(HMGA1 >= 0 & HSPB8 >= 0 & SRGN >= 0)

# Run the linear regression model
model <- lm(HMGA1 ~ HSPB8, data = filtered_meta)
summary(model)


# Create a scatter plot of HMGA1 expression vs IER2 expression

geneexpression <- ggplot(filtered_meta, aes(x = HMGA1, y = HSPB8)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(title = "HMGA1 Expression vs HSPB8 Expression",
       x = "HMGA1 Expression",
       y = "HSPB8 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

geneexpression

# Add the Rsquared value onto the plot
r_squared <- summary(model)$r.squared

geneexpression<- geneexpression + annotate("text", x = 2, y = 2.5, label = paste("R² = ", round(r_squared,3)), size=6)
geneexpression

# Colour the points by disease stage or tissue type
stagegeneexpression <- geneexpression + geom_point(aes(color=Disease_stage))
stagegeneexpression

tissuegeneexpression<- geneexpression + geom_point(aes(color=Tissue))
tissuegeneexpression


#########################################################################
# Correlate several genes against each other using a correlation matrix #
#########################################################################

# Remove values with a standard deviation of 0, as this will break the correlation matrix
secretoryexp <- secretoryexp[,apply(secretoryexp,2, sd) > 0]

# Correlate genes against each other
cordata <- cor(secretoryexp [colnames(secretoryexp)],
             secretoryexp[colnames(secretoryexp)])


# Plot a correlation matrix consisting of genes with a correlation of >6
upper_tri <- upper.tri(cordata)
high_corr<-which(upper_tri & abs(cordata)>0.90, arr.ind=TRUE)

corrplot <- cordata[rownames(cordata) %in% rownames(high_corr),]
corrplot <- corrplot[,colnames(corrplot) %in% rownames(high_corr)]

library(pheatmap)

pdf(file = "~/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/r_plots/correlation_matrix.pdf",
     height = 16,
     width = 16)
pheatmap(corrplot, dendrogram = "none", treeheight_row = 0, treeheight_col = 0,
         show_rownames = TRUE, show_colnames = TRUE,
         color = hcl.colors(50, "BluYl"))
dev.off()

####################################
# Correlation matrix per cell type #
####################################

cell_type_list <- unique(metadata$Celltype)
cell_type_name_list <- c("SE1", "SE2", "CE", "F", "I", "STIC", "HGC")

corr_mat_list <- vector("list", 7) 
i <- 0

for (cell_type in cell_type_list) {
  i <- i+1
  
  print(paste("Preparing correlation matrix for", cell_type))
  
  # Subset the data to High Grade Carcinoma
  secretorymeta<- subset(metadata, metadata$Celltype %in% c(cell_type))
  secretoryexp<- subset(exp_data, rownames(exp_data) %in% secretorymeta$CellID)
  
  print(paste("secretorymeta dimensions for", cell_type))
  print(paste(dim(secretorymeta)))
  print(paste("secretoryexp dimensions for", cell_type))
  print(paste(dim(secretoryexp)))
  
  # Remove values with a standard deviation of 0, as this will break the correlation matrix
  secretoryexp <- secretoryexp[,apply(secretoryexp,2, sd) > 0]
  
  # Correlate genes against each other
  cordata <- cor(secretoryexp [colnames(secretoryexp)],
                 secretoryexp[colnames(secretoryexp)])
  
  
  # Plot a correlation matrix consisting of genes with a correlation of >6
  upper_tri <- upper.tri(cordata)
  high_corr<-which(upper_tri & abs(cordata)>0.90, arr.ind=TRUE)
  
  corrplot <- cordata[rownames(cordata) %in% rownames(high_corr),]
  corrplot <- corrplot[,colnames(corrplot) %in% rownames(high_corr)]
  
  corr_mat <- pheatmap(corrplot, dendrogram = "none", treeheight_row = 0, treeheight_col = 0,
                       show_rownames = TRUE, show_colnames = TRUE,
                       color = hcl.colors(50, "BluYl"))
  
  corr_mat_list[[i]] <- corr_mat
}

names(corr_mat_list) <- cell_type_name_list

library(pheatmap)
library(gridExtra)

#grid.arrange(corr_mat_list$SE1$gtable, corr_mat_list$SE2$gtable, corr_mat_list$CE$gtable, corr_mat_list$F$gtable)


pdf(file = "~/Documents/STP/University_of_Manchester/Modules/BIOL72241 - Year 3/Applied Statistics and Data Science/Assignment/r_plots/correlation_matrix_per_celltype.pdf",
    height = 8,
    width = 20)
grid.arrange(corr_mat_list$SE1$gtable, corr_mat_list$SE2$gtable, ncol = 2)
dev.off()

