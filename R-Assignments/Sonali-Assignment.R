
# MSB1005_Assignment_2021.R                                                   #
# 										                                      								  #
# Version: 1.0   															                                #
# Date: Dec 15, 2021											                                    #
# Author: Sonali Mittal; Msc. Student at Systems Biology                      #
# University: Maastricht University                                           #
# Supervisor: Michiel Adriaens, PhD; MaCSBio
# History:																	                                  #
#  1.0: Creation															                                #

# This R script includes the assignment for experimental data and design course#
# We follow the data priniciples of Findable, Accessible, Interoperable,
# Reusable(FAIR) science.


# Resources for the data script.


# This dataset has Resources of left ventricular transcriptomes of participants.
# It is part of Myocardial Applied Genomics Network (MAGnet).
# This network is a data bank for human cardiac tissue data to be used for
# genetic research.
# It can reached at following url,
# <a href="www.med.upenn.edu/magnet">A link</a>
# Additionally, the raw data set of this specific file can be found at:
# Gene Expression Omnibus, Series GSE141910,url for same is:
#  <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910">A link</a>
# We use the gene id numbers as given by ensemble data set for homo sapiens.
# They can be accessed at url:
# <a href="https://www.ensembl.org/index.html">A link</a>


# Files needed for the analysis
#-----------------------------------------------------------------------------#

# This analysis uses 3 files.
# The first text file we need is MAGNET_GeneExpressionData_CPM_19112020.txt
# This file contains the gene id as the row names and their level of expression
# for each sampled participant as column names.
# The next text file we need is MAGNET_SampleData_19112020.txt.
# This contains ID numbers of participants as represented by SampleID.
# Additionally, it contains the age, sex, disease condition and ethnicity
# of every participants.
# The third file we need is MAGNET_exonLengths.txt
# This file contains the ensemble gene id and length of exons.






#-----------------------------------------------------------------------------#
# The pre-processing done on raw data.
#-----------------------------------------------------------------------------#
# The data in this analysis has been preprocessed and log transformed
# to facilitate comparison between samples.


#-----------------------------------------------------------------------------#
# Information about the data script.
#-----------------------------------------------------------------------------#
# In the data structure abbreviations are used, as follows:
# Donor stands for non-failing healthy donors.
# DCM stands for dilated cardiomyopathy.
# HCM stands for hypertrophic cardiomyopathy.
# PPCm stands for Peripartum cardiomyopathie.


#-----------------------------------------------------------------------------#
# Question 1: Importing the data and inspecting sample information
#-----------------------------------------------------------------------------#

# downloading and installing required libraries of tidyverse, ggplot, tidyr

#install.packages("tidyr")
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)

#to install package for formatting code.
# install.packages("styler")

# To load package for writing tab delimited files.
library(readr)

# 1.a. A new working directory is set. Then data for gene expression, sample data
# and gene expression data files are imported as separate object.

# code for windows
#gxData <-
#  read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1)
# gxData <- cbind(EnsembleGeneID = rownames(gxData), gxData)
# rownames(gxData) <- 1:nrow(gxData)
#sampleInfo <-
#  read.delim("MAGNET_SampleData_19112020.txt", row.names = 1)
#gxDataPCA <-
#  read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt")


# 1.a. This is an alternate version for importing files.

# code for macbook.

setwd("~/Desktop/MAGNET_GX")
gxData <-
  read.delim("~/Desktop/MAGNET_GX/MAGNET_GeneExpressionData_CPM_19112020.txt",
             row.names = 1)
sampleInfo <-
  read.delim("~/Desktop/MAGNET_GX/MAGNET_SampleData_19112020.txt",
             row.names = 1)
gxDataPCA <-
  read.delim("~/Desktop/MAGNET_GX/MAGNET_GeneExpressionData_CPM_19112020.txt")

# 1.b. Exporting a summary table and/or figure(s) of participant characteristics.

# figures of participants characteristics
##  bar plot for gender distribution
png("gender_distribution.png")
## generating bar plot for gender.
sex_distribution <- with(sampleInfo, table(Sex))
barplot(sex_distribution,
        main = "Gender distribution of participants",
        xlab = "Gender",
        ylab = "Counts")

dev.off()

## table for sex distribution.
table(sampleInfo$Sex)
## saving result as dataframe.
summary_sex = as.data.frame(table(sampleInfo$Sex))
## Exporting result in a file.
write_tsv(summary_sex, file = "~/Desktop/MAGNET_GX/gender.txt")

## table showing distribution of sex in each disease group.
with(sampleInfo, table(Disease, Sex, Ethnicity))
## saving result as dataframe.
summary_distribution = as.data.frame(with(sampleInfo, table(Disease, Sex, Ethnicity)))
# Exporting result in a file.
write_tsv(summary_distribution, file = "~/Desktop/MAGNET_GX/diseaseDistribution.txt")

## exporting bar plot for ethnicity distribution
png("ethnicity_distribution.png")
## creating bar plot ethnicity distribution
ethnicity_distribution <- with(sampleInfo, table(Ethnicity))
barplot(
  ethnicity_distribution,
  main = "Ethnicity distribution of participants",
  xlab = "Ethnicity",
  ylab = "Counts"
)

dev.off()

## table for ethnicity distribution
table(sampleInfo$Ethnicity)
## saving table as dataframe
summary_ethnicity = as.data.frame(table(sampleInfo$Ethnicity))
## exporting dataframe as a file.
write_tsv(summary_ethnicity, file = "~/Desktop/MAGNET_GX/ethnicity.txt")




## bar plot for disease conditions distribution
png("disease_distribution.png")

disease_distribution <- with(sampleInfo, table(Disease))
barplot(disease_distribution,
        main = "Disease distribution of participants",
        xlab = "Disease",
        ylab = "Counts")

dev.off()

## table for disease conditions
table(sampleInfo$Disease)
## saving table as dataframe
summary_disease = as.data.frame(table(sampleInfo$Disease))
## exporting dataframe as a file.
write_tsv(summary_disease, file = "~/Desktop/MAGNET_GX/diseaseConditions.txt")



## barplot for age distribution.
png("age_distribution.png")

barplot(sampleInfo$Age,
        main = "Age distribution of participants",
        xlab = "Age",
        ylab = "Counts")

dev.off()

## table for age distribution.
summary(sampleInfo$Age)
## saving table as a dataframe.
summary_age = as.data.frame(summary(sampleInfo$Age))
## exporting dataframe as a file.
write_tsv(summary_disease, file = "~/Desktop/MAGNET_GX/ageDistribution.txt")


# CPM expression values can not be directly used to compare across genes.
# So, we need to convert the CPM expression values to FPKM expression values.
# FPKM corrects for gene length.
#
# We load the file which contains information about exon length.
geneTotExonLengths <-
  read.delim("MAGNET_exonLengths.txt",
             as.is = T,
             row.names = 1)
# We make a validation check to make sure that all genes in main data frame are present in the new file.
# Proceed forward if result is 'TRUE'.
all(rownames(geneTotExonLengths) == rownames(gxData))

# We define a function which will convert cpm into fpkm, when called.
cpm2fpkm <- function(x) {
  .t <- 2 ^ (x) * 1E3 / geneTotExonLengths[, 1]
}
# We use the function and store the converted values into a appended version of gxdata dataframe. 
#We call this appended version gxData_fpkm.
gxData_fpkm <- cpm2fpkm(gxData)
# To export the dataframe.
write_tsv(gxData_fpkm, file = "~/Desktop/MAGNET_GX/gxData_fpkm.txt")
#-----------------------------------------------------------------------------#
# Question 2: Diagnostic Plots
#-----------------------------------------------------------------------------#
# 2.a. Boxplots of gene expression Values
# Finding 5 most expressed genes
# calculate means of each column to get the average CPM count of gene expression
# and store as a data frame.
R_means <- data.frame(avg_cpm = rowMeans(gxData))
# convert rownames to a column called "Gene_exp"
R_means <- R_means %>% rownames_to_column(var = "Gene_exp")
# sort the gene_exp in descending order based on avg CPM values
R_means <- R_means[order(-R_means$avg_cpm),]
# To export the dataframe showing average CPM count of gene expression.
write_tsv(R_means, file = "~/Desktop/MAGNET_GX/avg_CPMcount.txt")

# list of 5 most expressed genes Ensemble ID.
top_5GX <- strsplit(head(R_means$Gene_exp, n = 5), " ")
top_5GX
# subsetting the data for 5 most expressed genes Ensemble ID.
gxData_5 <- gxData %>% filter(row.names(gxData) %in% top_5GX)
# melt the gxData_5 data frame
melt_df <- gather(gxData_5, key = "Sample", value = "CPM")
# merge the melt_df with sampleInfo
merge_df <-
  merge(melt_df,
        sampleInfo,
        by.x = "Sample",
        by.y = 0,
        all.x = TRUE)
# make another data frame using gxData_5 by moving rownames to a column
gene_data <- gxData_5 %>% rownames_to_column(var = "Gene_id")
# melt the data frame
melt_df2 <-
  gather(gene_data, key = "Sample", value = "CPM",-Gene_id)


## Boxplot of 5 most expressed gene for DCM

# filtering data for covariates with DCM disease
DCM_data <- merge_df %>%
  filter(Disease == "DCM") 

p2 <- ggplot(DCM_data, aes(x = Sample, y = CPM, fill = Sample)) +
  geom_boxplot() +
  ggtitle("boxplots of the 5 most expressed gene for DCM") +
  ylab("logCPM") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
p2

ggsave("boxplot_DCM.png") # Export graph

## Boxplot of 5 most expressed gene for Donor

# filtering data for covariates with Donor disease
Donor_data <- merge_df %>%
  filter(Disease == "Donor")

p3 <- ggplot(Donor_data, aes(x = Sample, y = CPM, fill = Sample)) +
  geom_boxplot() +
  ggtitle("boxplots of the 5 most expressed gene for Donor") +
  ylab("logCPM") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
p3

ggsave("boxplot_Donor.png") # Export graph

## Boxplot of 5 most expressed gene for HCM

# filtering data for covariates with HCM disease
HCM_data <- merge_df %>%
  filter(Disease == "HCM") 

p4 <- ggplot(HCM_data, aes(x = Sample, y = CPM, fill = Sample)) +
  geom_boxplot() +
  ggtitle("boxplots of the 5 most expressed gene for HCM") +
  ylab("logCPM") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
p4

# Export graph
ggsave("boxplot_HCM.png") 

## Boxplot of 5 most expressed gene for PPCM

# filtering data for covariates with PPCM disease
PPCM_data <- merge_df %>%
  filter(Disease == "PPCM") 

p5 <- ggplot(PPCM_data, aes(x = Sample, y = CPM, fill = Sample)) +
  geom_boxplot() +
  ggtitle("boxplots of the 5 most expressed gene for PPCM") +
  ylab("logCPM") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
p5
# Export graph
ggsave("boxplot_PPCM.png") 

# 2.b PCA plot showing the sample clustering colored by all available covariates
## Calculating PCA

# perform PCA on gxData
pca_out <- prcomp(gxData, scale = TRUE) 
# shows the top results of rotation column in pca_out.
head(pca_out$rotation)
# Caluculating the contributions of individuals to the principal components:
pca_perc <- round(100 * pca_out$sdev ^ 2 / sum(pca_out$sdev ^ 2), 1)

pca_df <- data.frame(pca_out$rotation,
                     sample = colnames(gxData),
                     Disease = sampleInfo$Disease)
# Exporting PCA results.
write_tsv(pca_df, file = "~/Desktop/MAGNET_GX/pca_disease.txt")

# Plotting PCA
ggplot(pca_df, aes(PC1, PC2, color = sample)) +
  geom_point(size = 1) +
  labs(x = paste0("PC1 (", pca_perc[1], ")"),
       y = paste0("PC2 (", pca_perc[2], ")")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90))
# Export graph
ggsave("PCA.png")
# plotting PCA by covariate of disease
ggplot(pca_df, aes(PC1, PC2, color = Disease, group = Disease)) +
  geom_point(size = 1.5) +
  labs(x = paste0("PC1 (", pca_perc[1], ")"),
       y = paste0("PC2 (", pca_perc[2], ")")) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90))
# Export graph
ggsave("PCA_disease.png")
#-----------------------------------------------------------------------------#
# Question 3: Statistical Analysis
#-----------------------------------------------------------------------------#


# Code for installing required package
## if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("limma")

library(limma)
# Differential gene expression analysis comparing DCM and donors.

# create a list of all samples belonging to DCM and Donor category

# list of samples having DCM
DCM_samples <-
  unique(strsplit(DCM_data$Sample, " ")) 

# list of Donor samples
donor_samples <-
  unique(strsplit(Donor_data$Sample, " ")) 

# combined list of both samples
donor_dcm_samples <-
  append(DCM_samples, donor_samples) 


# keep CPM counts for samples belonging to DCM and Donor cataegory only
trans_gxData <- as.data.frame(t(gxData))
donor_dcm_data <- trans_gxData %>% filter(row.names(trans_gxData)
                                          %in% donor_dcm_samples)
donor_dcm_data <- as.data.frame(t(donor_dcm_data))

# export cpm counts for samples belonging to DCM and Donor cataegory only.
write_tsv(donor_dcm_data, file = "~/Desktop/MAGNET_GX/donor_dcm_data.txt")

# create a design matrix
disease <- factor(rep(c("DCM", "Donor"), c(
  length(DCM_samples),
  length(donor_samples)
)))
design <- model.matrix( ~ 0 + disease)

# perform DGE using limma
fit <- lmFit(donor_dcm_data, design)
# calculating statistical measures 
fit <- eBayes(fit, trend = TRUE)
topTable(fit, coef = ncol(design))
# converting into dataframe for exporting
fit1 = as.data.frame(fit)
# export results
write_tsv(fit1, file = "~/Desktop/MAGNET_GX/differentialGeneExpression.txt")
#-----------------------------------------------------------------------------#
# Question 4: Gene Annotation
#-----------------------------------------------------------------------------#

# Code for installing additional libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)
# Adding gene symbols, gene names and gene descriptions to the data based on
# the provided Ensembl gene identifiers.
# make a copy of original data
gxData_copy <- gxData

# select the ensemble ids
ensemble_ids <- gxDataPCA$EnsemblGeneID

# creart ensemble mart
ensembl <-
  useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# get the gene_id, symbol and description
gene_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  values = ensemble_ids,
  mart = ensembl
)


result_eid <- strsplit(gene_list$ensembl_gene_id, " ")
eid <- strsplit(gxDataPCA$EnsemblGeneID, " ")

print(paste("The no. of missing gene info using Biomart are:",
            length(eid[!eid %in% result_eid])))


# creating annotated data
dcm_data <-
  as.data.frame(t(trans_gxData %>% filter(row.names(trans_gxData)
                                          %in% DCM_samples)))
dcm_avg <- data.frame(DCM_avg = rowMeans(dcm_data))
annotated_data <- merge(gene_list,
                        dcm_avg,
                        by.x = "ensembl_gene_id",
                        by.y = 0,
                        all.y = TRUE)
annotated_data <- merge(
  annotated_data,
  R_means,
  by.x = "ensembl_gene_id",
  by.y = "Gene_exp",
  all.x = TRUE
)


# export the annotated data frame as tab delimited text file
write_tsv(annotated_data, file = "~/Desktop/MAGNET_GX/annotated_dataset.txt")

#-----------------------------------------------------------------------------#
# Question 5: Relative Expression Levels
#-----------------------------------------------------------------------------#


# Assess for each gene in the dataset whether it is expressed above noise level.


female_samples <-
  strsplit(rownames(sampleInfo %>% filter(Sex == "Female")), " ")
female_data <-
  trans_gxData %>% filter(row.names(trans_gxData) %in% female_samples)
female_data <- as.data.frame(t(female_data))

# calculate Baseline value of gene: ENSG00000177732 for the female sample data.
# This gene is an SRY and should only be detected in males. Presence of gene in
# females shows noise.

baseline_val <-
  rowMeans(female_data %>% filter(row.names(female_data) == "ENSG00000177732"))[[1]]
print(paste("The basline value ENSG00000177732 gene is", baseline_val))



# Assess for each gene in the data set whether it is expressed above background (noise) level.
# A 'True' entry  means the gene expression level is above noise level and 
#'False' value means the gene expression is below noise level.
above_noise = as.data.frame(gxData > baseline_val)


# Export the result

write_tsv(above_noise, file = "~/Desktop/MAGNET_GX/noiseThreshold.txt")
