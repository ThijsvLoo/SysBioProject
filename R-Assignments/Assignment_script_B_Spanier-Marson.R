#==========================================================================================#
#FINAL SKILLS ASSIGNMENT 
#
#Subject: Experimental Design and Data Management
#Done by: Brittany Spanier-Marson 
#Author's student ID: i6225516
#Institution: Maastricht University 
#Department: Systems Biology (Faculty of Science and Engineering)
#Position: Master's Student 
#==========================================================================================#
#This author strives to adhere to all FAIR principles. If anything is unclear or more information is 
#needed regarding this analysis please contact the author at: 
#b.spanier-marson@student.maastrichtuniversity.nl
#
#
#This file contains code that explores a RNA-sequence dataset of left ventricular 
#biopsies from patients with heart failure and then healthy controls. All participants come from the 
#MAGNET consortium 
#
#For more information regarding the collection of samples please see:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910
#
#The following abbreviations for the diseases are used throughout: 
#DCM- Dilated cardiomyopathy
#HCM- Hypertrophic cardiomyopathy
#PPCM- Peripartum cardiomyopathy
#and "Donor'refers to the healthy controls
#
#Some more information regarding the diseases
#
#Dilated cardiomyopathy is where the ventricle heart muscle is stretched and thin which prevents 
#sufficient pumping of the heart.
#
#Hypertrophic cardiomyopathy is where the heart muscle is abnormally thickened this also  prevents 
#the heart from pumping properly. 
#
#Peripartum cardiomyopathy is a weakening of the heart in the final stages of pregnancy 
#(therefore there are only female samples of PPCM)


---------------------------------------------------------------------------------------------------------------------

# load packages needed for analysis (**If any of these packages are already loaded on your system you can comment relevant code)
install.packages("ggplot2")
library("ggplot2")
install.packages("pcaMethods")
library("pcaMethods")
install.packages("tidyr")
library("tidyr")
install.packages("Biobase")
library("Biobase")
install.packages("limma")
library("limma")
install.packages("BiocManager")
library("BiocManager")
install.packages("biomaRt")
library("biomaRt")
install.packages("gridExtra")
library("gridExtra")

#**Please download and unzip file named MAGNET_GX 
#
#This folder contains 3 text files: 
#1. MAGNET_exonLengths-  which contains the exon lengths of the genes whose expression was measured.
#2. MAGNET_GeneExpressionData_CPM_19112020- which contains the expression data for each gene in each 
#sample (Columns= samples, Rows= Gene)
#3. MAGNET_SampleData_19112020- which contains information about the samples. This includes the known
#covariates: Age, Sex and Ethnicity.

# open file that contains dataset
#**PLEASE ADJUST THIS TO OPEN THE FILE ON YOUR SYSTEM AS WELL ALL FUTURE CODE WHERE AN EXPORT IS DONE
setwd("C:/Users/bspan/Desktop/Systems Biology/Experimental data and design/Assignment/MAGNET_GX")



#============================================PART ONE- Data Import================================================================

#import data--------------------------------------------------------------------------------------------------------- 

gxData.R<- read.table("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names=1, header=T)
sampleInfo.R<- read.table("MAGNET_SampleData_19112020.txt", row.names=1, header= T)

sampleInfoSplit.R<- split(sampleInfo.R, sampleInfo.R$Disease)
DCMsamples<- sampleInfoSplit.R[[1]]
Donorsamples<- sampleInfoSplit.R[[2]]
HCMsamples<- sampleInfoSplit.R[[3]]
PPCMsamples<- sampleInfoSplit.R[[4]]
sampleInfoSplit.R[c(1,2)]  

# Summary of participants characteristics----------------------------------------------------------------------------                     

disease_distribution<- table(sampleInfo.R$Disease)

#Histogram displaying the distribution of age throughout the sample 

hist(sampleInfo.R$Age, breaks=10,xlab= "Age", main= "Age distribution of Sample") 

#Bar graph showing the sex ratios within the diseases and the donor

sex_vs_disease <- table(sampleInfo.R$Sex,sampleInfo.R$Disease) 
print(sex_vs_disease)
ggplot(sampleInfo.R, aes(x=Disease, fill=Sex)) +geom_bar(position="dodge")+ggtitle("Bar graph showing sex distribution in different diseases")

#Bar graph showing the distribution of ethnicity within the diseases and the donor

ethnicity_vs_disease<- table(sampleInfo.R$Ethnicity,sampleInfo.R$Disease) 
print(ethnicity_vs_disease)
ggplot(sampleInfo.R, aes(x=Disease, fill=Ethnicity)) +geom_bar(position="dodge")+ggtitle("Bar graph showing ethnic distribution in different diseases")


#Transform to FPKM values---------------------------------------------------------------------------------------------

geneTotExonLengths <- read.delim("MAGNET_exonLengths.txt", as.is = T, 
                                 row.names = 1) 
all(rownames(geneTotExonLengths) == rownames(gxData.R)) # TRUE (just a check)
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
  #	.t <- 2^(x) * 1E3 / geneTotExonLengths[, 1] # this does the same, but shorter
}
gxData_fpkm <- cpm2fpkm(gxData.R)




#===========================================PART TWO- Diagnostic plots====================================================================

#Create boxplots--------------------------------------------------------------------------------------------------------
#These box plots describe the distribution of expression of genes in each sample. The samples have been divided 
#into disease states so each plot only shows the samples of that disease state. The expression levels are measured
#in CPM. 
#Please see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ for more 
#information regarding CPM and FPKM values. 

require(tidyr)
require(ggplot2)

#Boxplot for DCM data
plotData <- gather(gxData.R[, sampleInfo.R$Disease == "Donor"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()
ggsave("DonorBoxPlot.pdf")

#Boxplot for Donor data
plotData <- gather(gxData.R[, sampleInfo.R$Disease == "DCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()
ggsave("DCMBoxPlot.pdf")

#Boxplot for HCM data
plotData <- gather(gxData.R[, sampleInfo.R$Disease == "HCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()
ggsave("HCMBoxPlot.pdf")

#Boxplot for PPCM
plotData <- gather(gxData.R[, sampleInfo.R$Disease == "PPCM"], 
                   key = "SampleID", value = "CPM")
ggplot(plotData, aes(x = SampleID, y = CPM)) + 
  geom_boxplot() + theme_classic()
ggsave("PPCMBoxPlot.pdf")


#PCA saved as a pdf --------------------------------------------------------------------------------------------------------------------
#These PCAs have been created to analyse the effect of the covariants on gene expression. It can be seen that none 
#of the covariants result in significant differences in gene expression. The output of this section is a pdf named
#"pcaGrids2.pdf"

require(pcaMethods)


pcaGx <- pca(t(gxData.R), nPcs = 10)
plot(pcaGx)
plotPcs(pcaGx, c(1,2))


all(rownames(pcaGx@scores) == rownames(sampleInfo.R)) # TRUE  

# Create ggplot plotting data:
plotData <- cbind(data.frame(pcaGx@scores), sampleInfo.R)

# Create PCA plots:
g1<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = 1, col = Ethnicity))           #can use gridextra to put on the same page 

g2<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = 1, col = Sex))

g3<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = 1, col = Age))

g4<-ggplot(plotData, aes(x = PC1, y = PC2)) + 
  geom_point(aes(size = 1, col = Disease))

pcaPlots<-grid.arrange(g1, g2, g3,g4, ncol=2, nrow=2)

ggsave("pcaGrids2.pdf", plot=pcaPlots)




#================================================PART 3- Statistical analysis==================================================================================
#In this section a differential gene expression analysis is done which is comparing DCM patients and healthy donors.
#The output shows the number of genes which are up or down regulated and those with no significant change in expression 
#in Donors compared to DCM patients. 


require(biomaRt)
require(limma)

eset<- ExpressionSet(data.matrix(gxData.R), phenoData = AnnotatedDataFrame(sampleInfo.R))

design <- model.matrix(~0+Disease, data = sampleInfo.R)

colSums(design)

cm <- makeContrasts(status = DiseaseDonor - DiseaseDCM,
                    levels = design)

# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
diffGeneAnalysis<-data.frame(summary(results))

# Export Results
write.table(diffGeneAnalysis,file="C:/Users/bspan/Desktop/Systems Biology/Experimental data and design/Assignment/MAGNET_GX/DiffGxAnalysis.txt",row.names=FALSE,col.names = FALSE, sep="\t")

#Differential gene analysis is now done while taking sex as a covariant into account ------------------------------------------------------------------------------------------------------------

design2 <- model.matrix(~0+Disease+Sex, data = sampleInfo.R)

colSums(design2)

cm2 <- makeContrasts(status = DiseaseDonor - DiseaseDCM- SexMale,
                    levels = design2)

# Fit the model
fit3 <- lmFit(eset, design2)

# Fit the contrasts
fit4 <- contrasts.fit(fit3, contrasts = cm2)

# Calculate the t-statistics for the contrasts
fit4 <- eBayes(fit4)

# Summarize results
resultsPlusSex <- decideTests(fit4)
diffGeneAnalysisSex<-data.frame(summary(resultsPlusSex))

#Export Results 
write.table(diffGeneAnalysisSex,file="C:/Users/bspan/Desktop/Systems Biology/Experimental data and design/Assignment/MAGNET_GX/DiffGxAnalysis_withSex.txt",row.names=FALSE,col.names = FALSE, sep="\t")




#===============================================PART 4- Additional gene annotation====================================================================================

listEnsembl()
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
geneAnnotation <- getBM(c("external_gene_name", "hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", rownames(geneTotExonLengths), mart)

#Export results 
write.table(geneAnnotation,file="C:/Users/bspan/Desktop/Systems Biology/Experimental data and design/Assignment/MAGNET_GX/AdditionalGeneAnnotation.txt",row.names=FALSE,col.names = TRUE, sep="\t")




#==============================================PART 5- Relative expression levels=====================================================================================


sampleInfoSplitSex.R<- split(sampleInfo.R, sampleInfo.R$Sex)
Femalesamples<- sampleInfoSplitSex.R[[1]]

femaleGxData<-gxData_fpkm[, sampleInfo.R$Sex == "Female"]

Ychrom_GxData<- femaleGxData[c("ENSG00000012817","ENSG00000067048","ENSG00000114374",
"ENSG00000129824","ENSG00000198692","ENSG00000215414","ENSG00000229238","ENSG00000231341","ENSG00000235001","ENSG00000235175"),]

meanYchromEx<-mean(data.matrix(Ychrom_GxData))
meanGxData_fpkm<- data.frame(rowMeans(data.matrix(gxData_fpkm)))
meanGxData_fpkmDF<-data.frame(meanGxData_fpkm)


expressed_above_noise<- data.frame(meanGxData_fpkm>meanYchromEx) #True if on average the gene is expressed above 
                                                                 #background noise. False if on average the gene 
                                                                 #is expressed below background noise 


#Export data
write.table(expressed_above_noise,file="C:/Users/bspan/Desktop/Systems Biology/Experimental data and design/Assignment/MAGNET_GX/GeneExpressedAboveNoise.txt",row.names=TRUE,col.names = TRUE, sep="\t")
