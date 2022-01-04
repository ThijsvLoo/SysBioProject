###### R code: Assignment Experimental Design and Data Management ######
#                                                                      #
#     Author: Max L. Jacobs                                            #
#     Maastricht University                                            #
#     Student ID: i6296498                                             #
#     Date: 17-12-2021 (dd-mm-yyyy)                                    #
#                                                                      #
# didn't go exactly my way, I have to admit there are some concepts I  #
# don't understand to their full extend, or how to implement them...   #
####################### end of disclaimer ##############################
# load dependencies #####
# unload all non-base packages
# source: https://stackoverflow.com/a/57981202 , answer by user:7322615, 17-Sep-'19 at 19:53, accessed 10-Dec-'21
for (i in 1:2) {
  base::lapply(names(sessionInfo()$otherPkgs), function(pkgs) {
    detach(
      paste0("package:", pkgs),
      character.only = T,
      unload = T,
      force = T
    )
  })
}

# check if required, passive packages are installed, install if needed
# adapted from https://stackoverflow.com/a/4090208 , answer by user:163053, 03-Nov-'10 at 18:13, accessed 09-Dec-'21
requiredPackages <- c("cowplot", "ggbeeswarm", "viridisLite", "plyr", "data.table", "lares", "pcaMethods", "grid")
newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[, "Package"])]
if (length(newPackages) > 0) {
  cat("installing package:", packages[i])
  install.packages(newPackages)
}

rm(requiredPackages, newPackages)

# load active packages or install and load if needed
lapply(c("tidyverse", "biomaRt", "gridExtra"), require, character.only = T)
require("Biobase", character.only = T, attach.required = F)
# set R-session to interactive state for request on OS and file path
if (interactive() == F) {
  interactive <- T
}

# define functions & misc objects #####

## retrieve legend from plot and create a blank plot for gridExtra plot arrangements
## source: unknown
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
blankPlot <- ggplot() +
  geom_blank(aes(1, 1)) +
  cowplot::theme_nothing()

## convert CPM to FPKM for inter-experiment comparability
## source: kindly provided as course material by Michiel Adriaens, Maastricht University, MaCSBio
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
}
# get_legend() and blankPlot currently unused, remove if feasible ##################################################
# set working directory on different OS #####

# for me the working directory is set by the R_project path
# for everyone else this is an opportunity to run the script without changing the wd inside the script

# don't look at it too long, just run it as described in the disclaimer (row 8)

# require wd path

if (exists("path")) { # && path != "") {
  cat(
    "Object 'path' has been identified.\npath =", paste(path), "\nYour current working directory is:",
    getwd(), "\nDo you want to specify your path again?\n[1] Yes  [2] No \n > "
  )
  reset_wd <- readline()
  if (reset_wd == 1) {
    cat("\n\n##### set working directory #####\nPlease enter the path for your working directory or enter 'skip' or '0' to skip this step. Skipping might result in an error message and the script to be terminated.\nThe working directory should include a folder 'MAGNET_GX', which in turn contains the data.\nYou can type the desired path in the console. Copying and pasting works best, as it avoids typos.\nNo adjustments for OS differences in pathway syntax are needed. Confirm with the enter key.\n")
    path <- readline()
  } else if (reset_wd == 2) {
    cat("\nYour current working directory is:", getwd())
  } else {
    cat("Why did you type ", reset_wd, "? You know it is invalid", sep = "")
    for (i in 1:3) {
      Sys.sleep(1)
      cat(".")
    }
    Sys.sleep(3)
    cat("\nSure, go on without resetting the current working directory.\nYour current working
        directory is:", getwd())
    Sys.sleep(1.5)
  }
} else {
  cat("\n\n##### set working directory #####\nPlease enter the path for your working directory or enter 'skip' or '0' to skip this step. Skipping might result in an error message and the script to be terminated.\nThe working directory should include a folder 'MAGNET_GX', with the respective data inside it.\nYou can type the desired path in the console. Copying and pasting works best, as it avoids typos.\nNo adjustments for OS differences in pathway syntax are needed. Confirm with the enter key.\n")
  path <- readline()
  if (path == "") {
    stop("You didn't specify a working directory. The script will be terminated.")
  } else if (path == "skip" || path == "0") {
    cat("\nCongratulations! You've decided to skip setting your working directory.\nSkipping might result in an error message and the script to be terminated when trying to import the data.\nYour current working directory is:", getwd())
  } else {
    setwd(path)
    cat("\nYour working directory has been set to:", getwd())
  }
}

# require OS


if (exists("OS") == T) {
  cat(
    "Object 'OS' has been identified.\n[1] - Windows\n[2] - Mac\n[3] - Linux\nYour OS is set to: ", OS,
    "\nDo you want to choose your operating system again?\n [1] Yes  [2] No"
  )
  reset_OS <- readline()
  if (reset_OS == 1) {
    cat("\n\n##### ask for OS #####\nWhich operating system are you using?\n[1] - Windows\n[2] - Mac\n[3] - Linux \nPlease type the respective number in the console.\nConfirm with the enter key.\n")
    OS <- readline()
  } else if (reset_OS == 2) {
    cat("Alrighty then! Hope everything works for you!")
  } else {
    cat("Why did you type ", reset_OS, "? You know it is invalid", sep = "")
    for (i in 1:3) {
      Sys.sleep(1)
      cat(".")
    }
    Sys.sleep(3)
    cat("\nSure, go on without resetting the chosen operating system.")
    Sys.sleep(1.5)
  }
} else {
  cat("\n\n##### ask for OS #####\nWhich operating system are you using?\n[1] - Windows\n[2] - Mac\n[3] - Linux\nPlease type the respective number in the console.\nConfirm with the enter key.\n")
  OS <- readline()
}

if (OS == 1) {
  # path for Windows
  cat("\nFile paths are now optimized for Windows syntax.\n")
  path_SampleData <- "MAGNET_GX\\MAGNET_SampleData_19112020.txt"
  path_gxData <- "MAGNET_GX\\MAGNET_GeneExpressionData_CPM_19112020.txt"
  path_exonLengths <- "MAGNET_GX\\MAGNET_exonLengths.txt"
} else if (OS == 2 || OS == 3) {
  # path for Mac or Linux
  cat("\nPathways are now optimized for Mac and Linux syntax.\n")
  path_SampleData <- "MAGNET_GX/MAGNET_SampleData_19112020.txt"
  path_gxData <- "MAGNET_GX/MAGNET_GeneExpressionData_CPM_19112020.txt"
  path_exonLengths <- "MAGNET_GX/MAGNET_exonLengths.txt"
} else {
  stop("You didn't specify your operating system correctly. The script will be terminated.")
}

# import data #####
cat("\nData is being imported. Please wait.\n")

sampleInfo <- read.delim(paste(path_SampleData), header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
geneTotExonLengths <- read.delim(paste(path_exonLengths), as.is = T, row.names = 1)

# gxdata is quite large and takes time to load, so I skip it if it already exists in the GEnv
if (exists("gxData") == F) {
  gxData <- read.delim(paste(path_gxData), header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
}

# remove some of the smaller files we don't need anymore. Make memory space whenever possible!
rm(reset_wd, reset_OS, path_exonLengths, path_SampleData, path_gxData)

# check imported data structure for compatibility and completeness
# only useful if manually running the script
all(rownames(geneTotExonLengths) == rownames(gxData))
all(rownames(sampleInfo) == colnames(gxData))

## DATA EXPLORATION #######
# summarize sample info #####
cat("\nSummarizing sample Info...\n")

info_summary_fac_gr <- sampleInfo %>%
  group_by(Sex, Disease) %>%
  summarise(
    mean_age = mean(Age),
    n_Total = n(),
    n_Caucasian = sum(Ethnicity == "Caucasian"),
    n_African.American = sum(Ethnicity == "African.American")
  )

print(info_summary_fac_gr)

# visualize age distribution in box plot #####
cat("\nVisualizating sample info...\n")

info_to_vis <- sampleInfo %>%
  mutate(
    Disease_Sex = paste(Disease, Sex, sep = " "),
  )

hBox <- info_to_vis %>%
  ggplot(
    aes(Disease_Sex, Age)
  ) +
  geom_boxplot(outlier.shape = NA) + # outliers are hidden, as data points are visualized using the beeswarm plot
  ggbeeswarm::geom_beeswarm(
    aes(color = Ethnicity),
    alpha = 0.6,
    size = 2
  ) +
  scale_color_manual(values = viridisLite::viridis(2, option = "D", begin = 0.15, end = 0.85)) +
  theme_classic() +
  ylim(c(0, plyr::round_any(max(sampleInfo$Age), 5))) +
  labs(
    title = "Age distribution of particpants",
    subtitle = "by sex, ethnicity and disease state"
  ) +
  xlab("Sex and disease state") +
  ylab("Age") +
  cowplot::theme_cowplot(12) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = c(0.75, 0.1),
    legend.background = element_blank(),
    legend.key.size = unit(20, "pt"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# create table with number of participants per group #####

info_to_gtable <- info_summary_fac_gr %>%
  ungroup() %>%
  mutate(
    Disease_Sex = paste(Disease, Sex, sep = " "),
    .before = 1
  ) %>%
  dplyr::select(!c("Sex", "Disease", "mean_age")) %>%
  arrange(Disease_Sex) %>%
  column_to_rownames("Disease_Sex")

info_to_gtable_t <- info_to_gtable %>%
  mutate(
    n_Total = as.character(n_Total),
    n_Caucasian = as.character(n_Caucasian),
    n_African.American = as.character(n_African.American)
  ) %>%
  data.table::transpose()

rownames(info_to_gtable_t) <- info_to_gtable %>%
  colnames() %>%
  gsub(pattern = "^n_", replacement = "") %>%
  gsub(pattern = "\\.", replacement = " ")
colnames(info_to_gtable_t) <- rownames(info_to_gtable)

table_minimal <- gridExtra::ttheme_minimal(
  core = list(padding = unit(c(5, 5), "pt")),
  base_size = 12
)

info_gtable <- tableGrob(
  d = info_to_gtable_t,
  rows = rownames(info_to_gtable_t),
  cols = colnames(info_to_gtable_t),
  theme = table_minimal
)

table_title <- grid::grid.text(
  "Number of participants by group",
  gp = grid::gpar(fontsize = 14, fontface = "bold"),
  vjust = 1.5,
  hjust = 1.5
)

# combine and save the graph and table #####

ageDist_nPop <- grid.arrange(hBox, table_title, info_gtable,
  ncol = 1, nrow = 3,
  heights = c(10, 0.5, 2.5)
)

# save graph and table
lares::export_plot(
  ageDist_nPop,
  name = "gxData_summary_ageDistribution_nParticipants",
  width = 10,
  height = 10,
  res = 300,
  format = "png",
  dir = getwd(),
  quiet = T
)

rm(table_title, info_gtable, hBox, ageDist_nPop)

# convert cpm to fpkm #####
cat("\nConverting CPM to FPKM...\n")

gxData_fpkm <- cpm2fpkm(gxData)
gxData_fpkm <- log2(gxData_fpkm)
## DATA ANALYSIS #######
# box plot - gene expression #####
cat("\nVisualizing gene Expression...\n")

sampleInfo_ext <- sampleInfo %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    Disease = factor(Disease, levels = c(unique(sampleInfo$Disease)))
  )

Expr_plotData <- gxData %>%
  pivot_longer(cols = 1:length(gxData), names_to = "SampleID", values_to = "CPM") %>%
  right_join(sampleInfo_ext[, c(1, 3, 4)], by = "SampleID")

for (i in 1:length(levels(sampleInfo_ext$Disease))) {
  loopPlotData <- Expr_plotData %>%
    filter(Disease == levels(sampleInfo_ext$Disease)[i])

  plot_geneExp_byDis <- ggplot(data = loopPlotData, (aes(x = SampleID, y = CPM))) +
    geom_boxplot() +
    ggtitle(paste("gene expression in disease group", levels(sampleInfo_ext$Disease)[i])) +
    facet_wrap(~Sex, scales = "free_x") +
    theme_classic() +
    theme(
      axis.text = element_text(angle = 90)
    )

  lares::export_plot(
    plot_geneExp_byDis,
    name = paste("cpm_geneExpression_", levels(sampleInfo_ext$Disease)[i], sep = ""),
    width = 17,
    height = 10,
    res = 300,
    format = "png",
    dir = getwd(),
    quiet = T
  )
}

rm(loopPlotData, plot_geneExp_byDis, Expr_plotData)

# PCA #####
cat("\nCalculating and plotting PCAs...\n")

# How do we define the number of PCs to analyze? As far as I am concerned, we only have 4 different variables (Age, Sex, Disease, Ethnicity), 
# so we can only have 4 dimensions in the PCA
pcaRes <- pcaMethods::pca(t(gxData), nPcs = 4)
plot(pcaRes)
#  pcaMethods::plotPcs(pcaRes, c(1,2))

# retrieve residual scores
PCAplotData <- cbind(data.frame(pcaRes@scores), sampleInfo)

# how do I retrace a principal component to a variable?


# The size variation of the points we used in the skill session is horrible for interpretation from my point of view,
# so I made a second plot using colors for age as well and combined Disease and Age plots

################## PCA plots loop (experimental) ##################
# I would like to put this in a loop, but creating all possible combinations without duplication or combining a PC with itself seems difficult
# and renaming objects during a loop according to its content appears impossible
#
# PC_comb <- t(combn(x = (1:4), 2))
# PC_comb_names <- PC_comb %>% data.frame()
# base::colnames(PC_comb_names) <- c("c1", "c2")
# PC_comb_names <- PC_comb_names %>% transmute(comb = paste("PC", c1, ".", c2, "_", sep = ""))
#
# for (i in 1:nrow(PC_comb_names)) {
#   
# plot_PC1.2_Dis <- ggplot(PCAplotData, aes(x = PCAplotData[, PC_comb[i, 1]], y = PCAplotData[, PC_comb[i, 2]])
#                          ) +
#   geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
#   scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)
#                      ) +
#   xlab(paste("PC", PC_comb[i, 1], sep = " ")
#        ) +
#   ylab(paste("PC", PC_comb[i, 2], sep = " ")
#        ) +
#   theme(
#     panel.background = element_blank(),
#     plot.background = element_blank()
#   )
# 
# plot_PC1.2_Age <- ggplot(PCAplotData, aes(x = PCAplotData[, PC_comb[i, 1]], y = PCAplotData[, PC_comb[i, 2]])
#                          )+
#   geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
#   scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), 
#                                                       option = "D", begin = 0.1, end = 0.9)
#                         ) +
#   xlab(paste("PC", PC_comb[i, 1], sep = " ")
#        ) +
#   ylab(paste("PC", PC_comb[i, 2], sep = " ")
#        ) +
#   theme(
#     panel.background = element_blank(),
#     plot.background = element_blank()
#   )
# 
# }
################## PCA plots loop end ##################

# so I still have to plot all individually

# Plot PC1-PC2
plot_PC1.2_Dis <- ggplot(PCAplotData, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC1.2_Age <- ggplot(PCAplotData, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Plot PC1-PC3
plot_PC1.3_Dis <- ggplot(PCAplotData, aes(x = PC1, y = PC3)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC1.3_Age <- ggplot(PCAplotData, aes(x = PC1, y = PC3)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Plot PC1-PC4
plot_PC1.4_Dis <- ggplot(PCAplotData, aes(x = PC1, y = PC4)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC1.4_Age <- ggplot(PCAplotData, aes(x = PC1, y = PC4)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Plot PC2-PC3
plot_PC2.3_Dis <- ggplot(PCAplotData, aes(x = PC2, y = PC3)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC2.3_Age <- ggplot(PCAplotData, aes(x = PC2, y = PC3)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Plot PC2-PC4
plot_PC2.4_Dis <- ggplot(PCAplotData, aes(x = PC2, y = PC4)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC2.4_Age <- ggplot(PCAplotData, aes(x = PC2, y = PC4)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Plot PC3-PC4
plot_PC3.4_Dis <- ggplot(PCAplotData, aes(x = PC3, y = PC4)) +
  geom_point(aes(col = Disease), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_manual(values = viridisLite::viridis(4, option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot_PC3.4_Age <- ggplot(PCAplotData, aes(x = PC3, y = PC4)) +
  geom_point(aes(col = Age), size = 3, alpha = 0.5, stroke = 1.1) +
  scale_color_gradientn(colors = viridisLite::viridis(length(unique(sampleInfo$Age)), option = "D", begin = 0.1, end = 0.9)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

rm(pcaRes, PCAplotData)
# export PCA plots #####
# I couldn't export this using the same method as before, as the y axis appears to be not finite
# so I use the generic method

PCplots <- ls()[grep("plot_PC\\d{1}\\.\\d{1}_[Dis|Age]", ls())]
PCplots_export <- sub("_((Dis)|(Age))", "", PCplots) %>% unique()

for (i in PCplots_export) {
  png(
    file = paste(i, "_Dis_Age.png", sep = ""),
    width = 25, height = 10, units = "cm", res = 300
  )
  grid.arrange(get(paste(i, "_Dis", sep = ""), envir = .GlobalEnv), get(paste(i, "_Age", sep = ""), envir = .GlobalEnv),
               ncol = 2, nrow = 1
  )
  dev.off()
}

rm(
  plot_PC1.2_Age, plot_PC1.2_Dis, plot_PC1.3_Age, plot_PC1.3_Dis, plot_PC1.4_Age, plot_PC1.4_Dis,
  plot_PC2.3_Age, plot_PC2.3_Dis, plot_PC2.4_Age, plot_PC2.4_Dis, plot_PC3.4_Age, plot_PC3.4_Dis
)

# gene annotation #####
cat("\nRetrieving gene information from BioMart...\n")
# prepare biomart functions for ensembl human database
if (exists("geneInformation") == F) {
  ensembl_hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  geneInformation <- data.frame(rownames(geneTotExonLengths))
  base::colnames(geneInformation) <- "ensembl_gene_id"
  # check dataframe size
  length(unique(geneInformation$ensemble_gene_id)) == length(geneInformation$ensemble_gene_id)

  # I wanted to add "phenotype_description" and "entrezgene_accession" but it messed with the dataframe structure and length respectively, as it caused a lot of duplicates in ensebl_gene_id...
  retrieve_from_BM <- c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "chromosome_name", "description")

  BMattributes <- getBM(attributes = retrieve_from_BM, filters = "ensembl_gene_id", values = geneInformation$ensembl_gene_id, mart = ensembl_hs)
  geneInformation <- full_join(geneInformation, BMattributes, by = "ensembl_gene_id")


  # commands used for manual curation of gene objects
  # check number of rows compared to our data, then check for duplicates in the retrieved and joined attribute data
  nrow(geneInformation) == nrow(gxData_fpkm)
  any(duplicated(geneInformation$ensembl_gene_id) == T)
  which(duplicated(geneInformation$ensembl_gene_id) == T)

  # take a look at the gene indicated as duplicate, as well as its neighbors
  geneInformation[which(duplicated(geneInformation$ensembl_gene_id) == T) - 1, ]
  geneInformation[which(duplicated(geneInformation$ensembl_gene_id) == T), ]
  geneInformation[which(duplicated(geneInformation$ensembl_gene_id) == T) + 1, ]

  # manually remove the duplicate, here I chose to remove HGNC STRA6LP, as the description is "SUGT1P4-STRA6LP readthrough"
  geneInformation <- geneInformation[-which(duplicated(geneInformation$ensembl_gene_id) == T), ]
  any(duplicated(geneInformation$ensembl_gene_id) == T)

  # clean dataframe geneInformation
  geneInformation <- remove_rownames(geneInformation)
  geneInformation <- column_to_rownames(geneInformation, var = "ensembl_gene_id")
}

# remove unnecessary objects
rm(BMattributes, retrieve_from_BM, ensembl_hs)

# differential gene expression #####
cat("\nCalculating differential gene expression...\n")
# Perform at least a differential gene expression analysis comparing DCM patients and healthy donors.
#---- use cpm, not fpkm? ----#
gxMatrix <- as.matrix(gxData)

gx_eset <- ExpressionSet(
  assayData = gxMatrix,
  phenoData = AnnotatedDataFrame(sampleInfo),
  featureData = AnnotatedDataFrame(geneInformation)
)

rm(gxMatrix)

gx_subEset_Donor <- gx_eset[, gx_eset$Disease == "Donor"]
gx_subEset_PPCM <- gx_eset[, gx_eset$Disease == "PPCM"]
gx_subEset_HCM <- gx_eset[, gx_eset$Disease == "HCM"]
gx_subEset_DCM <- gx_eset[, gx_eset$Disease == "DCM"]

gx_subEset_DonorPPCM <- cbind(exprs(gx_subEset_Donor), exprs(gx_subEset_PPCM))
gx_subEset_DonorHCM <- cbind(exprs(gx_subEset_Donor), exprs(gx_subEset_HCM))
gx_subEset_DonorDCM <- cbind(exprs(gx_subEset_Donor), exprs(gx_subEset_DCM))



# my attempt to create a design matrix, not remembering to use model.matrix() ####
# as expected, they are the same, model.matrix() just adds information on contrast treatment?
# design_fpkm_DonorVsAll <- data.frame(Disease = sampleInfo$Disease, row.names = rownames(sampleInfo)) %>%
#   mutate(
#     Donor = as.numeric(Disease == "Donor"),
#     PPCM = as.numeric(Disease == "PPCM"),
#     HCM = as.numeric(Disease == "HCM"),
#     DCM = as.numeric(Disease == "DCM")
#   ) %>%
#   dplyr::select(-Disease)
# ####

library(limma)
# create design matrix the intended way
design_DonorVsAll <- model.matrix(~ 0 + Disease, data = pData(gx_eset))

# for all diseases vs control
fit <- lmFit(gx_eset, design_DonorVsAll)
cont.matrix <- makeContrasts(
  DonorVsPPCM = DiseaseDonor - DiseasePPCM, 
  DonorVsHCM = DiseaseDonor - DiseaseHCM, 
  DonorVsDCM = DiseaseDonor - DiseaseDCM,
  levels = design_DonorVsAll)
fit_DonorVsAll <- contrasts.fit(fit, contrasts = cont.matrix)
fit_DonorVsAll <- eBayes(fit_DonorVsAll)
diffExpr_DonorVsAll <- topTable(fit_DonorVsAll, adjust = "BH", number = nrow(gxData))
results_DonorVsAll <- decideTests(fit_DonorVsAll)
summary(results_DonorVsAll)

# I don't know how to retrieve the log-fold change from multiple comparisons

# DiseaseDonor vs DiseasePPCM

design_DonorVsPPCM <- model.matrix(~ 0 + Disease, data = pData(gx_eset)[which(sampleInfo$Disease == "Donor" | sampleInfo$Disease == "PPCM"),])

fit <- lmFit(gx_subEset_DonorPPCM, design_DonorVsPPCM)
cont.matrix <- makeContrasts(Donor_vs_PPCM = DiseaseDonor - DiseasePPCM, levels = design_DonorVsPPCM)
fit_DonorVsPPCM <- contrasts.fit(fit, cont.matrix)
fit_DonorVsPPCM <- eBayes(fit_DonorVsPPCM)
diffExpr_DonorVsPPCM <- topTable(fit_DonorVsPPCM, adjust = "BH", number = nrow(gxData)) %>%
  mutate(status = "DonorVsPPCM")

# DiseaseDonor vs DiseaseHCM
design_DonorVsHCM <- model.matrix(~ 0 + Disease, data = pData(gx_eset)[which(sampleInfo$Disease == "Donor" | sampleInfo$Disease == "HCM"),])

fit <- lmFit(gx_subEset_DonorHCM, design_DonorVsHCM)
cont.matrix <- makeContrasts(Donor_vs_HCM = DiseaseDonor - DiseaseHCM, levels = design_DonorVsHCM)
fit_DonorVsHCM <- contrasts.fit(fit, cont.matrix)
fit_DonorVsHCM <- eBayes(fit_DonorVsHCM)
diffExpr_DonorVsHCM <- topTable(fit_DonorVsHCM, adjust = "BH", number = nrow(gxData)) %>%
  mutate(status = "DonorVsHCM")

# DiseaseDonor vs DiseaseDCM
design_DonorVsDCM <- model.matrix(~ 0 + Disease, data = pData(gx_eset)[which(sampleInfo$Disease == "Donor" | sampleInfo$Disease == "DCM"),])


fit <- lmFit(gx_subEset_DonorDCM, design_DonorVsDCM)
cont.matrix <- makeContrasts(Donor_vs_DCM = DiseaseDonor - DiseaseDCM, levels = design_DonorVsDCM)
fit_DonorVsDCM <- contrasts.fit(fit, cont.matrix)
fit_DonorVsDCM <- eBayes(fit_DonorVsDCM)
diffExpr_DonorVsDCM <- topTable(fit_DonorVsDCM, adjust = "BH", number = nrow(gxData)) %>%
  mutate(status = "DonorVsDCM")

#rm(fit, cont.matrix, gx_subEset_DonorDCM, gx_subEset_DonorHCM, gx_subEset_DonorPPCM)

# Correct for relevant co-variates and add comments to the scripts explaining your decision
# population stratification by ethnicity?

# I do not understand how to determine possible co-variates. I sure can just build the model with each of the covariates and take a look how it changes correlations.
# But that seems like unnecessary effort, considering I already have done the PCA. But how do I retrace which Principal component links to which possible co-variate?

# relative expression #####
cat("\nCalculating background...\n")
# expressed above background?
# retrieve background from y chromosomal gene expression in females
yChromGenes <- rownames_to_column(geneInformation, var = "ensembl_gene_id") %>%
  dplyr::select(c("chromosome_name", "ensembl_gene_id"))

background <- gxData %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  left_join(yChromGenes, by = "ensembl_gene_id") %>%
  filter(
    chromosome_name == "Y"
  ) %>%
  dplyr::select(-c(chromosome_name, rownames(sampleInfo[which(sampleInfo$Sex == "Male"), ]))) %>%
  column_to_rownames(var = "ensembl_gene_id") 

# there is not really an optimal way to assess a specific value, 
# the arithmetic mean or the median over all samples and genes would be a possibility.
# As they are quite low I decided to use the threshold-value for the third and fourth quantile as baseline
background <- as.vector(t(background))
baseline <- max(background)#quantile(background)[4]

gxBackground_Donor <- (exprs(gx_subEset_Donor) - baseline) %>% data.frame() %>% filter(. > 0) 
gxBackground_PPCM <- (exprs(gx_subEset_PPCM) - baseline) %>% data.frame() %>% filter(. > 0) 
gxBackground_HCM <- (exprs(gx_subEset_HCM) - baseline) %>% data.frame() %>% filter(. > 0) 
gxBackground_DCM <- (exprs(gx_subEset_DCM) - baseline) %>% data.frame() %>% filter(. > 0) 

genesAboveBackground <- c(rownames(gxBackground_Donor),
  rownames(gxBackground_PPCM),
  rownames(gxBackground_HCM),
  rownames(gxBackground_DCM))
genesAboveBackground_All <- genesAboveBackground[duplicated(genesAboveBackground)]

# I don't know how to filter the data properly here, as we might have different genes expressed above noise level throughout the different diseases
# I probably could group them, but comparing multiple sets of gene IDs with differing length for each disease appears 

##### export results #####
cat("\nCreating summary and exporting data...\n")

ExprsResults <- diffExpr_DonorVsAll %>% 
  mutate(AboveBG_allSamples = rownames(diffExpr_DonorVsAll) %in% genesAboveBackground_All,
         AboveBG_Donor = rownames(diffExpr_DonorVsAll) %in% rownames(gxBackground_Donor),
         AboveBG_PPCM = rownames(diffExpr_DonorVsAll) %in% rownames(gxBackground_PPCM),
         AboveBG_HCM = rownames(diffExpr_DonorVsAll) %in% rownames(gxBackground_HCM),
         AboveBG_DCM = rownames(diffExpr_DonorVsAll) %in% rownames(gxBackground_DCM),
         logFC_DonorVsPPCM = diffExpr_DonorVsPPCM$logFC,
         logFC_DonorVsHCM = diffExpr_DonorVsHCM$logFC,
         logFC_DonorVsDCM = diffExpr_DonorVsDCM$logFC,
         AveExpr_Donor = rowMeans(exprs(gx_subEset_Donor)),
         AveExpr_PPCM = rowMeans(exprs(gx_subEset_PPCM)),
         AveExpr_HCM = rowMeans(exprs(gx_subEset_HCM)),
         AveExpr_DCM = rowMeans(exprs(gx_subEset_DCM))
  ) %>% 
  relocate("logFC_DonorVsPPCM", .after = "DonorVsPPCM") %>%
  relocate("logFC_DonorVsHCM", .after = "DonorVsHCM") %>%
  relocate("logFC_DonorVsDCM", .after = "DonorVsDCM") %>%
  relocate(c(AveExpr_Donor, AveExpr_PPCM, AveExpr_HCM, AveExpr_DCM), .after = "AveExpr")

library(openxlsx)

DiffExpr_summary_gxData <- createWorkbook(
  title = "DiffExpr_summary_gxData"
)
addWorksheet(DiffExpr_summary_gxData, sheetName = "Summary")
writeData(DiffExpr_summary_gxData, "Summary", ExprsResults)
saveWorkbook(DiffExpr_summary_gxData, file = "DiffExpr_summary_gxData.xlsx", overwrite = T)
detach("package:openxlsx", unload = T, character.only = T)
detach("package:limma", unload = T, character.only = T)

write.csv(ExprsResults, file = "DiffExpr_summary_gxData.csv")

