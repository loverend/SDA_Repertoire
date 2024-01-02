##Code to study tensors
## Read in functions
## lauren Overend 
## lauren.overend@oriel.ox.ac.uk
#source("/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/TensorDecomposition/post_processing_fns.R")
# Packages
module load R/4.1.0-foss-2021a 

library(VIM)
library(Hmisc)
library(missForest) 
library(data.table)
library(mice)
library(missCompare) 
library(xcms)
library(Amelia)
library(ggpubr)
library(missMDA)
library(reshape2)
library(ggplot2)
library(corrplot)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(umap)
library(psych)
library(lavaan)
library(cowplot)
library(ggrepel)
library(VennDiagram)
library(ggVennDiagram)
library(grid)
library(gridExtra)


source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/TensorDecomposition/ClusterComponents.R')
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/TensorDecomposition/Model_Components_HvsD.R')

## BCR
outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/"
feature_file <- "/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_FEATURES.txt"
data_file <-  "/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_DATA.txt"
sample_file <- "/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_SAMPLES.txt"
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
type_receptor <- "BCR"
## First lets get the robust components 
my_plots <- cluster_compontents(outputdir, sample_file, feature_file, data_file, eigenvectors)
##Now lets model in health and disease 
components <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/Tensor_Components_FINAL.txt'
p_hd_extended_hvsd <- correlate_components_hvd(components, metadata, outputdir, type_receptor, metahealth, NA)

## Want to plot the correlation for top pairs 
cors_M <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/Tensor_Components_vs_Modules_COR_ADJP.txt', sep="\t")
cors_R <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/TensorDecomposition/Tensor_Components_vs_Modules_COR_R.txt', sep="\t")

colnames(cors_M) <- gsub("Cluster_", "", colnames(cors_M))
colnames(cors_R) <- gsub("Cluster_", "", colnames(cors_R))
## We want to pull out the most positive and most negative significant correlations
cors_R2 <- cors_R
cors_P2 <- cors_M
for(i in 1:length(colnames(cors_R))){
	## Get maximum correlation
	max_cor <- max(cors_R2[,colnames(cors_R)[i]])
	min_cor <- min(cors_R2[,colnames(cors_R)[i]])
	## Make the rest 0 
	cors_R2[,colnames(cors_R)[i]][cors_R2[,colnames(cors_R)[i]] != max_cor & !cors_R2[,colnames(cors_R)[i]]==min_cor] <- NA
	rows <- rownames(cors_R2)[!is.na(cors_R2[,colnames(cors_R)[i]])]
	cors_P2[,colnames(cors_P2)[i]][!rownames(cors_P2) %in% rows] <- NA
}

cors_P2 <- as.matrix(cors_P2)
cors_R2 <- as.matrix(cors_R2)
pdf(paste0(outputdir, "/TensorDecomposition/Tensor_vs_WGCNA_TOPCORS.pdf"), height=12, width=7)
par(mfrow=c(2,1))
corrplot(as.matrix(cors_R), p.mat=as.matrix(cors_M), insig="label_sig",  sig.level=0.05, title =paste0("A. P Adj <0.05"), na.label = " ", method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
corrplot(cors_R2, p.mat=cors_P2, insig="label_sig",  sig.level=0.05, title =paste0("B. P Adj <0.05, Largest R (+ve and -ve)"), na.label = " ", method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
dev.off()

keeps <- colnames(cors_R)
keeps <- c("Component_10", "Component_120", "Component_2", "Component_29", "Component_39", "Component_42", "Component_47", "Component_69", "Component_8", "Component_95", "Component_25", "Component_26", "Component_54", "Component_58", "Component_84")
cors_P3 <- as.matrix(cors_M[,colnames(cors_M)[colnames(cors_M) %in% unique(keeps)]])
cors_R3 <- as.matrix(cors_R[,colnames(cors_R)[colnames(cors_R) %in% unique(keeps)]])
cors_P33 <- as.matrix(cors_P2[,colnames(cors_P2)[colnames(cors_P2) %in% unique(keeps)]])
cors_R33 <- as.matrix(cors_R2[,colnames(cors_R2)[colnames(cors_R2) %in% unique(keeps)]])

pdf(paste0(outputdir, "/TensorDecomposition/Tensor_vs_WGCNA_TOPCORS_significant.pdf"), height=7, width=8)
par(mfrow=c(1,2))
corrplot(cors_R3, p.mat=cors_P3, insig="label_sig",  sig.level=0.05, title =paste0("A. P Adj <0.05"), na.label = " ", method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
corrplot(as.matrix(cors_R33), p.mat=as.matrix(cors_P33),  sig.level=0.05, title =paste0("B. P Adj <0.05, Largest R (+ve and -ve)"), na.label = " ", method = "number",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 1,   number.cex = 0.3)
dev.off()

################################################################################
## Lets extract key relationships and use a ven diagram to plot the component vs module overlap 
## See if we are detecting the same biological signal 
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(RColorBrewer)
###########################
my_vens_positive <- list()
my_vens_negative <- list()
count_positive <- 1
count_negative <- 1
for(i in 1:length(keeps)){
	features <- read.delim(paste0(outputdir, "TensorDecomposition/Tensor_Components_Cluster_", keeps[i], "_KEYFEATURES.txt"))
	features2 <- features$Sample[!is.na(features$labels_use)]
	
	compare <- cors_R33[, colnames(cors_R33)==keeps[i]]
	compare <- compare[!is.na(compare)]
	compare_modules <- names(compare) 

	
	module_assignment <- read.delim('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Clustered_Features_assignment_BCR_PRODUCTIVE_NON_IMPUTED.txt', sep=" ")
	
	for(x in 1:length(compare_modules)){
		module_use <- unlist(str_split( compare_modules[x], "_", 2))[2]
		features_module <- module_assignment$Feature[module_assignment$Cluster==module_use]
		overlap <- intersect(features2, features_module)
		n <- list(A =features2, B=features_module)
		
		a <- gsub("omponent_", "",keeps[i])
		b<- gsub("odule_", "",compare_modules[x])
		names(n) <- c(a , b)
		
		nx <- 3
		color_palette <- brewer.pal(nx, "Pastel1")

		
		r_value <- cors_R2[compare_modules[x], keeps[i]]
		r_value <- round(r_value, digits = 2)
		pq <- ggVennDiagram(n, label_alpha = 0) + scale_fill_distiller(palette = "Pastel1") + scale_color_brewer(palette = "Pastel1")+guides(fill="none")+ggtitle(paste0("Rrmcorr=", r_value)) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
		
		if(r_value >0){
			my_vens_positive[[count_positive]] <- pq
			count_positive <- count_positive +1
		} else {
			my_vens_negative[[count_negative]] <- pq
			count_negative <- count_negative +1
		}
  }
}

theme_settings <- theme(
  plot.margin = unit(c(0, 0, 0, 0), "mm")
)

### Plot negative and positive correlations seperately !!!

##################
pdf(paste0(outputdir, "/TensorDecomposition/Vens.pdf"), height=15, width=15)
plot(do.call("grid.arrange", c(my_vens_positive, ncol=floor(sqrt(length(my_vens_positive))))), theme = theme_settings)
plot(do.call("grid.arrange", c(my_vens_negative, ncol=floor(sqrt(length(my_vens_positive))))), theme = theme_settings)
dev.off()

pdf(paste0(outputdir, "/TensorDecomposition/Composition.pdf"), height=15, width=12)
plot(do.call("grid.arrange", c(my_plots, ncol=5)))
dev.off()

## TCRAB 
outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/"
subsampleddepth <- ""
productivity <- "PRODUCTIVE"
type_use <- "TCRAB"
cluster_compontents(outputdir,subsampleddepth, productivity, type_use)

## TCRGD
outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/"
subsampleddepth <- ""
productivity <- "PRODUCTIVE"
type_use <- "TCRGD"
cluster_compontents(outputdir,subsampleddepth, productivity, type_use)
