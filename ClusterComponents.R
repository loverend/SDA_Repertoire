library(tidyverse)
library(purrr)
library(Hmisc)
library(dendextend)
library(dendextend)
library(cluster)
library(Rfast)
library(cluster)
library(psych)
library(corrplot)
library(ggrepel)
library(rmcorr)
library(RColorBrewer)
library(VIM)
library(missForest) 
library(data.table)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel) 
library(psych)
library(lavaan)
library(cowplot)
library(ggrepel)



#outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/"
#subsampleddepth <- 1200
#productivity <- "PRODUCTIVE"
#type_use <- "BCR"


cluster_compontents <- function(outputdir, sample_file, feature_file, data_file, eigenvectors){
		
		outputdirectory <- paste0(outputdir, "TensorDecomposition/")
		files_tensor <- list.files(outputdirectory, recursive=TRUE, full.name=TRUE)
		files_use <- grep( "it3000/A", files_tensor, value=TRUE)
		
		sampleids <- read.delim(paste0(sample_file), sep="\t", header=FALSE)		
		feature_ids <- read.delim(paste0(feature_file), sep="\t", header=FALSE)
		
		#-------------------------------------------------------------------------------
		## Using the loading scores!!!!
		data <- NULL
		data<- list()
		for (f in 1:length(files_use)){
			# open each file
			file <-read.delim(files_use[f], sep=" ", header=FALSE)
			colnames(file) <- paste0("Run_", f, "_", colnames(file))
			# append specific column of the file to the dataframe
			data[[f]] <- file
			}
		data <- data.frame(do.call(cbind,data))
		
		# Remove TENSORS that don't contain anything
		data <- data[colSums(data)!=0]
		data <- as.matrix(data)
		
		### Now lets correlate all the compontents 
		correlation <- cor(data, method="pearson", use="pairwise.complete.obs") 
		distance <- as.dist(1-abs(correlation))

		## What method to define structure 
		m <- c( "average", "single", "complete", "ward", "weighted")
		names(m) <- c( "average", "single", "complete", "ward", "weighted")
		print(paste0("Determining which cluster aglomeration method to use"))
		# function to compute coefficient
		ac <- function(x) {
		  agnes(distance, method = x)$ac
		}
		# Best ac score is best agglomerative coefficient 
		# Indicates amount of clustering structure found!
		ac_scores <- data.frame(map_dbl(m, ac))
		colnames(ac_scores) <- "agglomerative_coefficient"
		ac_scores$Method <- rownames(ac_scores)
		## Plotting cluster scores per method
		pdf(paste0(outputdirectory, "CLUSTER_AGLOMERATION_METHOD.pdf"), height=5, width=5)
		plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
		dev.off()
		## Select Method
		
		method_ac_use <- ac_scores$Method[ac_scores$agglomerative_coefficient ==max(ac_scores$agglomerative_coefficient)]
		##rename for using hclust
		if(method_ac_use=="ward"){
			method_ac_use="ward.D2"
		}
		print(paste0("Hericical Cluster Method: ",  method_ac_use))
	
		#### Lets actually do the clustering now 
		clust <- hclust(distance, method=method_ac_use)
		cut_avg <- cutree(clust, h=0.4)
		

		############ Lets colour branches by their run 
		unique_vars <- data.frame(names(distance))	
		order_dendro <- order.hclust(clust)
		unique_vars$run <- str_split_fixed(unique_vars[,1], "_", 3)[,2]
		color_count <- length(unique(unique_vars$run))
		get_palette <- colorRampPalette(brewer.pal(n = 8, name = "Set1"))
		palette <- get_palette(color_count) %>% as.data.frame() %>%rename("color" = ".") %>% rownames_to_column(var = "row_id")
		unique_varsx <- merge(unique_vars, palette, by.x="run", by.y="row_id", sort = FALSE)
		## Now we need to reorder 
		tryx <- unique_varsx[c(order_dendro),]
		run_col <- tryx$color
		names(run_col) <- tryx$names.distance.
		
		### Lets plot 
		pdf(paste0(outputdirectory, "ComponentDendrogram.pdf"), width=17, height=7)
		par(mfrow= c(1,1), mar = c(5,5,3,3))
		plot(clust, cex=0.8, hang = -1, labels=FALSE)
		plot(color_branches(clust,col=c(run_col)))
		abline(h = 0.4,  col="red", lwd=4)
		dev.off()

		## Now we want the PIP scores  
		files_tensor <- list.files(outputdirectory, recursive=TRUE, full.name=TRUE)
		files_use <- grep( "it3000/S", files_tensor, value=TRUE)
		
		pip <- NULL
		pip<- list()
		for (f in 1:length(files_use)){
		# open each file
		file <-read.delim(files_use[f], sep=" ", header=FALSE)
		file <- data.frame(file)
		rownames(file) <- paste0("Run_", f, "_V", rownames(file))
		file <- t(file)
		# append specific column of the file to the dataframe
		pip[[f]] <- file
		}
		pip <- data.frame(do.call(cbind,pip))
		
		## Now we want the Feature Scores
		files_tensor <- list.files(outputdirectory, recursive=TRUE, full.name=TRUE)
		files_use <- grep( "it3000/X", files_tensor, value=TRUE)		
		feat <- NULL
		feat<- list()
		for (f in 1:length(files_use)){
		# open each file
		file <-read.delim(files_use[f], sep=" ", header=FALSE)
		file <- data.frame(file)
		rownames(file) <- paste0("Run_", f, "_V", rownames(file))
		file <- t(file)
		# append specific column of the file to the dataframe
		feat[[f]] <- file
		}
		feat <- data.frame(do.call(cbind,feat))
		
		## Convert the sample scores back to a dataframe
		sample_scores_df <- data.frame(data)
		
		## Get the mean values
		cluster_component_sample_scores <- list()
		cluster_component_pip_scores <- list()
		cluster_component_feat_scores <- list()
		list_count <- 1 # we want to fill in the list sucessively 
		
		for(i in 1:length(unique(cut_avg))){
			cluster_component <- unique(cut_avg)[i]
			runs_in_cluster <- names(cut_avg)[cut_avg==cluster_component]
			## we only want clusters that have components in at least 5 out of the 10 runs 
			runs <- str_split_fixed(runs_in_cluster, "_", 3)
			runs <- runs[,2]

			if(length(unique(runs)) >= 5){
				print(paste0("Cluster ", i, " Component found in at least 5/10 runs"))
				## Sample Scores
				sample_scores <- sample_scores_df[, c(runs_in_cluster)]
				sample_means <- data.frame(rowMeans(sample_scores))
				colnames(sample_means) <- paste0("Cluster_Component_", cluster_component)
				cluster_component_sample_scores[[list_count]] <- sample_means
				
				## Pip Scores
				pip_scores <- pip[, c(runs_in_cluster)]
				pip_means <- data.frame(rowMeans(pip_scores))
				colnames(pip_means) <- paste0("Cluster_Component_", cluster_component)
				cluster_component_pip_scores[[list_count]] <- pip_means
				
				## Feature Scores
				feat_scores <- feat[, c(runs_in_cluster)]
				feat_means <- data.frame(rowMeans(feat_scores))
				colnames(feat_means) <- paste0("Cluster_Component_", cluster_component)
				cluster_component_feat_scores[[list_count]] <- feat_means
				
				## Update the list count
				list_count <- list_count+1
			} else {
				#print("Component found in <5/10 Runs")
			}
		} 
		cluster_component_sample_scores <- data.frame(do.call(cbind,cluster_component_sample_scores))
		cluster_component_pip_scores <- data.frame(do.call(cbind,cluster_component_pip_scores))
		cluster_component_feat_scores <- data.frame(do.call(cbind,cluster_component_feat_scores))
		
		rownames(cluster_component_sample_scores) <- sampleids[,1]
		rownames(cluster_component_pip_scores) <- feature_ids[,1]
		rownames(cluster_component_feat_scores) <- feature_ids[,1]
		
		### ########################
		## Get eigenvectors and plot correllogram 
		eigenvectors <- read.delim(eigenvectors)
		eigenvectors <- eigenvectors[, colnames(eigenvectors) %like% "Module|Component"]
		eigenvectors <-eigenvectors[rownames(eigenvectors) %in% rownames( cluster_component_sample_scores),]
		
		## Lets merge together 
		comps <-cluster_component_sample_scores[rownames(cluster_component_sample_scores) %in% rownames( eigenvectors),]
		
		eigenvectors2 <- eigenvectors
		eigenvectors2$Sample <- rownames(eigenvectors)
		cluster_component_sample_scores2 <- cluster_component_sample_scores
		cluster_component_sample_scores2$Sample <- rownames(cluster_component_sample_scores2)
		
		all_vectors <- merge(eigenvectors2, cluster_component_sample_scores2, by="Sample")

		#################################################################################################################
		#################################################################################################################

		## CLUSTER EIGENVECTORS VS TENSOR COMPONENTS 
		mat_eigenvectors <- all_vectors
		## Add barcode
		for(x in 1:length(mat_eigenvectors$Sample)){
				mat_eigenvectors$Barcode[x] <- str_split_fixed (mat_eigenvectors$Sample[x], "_", 2)[,1]
				if(mat_eigenvectors$Sample[x] %like% "HV"){
				 mat_eigenvectors$Barcode[x] <- paste0(str_split_fixed(mat_eigenvectors$Sample[x], "_", 3)[,1], "_", str_split_fixed(mat_eigenvectors$sample[x], "_", 3)[,2])
				}
			}
		users <- colnames(mat_eigenvectors)[colnames(mat_eigenvectors) %like% "Module|Cluster_Component"]
		rmcorrbigmatt <- rmcorr_mat(Barcode, users, mat_eigenvectors)
		
		type_1 <- "Module" 
		type_2 <- "Cluster_Component"
		
		to_plot <- rmcorrbigmatt$matrix 
		to_plot  <- to_plot[, colnames(to_plot)[colnames(to_plot) %like% type_2]]
		to_plot  <- to_plot[rownames(to_plot)[rownames(to_plot) %like% type_1],]
		
		## This is the comparison we want 
		#####################
		## get p valus 
		x <- rmcorrbigmatt$summary
		x <- x[, c("measure1", "measure2", "p.vals")]
		x <- x[x$measure1 %like% type_1 & x$measure2 %like% type_2 | x$measure1 %like% type_2 & x$measure2 %like% type_1,]
		x$p.vals <- as.numeric(x$p.vals)
		## Adjust p values
		x$p_adjust <- p.adjust(x$p.vals, method = "BH")
		
		## unajusted plot
		m <- matrix(NA, ncol = length(colnames(to_plot)), nrow = length(rownames(to_plot)))
		colnames(m) <- colnames(to_plot)
		rownames(m) <- rownames(to_plot)
		## fill in matrix will p values 
		for(i in 1:length(colnames(m))){
			measure2x <- colnames(m)[i]
				for(j in 1:length(rownames(m))){
					measure1x <- rownames(m)[j] 
					if(measure1x==measure2x){
						m[j,i] <- 1
					} else {
						pfil <- x$p.vals[x$measure1 == measure1x & x$measure2==measure2x]
						if(length(pfil)==0){
							pfil <- x$p.vals[x$measure1 == measure2x & x$measure2==measure1x]
						}
						m[j,i] <- pfil
					}
				}
			}
		
		## adjusted p values 
		madj <- matrix(NA, ncol = length(colnames(to_plot)), nrow = length(rownames(to_plot)))
		colnames(madj) <-  colnames(to_plot)
		rownames(madj) <- rownames(to_plot)
		## fill in matrix will p values 
		for(i in 1:length(colnames(madj))){
			measure2x <- colnames(madj)[i]
				for(j in 1:length(rownames(madj))){
					measure1x <- rownames(madj)[j] 
					if(measure1x==measure2x){
						madj[j,i] <- 1
					} else {
						pfil <- x$p_adjust[x$measure1 == measure1x & x$measure2==measure2x]
						if(length(pfil)==0){
							pfil <- x$p_adjust[x$measure1 == measure2x & x$measure2==measure1x]
						}
						madj[j,i] <- pfil
					}
				}
			}	
		pdf(paste0(outputdirectory, "/Plots/Tensor_vs_WGCNA.pdf"), height=7, width=7)
		corrplot(to_plot, p.mat=m, insig="label_sig",  sig.level=0.05, title =paste0("SDA Tensor Components vs WGCNA Modules\n RMCORR"), method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
		if(min(madj)<0.05){
			corrplot(to_plot, p.mat=madj, insig="label_sig",  sig.level=0.05, title =paste0("SDA Tensor Components vs WGCNA Modules\nRMCORR (BH Adj p)"), method = "circle",  order = 'original', tl.cex = 0.7, mar=c(0,0,2,0), tl.col = "black", pch.cex = 2)
		}
		dev.off()
	
		write.table(m, paste0(outputdirectory, "Tensor_Components_vs_Modules_COR_ADJP.txt"), sep="\t")
		write.table(madj, paste0(outputdirectory, "Tensor_Components_vs_Modules_COR_P.txt"), sep="\t")
		write.table(to_plot, paste0(outputdirectory, "Tensor_Components_vs_Modules_COR_R.txt"), sep="\t")

	
		#################################################################################################################
		#################################################################################################################
		########## Now lets look at what is in each component
		
		## plot the contribution of each component to the component!!!
		print("Plotting Component Contributions") 
		my_plots <- list()
		pdf(paste0(outputdirectory, "Plots/TENSOR_Components.pdf"), height=5, width=5)
		for(i in 1:(length(colnames(cluster_component_feat_scores)))){
				
				coluse <- colnames(cluster_component_feat_scores)[i]
				data_subset <- data.frame(cluster_component_feat_scores[,i])
				colnames(data_subset) <- "Feature_Loading_Score"
				data_subset$Sample <- rownames(cluster_component_feat_scores)	
				
				data_subset_probability <- data.frame(cluster_component_pip_scores[,i])
				colnames(data_subset_probability) <- "Posterior_Inclusion_Probability"
				data_subset_probability$Sample <- rownames(cluster_component_pip_scores)	
				
				## PIP and feature loading score 
				data_all <- merge(data_subset, data_subset_probability, by="Sample")
				data_all <- data_all[order(data_all$Sample),]
				
				## lets get mean and standard deviation of effect 
				meanx <- mean(data_all$Feature_Loading_Score)
				stdv <- sd(data_all$Feature_Loading_Score)
				up_thresh <-  meanx+stdv
				low_thresh <- meanx-stdv
				
				up_thresh2 <-  meanx+(stdv*4)
				low_thresh2 <- meanx-(stdv*4)
				
				## need to sort pior to plotting (cooincide with ggplot ordering!)
				data_all$labels_use <- NA
				data_all$labels_use[data_all$Posterior_Inclusion_Probability >= 0.5 & (data_all$Feature_Loading_Score >= up_thresh2 | data_all$Feature_Loading_Score <= low_thresh2)] <- data_all$Sample[data_all$Posterior_Inclusion_Probability >= 0.5 & (data_all$Feature_Loading_Score >= up_thresh2 | data_all$Feature_Loading_Score <= low_thresh2)]
				data_all$labels_use <- gsub("BCR_READS_", "", data_all$labels_use)
				data_all$labels_use <- gsub("__", ": ", data_all$labels_use)
				data_all$labels_use <- gsub("_", " ", data_all$labels_use)
				data_all$labels_use <- str_wrap(data_all$labels_use, width = 20)
				component_name <- gsub("_", " ", coluse)
				component_name <- gsub("Cluster", "SDA Tensor", component_name)
 				p1 <- ggplot(data_all, aes(y=Feature_Loading_Score, x=Posterior_Inclusion_Probability)) +geom_point(colour="red")+theme_classic() + xlab("Posterior Inclusion Probability") +ylab("Feature Loading Score")+ggtitle(component_name)+
				gghighlight::gghighlight(Posterior_Inclusion_Probability >= 0.5 & (Feature_Loading_Score >= up_thresh | Feature_Loading_Score <= low_thresh), use_direct_label = FALSE)+ 
				geom_hline(yintercept= up_thresh, col="blue") + geom_hline(yintercept =low_thresh, col="blue") + geom_vline(xintercept=0.5, col="blue")+
				geom_label_repel( aes(label = labels_use),  na.rm = TRUE, box.padding = 0.25, max.overlaps = Inf, col="black", size=2)
				plot(p1)
				### Lets save this for look_up tabel 
				component_name <- gsub("_", " ", coluse)
				component_name <- gsub("Cluster ", "", component_name)				
				p2 <- ggplot(data_all, aes(y=Feature_Loading_Score, x=Posterior_Inclusion_Probability)) +geom_point(colour="red")+theme_classic() + xlab("PIP") +ylab("Feature\nLoading")+ggtitle(component_name)+
				gghighlight::gghighlight(Posterior_Inclusion_Probability >= 0.5 & (Feature_Loading_Score >= up_thresh | Feature_Loading_Score <= low_thresh), use_direct_label = FALSE)+ 
				geom_hline(yintercept= up_thresh, col="blue") + geom_hline(yintercept =low_thresh, col="blue") + geom_vline(xintercept=0.5, col="blue")
				data_all$labels_use <- NA
				data_all$labels_use[data_all$Posterior_Inclusion_Probability >= 0.5 & (data_all$Feature_Loading_Score >= up_thresh | data_all$Feature_Loading_Score <= low_thresh)] <- data_all$Sample[data_all$Posterior_Inclusion_Probability >= 0.5 & (data_all$Feature_Loading_Score >= up_thresh | data_all$Feature_Loading_Score <= low_thresh)]
				write.table(data_all, paste0(outputdirectory, "Tensor_Components_", coluse, "_KEYFEATURES.txt"), sep="\t")
				my_plots[[i]] <- p2
				## Conditionally colour labels 
				#con <- ifelse(abs(data_subset$Feature_Loading_Score) >= 0.1, 'red', 'black')
				## Which labels to plot (makes easier to read)	
				#labels_use <- data_subset_probability$labels_use	
				#p1 <- ggplot(data_subset, aes(y=Feature_Loading_Score, x=Sample)) +geom_point()+xlab("Feature Loading Score") +theme()+ggtitle(coluse)  +geom_hline(yintercept=0.1, col="red") +geom_hline(yintercept=-0.1, col="red")+ theme_bw()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5, colour=con)) +scale_x_discrete(labels = labels_use)
				#p2 <- ggplot(data_subset_probability, aes(y=Posterior_Inclusion_Probability, x=Sample)) +xlab("Posterior Inclusion Probability")+geom_point() +theme()+ggtitle(coluse)  +geom_hline(yintercept=0.5, col="green")+ theme_bw()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5, colour=con)) +scale_x_discrete(labels = labels_use)
				#plot(plot_grid(p2, p1, ncol=1))
			}
		dev.off()

		#################################################################################################################
		#################################################################################################################
		print("Saving Files")
		## Now we want to save the output matrix so we can use it to do health vs disease analysis of components 
		cluster_component_sample_scores$sample <- rownames(cluster_component_sample_scores)
		cluster_component_sample_scores$DISEASE <- "SEPSIS"
		cluster_component_sample_scores$DISEASE[cluster_component_sample_scores$sample %like% "HV"] <- "HEALTH"
		cluster_component_sample_scores$DAY <- NA 
		cluster_component_sample_scores$DAY[grep("_1", cluster_component_sample_scores$sample)] <- "Day1"
		cluster_component_sample_scores$DAY[grep("_3", cluster_component_sample_scores$sample)] <- "Day3"
		cluster_component_sample_scores$DAY[grep("_5", cluster_component_sample_scores$sample)] <- "Day5"
		
		## Lets save 
		write.table(cluster_component_sample_scores, paste0(outputdirectory, "Tensor_Components_FINAL.txt"), sep="\t")
		write.table(cluster_component_pip_scores, paste0(outputdirectory, "Tensor_Components_PIP.txt"), sep="\t")
		write.table(cluster_component_feat_scores, paste0(outputdirectory, "Tensor_Components_Feature.txt"), sep="\t")
		
		print("Number of non-0 components")
		print(length(order_dendro))
		
		print("Number of Mean Components")
		print(length(colnames(cluster_component_feat_scores)))
		return(my_plots)
		
}

