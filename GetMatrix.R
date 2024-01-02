## Code to get matrix in tensor format from imputed data
## Uses the imputed data from wgcna pipeline (NOT SCALED!!!!! - we use regression here so not necessary). 


imputed_matrix <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt')
rownames(imputed_matrix) <- gsub("_productive", "", rownames(imputed_matrix))
## Remove JR and bad samples!!!

################################################################
################################################################
### ACTUALLY IN WGCNA we did this after ??
## Lets just leave it will be fine
print("Removing Samples where RNAseq suggests a sample mixup")
bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
imputed_matrix <- imputed_matrix[!rownames(imputed_matrix) %in% bad_ids,]

## remove technicals 	
imputed_matrix <- imputed_matrix[!rownames(imputed_matrix) %like% "JR1795_1003",]

### This is now the same as what we put into wgcna 
### lets save in tensor format 

sample_order <- rownames(imputed_matrix)
feature_order <- colnames(imputed_matrix)

write.table(imputed_matrix, paste0("/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_DATA.txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(sample_order, paste0("/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_SAMPLES.txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(feature_order, paste0("/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/BCR_MATRIX_FOR_TENSOR_PRODUCTIVE_FEATURES.txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

print("DONE PART 1")

## get sample and feature order which are needed for tensor analysis
