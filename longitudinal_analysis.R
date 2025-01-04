rm(list = ls())
library(Seurat) # Seurat_4.3.0 # SeuratObject_4.1.4
library(dplyr) # dplyr_1.1.4
set.seed(1)
options(width=160)

# Load 7 and 17 months data and make the cell type and genotype conditions conistent to merge both datasets

# 7 month data (t0)
load("/home/m.ali/Projects/UL/mice_snRNAseq_cortex/enrico_script/sctrans2.RData")
month7 <- sctrans2
month7
# An object of class Seurat
# 38561 features across 44582 samples within 2 assays
# Active assay: SCT (18682 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, umap, tsne

# 17 month data (t1)
load("/home/m.ali/Projects/UL/mice_snRNAseq_cortex_17months/data/sctrans2.RData")
month17 <- sctrans2
month17
# An object of class Seurat
# 36730 features across 28190 samples within 2 assays
# Active assay: SCT (17681 features, 3000 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
#  3 dimensional reductions calculated: pca, umap, tsne

table(Idents(month7))
table(month7$classint)
table(month7$conditions)
month7$classint <- gsub(" ", "_", month7$classint)
month7$conditions <- gsub(" ", "_", month7$conditions)
month7$conditions <- paste0(month7$conditions, "_M7", sep="")
Idents(month7) <- paste0(month7$classint, "_", month7$conditions, sep="")

table(Idents(month17))
table(month17$classint)
table(month17$conditions)
month17$classint <- gsub(" ", "_", month17$classint)
month17$conditions <- gsub(" ", "_", month17$conditions)
month17$conditions <- paste0(month17$conditions, "_M17", sep="")
Idents(month17) <- paste0(month17$classint, "_", month17$conditions, sep="")

table(month17$classint %in% month7$classint)
#  TRUE
# 28190

data_combined <- merge(month7, y = month17, add.cell.ids = c("M7", "M17"), project = "Longitudinal")
data_combined
# An object of class Seurat
# 40089 features across 72772 samples within 2 assays
# Active assay: SCT (19326 features, 0 variable features)
#  3 layers present: counts, data, scale.data
#  1 other assay present: RNA
 
as.data.frame(table(data_combined$orig.ident))
table(Idents(data_combined))
idents_table <- as.data.frame(table(Idents(data_combined)))

data_combined_copy <- data_combined
data_combined <- PrepSCTFindMarkers(object = data_combined, assay = "SCT", verbose = TRUE)
save(data_combined, file="ThyTau22_Combined_Month_7_17.RData")

# function to get unique DEGs that are FDR < 0.05 in TG condition and P > 0.05 in WT condition
unique_DEGs <- function(df1, df2, target_fdr = 0.05, exclude_pval = 0.5) {
 df1 <- df1[,c("gene", "avg_log2FC", "p_val", "p_val_adj")]
 colnames(df1) <- c("gene", "avg_log2FC_d1", "p_val_d1", "p_val_adj_d1")
 df2 <- df2[,c("gene", "avg_log2FC", "p_val", "p_val_adj")]
 colnames(df2) <- c("gene", "avg_log2FC_d2", "p_val_d2", "p_val_adj_d2")
 combinde_df <- inner_join(df1, df2, by="gene")
 dim(combinde_df) # 7178    7
 # get same direction DEGs, we will apply "fdr < 0.05" in dataset1 and "p > 0.5" in dataset2 filters on them
 same_direction_significant <- rbind(combinde_df[combinde_df$avg_log2FC_d1 > 0 & combinde_df$avg_log2FC_d2 > 0,], combinde_df[combinde_df$avg_log2FC_d1 < 0 & combinde_df$avg_log2FC_d2 < 0,])
 same_direction_significant <- same_direction_significant[same_direction_significant$p_val_adj_d1 < target_fdr & same_direction_significant$p_val_d2 > exclude_pval,]
 dim(same_direction_significant) # 59  7
 # get opposite direction DEGs, they just need to be "fdr < 0.05" in dataset1, no need to apply "p > 0.5" in dataset2 as they are in opposite direction. But following the above rationale - specified by Enrico 
 opposite_direction_significant <- setdiff(combinde_df, same_direction_significant)
# opposite_direction_significant <- opposite_direction_significant[opposite_direction_significant$p_val_adj_d1 < target_fdr & opposite_direction_significant$p_val_d2 > exclude_pval,]
  opposite_direction_significant <- opposite_direction_significant[opposite_direction_significant$p_val_adj_d1 < target_fdr,]
 dim(opposite_direction_significant) # 43  7
 significant_genes <- unique(c(same_direction_significant$gene, opposite_direction_significant$gene))
 return(significant_genes)
}

# Determining gender-specific and gender-dimorphic DEGs using adjusted significance and a nominal p-value specificity filter
# a minimum absolute logFC threshold for gender-dimorphic genes + min. abs. logFC for the target gender for gender-specific genes
# (still consider as gender-shared if significant in both genders with shared logFC, and abs logFC only above 0.25 in one gender)
gender_spec_genes = function(oligodendrocyte_M, oligodendrocyte_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50){

	# filter by target gene FDR
	male_canddeg = rownames(oligodendrocyte_M)[which(oligodendrocyte_M$p_val_adj < target_fdr)]
		
	# filter other gender non-significance (and not close to significance) nominal p-value threshold
	male_spec_genes = male_canddeg[which(oligodendrocyte_F[match(male_canddeg, rownames(oligodendrocyte_F)),]$p_val > exclude_pval)]
	
	# filter by abs. logFc threshold in target gender
	male_spec_genes = male_spec_genes[which(abs(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC) >minabslog)]
	
	print("Male-specific genes:")
	print(male_spec_genes)	
	female_canddeg = rownames(oligodendrocyte_F)[which(oligodendrocyte_F$p_val_adj < target_fdr)]
	female_spec_genes = female_canddeg[which(oligodendrocyte_M[match(female_canddeg, rownames(oligodendrocyte_M)),]$p_val > exclude_pval)]
	
	# filter by abs. logFc threshold in target gender
	female_spec_genes = female_spec_genes[which(abs(oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC) >minabslog)]	
	
	print("Female-specific genes:")
	print(female_spec_genes)
	
	# shared DEGs == gender-dimorphic	
	dimorphic_genes = NULL
	intdegs = intersect(male_canddeg, female_canddeg)
	
	if(length(intdegs)){
		# different sign of the fold-change
		logfcs_male = oligodendrocyte_M[match(intdegs, rownames(oligodendrocyte_M)),]$avg_log2FC
		logfcs_female = oligodendrocyte_F[match(intdegs, rownames(oligodendrocyte_F)),]$avg_log2FC	
		
		diff_fc = intersect(which(sign(logfcs_male) != sign(logfcs_female)), intersect(which(abs(logfcs_male)>minabslog), which(abs(logfcs_female)>minabslog)))
		dimorphic_genes = intdegs[diff_fc]	
		
		shared_fc = which(sign(logfcs_male) == sign(logfcs_female))
		# optionally add: minabslog fulfilled in at least one of the genders
		shared_genes = intdegs[shared_fc]		
	}
	
	print("Gender-dimorphic genes:")
	print(dimorphic_genes)
	
	print("Gender-shared genes:")
	print(shared_genes)
			
	dfres = data.frame("DEG type"=c(rep("male-specific",length(male_spec_genes)), rep("female-specific",length(female_spec_genes)), rep("gender-dimorphic",length(dimorphic_genes)), rep("gender-shared",length(shared_genes))), "Gene symbols"=c(male_spec_genes, female_spec_genes, dimorphic_genes, shared_genes), "Male avg. logFC"=c(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(female_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(dimorphic_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(shared_genes, rownames(oligodendrocyte_M)),]$avg_log2FC), "Female avg. logFC"=c(oligodendrocyte_F[match(male_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(dimorphic_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(shared_genes, rownames(oligodendrocyte_F)),]$avg_log2FC), "Male FDR"=c(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(female_spec_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(dimorphic_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(shared_genes, rownames(oligodendrocyte_M)),]$p_val_adj), "Female FDR"=c(oligodendrocyte_F[match(male_spec_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(dimorphic_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(shared_genes, rownames(oligodendrocyte_F)),]$p_val_adj))
	
	return(dfres)
	
}


## Microglial Cells

microglial_cell_TG_F <- FindMarkers(data_combined, ident.1 = "Microglial_cell_TG_F_M17", ident.2 = "Microglial_cell_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_TG_F) # 7716    5
nrow(microglial_cell_TG_F[microglial_cell_TG_F$p_val_adj < 0.05,]) # 1283
nrow(microglial_cell_TG_F[microglial_cell_TG_F$p_val < 0.05,]) # 4359
microglial_cell_TG_F$gene <- row.names(microglial_cell_TG_F)

microglial_cell_TG_M <- FindMarkers(data_combined, ident.1 = "Microglial_cell_TG_M_M17", ident.2 = "Microglial_cell_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_TG_M) # 7575    5
nrow(microglial_cell_TG_M[microglial_cell_TG_M$p_val_adj < 0.05,]) # 583
nrow(microglial_cell_TG_M[microglial_cell_TG_M$p_val < 0.05,]) # 3004
microglial_cell_TG_M$gene <- row.names(microglial_cell_TG_M)

microglial_cell_WT_F <- FindMarkers(data_combined, ident.1 = "Microglial_cell_WT_F_M17", ident.2 = "Microglial_cell_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_WT_F) # 7469    5
nrow(microglial_cell_WT_F[microglial_cell_WT_F$p_val_adj < 0.05,]) # 411
nrow(microglial_cell_WT_F[microglial_cell_WT_F$p_val < 0.05,]) # 2731
microglial_cell_WT_F$gene <- row.names(microglial_cell_WT_F)

microglial_cell_WT_M <- FindMarkers(data_combined, ident.1 = "Microglial_cell_WT_M_M17", ident.2 = "Microglial_cell_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_WT_M) # 7351    5
nrow(microglial_cell_WT_M[microglial_cell_WT_M$p_val_adj < 0.05,]) # 677
nrow(microglial_cell_WT_M[microglial_cell_WT_M$p_val < 0.05,]) # 3573
microglial_cell_WT_M$gene <- row.names(microglial_cell_WT_M)

save(microglial_cell_TG_F, microglial_cell_TG_M, microglial_cell_WT_F, microglial_cell_WT_M, file="Microglial_cell.RData")

microglial_F_DEGs <- unique_DEGs(microglial_cell_TG_F, microglial_cell_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(microglial_F_DEGs) # 1272
microglial_F_DEGs_Stats <- microglial_cell_TG_F[microglial_F_DEGs,]
dim(microglial_F_DEGs_Stats) # 1272    6

microglial_M_DEGs <- unique_DEGs(microglial_cell_TG_M, microglial_cell_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(microglial_M_DEGs) # 570
microglial_M_DEGs_Stats <- microglial_cell_TG_M[microglial_M_DEGs,]
dim(microglial_M_DEGs_Stats) # 570  6

table(row.names(microglial_M_DEGs_Stats) %in% row.names(microglial_F_DEGs_Stats))
# FALSE  TRUE
#   226   344

# identify gender-specific DEGs 
microglial_gender_spec_genes <- gender_spec_genes(microglial_M_DEGs_Stats, microglial_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(microglial_gender_spec_genes) # 344   6
table(microglial_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#               45              299


## Astrocyte

astrocyte_TG_F <- FindMarkers(data_combined, ident.1 = "Astrocyte_TG_F_M17", ident.2 = "Astrocyte_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_TG_F) # 7912    5
nrow(astrocyte_TG_F[astrocyte_TG_F$p_val_adj < 0.05,]) # 791
nrow(astrocyte_TG_F[astrocyte_TG_F$p_val < 0.05,]) # 3749
astrocyte_TG_F$gene <- row.names(astrocyte_TG_F)

astrocyte_TG_M <- FindMarkers(data_combined, ident.1 = "Astrocyte_TG_M_M17", ident.2 = "Astrocyte_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_TG_M) # 8038    5
nrow(astrocyte_TG_M[astrocyte_TG_M$p_val_adj < 0.05,]) # 219
nrow(astrocyte_TG_M[astrocyte_TG_M$p_val < 0.05,]) # 2197
astrocyte_TG_M$gene <- row.names(astrocyte_TG_M)

astrocyte_WT_F <- FindMarkers(data_combined, ident.1 = "Astrocyte_WT_F_M17", ident.2 = "Astrocyte_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_WT_F) # 8015    5
nrow(astrocyte_WT_F[astrocyte_WT_F$p_val_adj < 0.05,]) # 149
nrow(astrocyte_WT_F[astrocyte_WT_F$p_val < 0.05,]) # 1861
astrocyte_WT_F$gene <- row.names(astrocyte_WT_F)

astrocyte_WT_M <- FindMarkers(data_combined, ident.1 = "Astrocyte_WT_M_M17", ident.2 = "Astrocyte_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_WT_M) # 7843    5
nrow(astrocyte_WT_M[astrocyte_WT_M$p_val_adj < 0.05,]) # 321
nrow(astrocyte_WT_M[astrocyte_WT_M$p_val < 0.05,]) # 2627
astrocyte_WT_M$gene <- row.names(astrocyte_WT_M)

save(astrocyte_TG_F, astrocyte_TG_M, astrocyte_WT_F, astrocyte_WT_M, file="Astrocyte.RData")

astrocyte_F_DEGs <- unique_DEGs(astrocyte_TG_F, astrocyte_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(astrocyte_F_DEGs) # 791
astrocyte_F_DEGs_Stats <- astrocyte_TG_F[astrocyte_F_DEGs,]
dim(astrocyte_F_DEGs_Stats) # 791    6

astrocyte_M_DEGs <- unique_DEGs(astrocyte_TG_M, astrocyte_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(astrocyte_M_DEGs) # 219
astrocyte_M_DEGs_Stats <- astrocyte_TG_M[astrocyte_M_DEGs,]
dim(astrocyte_M_DEGs_Stats) # 219    6

table(row.names(astrocyte_M_DEGs_Stats) %in% row.names(astrocyte_F_DEGs_Stats))
# FALSE  TRUE
#    99   120

# identify gender-specific DEGs 
astrocyte_gender_spec_genes <- gender_spec_genes(astrocyte_M_DEGs_Stats, astrocyte_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(astrocyte_gender_spec_genes) # 120   6
table(astrocyte_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#               28               92


## Oligodendrocyte

oligodendrocyte_TG_F <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_TG_F_M17", ident.2 = "Oligodendrocyte_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_TG_F) # 7689    5
nrow(oligodendrocyte_TG_F[oligodendrocyte_TG_F$p_val_adj < 0.05,]) # 665
nrow(oligodendrocyte_TG_F[oligodendrocyte_TG_F$p_val < 0.05,]) # 3399
oligodendrocyte_TG_F$gene <- row.names(oligodendrocyte_TG_F)

oligodendrocyte_TG_M <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_TG_M_M17", ident.2 = "Oligodendrocyte_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_TG_M) # 7590    5
nrow(oligodendrocyte_TG_M[oligodendrocyte_TG_M$p_val_adj < 0.05,]) # 367
nrow(oligodendrocyte_TG_M[oligodendrocyte_TG_M$p_val < 0.05,]) # 2610
oligodendrocyte_TG_M$gene <- row.names(oligodendrocyte_TG_M)

oligodendrocyte_WT_F <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_WT_F_M17", ident.2 = "Oligodendrocyte_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_WT_F) # 7705    5
nrow(oligodendrocyte_WT_F[oligodendrocyte_WT_F$p_val_adj < 0.05,]) # 148
nrow(oligodendrocyte_WT_F[oligodendrocyte_WT_F$p_val < 0.05,]) # 1876
oligodendrocyte_WT_F$gene <- row.names(oligodendrocyte_WT_F)

oligodendrocyte_WT_M <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_WT_M_M17", ident.2 = "Oligodendrocyte_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_WT_M) # 7301    5
nrow(oligodendrocyte_WT_M[oligodendrocyte_WT_M$p_val_adj < 0.05,]) # 226
nrow(oligodendrocyte_WT_M[oligodendrocyte_WT_M$p_val < 0.05,]) # 2229
oligodendrocyte_WT_M$gene <- row.names(oligodendrocyte_WT_M)

save(oligodendrocyte_TG_F, oligodendrocyte_TG_M, oligodendrocyte_WT_F, oligodendrocyte_WT_M, file="Oligodendrocyte.RData")

oligodendrocyte_F_DEGs <- unique_DEGs(oligodendrocyte_TG_F, oligodendrocyte_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(oligodendrocyte_F_DEGs) # 665
oligodendrocyte_F_DEGs_Stats <- oligodendrocyte_TG_F[oligodendrocyte_F_DEGs,]
dim(oligodendrocyte_F_DEGs_Stats) # 665    6

oligodendrocyte_M_DEGs <- unique_DEGs(oligodendrocyte_TG_M, oligodendrocyte_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(oligodendrocyte_M_DEGs) # 367
oligodendrocyte_M_DEGs_Stats <- oligodendrocyte_TG_M[oligodendrocyte_M_DEGs,]
dim(oligodendrocyte_M_DEGs_Stats) # 367    6

table(row.names(oligodendrocyte_M_DEGs_Stats) %in% row.names(oligodendrocyte_F_DEGs_Stats))
# FALSE  TRUE
#   226   141

# identify gender-specific DEGs 
oligodendrocyte_gender_spec_genes <- gender_spec_genes(oligodendrocyte_M_DEGs_Stats, oligodendrocyte_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(oligodendrocyte_gender_spec_genes) # 141   6
table(oligodendrocyte_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#               24              117


## Endothelial cells

endothelial_cell_TG_F <- FindMarkers(data_combined, ident.1 = "Endothelial_cell_TG_F_M17", ident.2 = "Endothelial_cell_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_TG_F) # 8197    5
nrow(endothelial_cell_TG_F[endothelial_cell_TG_F$p_val_adj < 0.05,]) # 753
nrow(endothelial_cell_TG_F[endothelial_cell_TG_F$p_val < 0.05,]) # 3638
endothelial_cell_TG_F$gene <- row.names(endothelial_cell_TG_F)

endothelial_cell_TG_M <- FindMarkers(data_combined, ident.1 = "Endothelial_cell_TG_M_M17", ident.2 = "Endothelial_cell_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_TG_M) # 8018    5
nrow(endothelial_cell_TG_M[endothelial_cell_TG_M$p_val_adj < 0.05,]) # 242
nrow(endothelial_cell_TG_M[endothelial_cell_TG_M$p_val < 0.05,]) # 2331
endothelial_cell_TG_M$gene <- row.names(endothelial_cell_TG_M)

endothelial_cell_WT_F <- FindMarkers(data_combined, ident.1 = "Endothelial_cell_WT_F_M17", ident.2 = "Endothelial_cell_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_WT_F) # 7893    5
nrow(endothelial_cell_WT_F[endothelial_cell_WT_F$p_val_adj < 0.05,]) # 219
nrow(endothelial_cell_WT_F[endothelial_cell_WT_F$p_val < 0.05,]) # 2022
endothelial_cell_WT_F$gene <- row.names(endothelial_cell_WT_F)

endothelial_cell_WT_M <- FindMarkers(data_combined, ident.1 = "Endothelial_cell_WT_M_M17", ident.2 = "Endothelial_cell_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_WT_M) # 7709    5
nrow(endothelial_cell_WT_M[endothelial_cell_WT_M$p_val_adj < 0.05,]) # 357
nrow(endothelial_cell_WT_M[endothelial_cell_WT_M$p_val < 0.05,]) # 2717
endothelial_cell_WT_M$gene <- row.names(endothelial_cell_WT_M)

save(endothelial_cell_TG_F, endothelial_cell_TG_M, endothelial_cell_WT_F, endothelial_cell_WT_M, file="Endothelial_cell.RData")

endothelial_cell_F_DEGs <- unique_DEGs(endothelial_cell_TG_F, endothelial_cell_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(endothelial_cell_F_DEGs) # 751
endothelial_cell_F_DEGs_Stats <- endothelial_cell_TG_F[endothelial_cell_F_DEGs,]
dim(endothelial_cell_F_DEGs_Stats) # 751    6

endothelial_cell_M_DEGs <- unique_DEGs(endothelial_cell_TG_M, endothelial_cell_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(endothelial_cell_M_DEGs) # 241
endothelial_cell_M_DEGs_Stats <- endothelial_cell_TG_M[endothelial_cell_M_DEGs,]
dim(endothelial_cell_M_DEGs_Stats) # 241    6

table(row.names(endothelial_cell_M_DEGs_Stats) %in% row.names(endothelial_cell_F_DEGs_Stats))
# FALSE  TRUE
#   112   129

# identify gender-specific DEGs 
endothelial_cell_gender_spec_genes <- gender_spec_genes(endothelial_cell_M_DEGs_Stats, endothelial_cell_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(endothelial_cell_gender_spec_genes) # 129   6
table(endothelial_cell_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#               17              112


## Mural cell

mural_cell_TG_F <- FindMarkers(data_combined, ident.1 = "Mural_cell_TG_F_M17", ident.2 = "Mural_cell_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_TG_F) # 4542    5
nrow(mural_cell_TG_F[mural_cell_TG_F$p_val_adj < 0.05,]) # 19
nrow(mural_cell_TG_F[mural_cell_TG_F$p_val < 0.05,]) # 290
mural_cell_TG_F$gene <- row.names(mural_cell_TG_F)

mural_cell_TG_M <- FindMarkers(data_combined, ident.1 = "Mural_cell_TG_M_M17", ident.2 = "Mural_cell_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_TG_M) # 4728    5
nrow(mural_cell_TG_M[mural_cell_TG_M$p_val_adj < 0.05,]) # 8
nrow(mural_cell_TG_M[mural_cell_TG_M$p_val < 0.05,]) # 232
mural_cell_TG_M$gene <- row.names(mural_cell_TG_M)

mural_cell_WT_F <- FindMarkers(data_combined, ident.1 = "Mural_cell_WT_F_M17", ident.2 = "Mural_cell_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_WT_F) # 3845    5
nrow(mural_cell_WT_F[mural_cell_WT_F$p_val_adj < 0.05,]) # 7
nrow(mural_cell_WT_F[mural_cell_WT_F$p_val < 0.05,]) # 163
mural_cell_WT_F$gene <- row.names(mural_cell_WT_F)

mural_cell_WT_M <- FindMarkers(data_combined, ident.1 = "Mural_cell_WT_M_M17", ident.2 = "Mural_cell_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_WT_M) # 4882    5
nrow(mural_cell_WT_M[mural_cell_WT_M$p_val_adj < 0.05,]) # 10
nrow(mural_cell_WT_M[mural_cell_WT_M$p_val < 0.05,]) # 182
mural_cell_WT_M$gene <- row.names(mural_cell_WT_M)

save(mural_cell_TG_F, mural_cell_TG_M, mural_cell_WT_F, mural_cell_WT_M, file="Mural_cell.RData")

mural_cell_F_DEGs <- unique_DEGs(mural_cell_TG_F, mural_cell_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(mural_cell_F_DEGs) # 19
mural_cell_F_DEGs_Stats <- mural_cell_TG_F[mural_cell_F_DEGs,]
dim(mural_cell_F_DEGs_Stats) # 19    6

mural_cell_M_DEGs <- unique_DEGs(mural_cell_TG_M, mural_cell_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(mural_cell_M_DEGs) # 8
mural_cell_M_DEGs_Stats <- mural_cell_TG_M[mural_cell_M_DEGs,]
dim(mural_cell_M_DEGs_Stats) # 8    6

table(row.names(mural_cell_M_DEGs_Stats) %in% row.names(mural_cell_F_DEGs_Stats))
# FALSE  TRUE
#     1     7

# identify gender-specific DEGs 
mural_cell_gender_spec_genes <- gender_spec_genes(mural_cell_M_DEGs_Stats, mural_cell_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(mural_cell_gender_spec_genes) # 7 6
table(mural_cell_gender_spec_genes$DEG.type)
# gender-shared
#             7


## Ependymal cell

ependymal_cell_TG_F <- FindMarkers(data_combined, ident.1 = "Ependymal_cell_TG_F_M17", ident.2 = "Ependymal_cell_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_TG_F) # 4827    5
nrow(ependymal_cell_TG_F[ependymal_cell_TG_F$p_val_adj < 0.05,]) # 9
nrow(ependymal_cell_TG_F[ependymal_cell_TG_F$p_val < 0.05,]) # 165
ependymal_cell_TG_F$gene <- row.names(ependymal_cell_TG_F)

ependymal_cell_TG_M <- FindMarkers(data_combined, ident.1 = "Ependymal_cell_TG_M_M17", ident.2 = "Ependymal_cell_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_TG_M) # 6185    5
nrow(ependymal_cell_TG_M[ependymal_cell_TG_M$p_val_adj < 0.05,]) # 2
nrow(ependymal_cell_TG_M[ependymal_cell_TG_M$p_val < 0.05,]) # 174
ependymal_cell_TG_M$gene <- row.names(ependymal_cell_TG_M)

ependymal_cell_WT_F <- FindMarkers(data_combined, ident.1 = "Ependymal_cell_WT_F_M17", ident.2 = "Ependymal_cell_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_WT_F) # 5306    5
nrow(ependymal_cell_WT_F[ependymal_cell_WT_F$p_val_adj < 0.05,]) # 3
nrow(ependymal_cell_WT_F[ependymal_cell_WT_F$p_val < 0.05,]) # 176
ependymal_cell_WT_F$gene <- row.names(ependymal_cell_WT_F)

ependymal_cell_WT_M <- FindMarkers(data_combined, ident.1 = "Ependymal_cell_WT_M_M17", ident.2 = "Ependymal_cell_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_WT_M) # 3972    5
nrow(ependymal_cell_WT_M[ependymal_cell_WT_M$p_val_adj < 0.05,]) # 5
nrow(ependymal_cell_WT_M[ependymal_cell_WT_M$p_val < 0.05,]) # 160
ependymal_cell_WT_M$gene <- row.names(ependymal_cell_WT_M)

save(ependymal_cell_TG_F, ependymal_cell_TG_M, ependymal_cell_WT_F, ependymal_cell_WT_M, file="Ependymal_cell.RData")

ependymal_cell_F_DEGs <- unique_DEGs(ependymal_cell_TG_F, ependymal_cell_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(ependymal_cell_F_DEGs) # 9
ependymal_cell_F_DEGs_Stats <- ependymal_cell_TG_F[ependymal_cell_F_DEGs,]
dim(ependymal_cell_F_DEGs_Stats) # 9    6

ependymal_cell_M_DEGs <- unique_DEGs(ependymal_cell_TG_M, ependymal_cell_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(ependymal_cell_M_DEGs) # 2
ependymal_cell_M_DEGs_Stats <- ependymal_cell_TG_M[ependymal_cell_M_DEGs,]
dim(ependymal_cell_M_DEGs_Stats) # 2    6

table(row.names(ependymal_cell_M_DEGs_Stats) %in% row.names(ependymal_cell_F_DEGs_Stats))
# FALSE  TRUE
#     1     1

# identify gender-specific DEGs 
ependymal_cell_gender_spec_genes <- gender_spec_genes(ependymal_cell_M_DEGs_Stats, ependymal_cell_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(ependymal_cell_gender_spec_genes) # 1 6
table(ependymal_cell_gender_spec_genes$DEG.type)
# gender-shared
#             1


## OPCs

oligodendrocyte_precursor_cell_TG_F <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_precursor_cell_TG_F_M17", ident.2 = "Oligodendrocyte_precursor_cell_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_TG_F) # 10030     5
nrow(oligodendrocyte_precursor_cell_TG_F[oligodendrocyte_precursor_cell_TG_F$p_val_adj < 0.05,]) # 605
nrow(oligodendrocyte_precursor_cell_TG_F[oligodendrocyte_precursor_cell_TG_F$p_val < 0.05,]) # 3664
oligodendrocyte_precursor_cell_TG_F$gene <- row.names(oligodendrocyte_precursor_cell_TG_F)

oligodendrocyte_precursor_cell_TG_M <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_precursor_cell_TG_M_M17", ident.2 = "Oligodendrocyte_precursor_cell_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_TG_M) # 9641    5
nrow(oligodendrocyte_precursor_cell_TG_M[oligodendrocyte_precursor_cell_TG_M$p_val_adj < 0.05,]) # 324
nrow(oligodendrocyte_precursor_cell_TG_M[oligodendrocyte_precursor_cell_TG_M$p_val < 0.05,]) # 2171
oligodendrocyte_precursor_cell_TG_M$gene <- row.names(oligodendrocyte_precursor_cell_TG_M)

oligodendrocyte_precursor_cell_WT_F <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_precursor_cell_WT_F_M17", ident.2 = "Oligodendrocyte_precursor_cell_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_WT_F) # 7688    5
nrow(oligodendrocyte_precursor_cell_WT_F[oligodendrocyte_precursor_cell_WT_F$p_val_adj < 0.05,]) # 170
nrow(oligodendrocyte_precursor_cell_WT_F[oligodendrocyte_precursor_cell_WT_F$p_val < 0.05,]) # 1222
oligodendrocyte_precursor_cell_WT_F$gene <- row.names(oligodendrocyte_precursor_cell_WT_F)

oligodendrocyte_precursor_cell_WT_M <- FindMarkers(data_combined, ident.1 = "Oligodendrocyte_precursor_cell_WT_M_M17", ident.2 = "Oligodendrocyte_precursor_cell_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_WT_M) # 8437    5
nrow(oligodendrocyte_precursor_cell_WT_M[oligodendrocyte_precursor_cell_WT_M$p_val_adj < 0.05,]) # 195
nrow(oligodendrocyte_precursor_cell_WT_M[oligodendrocyte_precursor_cell_WT_M$p_val < 0.05,]) # 1250
oligodendrocyte_precursor_cell_WT_M$gene <- row.names(oligodendrocyte_precursor_cell_WT_M)

save(oligodendrocyte_precursor_cell_TG_F, oligodendrocyte_precursor_cell_TG_M, oligodendrocyte_precursor_cell_WT_F, oligodendrocyte_precursor_cell_WT_M, file="Oligodendrocyte_precursor_cell.RData")

oligodendrocyte_precursor_cell_F_DEGs <- unique_DEGs(oligodendrocyte_precursor_cell_TG_F, oligodendrocyte_precursor_cell_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(oligodendrocyte_precursor_cell_F_DEGs) # 602
oligodendrocyte_precursor_cell_F_DEGs_Stats <- oligodendrocyte_precursor_cell_TG_F[oligodendrocyte_precursor_cell_F_DEGs,]
dim(oligodendrocyte_precursor_cell_F_DEGs_Stats) # 602    6

oligodendrocyte_precursor_cell_M_DEGs <- unique_DEGs(oligodendrocyte_precursor_cell_TG_M, oligodendrocyte_precursor_cell_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(oligodendrocyte_precursor_cell_M_DEGs) # 324
oligodendrocyte_precursor_cell_M_DEGs_Stats <- oligodendrocyte_precursor_cell_TG_M[oligodendrocyte_precursor_cell_M_DEGs,]
dim(oligodendrocyte_precursor_cell_M_DEGs_Stats) # 324    6

table(row.names(oligodendrocyte_precursor_cell_M_DEGs_Stats) %in% row.names(oligodendrocyte_precursor_cell_F_DEGs_Stats))
# FALSE  TRUE
#   119   205

# identify gender-specific DEGs 
oligodendrocyte_precursor_cell_gender_spec_genes <- gender_spec_genes(oligodendrocyte_precursor_cell_M_DEGs_Stats, oligodendrocyte_precursor_cell_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(oligodendrocyte_precursor_cell_gender_spec_genes) # 205   6
table(oligodendrocyte_precursor_cell_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#                4              201


## Neuron

neuron_TG_F <- FindMarkers(data_combined, ident.1 = "Neuron_TG_F_M17", ident.2 = "Neuron_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_TG_F) # 9502    5
nrow(neuron_TG_F[neuron_TG_F$p_val_adj < 0.05,]) # 171
nrow(neuron_TG_F[neuron_TG_F$p_val < 0.05,]) # 2045
neuron_TG_F$gene <- row.names(neuron_TG_F)

neuron_TG_M <- FindMarkers(data_combined, ident.1 = "Neuron_TG_M_M17", ident.2 = "Neuron_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_TG_M) # 10138     5
nrow(neuron_TG_M[neuron_TG_M$p_val_adj < 0.05,]) # 225
nrow(neuron_TG_M[neuron_TG_M$p_val < 0.05,]) # 2488
neuron_TG_M$gene <- row.names(neuron_TG_M)

neuron_WT_F <- FindMarkers(data_combined, ident.1 = "Neuron_WT_F_M17", ident.2 = "Neuron_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_WT_F) # 9514    5
nrow(neuron_WT_F[neuron_WT_F$p_val_adj < 0.05,]) # 123
nrow(neuron_WT_F[neuron_WT_F$p_val < 0.05,]) # 1523
neuron_WT_F$gene <- row.names(neuron_WT_F)

neuron_WT_M <- FindMarkers(data_combined, ident.1 = "Neuron_WT_M_M17", ident.2 = "Neuron_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_WT_M) # 9730    5
nrow(neuron_WT_M[neuron_WT_M$p_val_adj < 0.05,]) # 155
nrow(neuron_WT_M[neuron_WT_M$p_val < 0.05,]) # 2074
neuron_WT_M$gene <- row.names(neuron_WT_M)

save(neuron_TG_F, neuron_TG_M, neuron_WT_F, neuron_WT_M, file="Neuron.RData")

neuron_F_DEGs <- unique_DEGs(neuron_TG_F, neuron_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(neuron_F_DEGs) # 171
neuron_F_DEGs_Stats <- neuron_TG_F[neuron_F_DEGs,]
dim(neuron_F_DEGs_Stats) # 171    6

neuron_M_DEGs <- unique_DEGs(neuron_TG_M, neuron_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(neuron_M_DEGs) # 224
neuron_M_DEGs_Stats <- neuron_TG_M[neuron_M_DEGs,]
dim(neuron_M_DEGs_Stats) # 224    6

table(row.names(neuron_M_DEGs_Stats) %in% row.names(neuron_F_DEGs_Stats))
# FALSE  TRUE
#   126    98

# identify gender-specific DEGs 
neuron_gender_spec_genes <- gender_spec_genes(neuron_M_DEGs_Stats, neuron_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(neuron_gender_spec_genes) # 98  6
table(neuron_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#                2               96


## Macrophage

macrophage_TG_F <- FindMarkers(data_combined, ident.1 = "Macrophage_TG_F_M17", ident.2 = "Macrophage_TG_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_TG_F) # 5615    5
nrow(macrophage_TG_F[macrophage_TG_F$p_val_adj < 0.05,]) # 35
nrow(macrophage_TG_F[macrophage_TG_F$p_val < 0.05,]) # 505
macrophage_TG_F$gene <- row.names(macrophage_TG_F)

macrophage_TG_M <- FindMarkers(data_combined, ident.1 = "Macrophage_TG_M_M17", ident.2 = "Macrophage_TG_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_TG_M) # 5669    5
nrow(macrophage_TG_M[macrophage_TG_M$p_val_adj < 0.05,]) # 34
nrow(macrophage_TG_M[macrophage_TG_M$p_val < 0.05,]) # 434
macrophage_TG_M$gene <- row.names(macrophage_TG_M)

macrophage_WT_F <- FindMarkers(data_combined, ident.1 = "Macrophage_WT_F_M17", ident.2 = "Macrophage_WT_F_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_WT_F) # 4937    5
nrow(macrophage_WT_F[macrophage_WT_F$p_val_adj < 0.05,]) # 13
nrow(macrophage_WT_F[macrophage_WT_F$p_val < 0.05,]) # 281
macrophage_WT_F$gene <- row.names(macrophage_WT_F)

macrophage_WT_M <- FindMarkers(data_combined, ident.1 = "Macrophage_WT_M_M17", ident.2 = "Macrophage_WT_M_M7", verbose = TRUE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_WT_M) # 4217    5
nrow(macrophage_WT_M[macrophage_WT_M$p_val_adj < 0.05,]) # 13
nrow(macrophage_WT_M[macrophage_WT_M$p_val < 0.05,]) # 240
macrophage_WT_M$gene <- row.names(macrophage_WT_M)

save(macrophage_TG_F, macrophage_TG_M, macrophage_WT_F, macrophage_WT_M, file="Macrophage.RData")

macrophage_F_DEGs <- unique_DEGs(macrophage_TG_F, macrophage_WT_F, target_fdr = 0.05, exclude_pval = 0.5)
length(macrophage_F_DEGs) # 35
macrophage_F_DEGs_Stats <- macrophage_TG_F[macrophage_F_DEGs,]
dim(macrophage_F_DEGs_Stats) # 35    6

macrophage_M_DEGs <- unique_DEGs(macrophage_TG_M, macrophage_WT_M, target_fdr = 0.05, exclude_pval = 0.5)
length(macrophage_M_DEGs) # 34
macrophage_M_DEGs_Stats <- macrophage_TG_M[macrophage_M_DEGs,]
dim(macrophage_M_DEGs_Stats) # 34    6

table(row.names(macrophage_M_DEGs_Stats) %in% row.names(macrophage_F_DEGs_Stats))
# FALSE  TRUE
#    22    12

# identify gender-specific DEGs 
macrophage_gender_spec_genes <- gender_spec_genes(macrophage_M_DEGs_Stats, macrophage_F_DEGs_Stats, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
dim(macrophage_gender_spec_genes) # 12  6
table(macrophage_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared
#                3                9

microglial_gender_spec_genes$Cell_Type <- "Microglial cells"
astrocyte_gender_spec_genes$Cell_Type <- "Astrocytes"
oligodendrocyte_gender_spec_genes$Cell_Type <- "Oligodendrocytes"
endothelial_cell_gender_spec_genes$Cell_Type <- "Endothelial cells"
mural_cell_gender_spec_genes$Cell_Type <- "Mural cells"
ependymal_cell_gender_spec_genes$Cell_Type <- "Ependymal cells"
oligodendrocyte_precursor_cell_gender_spec_genes$Cell_Type <- "OPCs"
neuron_gender_spec_genes$Cell_Type <- "Neurons"
macrophage_gender_spec_genes$Cell_Type <- "Macrophages"

longitudinal_DEGs <- rbind(microglial_gender_spec_genes, astrocyte_gender_spec_genes, oligodendrocyte_gender_spec_genes, endothelial_cell_gender_spec_genes,
	mural_cell_gender_spec_genes, ependymal_cell_gender_spec_genes, oligodendrocyte_precursor_cell_gender_spec_genes,
	neuron_gender_spec_genes, macrophage_gender_spec_genes)
dim(longitudinal_DEGs) # 1057    7
head(sort(unique(longitudinal_DEGs$Male.FDR)))
# [1]  0.000000e+00 1.580394e-259 1.841732e-252 9.277304e-245 9.850767e-245 1.544446e-236
longitudinal_DEGs$Male.FDR <- ifelse(longitudinal_DEGs$Male.FDR == 0, 1.580394e-259, longitudinal_DEGs$Male.FDR)
head(sort(unique(longitudinal_DEGs$Female.FDR)))
# [1] 0.000000e+00 2.109539e-280 2.741659e-278 5.067706e-260 1.076462e-229 4.782702e-224
longitudinal_DEGs$Female.FDR <- ifelse(longitudinal_DEGs$Female.FDR == 0, 2.109539e-280, longitudinal_DEGs$Female.FDR)
write.table(longitudinal_DEGs, file="ThyTau22_Longitudinal_DEGs.txt", sep="\t", row.names=F, quote=F)
# write.csv(longitudinal_DEGs, file="ThyTau22_Longitudinal_DEGs.csv", row.names=F, quote=F)

as.data.frame(table(longitudinal_DEGs$DEG.type, longitudinal_DEGs$Cell_Type))
#               Var1              Var2 Freq
# 1  gender-dimorphic        Astrocytes   28
# 2     gender-shared        Astrocytes   92
# 3  gender-dimorphic Endothelial cells   17
# 4     gender-shared Endothelial cells  112
# 5  gender-dimorphic   Ependymal cells    0
# 6     gender-shared   Ependymal cells    1
# 7  gender-dimorphic       Macrophages    3
# 8     gender-shared       Macrophages    9
# 9  gender-dimorphic  Microglial cells   45
# 10    gender-shared  Microglial cells  299
# 11 gender-dimorphic       Mural cells    0
# 12    gender-shared       Mural cells    7
# 13 gender-dimorphic           Neurons    2
# 14    gender-shared           Neurons   96
# 15 gender-dimorphic  Oligodendrocytes   24
# 16    gender-shared  Oligodendrocytes  117
# 17 gender-dimorphic              OPCs    4
# 18    gender-shared              OPCs  201
length(unique(longitudinal_DEGs[longitudinal_DEGs$DEG.type == "gender-dimorphic",]$Gene.symbols)) # 89
length(unique(longitudinal_DEGs[longitudinal_DEGs$DEG.type == "gender-shared",]$Gene.symbols)) # 700

save(microglial_gender_spec_genes, astrocyte_gender_spec_genes, oligodendrocyte_gender_spec_genes, endothelial_cell_gender_spec_genes,
	mural_cell_gender_spec_genes, ependymal_cell_gender_spec_genes, oligodendrocyte_precursor_cell_gender_spec_genes,
	neuron_gender_spec_genes, macrophage_gender_spec_genes, file="ThyTau22_Longitudinal_SexSpecific_DEGs.RData")

# Violinplot for key genes
features <- c("Malat1", "Gm42418", "Aldoc", "Apoe", "Ddx3y")

png("Mbp_Oligodendrocyte.png", width=16, height=8, units="in", res=300)
p <- VlnPlot(data_combined, features = c("Mbp"), assay="RNA", log = TRUE, idents=c("Oligodendrocyte_TG_F_M17", "Oligodendrocyte_TG_F_M7", "Oligodendrocyte_TG_M_M17", "Oligodendrocyte_TG_M_M7", "Oligodendrocyte_WT_F_M17", "Oligodendrocyte_WT_F_M7", "Oligodendrocyte_WT_M_M17", "Oligodendrocyte_WT_M_M7") , pt.size = 0, combine = FALSE)
lapply(p, function(x){x + labs(title = "Mbp in Oligodendrocytes (sex-shared)")})
dev.off()

png("Malat1_Oligodendrocyte.png", width=16, height=8, units="in", res=300)
p <- VlnPlot(data_combined, features = c("Malat1"), assay="RNA", log = TRUE, idents=c("Oligodendrocyte_TG_F_M17", "Oligodendrocyte_TG_F_M7", "Oligodendrocyte_TG_M_M17", "Oligodendrocyte_TG_M_M7", "Oligodendrocyte_WT_F_M17", "Oligodendrocyte_WT_F_M7", "Oligodendrocyte_WT_M_M17", "Oligodendrocyte_WT_M_M7") , pt.size = 0, combine = FALSE)
lapply(p, function(x){x + labs(title = "Malat1 in Oligodendrocytes (sex-dimorphic)")})
dev.off()

png("Gm42418_Astrocytes.png", width=16, height=8, units="in", res=300)
p <- VlnPlot(data_combined, features = c("Gm42418"), assay="RNA", log = TRUE, idents=c("Astrocyte_TG_F_M17", "Astrocyte_TG_F_M7", "Astrocyte_TG_M_M17", "Astrocyte_TG_M_M7", "Astrocyte_WT_F_M17", "Astrocyte_WT_F_M7", "Astrocyte_WT_M_M17", "Astrocyte_WT_M_M7") , pt.size = 0, combine = FALSE)
lapply(p, function(x){x + labs(title = "Gm42418 in Astrocytes (sex-dimorphic)")})
dev.off()

png("C1qc_Astrocytes.png", width=16, height=8, units="in", res=300)
p <- VlnPlot(data_combined, features = c("C1qc"), assay="RNA", log = TRUE, idents=c("Astrocyte_TG_F_M17", "Astrocyte_TG_F_M7", "Astrocyte_TG_M_M17", "Astrocyte_TG_M_M7", "Astrocyte_WT_F_M17", "Astrocyte_WT_F_M7", "Astrocyte_WT_M_M17", "Astrocyte_WT_M_M7") , pt.size = 0, combine = FALSE)
lapply(p, function(x){x + labs(title = "C1qc in Astrocytes (sex-shared)")})
dev.off()
