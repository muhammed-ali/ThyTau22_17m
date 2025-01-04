##
## Global Differential expression analysis (DEA) ##
##

library(Seurat) # SeuratObject_4.1.4 Seurat_4.3.0
library(openxlsx) # openxlsx_4.2.5.2
library(ggplot2) # ggplot2_3.5.1
set.seed(1)

load("./data/sctrans2.RData")

# Set sample identities to celltype/condition combinations
orig.ident <- sctrans2$orig.ident
sctrans2$orig.ident <- paste(sctrans2$classint, sctrans2$conditions)
Idents(sctrans2) <- paste(sctrans2$classint, sctrans2$conditions)
table(sctrans2@meta.data$conditions)
# TG F TG M WT F WT M
# 8838 6988 4821 7543
# sctrans2$cond <- conditions_cells
sctrans2$cond <- sctrans2@meta.data$conditions
Idents(sctrans2) <- "cond"

# don't use the following parameters: min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 0, min.cells.group = 0, 
DEG_Global_Poisson_F <- FindMarkers(sctrans2, ident.1 = "TG F", ident.2 = "WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos = FALSE, logfc.threshold = -Inf)
nrow(DEG_Global_Poisson_F[DEG_Global_Poisson_F$p_val_adj < 0.05,]) # 61
nrow(DEG_Global_Poisson_F[DEG_Global_Poisson_F$p_val < 0.05,]) # 330
dim(DEG_Global_Poisson_F) # 705   5
write.csv(DEG_Global_Poisson_F, file="./results/files/DEG_Global_Poisson_Female.csv", quote=F)

DEG_Global_Poisson_M <- FindMarkers(sctrans2, ident.1 = "TG M", ident.2 = "WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos = FALSE, logfc.threshold = -Inf)
nrow(DEG_Global_Poisson_M[DEG_Global_Poisson_M$p_val_adj < 0.05,]) # 367
nrow(DEG_Global_Poisson_M[DEG_Global_Poisson_M$p_val < 0.05,]) # 605
dim(DEG_Global_Poisson_M) # 764   5
write.csv(DEG_Global_Poisson_M, file="./results/files/DEG_Global_Poisson_Male.csv", quote=F)

# Calculate and print proportion of genes with higher variation in females
prop_higher_in_females <- mean(female_sd > male_sd, na.rm=TRUE)
print(paste("Proportion of genes with higher SD in females:", prop_higher_in_females))


# Plot top DEGs expression dotplot
table(Idents(sctrans2))
# WT M TG F TG M WT F
# 7543 8838 6988 4821
myLevels <- c("WT F", "TG F", "WT M", "TG M")
Idents(sctrans2) <- factor(Idents(sctrans2), levels= myLevels)
# DEGs <- c("Hsp90aa1", "Mag", "Lamp1", "Slc25a5", "Cst3", "Mbp", "Apod", "Fth1", "Malat1", "Plp1")
DEGs <- c("Hsp90aa1", "Scd2", "Trf", "Zc3h13", "Cst3", "Mbp", "Apod", "Fth1", "Malat1", "Plp1")
# DotPlot(object = sctrans2, features = DEGs, split.by="conditions", cols="RdBu", dot.scale=10)+RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=12), axis.title = element_blank())+scale_y_discrete(labels=c('WT_Female', 'TG_Female', 'WT_Male', 'TG_Male'))
png("./results/figures/DotPlot_Top3_Global_DEGs.png", width=8, height=4, units="in", res=300)
DotPlot(object = sctrans2, features = DEGs, split.by="conditions", cols="RdBu", dot.scale=10)+theme(axis.text.x=element_text(size=14, angle = 60, vjust = 1, hjust=1), axis.text.y=element_text(size=14), axis.title = element_blank())+scale_y_discrete(labels=c('WT_Female', 'TG_Female', 'WT_Male', 'TG_Male'))
dev.off()


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

global_gender_spec_genes = gender_spec_genes(DEG_Global_Poisson_M, DEG_Global_Poisson_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
save(DEG_Global_Poisson_M, DEG_Global_Poisson_F, global_gender_spec_genes, file="./results/files/Global_DEGs.RData")
# write.table(global_gender_spec_genes, file="./results/files/global_gender_spec_genes.txt", sep="\t", row.names=F, quote=F)#
write.table(global_gender_spec_genes, file="./results/files/global_gender_spec_genes_Nov2024.txt", sep="\t", row.names=F, quote=F)

# save one version of min abs. logFC threshold = 0.50
head(sort(global_gender_spec_genes$Male.FDR))
# [1]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 2.498074e-274
head(sort(global_gender_spec_genes$Female.FDR))
# [1] 5.175763e-66 1.138992e-45 1.370410e-44 1.021988e-26 9.354951e-26 2.534532e-23
global_gender_spec_genes$Male.FDR <- ifelse(global_gender_spec_genes$Male.FDR == 0, 2.498074e-274, global_gender_spec_genes$Male.FDR)
dim(global_gender_spec_genes) # 47  6
write.table(global_gender_spec_genes, file="./results/files/global_gender_spec_genes_Nov2024.txt", sep="\t", row.names=F, quote=F)


##
## Cell-type specific DEA ##
##

#
# Microglial cells
#

Idents(sctrans2) <- paste(sctrans2$classint, sctrans2$conditions)

#microglial_cell_F_logfcfilt <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG F", ident.2 = "Microglial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE)
#dim(microglial_cell_F_logfcfilt) # 20  5

microglial_cell_F <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG F", ident.2 = "Microglial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_F) # 770   5
nrow(microglial_cell_F[microglial_cell_F$p_val_adj < 0.05,]) # 27
nrow(microglial_cell_F[microglial_cell_F$p_val < 0.05,]) # 177

microglial_cell_M <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG M", ident.2 = "Microglial cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_M) # 823   5
nrow(microglial_cell_M[microglial_cell_M$p_val_adj < 0.05,]) # 227
nrow(microglial_cell_M[microglial_cell_M$p_val < 0.05,]) # 513

microglia_gender_spec_genes = gender_spec_genes(microglial_cell_M, microglial_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(microglia_gender_spec_genes) # 30  6
table(microglia_gender_spec_genes$DEG.type)
#  female-specific gender-dimorphic    gender-shared    male-specific
#                1                6                9               14


#
# Oligodendrocytes
#

oligodendrocyte_cell_F <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte TG F", ident.2 = "Oligodendrocyte WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_cell_F) # 827   5
nrow(oligodendrocyte_cell_F[oligodendrocyte_cell_F$p_val_adj < 0.05,]) # 14
nrow(oligodendrocyte_cell_F[oligodendrocyte_cell_F$p_val < 0.05,]) # 189

oligodendrocyte_cell_M <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte TG M", ident.2 = "Oligodendrocyte WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_cell_M) # 845   5
nrow(oligodendrocyte_cell_M[oligodendrocyte_cell_M$p_val_adj < 0.05,]) # 124
nrow(oligodendrocyte_cell_M[oligodendrocyte_cell_M$p_val < 0.05,]) # 458

oligodendrocyte_gender_spec_genes = gender_spec_genes(oligodendrocyte_cell_M, oligodendrocyte_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(oligodendrocyte_gender_spec_genes) # 21  6
table(oligodendrocyte_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#                7                3               11


#
# Neuron
# 

neuron_cell_F <- FindMarkers(sctrans2, ident.1 = "Neuron TG F", ident.2 = "Neuron WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_cell_F) # 812   5
nrow(neuron_cell_F[neuron_cell_F$p_val_adj < 0.05,]) # 26
nrow(neuron_cell_F[neuron_cell_F$p_val < 0.05,]) # 277

neuron_cell_M <- FindMarkers(sctrans2, ident.1 = "Neuron TG M", ident.2 = "Neuron WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_cell_M) # 812   5
nrow(neuron_cell_M[neuron_cell_M$p_val_adj < 0.05,]) # 163
nrow(neuron_cell_M[neuron_cell_M$p_val < 0.05,]) # 568

neuron_gender_spec_genes = gender_spec_genes(neuron_cell_M, neuron_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(neuron_gender_spec_genes) # 22  6
table(neuron_gender_spec_genes$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#               11                4                7


#
# Astrocytes
#

astrocyte_cell_F <- FindMarkers(sctrans2, ident.1 = "Astrocyte TG F", ident.2 = "Astrocyte WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_cell_F) # 669   5
nrow(astrocyte_cell_F[astrocyte_cell_F$p_val_adj < 0.05,]) # 19
nrow(astrocyte_cell_F[astrocyte_cell_F$p_val < 0.05,]) # 189

astrocyte_cell_M <- FindMarkers(sctrans2, ident.1 = "Astrocyte TG M", ident.2 = "Astrocyte WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_cell_M) # 709   5
nrow(astrocyte_cell_M[astrocyte_cell_M$p_val_adj < 0.05,]) # 107
nrow(astrocyte_cell_M[astrocyte_cell_M$p_val < 0.05,]) # 373

astrocyte_gender_spec_genes = gender_spec_genes(astrocyte_cell_M, astrocyte_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(astrocyte_gender_spec_genes) # 24  6
table(astrocyte_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1               12                2                9


#
# Endothelial cells
#

endothelial_cell_F <- FindMarkers(sctrans2, ident.1 = "Endothelial cell TG F", ident.2 = "Endothelial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_F) # 877   5
nrow(endothelial_cell_F[endothelial_cell_F$p_val_adj < 0.05,]) # 22
nrow(endothelial_cell_F[endothelial_cell_F$p_val < 0.05,]) # 218

endothelial_cell_M <- FindMarkers(sctrans2, ident.1 = "Endothelial cell TG M", ident.2 = "Endothelial cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_M) # 913   5
nrow(endothelial_cell_M[endothelial_cell_M$p_val_adj < 0.05,]) # 132
nrow(endothelial_cell_M[endothelial_cell_M$p_val < 0.05,]) # 447

endothelial_gender_spec_genes = gender_spec_genes(endothelial_cell_M, endothelial_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(endothelial_gender_spec_genes) # 19  6
table(endothelial_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               2               10                2                5


#
# Oligodendrocyte Precursor Cells
#

oligodendrocyte_precursor_cell_F <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte precursor cell TG F", ident.2 = "Oligodendrocyte precursor cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_F) # 941   5
nrow(oligodendrocyte_precursor_cell_F[oligodendrocyte_precursor_cell_F$p_val_adj < 0.05,]) # 3
nrow(oligodendrocyte_precursor_cell_F[oligodendrocyte_precursor_cell_F$p_val < 0.05,]) # 100

oligodendrocyte_precursor_cell_M <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte precursor cell TG M", ident.2 = "Oligodendrocyte precursor cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_M) # 1015    5
nrow(oligodendrocyte_precursor_cell_M[oligodendrocyte_precursor_cell_M$p_val_adj < 0.05,]) # 9
nrow(oligodendrocyte_precursor_cell_M[oligodendrocyte_precursor_cell_M$p_val < 0.05,]) # 135

oligodendrocyte_precursor_gender_spec_genes = gender_spec_genes(oligodendrocyte_precursor_cell_M, oligodendrocyte_precursor_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(oligodendrocyte_precursor_gender_spec_genes) # 4 6
table(oligodendrocyte_precursor_gender_spec_genes$DEG.type)
# gender-shared male-specific
#             2             2


#
# Macrophage
#

macrophage_cell_F <- FindMarkers(sctrans2, ident.1 = "Macrophage TG F", ident.2 = "Macrophage WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_cell_F) # 1004    5
nrow(macrophage_cell_F[macrophage_cell_F$p_val_adj < 0.05,]) # 2
nrow(macrophage_cell_F[macrophage_cell_F$p_val < 0.05,]) # 119

macrophage_cell_M <- FindMarkers(sctrans2, ident.1 = "Macrophage TG M", ident.2 = "Macrophage WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_cell_M) # 1041    5
nrow(macrophage_cell_M[macrophage_cell_M$p_val_adj < 0.05,]) # 10
nrow(macrophage_cell_M[macrophage_cell_M$p_val < 0.05,]) # 166

macrophage_gender_spec_genes = gender_spec_genes(macrophage_cell_M, macrophage_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(macrophage_gender_spec_genes) # 3 6
table(macrophage_gender_spec_genes$DEG.type)
# gender-dimorphic    male-specific
#                2                1


#
# Mural cells
#

mural_cell_F <- FindMarkers(sctrans2, ident.1 = "Mural cell TG F", ident.2 = "Mural cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_F) # 1413    5
nrow(mural_cell_F[mural_cell_F$p_val_adj < 0.05,]) # 3
nrow(mural_cell_F[mural_cell_F$p_val < 0.05,]) # 92

mural_cell_M <- FindMarkers(sctrans2, ident.1 = "Mural cell TG M", ident.2 = "Mural cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_M) # 1186    5
nrow(mural_cell_M[mural_cell_M$p_val_adj < 0.05,]) # 2
nrow(mural_cell_M[mural_cell_M$p_val < 0.05,]) # 63

mural_gender_spec_genes = gender_spec_genes(mural_cell_M, mural_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(mural_gender_spec_genes) # 2 6
table(mural_gender_spec_genes$DEG.type)
# gender-shared male-specific
#             1             1


#
# Ependymal cells
#

ependymal_cell_F <- FindMarkers(sctrans2, ident.1 = "Ependymal cell TG F", ident.2 = "Ependymal cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_F) # 1267    5
nrow(ependymal_cell_F[ependymal_cell_F$p_val_adj < 0.05,]) # 1
nrow(ependymal_cell_F[ependymal_cell_F$p_val < 0.05,]) # 60

ependymal_cell_M <- FindMarkers(sctrans2, ident.1 = "Ependymal cell TG M", ident.2 = "Ependymal cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_M) # 1103    5
nrow(ependymal_cell_M[ependymal_cell_M$p_val_adj < 0.05,]) # 1
nrow(ependymal_cell_M[ependymal_cell_M$p_val < 0.05,]) # 39

ependymal_gender_spec_genes = gender_spec_genes(ependymal_cell_M, ependymal_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.50)
dim(ependymal_gender_spec_genes) # 1 6
table(ependymal_gender_spec_genes$DEG.type)
# gender-shared
#             1

save(microglial_cell_M, microglial_cell_F, oligodendrocyte_cell_F, oligodendrocyte_cell_M, 
 neuron_cell_F, neuron_cell_M, astrocyte_cell_F, astrocyte_cell_M,
 endothelial_cell_F, endothelial_cell_M, oligodendrocyte_precursor_cell_F, oligodendrocyte_precursor_cell_M,
 macrophage_cell_F, macrophage_cell_M, mural_cell_F, mural_cell_M, ependymal_cell_F, ependymal_cell_M,
 file="./results/files/TauCortex17m_DEGs_9CT_Poisson_SCT.RData")

# Saving cell type specific DEGs object (minabslog=0.50)
microglia_gender_spec_genes$cell_type <- "microglialCells"
oligodendrocyte_gender_spec_genes$cell_type <- "oligodendrocytes"
neuron_gender_spec_genes$cell_type <- "neurons"
astrocyte_gender_spec_genes$cell_type <- "astrocytes"

astrocytes_res_tauCortex_17m <- astrocyte_gender_spec_genes
microglial_res_tauCortex_17m <- microglia_gender_spec_genes
neurons_res_tauCortex_17m <- neuron_gender_spec_genes
oligodendrocyte_res_tauCortex_17m <- oligodendrocyte_gender_spec_genes
save(oligodendrocyte_res_tauCortex_17m, microglial_res_tauCortex_17m, astrocytes_res_tauCortex_17m, neurons_res_tauCortex_17m, 
  file="./results/files/TauCortex17m_DEGs_4CT_astro_oligo_microg_neuron_50percent_minabslog.RData")

endothelial_gender_spec_genes$cell_type <- "endothelial"
oligodendrocyte_precursor_gender_spec_genes$cell_type <- "OPCs"
mural_gender_spec_genes$cell_type <- "mural"
macrophage_gender_spec_genes$cell_type <- "macrophages"
ependymal_gender_spec_genes$cell_type <- "ependymals"

endothelial_res_tauCortex_17m <- endothelial_gender_spec_genes
OPC_res_tauCortex_17m <- oligodendrocyte_precursor_gender_spec_genes
mural_res_tauCortex_17m <- mural_gender_spec_genes
macrophage_res_tauCortex_17m <- macrophage_gender_spec_genes
ependymal_res_tauCortex_17m <- ependymal_gender_spec_genes

save(oligodendrocyte_res_tauCortex_17m, microglial_res_tauCortex_17m, astrocytes_res_tauCortex_17m, neurons_res_tauCortex_17m, endothelial_res_tauCortex_17m, OPC_res_tauCortex_17m, mural_res_tauCortex_17m, macrophage_res_tauCortex_17m, ependymal_res_tauCortex_17m, file="./results/files/TauCortex17m_DEGs_All9CTs_50percent_minabslog.RData")

CT_specific_DEGs <- rbind(oligodendrocyte_res_tauCortex_17m, microglial_res_tauCortex_17m, astrocytes_res_tauCortex_17m, neurons_res_tauCortex_17m, endothelial_res_tauCortex_17m, OPC_res_tauCortex_17m, mural_res_tauCortex_17m, macrophage_res_tauCortex_17m, ependymal_res_tauCortex_17m)
dim(CT_specific_DEGs) # 126   7
table(CT_specific_DEGs$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               4               48               24               50
               
# Sanity Check: if the there is any gender-dimorphic DEG with abs(LogFC) < 0.5
unique(CT_specific_DEGs[abs(CT_specific_DEGs$Male.avg..logFC) < 0.5,]$DEG.type) # "gender-shared"   "female-specific"
unique(CT_specific_DEGs[abs(CT_specific_DEGs$Female.avg..logFC) < 0.5,]$DEG.type) # "male-specific" "gender-shared"

nrow(CT_specific_DEGs[abs(CT_specific_DEGs$Male.avg..logFC) < 0.5,]) # 2
nrow(CT_specific_DEGs[abs(CT_specific_DEGs$Female.avg..logFC) < 0.5,]) # 28

head(sort(CT_specific_DEGs$Male.FDR), 10)
# [1]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [9]  0.000000e+00 1.823309e-266
head(sort(CT_specific_DEGs$Female.FDR))
# [1] 7.556180e-142 4.890386e-119 4.298870e-110  1.951927e-65  2.360167e-49  5.508267e-42
CT_specific_DEGs$Male.FDR <- ifelse(CT_specific_DEGs$Male.FDR == 0, 1.823309e-266, CT_specific_DEGs$Male.FDR)

write.csv(CT_specific_DEGs, file="./results/files/TauCortex17m_DEGs_All9CTs_50percent_minabslog.csv", row.names=F, quote=F)

# Saving complete cell type specific DEGs object
microglial_cell_F$Sex <- "Female-specific"
microglial_cell_F$Cell_Type <- "microglial_cells"
microglial_cell_M$Sex <- "Male-specific"
microglial_cell_M$Cell_Type <- "Microglial_cells"
microglial_cell <- rbind(microglial_cell_F, microglial_cell_M)
dim(microglial_cell) # 1593    7

oligodendrocyte_cell_F$Sex <- "Female-specific"
oligodendrocyte_cell_F$Cell_Type <- "Oligodendrocytes"
oligodendrocyte_cell_M$Sex <- "Male-specific"
oligodendrocyte_cell_M$Cell_Type <- "Oligodendrocytes"
oligodendrocyte_cell <- rbind(oligodendrocyte_cell_F, oligodendrocyte_cell_M)
dim(oligodendrocyte_cell) # 1672    7

neuron_cell_F$Sex <- "Female-specific"
neuron_cell_F$Cell_Type <- "Neurons"
neuron_cell_M$Sex <- "Male-specific"
neuron_cell_M$Cell_Type <- "Neurons"
neuron_cell <- rbind(neuron_cell_F, neuron_cell_M)
dim(neuron_cell) # 1624    7

astrocyte_cell_F$Sex <- "Female-specific"
astrocyte_cell_F$Cell_Type <- "Astrocytes"
astrocyte_cell_M$Sex <- "Male-specific"
astrocyte_cell_M$Cell_Type <- "Astrocytes"
astrocyte_cell <- rbind(astrocyte_cell_F, astrocyte_cell_M)
dim(astrocyte_cell) # 1378    7

endothelial_cell_F$Sex <- "Female-specific"
endothelial_cell_F$Cell_Type <- "Endothelial_cells"
endothelial_cell_M$Sex <- "Male-specific"
endothelial_cell_M$Cell_Type <- "Endothelial_cells"
endothelial_cell <- rbind(endothelial_cell_F, endothelial_cell_M)
dim(endothelial_cell) # 1790    7

oligodendrocyte_precursor_cell_F$Sex <- "Female-specific"
oligodendrocyte_precursor_cell_F$Cell_Type <- "Oligodendrocyte_precursor_cells"
oligodendrocyte_precursor_cell_M$Sex <- "Male-specific"
oligodendrocyte_precursor_cell_M$Cell_Type <- "Oligodendrocyte_precursor_cells"
oligodendrocyte_precursor_cell <- rbind(oligodendrocyte_precursor_cell_F, oligodendrocyte_precursor_cell_M)
dim(oligodendrocyte_precursor_cell) # 1956    7

macrophage_cell_F$Sex <- "Female-specific"
macrophage_cell_F$Cell_Type <- "Macrophages"
macrophage_cell_M$Sex <- "Male-specific"
macrophage_cell_M$Cell_Type <- "Macrophages"
macrophage_cell <- rbind(macrophage_cell_F, macrophage_cell_M)
dim(macrophage_cell) # 2045    7

mural_cell_F$Sex <- "Female-specific"
mural_cell_F$Cell_Type <- "Mural_cells"
mural_cell_M$Sex <- "Male-specific"
mural_cell_M$Cell_Type <- "Mural_cells"
mural_cell <- rbind(mural_cell_F, mural_cell_M)
dim(mural_cell) # 2599    7

ependymal_cell_F$Sex <- "Female-specific"
ependymal_cell_F$Cell_Type <- "Ependymal_cells"
ependymal_cell_M$Sex <- "Male-specific"
ependymal_cell_M$Cell_Type <- "Ependymal_cells"
ependymal_cell <- rbind(ependymal_cell_F, ependymal_cell_M)
dim(ependymal_cell) # 2370    7

save(microglial_cell_F, microglial_cell_M, oligodendrocyte_cell_F, oligodendrocyte_cell_M, neuron_cell_F, neuron_cell_M, astrocyte_cell_F, astrocyte_cell_M, endothelial_cell_F, endothelial_cell_M, oligodendrocyte_precursor_cell_F, oligodendrocyte_precursor_cell_M, macrophage_cell_F, macrophage_cell_M, mural_cell_F, mural_cell_M,  ependymal_cell_F, ependymal_cell_M, file="./results/files/TauCortex17m_DEGs_All9CTs_Complete_Sumstats.RData")

CT_specific_DEGs_Complete <- rbind(microglial_cell, oligodendrocyte_cell, neuron_cell, astrocyte_cell, endothelial_cell, oligodendrocyte_precursor_cell, macrophage_cell, mural_cell, ependymal_cell)
dim(CT_specific_DEGs_Complete) # 17027     7
head(unique(sort(CT_specific_DEGs_Complete$p_val_adj)))
# [1]  0.000000e+00 1.823309e-266 3.880433e-255 1.614712e-202 7.390897e-191 1.314476e-186
head(unique(sort(CT_specific_DEGs_Complete$p_val)))
# [1]  0.000000e+00 1.031225e-270 2.194691e-259 9.132469e-207 4.180135e-195 7.434399e-191
CT_specific_DEGs_Complete$p_val_adj <- ifelse(CT_specific_DEGs_Complete$p_val_adj == 0, 1.823309e-266, CT_specific_DEGs_Complete$p_val_adj)
CT_specific_DEGs_Complete$p_val <- ifelse(CT_specific_DEGs_Complete$p_val == 0, 1.031225e-270, CT_specific_DEGs_Complete$p_val)
write.table(CT_specific_DEGs_Complete, file="./results/files/TauCortex17m_DEGs_All9CTs_Complete_Sumstats.txt", sep="\t", quote=F)