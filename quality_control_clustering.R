## Quality control and pre-processing ##

R # R version 4.4.0

library(Seurat) # Seurat_5.1.0
library(glmGamPoi) # glmGamPoi_1.16.0
library(cluster) # cluster_2.1.6
library(ggplot2) # ggplot2_3.5.1

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1)
options(width=180)

load("./data/thy21_tau_mouse_cortex_scrnaseq.Rdata")

dim(merged) # 19049 28530
merged
# An object of class Seurat
# 19049 features across 28530 samples within 1 assay
# Active assay: RNA (19049 features, 0 variable features)

table(merged@meta.data$conditions)
# TG F  TG M  WT F  WT M
# 8940  7050  4890  7650


# Percentage mitochondrial genes
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

summary(merged@meta.data$percent.mt)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#        0       0       0       0       0       0

# no mitochondrial contamination

pdf("./results/figures/violinplot.pdf")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# We filter cells that have unique feature counts over 7000 or less than 200
# We filter cells that have >5% mitochondrial counts
dat_filt <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)
dat_filt
# An object of class Seurat
# 19049 features across 28190 samples within 1 assay
# Active assay: RNA (19049 features, 0 variable features)

pdf("./results/figures/violinplot_PostCellFiltration.pdf")
VlnPlot(dat_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# Warning message:
# In SingleExIPlot(type = type, data = data[, x, drop = FALSE], idents = idents,  :
#   All cells have the same value of percent.mt.

# use new pre-processing approach
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
sctrans2 = SCTransform(dat_filt, vst.flavor = "v2")
# An object of class Seurat
# 36730 features across 28190 samples within 2 assays
# Active assay: SCT (17681 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
#  1 other assay present: RNA
dim(sctrans2@assays$SCT) # 17681 28190


## Clustering ##

sctrans2 <- RunPCA(sctrans2, npcs = 30, verbose = FALSE) # only 30 PCs used in next analysis

sctrans2 <- RunUMAP(sctrans2, reduction = "pca", dims = 1:30, verbose = FALSE) # reduction = "pca" is default parameter
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'

sctrans2 <- FindNeighbors(sctrans2, reduction = "pca", dims = 1:30, verbose = FALSE)
sctrans2 <- FindClusters(sctrans2, resolution = 0.8, verbose = FALSE) 
# resolution, Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. default = 0.8
# consider adjusting the resolution parameter later for the clustering to obtain a smaller number of communities

pdf("./results/figures/dimplot_umap_v2.pdf")
DimPlot(sctrans2, label = TRUE) + NoLegend()
dev.off()


## Evaluate clustering - choose optimal number of clusters automatically ##

# Elbow plot: a ranking of principle components based on the percentage of variance explained by each one
pdf("./results/figures/ElbowPlot_PC30_mm_Cortex.pdf")
ElbowPlot(sctrans2)
dev.off()

# "A more quantitative approach may be a bit more reliable. We can calculate where the principal components start to elbow by taking the larger value of:
# The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# The point where the percent change in variation between the consecutive PCs is less than 0.1%.
# We will start by calculating the first metric:"

# Determine percent of variation associated with each PC
pct <- sctrans2[["pca"]]@stdev / sum(sctrans2[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
# 25

#The first metric returns PC25 as the PC matching these requirements. Lets check the second metric, which identifies the PC where the percent change in variation between consecutive PCs is less than 0.1%:

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
#19

#This second metric returns PC19. Usually, we would choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.

# Minimum of the two calculations
pcs <- min(co1, co2)

pcs
# 19


# Create a dataframe with values
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize
pdf("./results/figures/colored_elbow_plot.pdf")
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()


#
# Alternative: Silhouette width approach
#

# Make silhoutte plots at different clustering resolutions (of seurat) to identify best cluster number/size

reduction <- "pca"
dims <- 1:pcs
dist.matrix <- dist(x = Embeddings(object = sctrans2[[reduction]])[, dims])

resolution <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# resolution <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

silhoutte_width <- data.frame()
# par(mfrow=c(4,5))
# png("silhoutte_width_plots.png", width=8, height=6, units="in", res=300)
for (i in 1:length(resolution)) {
	print(resolution[i])
	sctrans2 <- FindClusters(sctrans2, resolution = resolution[i])
	clusters <- eval(parse(text = paste0("sctrans2$SCT_snn_res.",resolution[i])))
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	#plot(sil, main=resolution[i]) # silhoutte_5_cluster_res0.01 # Average silhoutte width: 0.63
	
	if(length(table(clusters)) < 2)
		next
	
	print(mean(sil[, "sil_width"])) # 0.5256185
	result <- cbind(resolution[i], length(unique(clusters)), mean(sil[, "sil_width"]))
	colnames(result) <- c("resolution", "clusters", "sil_width")
	silhoutte_width <- rbind(silhoutte_width, result)
}
# dev.off

silhoutte_width
#    resolution clusters sil_width
# 1       0.005        4 0.3844026
# 2       0.010        5 0.3708660
# 3       0.020        7 0.4536703
# 4       0.030        8 0.4574719
# 5       0.040        8 0.4597888
# 6       0.050        9 0.4607528
# 7       0.060        9 0.4627498 # Max. value ==> 9 clusters
# 8       0.070        9 0.3718682
# 9       0.080       10 0.3375138
# 10      0.100       10 0.3375459
# 11      0.200       13 0.2605275
# 12      0.300       15 0.2587615
# 13      0.400       19 0.2010725
# 14      0.500       19 0.2182632
# 15      0.600       19 0.2237676
# 16      0.700       22 0.2196651
# 17      0.800       23 0.1920736
# 18      0.900       26 0.1760318
# 19      1.000       27 0.1759401
pdf("./results/figures/silhoutte_width_mm_Cortex.pdf")
ggplot(silhoutte_width, aes(x = clusters, y = sil_width)) +
  geom_line() + geom_point() + scale_x_continuous(breaks = 2:31) + theme_bw()
dev.off()

# NOTE: Best result are at "0.06" resolution: 9 clusters with avg. sil_width = 0.462 (highest)
bestres = 0.06

sctrans2 <- FindClusters(sctrans2, resolution = bestres)
#Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#Maximum modularity in 10 random starts: 0.9850
#Number of communities: 9
#Elapsed time: 13 seconds

# Look at the number of cells in each cluster
table(sctrans2$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
# 8263 5505 5299 5086 2183 1254  307  188  105


# Look at cluster IDs of the first 5 cells
head(Idents(sctrans2), 5)

head(sctrans2@meta.data$conditions)
sctrans2$stim <- sctrans2@meta.data$conditions

pdf("./results/figures/DimPlot_umap_mm_Cortex.pdf")
DimPlot(sctrans2, reduction = "umap", split.by = "stim")
dev.off()

Idents(sctrans2) = sctrans2$seurat_clusters

png("./results/figures/DimPlot_umap_mm_Cortex.png", width=10, height=6, units="in", res=300)
DimPlot(sctrans2, reduction = "umap", split.by = "stim")
dev.off()

sctrans2 <- RunTSNE(sctrans2, reduction = "pca", dims = 1:pcs)
pdf("./results/figures/DimPlot_tSNE_mm_Cortex.pdf")
DimPlot(sctrans2, reduction = "tsne", split.by = "stim")
dev.off()

pdf("./results/figures/DimPlot_umap_mm_Cortex_nosplit.pdf")
DimPlot(sctrans2, reduction = "umap", label = T, repel = T) + ggtitle("Unsupervised clustering")
dev.off()

save(sctrans2, file="./data/sctrans2.RData")


## Clustering Annotation ##

rm(list = ls())
library(Seurat) # Seurat_5.1.0
library(ggplot2) # ggplot2_3.5.1
library(cluster) # cluster_2.1.6
library(tidyverse) # tidyverse_2.0.0
library(dplyr) # dplyr_1.1.4
library(HGNChelper) # HGNChelper_0.8.15
library(openxlsx) # openxlsx_4.2.5.2
set.seed(1)

load("./data/sctrans2.RData")


DefaultAssay(sctrans2) <- "SCT"

cluster_markers <- FindAllMarkers(sctrans2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)

table(cluster_markers$cluster)
#   0   1   2   3   4   5   6   7   8
# 152 120 113 147  73 117 134 168 137

cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC) %>% print(n = 30)
# Groups:   cluster [9]
#       p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene
#       <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
#  1 0              2.35 0.957 0.078  0        0       Ctss
#  2 0              2.17 0.989 0.461  0        0       Cst3
#  3 0              1.94 0.894 0.066  0        0       Cx3cr1
#  4 0              5.04 0.999 0.291  0        1       Plp1
#  5 0              1.99 0.888 0.029  0        1       Ermn
#  6 0              1.90 0.953 0.42   0        1       Fth1
#  7 0              2.52 0.98  0.177  0        2       Slc1a2
#  8 0              2.07 0.964 0.202  0        2       Atp1a2
#  9 0              1.84 0.945 0.312  0        2       Sparcl1
# 10 0              2.02 0.89  0.047  0        3       Flt1
# 11 0              1.55 0.775 0.087  0        3       Bsg
# 12 0              1.41 0.814 0.168  0        3       Sptbn1
# 13 0              3.14 0.71  0.054  0        4       Meg3
# 14 0              1.45 0.546 0.02   0        4       Snhg11
# 15 0              1.07 0.48  0.014  0        4       Syt1
# 16 0              1.46 0.861 0.156  0        5       Ptprz1
# 17 0              1.44 0.753 0.015  0        5       Cacng4
# 18 0              1.23 0.688 0.012  0        5       Tnr
# 19 3.77e-61       1.76 0.694 0.389  6.67e-57 6       Apoe
# 20 0              1.40 0.267 0.01   0        6       Cd74
# 21 0              1.39 0.577 0.045  0        6       Lyz2
# 22 0              2.42 0.915 0.025  0        7       Mgp
# 23 0              2.28 0.793 0.002  0        7       Prg4
# 24 0              1.80 0.718 0.019  0        7       Igfbp5
# 25 2.00e-96       2.27 0.952 0.256  3.54e-92 8       Dbi
# 26 1.18e-75       1.72 0.924 0.3    2.09e-71 8       Hsp90aa1
# 27 1.46e-74       1.58 0.905 0.262  2.58e-70 8       Gm12346

cluster_markers_top10 <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) %>% print(n = 100)
write.csv(as.data.frame(cluster_markers_top10), file="./results/files/cluster_markers_top10.txt", quote=F)


cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 4, wt = avg_log2FC) -> top4

png("./results/figures/Featureplot_top4_marker_genes.png", width=12, height=14, units="in", res=300)
FeaturePlot(sctrans2, features = top4$gene, 
  label = TRUE,
  label.size = 4,
  label.color = "red",)
dev.off()


# remove genes that do not exist in scale.data as it will give error gene not found
gene_panel_intersect <- intersect(top4$gene, rownames(sctrans2@assays$SCT@scale.data))
length(gene_panel_intersect) # 35

# Down-sample the cells because with all cells, heatmap does not work
# https://github.com/satijalab/seurat/issues/2724
png("./results/figures/Heatmap_top_marker_genes.png", width=7, height=7, units="in", res=300)
DoHeatmap(subset(sctrans2, downsample = 100), 
          features = gene_panel_intersect, 
                  assay = "SCT")
dev.off()


## Automated cell-type mapping and annotation ##

# download SCType source code
#system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R')
#system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R')

# load downloaded source
source("./data/gene_sets_prepare.R");
source("./data/sctype_score_.R")

# DB file
# system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx')

db_ = "./data/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain


# Load prepared mouse annotations from Cell Marker database
load(file="./data/cellmarker_scrnaseq.Rdata")


#
# Adjusted version of the sctype_score function (avoiding human gene format in upper case letter format for compatibility with mouse MGI symbols)
#

sctype_score_noupper = function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T);
 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)
	# higher number of occurrence = lower sensitivity


  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

  # combine scores (across marker genes per cell type)
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  es.max
}


gspositive_lst = sapply(cellmarker_scrnaseq$Cell.Marker, function(x) strsplit(x,", ")[[1]])
names(gspositive_lst) = cellmarker_scrnaseq$Cell.Type

# r merge named vectors with duplicate names
# https://stackoverflow.com/questions/57020599/r-list-combine-elements-with-same-name
gspositive_combined = tapply(unlist(gspositive_lst, use.names = FALSE), rep(names(gspositive_lst), lengths(gspositive_lst)), FUN = c)

# make unique
gspositive_combined  = sapply(gspositive_combined,function(x) unique(x))

es.max_new = sctype_score_noupper(scRNAseqData = sctrans2[["SCT"]]@scale.data, scaled = TRUE, gs = gspositive_combined, gs2 = NULL)
dim(es.max_new)
# 36 28190

cL_resutls_new = do.call("rbind", lapply(unique(sctrans2@meta.data$seurat_clusters), function(cl){
		#print(cl)
    es.max.cl = sort(rowSums(es.max_new[ ,rownames(sctrans2@meta.data[sctrans2@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sctrans2@meta.data$seurat_clusters==cl)), 10)
}))

dim(cL_resutls_new)
# 90  4

# show best annotations and corresponding scores for cluster 0
head(cL_resutls_new[which(cL_resutls_new$cluster == 0),], 10)


# filter clusters with too low scores as "unknown"
sctype_scores_new = cL_resutls_new %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

sctype_scores_new$type[as.numeric(as.character(sctype_scores_new$scores)) < sctype_scores_new$ncells/4] = "Unknown"


# show overview of cluster assignments
print(sctype_scores_new[,1:4])
scType_annotation <- as.data.frame(sctype_scores_new[,1:4])
scType_annotation[order(scType_annotation$cluster),]
#   cluster                           type    scores ncells
# 6       0                Microglial cell 49445.513   8263
# 8       1                Oligodendrocyte 52785.966   5505
# 7       2                      Astrocyte 31165.083   5299
# 1       3               Endothelial cell 41353.092   5086
# 2       4                         Neuron 13885.668   2183
# 4       5 Oligodendrocyte precursor cell 11347.552   1254
# 9       6                     Macrophage  1666.640    307
# 5       7                     Mural cell   805.181    188
# 3       8                 Ependymal cell  1105.206    105


# Look at the number of cells in each cluster
table(sctrans2$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
# 8263 5505 5299 5086 2183 1254  307  188  105

sctrans2@meta.data$classint = ""
for(j in unique(sctype_scores_new$cluster)){
  cl_type = sctype_scores_new[sctype_scores_new$cluster==j,]; 
  sctrans2@meta.data$classint[sctrans2@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Dimplot with cluster annotation
pdf("./results/figures/umap_plot_celltype_annot_2024-01-24.pdf")
DimPlot(sctrans2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'classint')  
dev.off()

# Dimplot with cluster annotation
png("./results/figures/Dimplot_SCtype_cluster_labelling.png", width=7, height=7, units="in", res=300)
DimPlot(sctrans2, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

# Dimplot with cluster annotation
png("./results/figures/Dimplot_SCtype_cluster_labelling_v2.png", width=7, height=5, units="in", res=300)
DimPlot(sctrans2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'classint') + ggtitle(NULL) +
  xlab("Dimension 1") + ylab("Dimension 2")
dev.off()

save(sctrans2, file="./data/sctrans2.RData")


# Regenerate Marker Heatmap figure
Idents(sctrans2) <- sctrans2$classint
table(Idents(sctrans2))
cluster_markers <- FindAllMarkers(sctrans2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)

cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 4, wt = avg_log2FC) -> top4

# remove genes that do not exist in scale.data as it will give error gene not found
gene_panel_intersect <- intersect(top4$gene, rownames(sctrans2@assays$SCT@scale.data))
length(gene_panel_intersect) # 35

# Down-sample the cells because with all cells, heatmap does not work
# https://github.com/satijalab/seurat/issues/2724
png("./results/figures/Heatmap_top_marker_genes_labelled.png", width=7, height=7, units="in", res=300)
DoHeatmap(subset(sctrans2, downsample = 100), 
          features = gene_panel_intersect, 
                  assay = "SCT", label = FALSE)
dev.off()



# Automatically detect the tissue type of the dataset

# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
png("./results/figures/Tissue_labelling.png", width=12, height=7, units="in", res=300)
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = sctrans2, scaled = TRUE, assay = "SCT")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used
dev.off()
dim(tissue_guess) # 15  2
print(tissue_guess)
#           tissue     score
# 14       Stomach 10081.338
# 6          Brain  9553.198
# 10     Intestine  7712.973
# 4            Eye  5784.283
# 7           Lung  4674.203
# 3          Liver  4442.168
# 1  Immune system  4043.358
# 8        Adrenal  3900.773
# 9          Heart  3394.297
# 5         Kidney  3303.271
# 12      Placenta  3017.816
# 11        Muscle  3002.584
# 2       Pancreas  2849.521
# 15        Thymus  2491.279
# 13        Spleen  1589.338