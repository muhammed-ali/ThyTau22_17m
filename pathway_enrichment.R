## Global Gene set enrichment analysis ##

library(clusterProfiler) # clusterProfiler_4.12.0
library(org.Mm.eg.db) # org.Mm.eg.db_3.19.1
library(enrichplot) # enrichplot_1.24.0 {http://yulab-smu.top/clusterProfiler-book/chapter12.html}
library(tidyverse) # tidyverse_2.0.0
library(cowplot) # cowplot_1.1.3
library(venn) # venn_1.12
library(ggplot2) # ggplot2_3.5.1
library(rrvgo) # rrvgo_1.16.0
library(openxlsx) # openxlsx_4.2.5.2
set.seed(1)

load("./results/files/Global_DEGs.RData")

DEG_Global_Poisson <- rbind(DEG_Global_Poisson_M, DEG_Global_Poisson_F)
DEG_Global_Poisson$SYMBOL <- rownames(DEG_Global_Poisson)

table(row.names(DEG_Global_Poisson_F) %in% row.names(DEG_Global_Poisson_M))
# FALSE  TRUE
#    18   687

table(global_gender_spec_genes$DEG.type)
#  female-specific gender-dimorphic    gender-shared    male-specific
#                1               21               10               15


#
# Sex-dimorphic
#

gene.df <- bitr(global_gender_spec_genes$Gene.symbols[which(global_gender_spec_genes$DEG.type=="gender-dimorphic")], fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 5.88% of input gene IDs are fail to map...
dim(gene.df) # 19  2
GOI <- left_join(gene.df, DEG_Global_Poisson, by="SYMBOL")
dim(GOI) # 19  7

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 60  9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 16  9

BP_GenderDimorphic <- BioProc
MF_GenderDimorphic <- MolFun


# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

DEG_GenderDimorphic <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-dimorphic",]
dim(DEG_GenderDimorphic) # 21  6

# Gender-Shared
MF_GenderDimorphic@result$up <- NA
MF_GenderDimorphic@result$down <- NA
BP_GenderDimorphic@result$up <- NA
BP_GenderDimorphic@result$down <- NA
get_up_down_genes_GS <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0 & deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0 & deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_GenderDimorphic <- get_up_down_genes_GS(MF_GenderDimorphic, DEG_GenderDimorphic)
BP_GenderDimorphic <- get_up_down_genes_GS(BP_GenderDimorphic, DEG_GenderDimorphic)

png("./results/figures/Global_gender_dimorphic_DEGs_enrichment_BP_MF.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_GenderDimorphic, showCategory=18) + ggtitle("Sex-dimorphic: BP Enrichment")
p2 <- dotplot(MF_GenderDimorphic, showCategory=18) + ggtitle("Sex-dimorphic: MF Enrichment")
plot_grid(p1, p2, ncol=2)
dev.off()

# Treeplot
BP_GD_simMatrix <- calculateSimMatrix(BP_GenderDimorphic$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_GD_scores <- setNames(-log10(BP_GenderDimorphic$qvalue), BP_GenderDimorphic$ID)
BP_GD_reducedTerms <- reduceSimMatrix(BP_GD_simMatrix,
                                BP_GD_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Global_gender_dimorphic_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_GD_reducedTerms)
dev.off()

MF_GD_simMatrix <- calculateSimMatrix(MF_GenderDimorphic$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_GD_scores <- setNames(-log10(MF_GenderDimorphic$qvalue), MF_GenderDimorphic$ID)
MF_GD_reducedTerms <- reduceSimMatrix(MF_GD_simMatrix,
                                MF_GD_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Global_gender_dimorphic_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_GD_reducedTerms)
dev.off()


#
# male-specific
#

gene.df <- bitr(global_gender_spec_genes$Gene.symbols[which(global_gender_spec_genes$DEG.type=="male-specific")], fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 3.23% of input gene IDs are fail to map...
dim(gene.df) # 15  2

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 4 9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 26  9

BP_Male <- BioProc
MF_Male <- MolFun


# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

DEG_Male <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
dim(DEG_Male) # 15  6


# Male
MF_Male@result$up <- NA
MF_Male@result$down <- NA
BP_Male@result$up <- NA
BP_Male@result$down <- NA
get_up_down_genes_Male <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Male <- get_up_down_genes_Male(MF_Male, DEG_Male)
BP_Male <- get_up_down_genes_Male(BP_Male, DEG_Male)

png("./results/figures/Globa_Male_DEGs_enrichment_BP_MF_Splitted_upDown.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Male, showCategory=12) + ggtitle("Male-specific: BP Enrichment") + facet_grid(.~upDown)
p2 <- dotplot(MF_Male, showCategory=8) + ggtitle("Male-specific: MF Enrichment") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()


# Treeplot
BP_Male_simMatrix <- calculateSimMatrix(BP_Male$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_Male_scores <- setNames(-log10(BP_Male$qvalue), BP_Male$ID)
BP_Male_reducedTerms <- reduceSimMatrix(BP_Male_simMatrix,
                                BP_Male_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Global_Male_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_Male_reducedTerms)
dev.off()

MF_Male_simMatrix <- calculateSimMatrix(MF_Male$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_Male_scores <- setNames(-log10(MF_Male$qvalue), MF_Male$ID)
MF_Male_reducedTerms <- reduceSimMatrix(MF_Male_simMatrix,
                                MF_Male_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Global_Male_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_Male_reducedTerms)
dev.off()

save(BP_GenderDimorphic, MF_GenderDimorphic, BP_Male, MF_Male, file="./results/files/Global_Enrichment_Male_GenderDimorphic_Dec2024.RData")
BP_Male@result$Category <- "BP"
MF_Male@result$Category <- "MF"
Male_Pathways <- rbind(BP_Male@result, MF_Male@result)
Male_Pathways <- Male_Pathways[Male_Pathways$p.adjust < 0.05,]
dim(Male_Pathways) # 30 13
write.xlsx(Male_Pathways, file="./results/files/Male_Pathways_Global.xlsx")

BP_GenderDimorphic@result$Category <- "BP"
MF_GenderDimorphic@result$Category <- "MF"
GD_Pathways <- rbind(BP_GenderDimorphic@result, MF_GenderDimorphic@result)
GD_Pathways <- GD_Pathways[GD_Pathways$p.adjust < 0.05,]
dim(GD_Pathways) # 76 13
write.xlsx(GD_Pathways, file="./results/files/GenderDimorphic_Pathways_Global.xlsx")


## Cell-type-specific Gene set enrichment analysis ##

#
# Microglial Cells
#

load("./results/files/TauCortex17m_DEGs_4CT_astro_oligo_microg_neuron_50percent_minabslog.RData")

dim(microglial_res_tauCortex_17m) # 30  7
table(microglial_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1                6                9               14

gene.df <- bitr(microglial_res_tauCortex_17m$Gene.symbols, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 10% of input gene IDs are fail to map...
dim(gene.df) # 27  2
GOI <- left_join(gene.df, microglial_res_tauCortex_17m, by=c("SYMBOL"="Gene.symbols"))
dim(GOI) # 27  8

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 91   9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 12  9

BP_Microglial <- BioProc
MF_Microglial <- MolFun

# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

# Gender-Shared
MF_Microglial@result$up <- NA
MF_Microglial@result$down <- NA
BP_Microglial@result$up <- NA
BP_Microglial@result$down <- NA
get_up_down_genes_GS <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0 & deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0 & deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Microglial <- get_up_down_genes_GS(MF_Microglial, microglial_res_tauCortex_17m)
BP_Microglial <- get_up_down_genes_GS(BP_Microglial, microglial_res_tauCortex_17m)

png("./results/figures/Microglial_DEGs_enrichment_BP_MF.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Microglial, showCategory=16) + ggtitle("Microglial cells: BP Enrichment") + facet_grid(.~upDown)
p2 <- dotplot(MF_Microglial, showCategory=18) + ggtitle("Microglial cells: MF Enrichment") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2, rel_widths = c(1.5, 1))
dev.off()

# Treeplot
BP_simMatrix <- calculateSimMatrix(BP_Microglial$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_scores <- setNames(-log10(BP_Microglial$qvalue), BP_Microglial$ID)
BP_reducedTerms <- reduceSimMatrix(BP_simMatrix,
                                BP_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Microglial_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_reducedTerms)
dev.off()

MF_simMatrix <- calculateSimMatrix(MF_Microglial$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_scores <- setNames(-log10(MF_Microglial$qvalue), MF_Microglial$ID)
MF_reducedTerms <- reduceSimMatrix(MF_simMatrix,
                                MF_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Microglial_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_reducedTerms)
dev.off()


#
# Oligodendrocyte
#

dim(oligodendrocyte_res_tauCortex_17m) # 21  7
table(oligodendrocyte_res_tauCortex_17m$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#                7                3               11

gene.df <- bitr(oligodendrocyte_res_tauCortex_17m$Gene.symbols, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 23.81% of input gene IDs are fail to map...
dim(gene.df) # 16  2
GOI <- left_join(gene.df, oligodendrocyte_res_tauCortex_17m, by=c("SYMBOL"="Gene.symbols"))
dim(GOI) # 16  8

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 109   9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 5 9

BP_Oligo <- BioProc
MF_Oligo <- MolFun

# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

# Gender-Shared
MF_Oligo@result$up <- NA
MF_Oligo@result$down <- NA
BP_Oligo@result$up <- NA
BP_Oligo@result$down <- NA

MF_Oligo <- get_up_down_genes_GS(MF_Oligo, oligodendrocyte_res_tauCortex_17m)
BP_Oligo <- get_up_down_genes_GS(BP_Oligo, oligodendrocyte_res_tauCortex_17m)

png("./results/figures/Oligodendrocytes_DEGs_enrichment_BP_MF.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Oligo, showCategory=16) + ggtitle("Oligodendrocytes: BP Enrichment") + facet_grid(.~upDown)
p2 <- dotplot(MF_Oligo, showCategory=18) + ggtitle("Oligodendrocytes: MF Enrichment") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2, rel_widths = c(1, 1.5))
dev.off()

# Treeplot
BP_simMatrix <- calculateSimMatrix(BP_Oligo$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_scores <- setNames(-log10(BP_Oligo$qvalue), BP_Oligo$ID)
BP_reducedTerms <- reduceSimMatrix(BP_simMatrix,
                                BP_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Oligodendrocytes_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_reducedTerms)
dev.off()

MF_simMatrix <- calculateSimMatrix(MF_Oligo$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_scores <- setNames(-log10(MF_Oligo$qvalue), MF_Oligo$ID)
MF_reducedTerms <- reduceSimMatrix(MF_simMatrix,
                                MF_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Oligodendrocytes_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_reducedTerms)
dev.off()


#
# Astrocyte
#

dim(astrocytes_res_tauCortex_17m) # 24  7
table(astrocytes_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1               12                2                9

gene.df <- bitr(astrocytes_res_tauCortex_17m$Gene.symbols, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 8.33% of input gene IDs are fail to map...
dim(gene.df) # 22  2
GOI <- left_join(gene.df, astrocytes_res_tauCortex_17m, by=c("SYMBOL"="Gene.symbols"))
dim(GOI) # 22  8

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 99  9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 22  9

BP_Astro <- BioProc
MF_Astro <- MolFun

# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

# Gender-Shared
MF_Astro@result$up <- NA
MF_Astro@result$down <- NA
BP_Astro@result$up <- NA
BP_Astro@result$down <- NA
get_up_down_genes_GS <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0 & deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0 & deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Astro <- get_up_down_genes_GS(MF_Astro, astrocytes_res_tauCortex_17m)
BP_Astro <- get_up_down_genes_GS(BP_Astro, astrocytes_res_tauCortex_17m)

png("./results/figures/Astrocytes_DEGs_enrichment_BP_MF.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Astro, showCategory=16) + ggtitle("Astrocytes: BP Enrichment") + facet_grid(.~upDown)
p2 <- dotplot(MF_Astro, showCategory=18) + ggtitle("Astrocytes: MF Enrichment") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()

# Treeplot
BP_simMatrix <- calculateSimMatrix(BP_Astro$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_scores <- setNames(-log10(BP_Astro$qvalue), BP_Astro$ID)
BP_reducedTerms <- reduceSimMatrix(BP_simMatrix,
                                BP_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Astrocytes_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_reducedTerms)
dev.off()

MF_simMatrix <- calculateSimMatrix(MF_Astro$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_scores <- setNames(-log10(MF_Astro$qvalue), MF_Astro$ID)
MF_reducedTerms <- reduceSimMatrix(MF_simMatrix,
                                MF_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Astrocytes_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_reducedTerms)
dev.off()

save(BP_Microglial, MF_Microglial, BP_Oligo, MF_Oligo, BP_Astro, MF_Astro, 
    file="./results/files/Cell_Type_Specific_Enrichment_Microg_Oligo_Astro_Dec2024.RData")


#
# Neurons
#

dim(neurons_res_tauCortex_17m) # 22  7
table(neurons_res_tauCortex_17m$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#               11                4                7

gene.df <- bitr(neurons_res_tauCortex_17m$Gene.symbols, fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 4.55% of input gene IDs are fail to map...
dim(gene.df) # 21  2
GOI <- left_join(gene.df, neurons_res_tauCortex_17m, by=c("SYMBOL"="Gene.symbols"))
dim(GOI) # 21  8

BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 81  9

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 59  9

BP_Neuron <- BioProc
MF_Neuron <- MolFun

# Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs

# Gender-Shared
MF_Neuron@result$up <- NA
MF_Neuron@result$down <- NA
BP_Neuron@result$up <- NA
BP_Neuron@result$down <- NA
get_up_down_genes_GS <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0 & deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0 & deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Neuron <- get_up_down_genes_GS(MF_Neuron, neurons_res_tauCortex_17m)
BP_Neuron <- get_up_down_genes_GS(BP_Neuron, neurons_res_tauCortex_17m)

png("./results/figures/Neurons_DEGs_enrichment_BP_MF.png", width=18, height=8, units="in", res=300)
p1 <- dotplot(BP_Neuron, showCategory=16) + ggtitle("Neurons: BP Enrichment") + facet_grid(.~upDown)
p2 <- dotplot(MF_Neuron, showCategory=16) + ggtitle("Neurons: MF Enrichment") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()

# Treeplot
BP_simMatrix <- calculateSimMatrix(BP_Neuron$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")
BP_scores <- setNames(-log10(BP_Neuron$qvalue), BP_Neuron$ID)
BP_reducedTerms <- reduceSimMatrix(BP_simMatrix,
                                BP_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Neurons_DEGs_enrichment_BP_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(BP_reducedTerms)
dev.off()

MF_simMatrix <- calculateSimMatrix(MF_Neuron$ID,
                                orgdb="org.Mm.eg.db",
                                ont="MF",
                                method="Rel")
MF_scores <- setNames(-log10(MF_Neuron$qvalue), MF_Neuron$ID)
MF_reducedTerms <- reduceSimMatrix(MF_simMatrix,
                                MF_scores,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
png("./results/figures/Neurons_DEGs_enrichment_MF_Treeplot.png", width=10, height=8, units="in", res=300)
treemapPlot(MF_reducedTerms)
dev.off()

save(BP_Microglial, MF_Microglial, BP_Oligo, MF_Oligo, BP_Astro, MF_Astro, BP_Neuron, MF_Neuron, 
    file="./results/files/Cell_Type_Specific_Enrichment_Microg_Oligo_Astro_Neurons_Dec2024.RData")

# Saving Pathway Results in Excel files

BP_Microglial@result$Category <- "BP"
MF_Microglial@result$Category <- "MF"
Microglial_Pathways <- rbind(BP_Microglial@result, MF_Microglial@result)
Microglial_Pathways <- Microglial_Pathways[Microglial_Pathways$p.adjust < 0.05,]
dim(Microglial_Pathways) # 103  13
write.xlsx(Microglial_Pathways, file="./results/files/Microglial_Pathways.xlsx")

BP_Oligo@result$Category <- "BP"
MF_Oligo@result$Category <- "MF"
Oligo_Pathways <- rbind(BP_Oligo@result, MF_Oligo@result)
Oligo_Pathways <- Oligo_Pathways[Oligo_Pathways$p.adjust < 0.05,]
dim(Oligo_Pathways) # 114  13
write.xlsx(Oligo_Pathways, file="./results/files/Oligo_Pathways.xlsx")

BP_Astro@result$Category <- "BP"
MF_Astro@result$Category <- "MF"
Astro_Pathways <- rbind(BP_Astro@result, MF_Astro@result)
Astro_Pathways <- Astro_Pathways[Astro_Pathways$p.adjust < 0.05,]
dim(Astro_Pathways) # 121  13
write.xlsx(Astro_Pathways, file="./results/files/Astro_Pathways.xlsx")

BP_Neuron@result$Category <- "BP"
MF_Neuron@result$Category <- "MF"
Neuron_Pathways <- rbind(BP_Neuron@result, MF_Neuron@result)
Neuron_Pathways <- Neuron_Pathways[Neuron_Pathways$p.adjust < 0.05,]
dim(Neuron_Pathways) # 140  13
write.xlsx(Neuron_Pathways, file="./results/files/Neuron_Pathways.xlsx")

##
## Pathway analyses for Longitudinal DEGs ##
##

# Function to perform pathway analysis for a specific cell type and gene category
perform_pathway_analysis <- function(data, cell_type, category = NULL) {
  # Filter data for cell type
  cell_data <- data[data$Cell_Type == cell_type,]
  
  # If category is provided (sex-specific/dimorphic), filter accordingly
  if (!is.null(category)) {
    cell_data <- cell_data[cell_data$DEG.type == category,]
  }
  
  # Get gene list
  genes <- unique(cell_data$Gene.symbols)
  
  # Convert gene symbols to ENTREZ IDs
  gene_ids <- mapIds(org.Mm.eg.db,
                     keys = genes,
                     keytype = "SYMBOL",
                     column = "ENTREZID",
                     multiVals = "first")
  
  # Remove NAs
  gene_ids <- gene_ids[!is.na(gene_ids)]
  
  if(length(gene_ids) < 2) {
    return(NULL)
  }
  
  # Perform GO enrichment analysis
  ego_BP <- enrichGO(gene = gene_ids,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
  
  ego_MF <- enrichGO(gene = gene_ids,
                     OrgDb = org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
  
  return(list(BP = ego_BP, MF = ego_MF))
}

# Function to create visualization
create_enrichment_plot <- function(ego_result, title) {
  if(is.null(ego_result) || nrow(ego_result@result) == 0) {
    return(NULL)
  }
  
  # Get top 15 terms
  plotData <- head(ego_result@result, 15)
  
  # Create dot plot
  p <- ggplot(plotData, 
              aes(x = Count/ego_result@geneSets %>% lengths(), 
                  y = reorder(Description, Count/ego_result@geneSets %>% lengths()),
                  size = Count,
                  color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "purple", high = "lightpink") +
    theme_minimal() +
    labs(x = "Gene Ratio",
         y = NULL,
         title = title,
         color = "Adjusted p-value",
         size = "Gene count") +
    theme(axis.text.y = element_text(size = 8))
  
  return(p)
}

# Read and prepare data
data <- read.table("./data/longitudinal_degs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# List of cell types to analyze
cell_types <- c("Microglial cells", "Astrocytes", "Neurons", "Oligodendrocytes")
categories <- c("sex-dimorphic", "sex-neutral")

# Perform analysis for each cell type and category
results <- list()
plots <- list()

for(cell in cell_types) {
  for(cat in categories) {
    key <- paste(cell, cat, sep = "_")
    results[[key]] <- perform_pathway_analysis(data, cell, cat)
    
    if(!is.null(results[[key]])) {
      # Create plots for BP and MF
      plots[[paste0(key, "_BP")]] <- create_enrichment_plot(
        results[[key]]$BP,
        paste0(cell, " - ", cat, " - Biological Process"))
      
      plots[[paste0(key, "_MF")]] <- create_enrichment_plot(
        results[[key]]$MF,
        paste0(cell, " - ", cat, " - Molecular Function"))
    }
  }
}


# Function to create valid and unique sheet name
create_sheet_name <- function(cell, category, ontology) {
    # Shorten cell type names
    cell_short <- switch(cell,
        "Microglial cells" = "MG",
        "Astrocytes" = "AST",
        "Neurons" = "NEU",
        "Oligodendrocytes" = "OLG",
        gsub(" ", "", cell)
    )
    
    # Shorten category
    cat_short <- switch(category,
        "sex-dimorphic" = "DIM",
        "sex-neutral" = "NEU",
        category
    )
    
    # Create base name
    sheet_name <- paste0(cell_short, "_", cat_short, "_", ontology)
    
    # Ensure it's a valid Excel sheet name
    sheet_name <- gsub("[\\[\\]\\*\\?/\\\\]", "", sheet_name)
    
    return(sheet_name)
}

# Helper function to get direction information for genes
get_direction_info <- function(data, cell_type, category) {
    # Filter data for cell type and category
    cell_data <- data[data$Cell_Type == cell_type & data$DEG.type == category,]
    
    # Calculate direction based on log fold change
    # Using male logFC as reference for consistency with the paper
    direction <- sign(cell_data$Male.avg..logFC)
    
    # Name the vector with gene symbols for easy lookup
    names(direction) <- cell_data$Gene.symbols
    
    return(direction)
}

format_enrichment_results <- function(ego_result, direction_info = NULL) {
    if(is.null(ego_result) || nrow(ego_result@result) == 0) {
        return(NULL)
    }
    
    # Get basic results and convert to data frame
    results_df <- data.frame(
        ID = ego_result@result$ID,
        Description = ego_result@result$Description,
        GeneRatio = ego_result@result$GeneRatio,
        BgRatio = ego_result@result$BgRatio,
        pvalue = ego_result@result$pvalue,
        p.adjust = ego_result@result$p.adjust,
        qvalue = ego_result@result$qvalue,
        geneID = ego_result@result$geneID,
        Count = sapply(strsplit(ego_result@result$GeneRatio, "/"), 
                       function(x) as.numeric(x[1])),
        stringsAsFactors = FALSE
    )
    
    # Calculate direction counts
    if(!is.null(direction_info)) {
        results_df$up <- 0
        results_df$down <- 0
        results_df$upDown <- "equal"
        
        # Process each term
        for(i in 1:nrow(results_df)) {
            genes <- unlist(strsplit(results_df$geneID[i], "/"))
            gene_directions <- direction_info[genes]
            up_count <- sum(gene_directions > 0, na.rm = TRUE)
            down_count <- sum(gene_directions < 0, na.rm = TRUE)
            
            results_df$up[i] <- up_count
            results_df$down[i] <- down_count
            
            if(up_count > down_count) {
                results_df$upDown[i] <- "up"
            } else if(down_count > up_count) {
                results_df$upDown[i] <- "down"
            }
        }
    } else {
        results_df$up <- 0
        results_df$down <- 0
        results_df$upDown <- "equal"
    }
    
    # Add Category
    results_df$Category <- ifelse(grepl("^GO:BP", results_df$ID), "BP", "MF")
    
    return(results_df)
}

# Export the results to an Excel table with one sheet per cell type
export_to_excel <- function(results, cell_types, categories, filename) {
    # Create workbook
    wb <- createWorkbook()
    
    # Create styles
    headerStyle <- createStyle(
        fontName = "Arial",
        fontSize = 10,
        textDecoration = "bold",
        border = "bottom",
        borderStyle = "thin",
        halign = "left"
    )
    
    dataStyle <- createStyle(
        fontName = "Arial",
        fontSize = 10,
        halign = "left"
    )
    
    numStyle <- createStyle(
        fontName = "Arial",
        fontSize = 10,
        numFmt = "0.00E+00"
    )
    
    # Keep track of created sheet names
    used_sheet_names <- c()
    
    # Process each cell type and category
    for(cell in cell_types) {
        for(cat in categories) {
            key <- paste(cell, cat, sep = "_")
            
            if(!is.null(results[[key]])) {
                # Process BP results
                if(!is.null(results[[key]]$BP) && 
                   nrow(results[[key]]$BP@result) > 0) {
                    
                    bp_results <- format_enrichment_results(
                        results[[key]]$BP,
                        direction_info = get_direction_info(data, cell, cat)
                    )
                    
                    if(!is.null(bp_results) && nrow(bp_results) > 0) {
                        sheet_name <- create_sheet_name(cell, cat, "BP")
                        
                        # Verify sheet name is unique
                        if(sheet_name %in% used_sheet_names) {
                            stop(paste("Duplicate sheet name generated:", sheet_name))
                        }
                        used_sheet_names <- c(used_sheet_names, sheet_name)
                        
                        addWorksheet(wb, sheet_name)
                        writeData(wb, sheet_name, bp_results)
                        
                        # Apply styles
                        addStyle(wb, sheet_name, headerStyle, rows = 1, 
                               cols = 1:ncol(bp_results))
                        
                        # Apply data styles
                        for(i in 2:(nrow(bp_results) + 1)) {
                            addStyle(wb, sheet_name, dataStyle, rows = i, 
                                   cols = 1:ncol(bp_results))
                        }
                        
                        # Apply numeric style to p-value columns
                        pvalue_cols <- which(colnames(bp_results) %in% 
                                          c("pvalue", "p.adjust", "qvalue"))
                        for(col in pvalue_cols) {
                            addStyle(wb, sheet_name, numStyle, 
                                   rows = 2:(nrow(bp_results) + 1), 
                                   cols = col, 
                                   stack = TRUE)
                        }
                        
                        setColWidths(wb, sheet_name, cols = 1:ncol(bp_results), 
                                   widths = "auto")
                    }
                }
                
                # Process MF results
                if(!is.null(results[[key]]$MF) && 
                   nrow(results[[key]]$MF@result) > 0) {
                    
                    mf_results <- format_enrichment_results(
                        results[[key]]$MF,
                        direction_info = get_direction_info(data, cell, cat)
                    )
                    
                    if(!is.null(mf_results) && nrow(mf_results) > 0) {
                        sheet_name <- create_sheet_name(cell, cat, "MF")
                        
                        # Verify sheet name is unique
                        if(sheet_name %in% used_sheet_names) {
                            stop(paste("Duplicate sheet name generated:", sheet_name))
                        }
                        used_sheet_names <- c(used_sheet_names, sheet_name)
                        
                        addWorksheet(wb, sheet_name)
                        writeData(wb, sheet_name, mf_results)
                        
                        # Apply styles
                        addStyle(wb, sheet_name, headerStyle, rows = 1, 
                               cols = 1:ncol(mf_results))
                        
                        # Apply data styles
                        for(i in 2:(nrow(mf_results) + 1)) {
                            addStyle(wb, sheet_name, dataStyle, rows = i, 
                                   cols = 1:ncol(mf_results))
                        }
                        
                        # Apply numeric style to p-value columns
                        pvalue_cols <- which(colnames(mf_results) %in% 
                                          c("pvalue", "p.adjust", "qvalue"))
                        for(col in pvalue_cols) {
                            addStyle(wb, sheet_name, numStyle, 
                                   rows = 2:(nrow(mf_results) + 1), 
                                   cols = col, 
                                   stack = TRUE)
                        }
                        
                        setColWidths(wb, sheet_name, cols = 1:ncol(mf_results), 
                                   widths = "auto")
                    }
                }
            }
        }
    }
    
    # Only save if we have any worksheets
    if(length(wb$worksheets) > 0) {
        saveWorkbook(wb, filename, overwrite = TRUE)
        cat("Results exported successfully to", filename, "\n")
    } else {
        cat("No results to export\n")
    }
}

# Export results
export_to_excel(results, 
                cell_types, 
                categories, 
                "./results/files/longitudinal_pathway_analysis.xlsx")



check_all_analyses <- function(data, cell_types, categories) {
    results_summary <- data.frame(
        cell_type = character(),
        category = character(),
        n_genes = numeric(),
        n_BP_terms = numeric(),
        n_MF_terms = numeric(),
        stringsAsFactors = FALSE
    )
    
    for(cell in cell_types) {
        for(cat in categories) {
            cat("\nAnalyzing", cell, "-", cat, "\n")
            
            # Get results
            res <- perform_pathway_analysis(data, cell, cat)
            
            # Get gene count
            cell_data <- data[data$Cell_Type == cell & data$DEG.type == cat,]
            n_genes <- length(unique(cell_data$Gene.symbols))
            
            # Add to summary
            results_summary <- rbind(results_summary, data.frame(
                cell_type = cell,
                category = cat,
                n_genes = n_genes,
                n_BP_terms = if(!is.null(res)) nrow(res$BP@result) else 0,
                n_MF_terms = if(!is.null(res)) nrow(res$MF@result) else 0,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    return(results_summary)
}


# Run summary analysis
summary_results <- check_all_analyses(data, cell_types, categories)



# Print summary
cat("\nOverall Analysis Summary:\n")
print(summary_results)
#         cell_type      category n_genes n_BP_terms n_MF_terms
#1 Microglial cells sex-dimorphic      45       1075        145
#2 Microglial cells   sex-neutral     299       3719        581
#3       Astrocytes sex-dimorphic      28        894         90
#4       Astrocytes   sex-neutral      92       1968        288
#5          Neurons sex-dimorphic       2         55         10
#6          Neurons   sex-neutral      96       2270        264
#7 Oligodendrocytes sex-dimorphic      24        869        101
#8 Oligodendrocytes   sex-neutral     117       2361        311



process_enrichment_plots <- function(results, cell_types, categories) {
    # Initialize plot list
    plots <- list()
    
    # Create plots for each combination
    for(cell in cell_types) {
        for(cat in categories) {
            key <- paste(cell, cat, sep = "_")
            if(!is.null(results[[key]])) {
                # Process BP results
                if(!is.null(results[[key]]$BP) && 
                   nrow(results[[key]]$BP@result) > 0) {
                    plot_name <- paste0(gsub(" ", "_", cell), "_", cat, "_BP")
                    plots[[plot_name]] <- create_enrichment_plot(
                        results[[key]]$BP,
                        paste(cell, "-", cat, "- Biological Process")
                    )
                }
                
                # Process MF results
                if(!is.null(results[[key]]$MF) && 
                   nrow(results[[key]]$MF@result) > 0) {
                    plot_name <- paste0(gsub(" ", "_", cell), "_", cat, "_MF")
                    plots[[plot_name]] <- create_enrichment_plot(
                        results[[key]]$MF,
                        paste(cell, "-", cat, "- Molecular Function")
                    )
                }
            }
        }
    }
    
    return(plots)
}



display_enrichment_plots = function(plots, cell_type) {
    # Filter plots for the specific cell type
    cell_plots <- plots[grep(paste0("^", gsub(" ", "_", cell_type)), names(plots))]
    
    # Remove NULL plots
    cell_plots <- cell_plots[!sapply(cell_plots, is.null)]
    
    if(length(cell_plots) == 0) {
        message(paste("No significant enrichment found for", cell_type))
        return(NULL)
    }
    
    # Calculate grid dimensions
    n_plots <- length(cell_plots)
    n_cols <- min(2, n_plots)
    n_rows <- ceiling(n_plots / n_cols)
    
    # Create layout
    if(n_plots == 1) {
        return(cell_plots[[1]])
    } else {
        # Arrange multiple plots
        combined_plot <- gridExtra::arrangeGrob(
            grobs = cell_plots,
            ncol = n_cols,
            top = textGrob(paste("Pathway Enrichment -", cell_type),
                           gp = gpar(fontsize = 12, fontface = "bold"))
        )
        return(combined_plot)
    }
}


# Function to save all plots with proper pagination
save_all_plots <- function(results, cell_types, categories, filename = "./results/figures/pathway_enrichment_plots.pdf") {
    # Process all plots
    plots <- process_enrichment_plots(results, cell_types, categories)
    
    # Save plots to PDF
    pdf(filename, width = 12, height = 8)
    
    # Display plots for each cell type
    for(cell_type in cell_types) {
        # Start a new page for each cell type
        grid.newpage()
        
        # Filter and process plots for this cell type
        plot_obj <- display_enrichment_plots(plots, cell_type)
        if(!is.null(plot_obj)) {
            # Draw the plot
            grid::grid.draw(plot_obj)
        }
    }
    
    dev.off()
    
    message("Plots saved to ", filename)
}

save_all_plots(results, cell_types, categories)