## Differential Expression with Sex*Genotype Interaction ##

 library(Seurat) # Seurat_4.3.0 # SeuratObject_4.1.4
library(openxlsx) # openxlsx_4.2.5.2
library(ggplot2)# ggplot2_3.5.1
library(edgeR) # edgeR_4.2.0
library(limma) # limma_3.60.3
library(scran) # scran_1.32.0
library(readxl) # readxl_1.4.3
library(SingleCellExperiment) # SingleCellExperiment_1.26.0
set.seed(1)
options(width=160)

load("./data/sctrans2.RData")

# Set sample identities to celltype/condition combinations
orig.ident = sctrans2$orig.ident
table(orig.ident)
# orig.ident
# 01WTM 02TGF 03WTM 04WTM 05WTM 06TGF 07TGM 08TGM 09TGF 10TGM 11TGM 12WTF 13WTF
#  2509  1946  1770  1870  1394  1347  1613  1956  2013  1375  2044  1359  1109
# 14TGF 15TGF 16WTF
#  1405  2127  2353

sctrans2$orig.ident = paste(sctrans2$classint, sctrans2$conditions)

Idents(sctrans2) <- paste(sctrans2$classint, sctrans2$conditions)
table(Idents(sctrans2))
sctrans2$cond <- sctrans2@meta.data$conditions


grub <- sctrans2
str(grub@meta.data)
grub$cond <- sub(" ", "_", grub$cond)
table(grub$cond)
# TG_F TG_M WT_F WT_M
# 8838 6988 4821 7543

grub$classint <- sub(" ", "_", grub$classint)
grub$classint <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", grub$classint)
grub$classint <- as.factor(grub$classint)
table(grub$classint)
#                     Astrocyte               Endothelial_cell                 Ependymal_cell                     Macrophage                Microglial_cell
#                          5299                           5086                            105                            307                           8263
#                    Mural_cell                         Neuron                Oligodendrocyte Oligodendrocyte_precursor_cell
#                           188                           2183                           5505                           1254

grub$diag_cell_sex <- paste0(grub$cond, "_", grub$classint, "_17m")
length(unique(grub$diag_cell_sex)) # 36
Idents(grub) <- "diag_cell_sex"
table(Idents(grub))
grub$sex <- sapply(strsplit(as.character(grub$stim), " "), '[', 2)
grub$sex <- as.factor(grub$sex)
table(grub$sex)
#     F     M
# 13659 14531
grub$status <- sapply(strsplit(as.character(grub$cond), "_"), '[', 1)
grub$status <- as.factor(grub$status)
table(grub$status)
#    TG    WT
# 15826 12364

DefaultAssay(grub) <- "SCT" # Previously, "RNA"
grub$diag_cell_sex <- as.factor(grub$diag_cell_sex)
str(grub@meta.data)

design <- model.matrix(~0 + sex:diag_cell_sex, data = grub@meta.data) # https://support.bioconductor.org/p/92225/
dim(design) # 28190    72

colnames(design) <- gsub("diag_cell_sex", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design)

dge <-  as.matrix(GetAssayData(grub, slot="counts"))
# Warning message: In asMethod(object) : sparse->dense coercion: allocating vector of size 3.7 GiB
dim(dge) # 17681 28190
dge <- DGEList(counts= dge)
keep <- filterByExpr(dge, design)
summary(keep)
#    Mode   FALSE    TRUE
# logical   17621      60
dge <- dge[keep,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMMwsp")

#Run Voom and make contrasts
vm <- voom(dge, design, plot = TRUE)
# Coefficients not estimable: TG M WT F WT M
# Partial NA coefficients for 46 probe(s)
fit <- lmFit(vm, design)
# Coefficients not estimable: M
head(coef(fit))
head(coef(fit))[grep("Neuron", head(coef(fit)))]

# create design matrix
colnames(design) <- gsub(":", "_", colnames(design))
colnames(design) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", colnames(design))
# fit_copy <- fit
fit <- fit_copy
colnames(fit$design) <- gsub(":", "_", colnames(fit$design))
colnames(fit$design) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", colnames(fit$design))
colnames(fit$coefficients) <- gsub(":", "_", colnames(fit$coefficients))
colnames(fit$coefficients) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", colnames(fit$coefficients))
colnames(fit$stdev.unscaled) <- gsub(":", "_", colnames(fit$stdev.unscaled))
colnames(fit$stdev.unscaled) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", colnames(fit$stdev.unscaled))
row.names(fit$cov.coefficients) <- gsub(":", "_", row.names(fit$cov.coefficients))
row.names(fit$cov.coefficients) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", row.names(fit$cov.coefficients))
colnames(fit$cov.coefficients) <- gsub(":", "_", colnames(fit$cov.coefficients))
colnames(fit$cov.coefficients) <- sub("Oligodendrocyte_precursor cell", "Oligodendrocyte_precursor_cell", colnames(fit$cov.coefficients))

# create contrast matrix
contrasts.matrix <- makeContrasts(AstM= (M_TG_M_Astrocyte_17m - M_WT_M_Astrocyte_17m) - (F_TG_F_Astrocyte_17m - F_WT_F_Astrocyte_17m),
                                  AstF= (F_TG_F_Astrocyte_17m - F_WT_F_Astrocyte_17m) - (M_TG_M_Astrocyte_17m - M_WT_M_Astrocyte_17m),
                                  EndM= (M_TG_M_Endothelial_cell_17m - M_WT_M_Endothelial_cell_17m) - (F_TG_F_Endothelial_cell_17m - F_WT_F_Endothelial_cell_17m),
                                  EndF = (F_TG_F_Endothelial_cell_17m - F_WT_F_Endothelial_cell_17m) - (M_TG_M_Endothelial_cell_17m - M_WT_M_Endothelial_cell_17m),
                                  MicM= (M_TG_M_Microglial_cell_17m - M_WT_M_Microglial_cell_17m) - (F_TG_F_Microglial_cell_17m - F_WT_F_Microglial_cell_17m),
                                  MicF = (F_TG_F_Microglial_cell_17m - F_WT_F_Microglial_cell_17m) - (M_TG_M_Microglial_cell_17m - M_WT_M_Microglial_cell_17m),
                                  OliM= (M_TG_M_Oligodendrocyte_17m - M_WT_M_Oligodendrocyte_17m) - (F_TG_F_Oligodendrocyte_17m - F_WT_F_Oligodendrocyte_17m),
                                  OliF= (F_TG_F_Oligodendrocyte_17m - F_WT_F_Oligodendrocyte_17m) - (M_TG_M_Oligodendrocyte_17m - M_WT_M_Oligodendrocyte_17m),
                                  OpcM= (M_TG_M_Oligodendrocyte_precursor_cell_17m - M_WT_M_Oligodendrocyte_precursor_cell_17m) - (F_TG_F_Oligodendrocyte_precursor_cell_17m - F_WT_F_Oligodendrocyte_precursor_cell_17m),
                                  OpcF= (F_TG_F_Oligodendrocyte_precursor_cell_17m - F_WT_F_Oligodendrocyte_precursor_cell_17m) - (M_TG_M_Oligodendrocyte_precursor_cell_17m - M_WT_M_Oligodendrocyte_precursor_cell_17m),
                                  EpeM= (M_TG_M_Ependymal_cell_17m - M_WT_M_Ependymal_cell_17m) - (F_TG_F_Ependymal_cell_17m - F_WT_F_Ependymal_cell_17m),
                                  EpeF= (F_TG_F_Ependymal_cell_17m - F_WT_F_Ependymal_cell_17m) - (M_TG_M_Ependymal_cell_17m - M_WT_M_Ependymal_cell_17m),
                                  MacM= (M_TG_M_Macrophage_17m - M_WT_M_Macrophage_17m) - (F_TG_F_Macrophage_17m - F_WT_F_Macrophage_17m),
                                  MacF= (F_TG_F_Macrophage_17m - F_WT_F_Macrophage_17m) - (M_TG_M_Macrophage_17m - M_WT_M_Macrophage_17m),
                                  MrlM= (M_TG_M_Mural_cell_17m - M_WT_M_Mural_cell_17m) - (F_TG_F_Mural_cell_17m - F_WT_F_Mural_cell_17m),
                                  MrlF= (F_TG_F_Mural_cell_17m - F_WT_F_Mural_cell_17m) - (M_TG_M_Mural_cell_17m - M_WT_M_Mural_cell_17m),
                                  NeuM= (M_TG_M_Neuron_17m - M_WT_M_Neuron_17m) - (F_TG_F_Neuron_17m - F_WT_F_Neuron_17m),
                                  NeuF= (F_TG_F_Neuron_17m - F_WT_F_Neuron_17m) - (M_TG_M_Neuron_17m - M_WT_M_Neuron_17m),
                                  levels = colnames(design))
 
fit <- contrasts.fit(fit, contrasts = contrasts.matrix) 
fit <- eBayes(fit)

#Run loop to get results 
delist_key <- c("AstM", "EndM",  "MicM", "MicM", "OliM" , "OpcM", "EpeM", "MacM", "MrlM", "NeuM", 
               "AstF", "EndF", "MicF", "MicF", "OliF", "OpcF", "EpeF", "MacF", "MrlF", "NeuF")
library(stringr)
markers2 <- NULL
markers3 <- NULL
for(key in delist_key){
  print(key)
  key1 <- str_sub(key,-4,-2) #Cell_type
  key2 <- gsub(key1,"",key) #Sex
  markers <- topTable(fit, coef= key, sort.by = "logFC", number = Inf, adjust.method = "BH")  
  markers$group <- key
  markers$Sex <- key2
  markers$cell_type <- key1
  markers$gene <- rownames(markers)
  markers$dir <- ifelse(markers$logFC < 0, "neg","pos")
  colnames(markers)[c(1,4,5)] <- c("avg_logFC", "p_val", "p_val_adj")
  markers2 <- rbind(markers2, markers) #no thresholds
  markers <- subset(markers, p_val_adj < 0.05 & abs(avg_logFC) > 0.25)  
  markers3 <- rbind(markers3, markers) 
}
table(markers3$group, markers3$dir)
#       neg pos
#  AstF  10   8
#  AstM   8  10
#  EndF  20   4
#  EndM   4  20
#  MacF   6   4
#  MacM   4   6
#  MicF  38  10
#  MicM  10  38
#  MrlF   0   5
#  MrlM   5   0
#  NeuF  12  12
#  NeuM  12  12
#  OliF  17   5
#  OliM   5  17
#  OpcF   2   1
#  OpcM   1   2
  
male.markersOG <- markers2[markers2$Sex == "M",]
dim(male.markersOG) # 600  11
female.markersOG <- markers2[markers2$Sex != "M",]
dim(female.markersOG) # 600  11
write.csv(male.markersOG, file="./results/files/ThyTau22_17m_limma_DEGs_SexInteraction.csv") #lnf= limma no filter
write.csv(female.markersOG, file="./results/files/ThyTau22_17m_limma_DEGs_SexInteraction_Female.csv") #lnf= limma no filter

nrow(male.markersOG[male.markersOG$p_val_adj < 0.05,]) # 288
write.csv(male.markersOG[male.markersOG$p_val_adj < 0.05,], file="./results/files/ThyTau22_limma_DEGs_SexInteraction_Significant.csv") #lnf= limma no filter
table(male.markersOG$group, male.markersOG$dir)
#       neg pos
#  AstM  35  25
#  EndM  11  49
#  EpeM  17  43
#  MacM  19  41
#  MicM  20 100
#  MrlM  43  17
#  NeuM  27  33
#  OliM   9  51
#  OpcM  17  43
save(fit, fit_copy, file="./results/files/ThyTau22_17m_limma_DEGs_SexInteraction.RData")

#
# Overlap of Seurat with edgeR (Sex*Genotype Interaction)
#

male_int <- read.csv("./results/files/ThyTau22_17m_limma_DEGs_SexInteraction.csv", header=T, stringsAsFactors=F)
dim(male_int) # 600  12
length(unique(male_int$gene)) # 60
male_int <- male_int[male_int$p_val_adj < 0.05,]
dim(male_int) # 288  12
male_int$X <- NULL
male_int$group <- NULL
male_int$Sex <- NULL

table(male_int$cell_type) # min. FDR in EpeM is 0.1294046
# Ast End Mac Mic Mrl Neu Oli Opc
#  27  49  10 108   5  37  49   3

load("./results/files/TauCortex17m_DEGs_All9CTs_50percent_minabslog.RData")

table(astrocytes_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1               12                2                9
table(microglial_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1                6                9               14
table(neurons_res_tauCortex_17m$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#               11                4                7
table(oligodendrocyte_res_tauCortex_17m$DEG.type)
# gender-dimorphic    gender-shared    male-specific
#                7                3               11
table(endothelial_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               2               10                2                5
table(OPC_res_tauCortex_17m$DEG.type)
# gender-shared male-specific
#             2             2
table(mural_res_tauCortex_17m$DEG.type)
# gender-shared male-specific
#             1             1
table(macrophage_res_tauCortex_17m$DEG.type)
# gender-dimorphic    male-specific
#                2                1
table(ependymal_res_tauCortex_17m$DEG.type)
# gender-shared
#             1

Astro <- male_int[male_int$cell_type == "Ast",]
Astro$cell_type <- NULL
Astro <- left_join(Astro, astrocytes_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Astro$Status <- ifelse(is.na(Astro$DEG.type), "unique", "common")
Astro <- Astro[order(Astro$Status, decreasing = T),]
Astro$cell_type <- "Astrocytes"
dim(Astro) # 27 15

Microg <- male_int[male_int$cell_type == "Mic",]
Microg$cell_type <- NULL
Microg <- left_join(Microg, microglial_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Microg$Status <- ifelse(is.na(Microg$DEG.type), "unique", "common")
Microg <- Microg[order(Microg$Status, decreasing = T),]
Microg$cell_type <- "Microglial cells"
dim(Microg) # 108  15

Neurons <- male_int[male_int$cell_type == "Neu",]
Neurons$cell_type <- NULL
Neurons <- left_join(Neurons, neurons_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Neurons$Status <- ifelse(is.na(Neurons$DEG.type), "unique", "common")
Neurons <- Neurons[order(Neurons$Status, decreasing = T),]
Neurons$cell_type <- "Neurons"
dim(Neurons) # 37 15

Oligo <- male_int[male_int$cell_type == "Oli",]
Oligo$cell_type <- NULL
Oligo <- left_join(Oligo, oligodendrocyte_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Oligo$Status <- ifelse(is.na(Oligo$DEG.type), "unique", "common")
Oligo <- Oligo[order(Oligo$Status, decreasing = T),]
Oligo$cell_type <- "Oligodendrocytes"
dim(Oligo) # 49 15

Endo <- male_int[male_int$cell_type == "End",]
Endo$cell_type <- NULL
Endo <- left_join(Endo, endothelial_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Endo$Status <- ifelse(is.na(Endo$DEG.type), "unique", "common")
Endo <- Endo[order(Endo$Status, decreasing = T),]
Endo$cell_type <- "Endothelial cells"
dim(Endo) # 49 15

OPCs <- male_int[male_int$cell_type == "Opc",]
OPCs$cell_type <- NULL
OPCs <- left_join(OPCs, OPC_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
OPCs$Status <- ifelse(is.na(OPCs$DEG.type), "unique", "common")
OPCs <- OPCs[order(OPCs$Status, decreasing = T),]
OPCs$cell_type <- "Oligodendrocyte precursor cells"
dim(OPCs) # 3 15

Mural <- male_int[male_int$cell_type == "Mrl",]
Mural$cell_type <- NULL
Mural <- left_join(Mural, mural_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Mural$Status <- ifelse(is.na(Mural$DEG.type), "unique", "common")
Mural <- Mural[order(Mural$Status, decreasing = T),]
Mural$cell_type <- "Mural cells"
dim(Mural) # 5 15

Macro <- male_int[male_int$cell_type == "Mac",]
Macro$cell_type <- NULL
Macro <- left_join(Macro, macrophage_res_tauCortex_17m, by=c("gene"="Gene.symbols"))
Macro$Status <- ifelse(is.na(Macro$DEG.type), "unique", "common")
Macro <- Macro[order(Macro$Status, decreasing = T),]
Macro$cell_type <- "Macrophages"
dim(Macro) # 10 15

seurat_edgeR <- rbind(Astro, Microg, Neurons, Oligo, Endo, OPCs, Mural, Macro)
dim(seurat_edgeR) # 288  15
seurat_edgeR <- seurat_edgeR %>% select(gene, avg_logFC, everything())
colnames(seurat_edgeR)[8] <- "direction.males"
seurat_edgeR <- unique(seurat_edgeR)
dim(seurat_edgeR) # 234  15
head(sort(unique(seurat_edgeR$Female.FDR)))
# 7.556180e-142 4.890386e-119  1.951927e-65  2.360167e-49  5.508267e-42  7.740723e-36
head(sort(unique(seurat_edgeR$Male.FDR)))
# 0.000000e+00 1.823309e-266 3.880433e-255 1.614712e-202 1.314476e-186 6.591949e-176
seurat_edgeR$Male.FDR <- ifelse(seurat_edgeR$Male.FDR == 0, 1.823309e-266, seurat_edgeR$Male.FDR)
dim(seurat_edgeR) # 234  15
write.table(seurat_edgeR, file="./results/files/ThyTau22_17m_limma_SexInteraction_and_Seurat_DEGs.txt", sep="\t", row.names=F, quote=F)