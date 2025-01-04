## Overlap between Tg2576, TauCortex, and Grubman data (4 CTs) ##

rm(list = ls())
library(dplyr) # dplyr_1.1.4
library(ggplot2) # ggplot2_3.5.1
library(ggVennDiagram) # ggVennDiagram_1.5.2
options(stringsAsFactors = FALSE)
options(width=120)
set.seed(1)


# Thy-Tau22 (17m)
# load("./results/files/TauCortex17m_DEGs_4CT_astro_oligo_microg_neuron.RData")
load("./results/files/TauCortex17m_DEGs_4CT_astro_oligo_microg_neuron_50percent_minabslog.RData")
ls()
# [1] "astrocytes_res_tauCortex_17m"      "microglial_res_tauCortex_17m"      "neurons_res_tauCortex_17m"
# [4] "oligodendrocyte_res_tauCortex_17m"
DEGs_TauCortex_17m <- rbind(astrocytes_res_tauCortex_17m, microglial_res_tauCortex_17m, neurons_res_tauCortex_17m, oligodendrocyte_res_tauCortex_17m)
dim(DEGs_TauCortex_17m) # 97  7
DEGs_TauCortex_17m$Gene.symbols <- toupper(DEGs_TauCortex_17m$Gene.symbols)
gender_specific_TauCortex_17m <- DEGs_TauCortex_17m[DEGs_TauCortex_17m$DEG.type == "female-specific" | DEGs_TauCortex_17m$DEG.type == "male-specific",]
dim(gender_specific_TauCortex_17m) # 43  7
gender_neutral_TauCortex_17m <- DEGs_TauCortex_17m[DEGs_TauCortex_17m$DEG.type == "gender-shared",]
dim(gender_neutral_TauCortex_17m) # 18  7
gender_dimorphic_TauCortex_17m <- DEGs_TauCortex_17m[DEGs_TauCortex_17m$DEG.type == "gender-dimorphic",]
dim(gender_dimorphic_TauCortex_17m) # 36 7


# Thy-Tau22
load("/home/m.ali/Projects/UL/mice_snRNAseq_cortex/enrico_script/TauCortex_DEGs_4CT_astro_oligo_microg_neuron.RData")
ls()
DEGs_TauCortex <- rbind(astrocytes_res_tauCortex, microglial_res_tauCortex, neurons_res_tauCortex, oligodendrocyte_res_tauCortex)
dim(DEGs_TauCortex) # 145  7
DEGs_TauCortex$Gene.symbols <- toupper(DEGs_TauCortex$Gene.symbols)
gender_specific_TauCortex <- DEGs_TauCortex[DEGs_TauCortex$DEG.type == "female-specific" | DEGs_TauCortex$DEG.type == "male-specific",]
dim(gender_specific_TauCortex) # 109  7
gender_neutral_TauCortex <- DEGs_TauCortex[DEGs_TauCortex$DEG.type == "gender-shared",]
dim(gender_neutral_TauCortex) # 35  7
gender_dimorphic_TauCortex <- DEGs_TauCortex[DEGs_TauCortex$DEG.type == "gender-dimorphic",]
dim(gender_dimorphic_TauCortex) # 1 7


# Tg2576
load("/home/m.ali/Projects/UL/published_snRNAseq/Grubman/mm_Cortex_v3_DEGs_Poisson_logFC0.RData")
ls()
DEGs_Tg2576 <- rbind(astrocytes_res, microglial_res, neurons_res, oligodendrocyte_res)
dim(DEGs_Tg2576) # 141   7
gender_specific_Tg2576 <- DEGs_Tg2576[DEGs_Tg2576$DEG.type == "female-specific" | DEGs_Tg2576$DEG.type == "male-specific",]
dim(gender_specific_Tg2576) # 48  7
gender_neutral_Tg2576 <- DEGs_Tg2576[DEGs_Tg2576$DEG.type == "gender-shared",]
dim(gender_neutral_Tg2576) # 90  7
gender_dimorphic_Tg2576 <- DEGs_Tg2576[DEGs_Tg2576$DEG.type == "gender-dimorphic",]
dim(gender_dimorphic_Tg2576) # 3 7


# Grubman et al
load("/home/m.ali/Projects/UL/published_snRNAseq/Grubman/Grubman_DEGs_4CT_astro_oligo_microg_neuron.RData")
ls()
DEGs_Grubman <- rbind(astrocytes_res_grubman, microglial_res_grubman, neurons_res_grubman, oligodendrocyte_res_grubman)
dim(DEGs_Grubman) # 364   7
gender_specific_Grubman <- DEGs_Grubman[DEGs_Grubman$DEG.type == "female-specific" | DEGs_Grubman$DEG.type == "male-specific",]
dim(gender_specific_Grubman) # 234   7
gender_neutral_Grubman <- DEGs_Grubman[DEGs_Grubman$DEG.type == "gender-shared",]
dim(gender_neutral_Grubman) # 111   7
gender_dimorphic_Grubman <- DEGs_Grubman[DEGs_Grubman$DEG.type == "gender-dimorphic",]
dim(gender_dimorphic_Grubman) # 19  7


library(UpSetR)
overlap = list(Human_AD = as.character(unique(DEGs_Grubman$Gene.symbols)), Tg2576 = unique(DEGs_Tg2576$Gene.symbols), ThyTau22 = unique(DEGs_TauCortex$Gene.symbols), ThyTau22_17m = unique(DEGs_TauCortex_17m$Gene.symbols))
png("./results/figures/Upset_Overlapping_DEGs_4datasets_4CTs_v3.png", width=12, height=10, units="in", res=300)
upset(fromList(overlap), sets=c("Human_AD", "Tg2576", "ThyTau22", "ThyTau22_17m"), order.by = "freq", text.scale = 3.5, point.size = 4) # order by "degree"
dev.off()

# 4 separate venn diagrams of cell-type-sepecific DEGs overlap between 4 different datasets
microglial_res_tauCortex_17m$Gene.symbols <- toupper(microglial_res_tauCortex_17m$Gene.symbols)
microglial_res_tauCortex$Gene.symbols <- toupper(microglial_res_tauCortex$Gene.symbols)

# library(venn)
overlap = list(Human_AD = as.character(unique(microglial_res_grubman$Gene.symbols)), Tg2576 = unique(microglial_res$Gene.symbols), ThyTau22 = unique(microglial_res_tauCortex$Gene.symbols), ThyTau22_17m = unique(microglial_res_tauCortex_17m$Gene.symbols))
png("./results/figures/Venn_Overlapping_DEGs_4datasets_Microglial_Cells.png", width=8, height=6, units="in", res=300)
ggVennDiagram(overlap,
        category.names = c("Human_AD","Tg2576","ThyTau22_7m", "ThyTau22_17m"), label = "count", label_size = 6, set_size = 6) + scale_x_continuous(expand = expansion(mult = .2)) + theme(legend.position = "none") + scale_fill_distiller(palette = "RdBu")
dev.off()

microg_common <- inner_join(microglial_res_grubman, microglial_res, by="Gene.symbols")
microg_common <- inner_join(microg_common, microglial_res_tauCortex_17m, by="Gene.symbols")
# HSP90AB1

oligodendrocyte_res_tauCortex_17m$Gene.symbols <- toupper(oligodendrocyte_res_tauCortex_17m$Gene.symbols)
oligodendrocyte_res_tauCortex$Gene.symbols <- toupper(oligodendrocyte_res_tauCortex$Gene.symbols)

overlap = list(Human_AD = as.character(unique(oligodendrocyte_res_grubman$Gene.symbols)), Tg2576 = unique(oligodendrocyte_res$Gene.symbols), ThyTau22 = unique(oligodendrocyte_res_tauCortex$Gene.symbols), ThyTau22_17m = unique(oligodendrocyte_res_tauCortex_17m$Gene.symbols))
png("./results/figures/Venn_Overlapping_DEGs_4datasets_Oligodendrocytes.png", width=8, height=6, units="in", res=300)
ggVennDiagram(overlap,
        category.names = c("Human_AD","Tg2576","ThyTau22_7m", "ThyTau22_17m"), label = "count", label_size = 6, set_size = 6) + scale_x_continuous(expand = expansion(mult = .2)) + theme(legend.position = "none") + scale_fill_distiller(palette = "RdBu")
dev.off()

oligo_common <- inner_join(oligodendrocyte_res_grubman, oligodendrocyte_res, by="Gene.symbols")
oligo_common <- inner_join(oligo_common, oligodendrocyte_res_tauCortex, by="Gene.symbols")
oligo_common <- inner_join(oligo_common, oligodendrocyte_res_tauCortex_17m, by="Gene.symbols")
# MALAT1 and MBP


neurons_res_tauCortex_17m$Gene.symbols <- toupper(neurons_res_tauCortex_17m$Gene.symbols)
neurons_res_tauCortex$Gene.symbols <- toupper(neurons_res_tauCortex$Gene.symbols)

overlap = list(Human_AD = as.character(unique(neurons_res_grubman$Gene.symbols)), Tg2576 = unique(neurons_res$Gene.symbols), ThyTau22 = unique(neurons_res_tauCortex$Gene.symbols), ThyTau22_17m = unique(neurons_res_tauCortex_17m$Gene.symbols))
png("./results/figures/Venn_Overlapping_DEGs_4datasets_Neurons.png", width=8, height=6, units="in", res=300)
ggVennDiagram(overlap,
        category.names = c("Human_AD","Tg2576","ThyTau22_7m", "ThyTau22_17m"), label = "count", label_size = 6, set_size = 6) + scale_x_continuous(expand = expansion(mult = .2)) + theme(legend.position = "none") + scale_fill_distiller(palette = "RdBu")
dev.off()

neuron_common <- inner_join(neurons_res_grubman, neurons_res_tauCortex, by="Gene.symbols")
neuron_common <- inner_join(neuron_common, neurons_res_tauCortex_17m, by="Gene.symbols")
# ATP1B1


astrocytes_res_tauCortex_17m$Gene.symbols <- toupper(astrocytes_res_tauCortex_17m$Gene.symbols)
astrocytes_res_tauCortex$Gene.symbols <- toupper(astrocytes_res_tauCortex$Gene.symbols)

overlap = list(Human_AD = as.character(unique(astrocytes_res_grubman$Gene.symbols)), Tg2576 = unique(astrocytes_res$Gene.symbols), ThyTau22 = unique(astrocytes_res_tauCortex$Gene.symbols), ThyTau22_17m = unique(astrocytes_res_tauCortex_17m$Gene.symbols))
png("./results/figures/Venn_Overlapping_DEGs_4datasets_Astrocytes.png", width=8, height=6, units="in", res=300)
ggVennDiagram(overlap,
        category.names = c("Human_AD","Tg2576","ThyTau22_7m", "ThyTau22_17m"), label = "count", label_size = 6, set_size = 6) + scale_x_continuous(expand = expansion(mult = .2)) + theme(legend.position = "none") + scale_fill_distiller(palette = "RdBu")
dev.off()