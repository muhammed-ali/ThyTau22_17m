## Global Male-specific DEG's GRN analysis ##

rm(list = ls())
set.seed(1)

load("./results/files/Global_DEGs.RData")

dim(global_gender_spec_genes) # 107   6
table(global_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1               34               10               62

#
# Male
#

global_male <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
dim(global_male) # 62  6
max(global_male$Male.FDR) # 0.02708744

global_male$Expression <- ifelse(global_male$Male.avg..logFC > 0, 1, 0)
table(global_male$Expression)
#  0  1
# 26 36
global_male$Gene <- toupper(global_male$Gene.symbols)
Expression <- global_male[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 62  2
head(Expression, 2)
#      Gene Expression
# 1 GM14303          1
# 2     MAG          1


# write GRN input files
write.table(Expression, file="./results/GRN/Male/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="./results/GRN/Male/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="./results/GRN/Male/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="./results/GRN/Male/N2_Mapping_Expr.txt", sep="\t", row.names=F, quote=F, col.names=F)

Mapping <- global_male[,c("Gene", "Male.avg..logFC")]
Mapping <- unique(Mapping)
Mapping <- na.omit(Mapping)
dim(Mapping) # 62  2
write.table(Mapping, file="./results/GRN/Male/N1_Mapping_MaleLogFC.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).


# Shell commands for network analysis (to be run in terminal)
cd ./results/GRN/Male/

cat('java -jar ~/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 31 interactions, 16 nodes
wc -l NetworkPhenotype2.txt # 21 interactions, 13 nodes

java -jar ~/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 6

java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 7


java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt

sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt
' > run_network_analysis.sh)


## Global Sex-dimorphic DEG's GRN analysis ##


rm(list = ls())
set.seed(1)

load("./results/files/Global_DEGs.RData")

dim(global_gender_spec_genes) # 107   6
table(global_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1               34               10               62

gender_dimorphic <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-dimorphic",]
dim(gender_dimorphic) # 34  6
gender_dimorphic$Expression <- ifelse(gender_dimorphic$Male.avg..logFC > 0, 1, 0)
table(gender_dimorphic$Expression)
#  0  1
# 10 24
gender_dimorphic$Gene <- toupper(gender_dimorphic$Gene.symbols)
Expression <- gender_dimorphic[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 34  2
head(Expression, 2)
#       Gene Expression
# 64    CST3          1
# 65 GM42418          0


# write GRN input files
write.table(Expression, file="./results/GRN/Gender_Dimorphic/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="./results/GRN/Gender_Dimorphic/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="./results/GRN/Gender_Dimorphic/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="./results/GRN/Gender_Dimorphic/N2_Mapping_Expr.txt", sep="\t", row.names=F, quote=F, col.names=F)

Mapping <- gender_dimorphic[,c("Gene", "Male.avg..logFC")]
Mapping <- unique(Mapping)
Mapping <- na.omit(Mapping)
dim(Mapping) # 34  2
write.table(Mapping, file="./results/GRN/Gender_Dimorphic/N1_Mapping_MaleLogFC.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).

# Shell commands for network analysis (to be run in terminal)
cd ./results/GRN/Gender_Dimorphic/

cat('java -jar ~/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 31 interactions, 16 nodes
wc -l NetworkPhenotype2.txt # 21 interactions, 13 nodes

java -jar ~/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 6

java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 7


java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt

sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt
' > run_network_analysis.sh)


## Male-specific Microglial DEG's GRN Analysis

rm(list = ls())
set.seed(1)

load("./results/files/TauCortex17m_DEGs_4CT_astro_oligo_microg_neuron.RData")

#
# Microglial GSEA
#

dim(microglial_res_tauCortex_17m) # 55  7
table(microglial_res_tauCortex_17m$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               1                8                9               37

global_male <- microglial_res_tauCortex_17m[microglial_res_tauCortex_17m$DEG.type == "male-specific",]
dim(global_male) # 37  7
max(global_male$Male.FDR) # 0.04975419

global_male$Expression <- ifelse(global_male$Male.avg..logFC > 0, 1, 0)
table(global_male$Expression)
#  0  1
#  9 28
global_male$Gene <- toupper(global_male$Gene.symbols)
Expression <- global_male[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 37  2
head(Expression, 2)
#     Gene Expression
# 1 MALAT1          1
# 2  H2-D1          1

# write GRN input files
write.table(Expression, file="./results/GRN/Microglial/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="./results/GRN/Microglial/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="./results/GRN/Microglial/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="./results/GRN/Microglial/N2_Mapping_Expr.txt", sep="\t", row.names=F, quote=F, col.names=F)

Mapping <- global_male[,c("Gene", "Male.avg..logFC")]
Mapping <- unique(Mapping)
Mapping <- na.omit(Mapping)
dim(Mapping) # 37  2
write.table(Mapping, file="./results/GRN/Microglial/N1_Mapping_MaleLogFC.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).



# Shell commands for network analysis (to be run in terminal)
cd ./results/GRN/Microglial/

cat('java -jar ~/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 31 interactions, 16 nodes
wc -l NetworkPhenotype2.txt # 21 interactions, 13 nodes

java -jar ~/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 6

java -jar ~/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 7


java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt

sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt
' > run_network_analysis.sh)


## GRN analysis of longitudinal DEGs (THY-Tau22 17-months vs. 7-months age) ##

# Clear workspace and set random seed
rm(list = ls())
set.seed(1)

# Read longitudinal DEGs from the input file
# Adjust file path as needed
long_degs <- read.table("longitudinal_degs.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Separate DEGs by type
sex_dimorphic <- long_degs[long_degs$DEG.type == "sex-dimorphic",]
sex_neutral <- long_degs[long_degs$DEG.type == "sex-neutral",]

# Function to prepare expression data for network analysis (binarize logFCs)
prepare_network_data <- function(genes, timepoint) {
  # Convert expression to binary based on logFC
  Expression <- data.frame(
    Gene = genes$Gene.symbols,
    Expression = ifelse(genes[[paste0(timepoint, ".avg..logFC")]] > 0, 1, 0)
  )
  Expression <- unique(Expression)
  Expression <- na.omit(Expression)
  return(Expression)
}

# Prepare data for sex-dimorphic DEGs at both timepoints
dim_expr_17m <- prepare_network_data(sex_dimorphic, "Male")
dim_expr_7m <- prepare_network_data(sex_dimorphic, "Female")

# Prepare data for sex-neutral DEGs at both timepoints
neut_expr_17m <- prepare_network_data(sex_neutral, "Male")
neut_expr_7m <- prepare_network_data(sex_neutral, "Female")


# Write network input files for sex-dimorphic analysis
write.table(dim_expr_17m, file="dimorphic_17m_expression.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
write.table(dim_expr_7m, file="dimorphic_7m_expression.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
write.table(dim_expr_17m$Gene, file="dimorphic_geneList.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

# Write network input files for sex-neutral analysis
write.table(neut_expr_17m, file="neutral_17m_expression.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
write.table(neut_expr_7m, file="neutral_7m_expression.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
write.table(neut_expr_17m$Gene, file="neutral_geneList.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

# Create mapping files for visualization
write.table(dim_expr_17m, file="dimorphic_N1_nodeColorMapping.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
dim_expr_17m_N2 <- dim_expr_17m
dim_expr_17m_N2$Expression <- ifelse(dim_expr_17m_N2$Expression < 1, 1, 0)
write.table(dim_expr_17m_N2, file="dimorphic_N2_Mapping_Expr.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

write.table(neut_expr_17m, file="neutral_N1_nodeColorMapping.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)
neut_expr_17m_N2 <- neut_expr_17m
neut_expr_17m_N2$Expression <- ifelse(neut_expr_17m_N2$Expression < 1, 1, 0)
write.table(neut_expr_17m_N2, file="neutral_N2_Mapping_Expr.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

# Create logFC mapping files
dim_mapping <- sex_dimorphic[,c("Gene.symbols", "Male.avg..logFC")]
dim_mapping <- unique(dim_mapping)
dim_mapping <- na.omit(dim_mapping)
write.table(dim_mapping, file="dimorphic_N1_Mapping_MaleLogFC.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

neut_mapping <- sex_neutral[,c("Gene.symbols", "Male.avg..logFC")]
neut_mapping <- unique(neut_mapping)
neut_mapping <- na.omit(neut_mapping)
write.table(neut_mapping, file="neutral_N1_Mapping_MaleLogFC.txt", 
            sep="\t", row.names=F, quote=F, col.names=F)

#
#1.) Upload genes from dimorphic_geneList.txt (later repeat with neutral_geneList.txt) on https://portal.genego.com/cgi/data_manager.cgi
#2.) GeneGo - Build Network - Copy list - Gene symbol (official only), deactivate: Group of proteins, complex of proteins, group of complexes
#Summary
#IDs in input file:	89
#Mapped on maps:	49
#Mapped in networks:	78
#3.) Network algorithm: Direct interactions
#4.) Prefilters: Species: (do not restrict to humans)
#5.) Object types: select all, remove DNA, RNA, Compound, Inorganic ion, Predicted metabolite
#6.) interaction types: B, TR, IE, cRT, Rg (Binding, Transcription regulation, Influence on expression, co-regulation of transcription (don't limit: bottom: effects: only activating/inhibiting interactions
#7.) File - Export: both nodemap (gene symbol/object mappings, select: mus musculus for the nodemap: twice --> also for "Through") + interactions, Excel format
#thytau22_longitudinal_network_nodes.xls
#thytau22_longitudinal_network_interactions.xls
#
#8.) for nodes.xls: need to remove hyperlinks first: mark all (C andD columns): copy - paste as values (works, blue hyperlink format may remain, but hyperlinks are removed)
#- then save as: xlsx with the same name
#for interactions.xls: just save as xlsx


library(xlsx)

# Function to preprocess gene lists and create network input files
# @param netobj DataFrame containing network node information
# @param net DataFrame containing network interaction information 
# @param nodefile Output filename for node mappings
# @param intfile Output filename for interactions
# @param genefile Output filename for gene list
preprocess_genelist <- function(netobj, net, nodefile="nodemap.txt", 
                              intfile="interactions.txt", genefile="genelist.txt") {
    
    # Create and write node mapping
    # Maps network object names to gene symbols
    nodemap <- cbind(as.matrix(netobj$Network.Object.Name), 
                    as.matrix(netobj$Gene.Symbol))
    nodemap <- unique(nodemap)
    write.table(nodemap, file=nodefile, sep="\t", quote=FALSE, 
                col.names=FALSE, row.names=FALSE)
    
    # Create network matrix with source, target, and interaction type
    network <- cbind(
        as.matrix(nodemap[match(net[,2], nodemap[,1]), 2]),  # Source genes
        as.matrix(nodemap[match(net[,4], nodemap[,1]), 2]),  # Target genes
        as.matrix(net[,6])  # Interaction types
    )
    
    # Write interactions in network format
    # Converts interaction types to arrow notation:
    # Unspecified -> --?
    # Activation -> -->
    # Inhibition -> --|
    cat("", file=intfile, append=FALSE)
    for(i in 1:nrow(network)) {
        interaction_symbol <- ifelse(network[i,3] == "Unspecified", "--?",
                                   ifelse(network[i,3] == "Activation", "-->", "--|"))
        write(paste(network[i,1], "\t", interaction_symbol, "\t", network[i,2], sep=""),
              file=intfile, append=TRUE)
    }
    
    # Write unique gene symbols to gene list file
    write.table(unique(as.matrix(netobj$Gene.Symbol)), file=genefile,
                quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# Function to process expression data and create binary activity states
# @param datmat Matrix with gene symbols (col 1) and log fold changes (col 2)
# @param outfile Output filename for expression data
expprocess <- function(datmat, outfile="expression.txt") {
    # Convert log fold changes to binary states (0 for negative, 1 for positive)
    genesymbol <- as.matrix(datmat[,1])
    logfcs <- datmat[,2]
    
    # Write binary expression states
    write.table(cbind(genesymbol, ifelse(logfcs < 0, 0, 1)), 
                sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE,
                file=outfile)
}


# Read network node and interaction data
netobj <- read.xlsx("thytau22_longitudinal_network_nodes.xlsx", startRow=3)
net <- read.xlsx("thytau22_longitudinal_network_interactions.xlsx", startRow=3)

# Process network files
preprocess_genelist(netobj, net, "nodemap.txt", "interactions.txt", "geneList.txt")

# Shell commands for network analysis (to be run in terminal)
cat('
# For sex-dimorphic network analysis:
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/Preprocessor.jar . nodemap.txt interactions.txt dimorphic_geneList.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/DifferentialNetworkAnalysis.jar dimorphic_17m_expression.txt adjacency.txt dimorphic_GAResult.txt 0 true 1000 100 .

# Check network sizes
wc -l dimorphic_NetworkPhenotype1.txt  # 17m network
wc -l dimorphic_NetworkPhenotype2.txt  # 7m network

java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/CommonNetworkGenerator.jar dimorphic_NetworkPhenotype1.txt dimorphic_NetworkPhenotype2.txt dimorphic_CommonNetwork.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/DifferentialNetworkGenerator.jar dimorphic_NetworkPhenotype1.txt dimorphic_NetworkPhenotype2.txt dimorphic_DiffNetwork.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/ComputeCycles.jar dimorphic_CommonNetwork.txt dimorphic_17m_expression.txt dimorphic_pos.txt dimorphic_neg.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/SteadyStateCalculator.jar dimorphic_17m_expression.txt dimorphic_NetworkPhenotype1.txt 1 dimorphic_SteadyState1.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/SteadyStateCalculator.jar dimorphic_17m_expression.txt dimorphic_NetworkPhenotype2.txt 2 dimorphic_SteadyState2.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/PerturbagenListGenerator.jar dimorphic_pos.txt dimorphic_neg.txt dimorphic_DiffNetwork.txt dimorphic_SteadyState1.txt dimorphic_Perturbagens1.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/PerturbagenListGenerator.jar dimorphic_pos.txt dimorphic_neg.txt dimorphic_DiffNetwork.txt dimorphic_SteadyState2.txt dimorphic_Perturbagens2.txt

# For sex-neutral network analysis:
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/Preprocessor.jar . nodemap.txt interactions.txt neutral_geneList.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/DifferentialNetworkAnalysis.jar neutral_17m_expression.txt adjacency.txt neutral_GAResult.txt 0 true 1000 100 .

wc -l neutral_NetworkPhenotype1.txt
wc -l neutral_NetworkPhenotype2.txt

java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/CommonNetworkGenerator.jar neutral_NetworkPhenotype1.txt neutral_NetworkPhenotype2.txt neutral_CommonNetwork.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/DifferentialNetworkGenerator.jar neutral_NetworkPhenotype1.txt neutral_NetworkPhenotype2.txt neutral_DiffNetwork.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/ComputeCycles.jar neutral_CommonNetwork.txt neutral_17m_expression.txt neutral_pos.txt neutral_neg.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/SteadyStateCalculator.jar neutral_17m_expression.txt neutral_NetworkPhenotype1.txt 1 neutral_SteadyState1.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/SteadyStateCalculator.jar neutral_17m_expression.txt neutral_NetworkPhenotype2.txt 2 neutral_SteadyState2.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/PerturbagenListGenerator.jar neutral_pos.txt neutral_neg.txt neutral_DiffNetwork.txt neutral_SteadyState1.txt neutral_Perturbagens1.txt
java -jar ~/ownCloud/Documents/NetworkAnalysis_Jars/PerturbagenListGenerator.jar neutral_pos.txt neutral_neg.txt neutral_DiffNetwork.txt neutral_SteadyState2.txt neutral_Perturbagens2.txt
' > run_network_analysis.sh)



