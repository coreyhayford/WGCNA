#####################################################################################
### SETTING UP THE ENVIRONMENT
#####################################################################################

# Installing necessary packages for RNASeq data and WGCNA
install.packages("rnaseqWrapper")
install.packages("WGCNA")
install.packages("car")

#loading libraries
library(rnaseqWrapper) # For Differential Expression analysis and variant analysis
library(WGCNA) # Includes variety of functions associated with WGCNA
library(car) # Fro regression analysis

# Setting the working directory to the file path - i.e. where TCGA data is

getwd()
setwd('/Users/Corey/Documents/QuarantaLab/WGCNA/GE_SKCM/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/')

## Determining if package is working properly with WGCNA tutorial data
#setwd('/Users/Corey/Desktop/R_class/Tutorial1')
#fem_data = read.csv("LiverFemale3600.csv")
#dim(fem_data)
#names(fem_data)


##Looking at the file and its properties
#setwd('/Users/Corey/Documents/WeaverLab/RNASeq_Jan/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3')
#file = read.csv("unc.edu.0ada6aee-0434-44fa-a452-c9a9fb69ef90.2355628.rsem.genes.results", sep="\t")
#file
#dim(file)
#names(file)

options(stringsAsFactors = FALSE)
enableWGCNAThreads() # Using multiple cores

#####################################################################################
### GETTING RNASeq DATA SET IN CORRECT FORMAT
#####################################################################################

# Make data file with level 3 analysis (TCGA metric, has high-level RNASeq data)
# only using the normalizes files (all the fileIDs with genes.normalized.results file
# type - they are tab-delimited. The first column (seqIDcol) has the gene names, and 
# we want to keep the the column named "normalized_count (column 2) - ignore idCols
# - and minMatchToMerge is the portion of gene names per file that need to match to 
# merge multiple datasets (there are 400+ genes_normalized_results files)
countData = mergeCountFiles ('Level_3/', 
                             fileID="*.genes.normalized_results$", 
                             fileSep="\t", 
                             seqIDcol=1,
                             colsToKeep=c("normalized_count"),
                             idCols=NULL,
                             minMatchToMerge=0.5)

# Visualizing dataframe and looking at dimensions of matrix
View(countData)
dim(countData)

# Make a variable that has the transposed countData matrix -- patients on x axis 
# and genes on y axis (typical way to look at data) and ensure the header still
# is present (or else it cuts off the data from the first patient)
myCountData = as.data.frame(t(countData), header=TRUE)

View(myCountData)
dim(myCountData)

# Make a variable that holds the dataframe (matrix) of a csv (actually txt) that has the 
# RNAseq data (comes in download from TCGA) - want to take patient barcode - 
# essentially a way they track patients - so you can match them to the RNASeq
# patients later
barcode = as.data.frame(read.csv("/Users/Corey/Documents/QuarantaLab/WGCNA/GE_SKCM/METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_SKCM.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt", header = TRUE, sep = '\t'))
View(barcode)
dim(barcode)

# Subset the barcode matrix for the row values associated with normalized gene
# expression values and assign it to a variable
UUID_to_barcode=barcode[barcode$Comment..TCGA.Data.Type..1=='RSEM_genes_normalized',]
View(UUID_to_barcode)
dim(UUID_to_barcode)

# Get extract name (aka barcode) and data-derived file (columns 2 and 22 of matrix)
# and assign it to a variable
samples_barcode=UUID_to_barcode[,c(2,22)]
View(samples_barcode)
dim(samples_barcode)

# Change the row names of the myCountData matrix so that the .normalized_count
# patients/barcodes are replaced with the .genes_normalized_results patients/barcodes
# - dimensions stay same
rownames(myCountData)=sub(".normalized_count*", ".genes.normalized_results", rownames(myCountData))
View(myCountData)
dim(myCountData)

# Make a variable with the row names of the count data matrix -- match those with
# the data file name in the barcode matrix (same name, why they can be matched) - 
# make a new count data file rearranged by these rows and change the row names of 
# the matrix from barcode to the TCGA barcode - write it to a tab delimited txt file
samples_uuid = rownames(myCountData)
head(samples_uuid)
samplesRows = match(samples_barcode$Derived.Data.File, samples_uuid)
View(samplesRows)
myCountDataRearranged = myCountData[samplesRows,]
rownames(myCountDataRearranged) = samples_barcode$Comment..TCGA.Barcode.

View(myCountDataRearranged)
dim(myCountDataRearranged)

write.table(myCountDataRearranged, file = "myCountDataRearranged.txt", sep="\t", col.names = NA)

final.data = myCountDataRearranged


#####################################################################################
### START WGCNA ANALYSIS
#####################################################################################

# Cluster gene expression data by distance (default = euclidean?) using complete
# linkage and plot the dendrogram (tree). Choose a signed network (so can see
# upregulation and downregulation) and setup graphical parameters. Change the 
# correlation (numbers) to a color heatmap (grey - no data, white - not 
# correlated, red - very correlated), and plot the dendrogram with trait
# heatmap for each trait -- summary of relatedness


# flashClust library does not have bugs associated with hclust - and is faster
# Sample Tree shows relatedness of patients according to gene expression data
# Depending on how closely related patients are, may not need to remove patients

library(flashClust)

#####################################################################################
### WGCNA ANALYSIS - Network Construction and Module Detection
#####################################################################################


# Using WGCNA naming convention
datExpr = final.data 

# Set the threshold for merging closely related modules
mergingThresh = 0.90

# Automatic network construction and module detection on large expression datasets 
# in a block-wise manner (using blockwisemodules function). Input the cleaned
# expression dataset and group by pearson correlation. Create the correlation matrix
# using the signed network equation - beta power = 12. Minimum and maximum module sizes
# (correlations) set, and merge closely associated modules. Module colors output 
# in numerical form.

### THIS IS A FIRST PASS AT THE ANALYSIS BY USING THE MOST BASIC, DEFAULT SETTINGS ###

net = blockwiseModules(datExpr, power = 12, TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                       pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "SKCM_TOM_2",
                       verbose = 3)
table(net$colors)
sizeGrWindow(12,9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "SKCM-networkconstruction.RData")

#####################################################################################
### WGCNA ANALYSIS - Choosing a soft threshold to match scale independence and mean connectivity
#####################################################################################

### i.e. Optimizing the settings depending on your dataset ###
# SFT = scale free topology

# Choose a set of soft-thresholding powers (1:10 by 1, 12:20 by 2) and use the data
# to pick the soft threshold where SFT > 0.90 and connectivty seems to flatten
# Default Beta for Signed Network (Soft Threshold Power = 12)
# Choose a Beta power - specific to each analysis
# If using R through Graphical User Interface (i.e. RStudio), you may need
# to disable parallel processing with disableWGCNAThreads(). But, make sure it
# is enabled for all further processes

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# Pick soft threshold requires you to disableWGCNAtheads() because it 
# does not work with parallel processing enabled. Enable immediately afterwards!
disableWGCNAThreads()
sft = pickSoftThreshold(final.data, powerVector = powers, verbose = 5, networkType = "signed")
enableWGCNAThreads()
# Plot the results in a certain plot size and graphical parameters
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

par(mfrow=c(1,2))
# SFT index as a function of different powers
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
# this line corresponds to using an R^2 cut-off of tree height
abline(h=0.90,col="red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")

#####################################################################################
### WGCNA ANALYSIS - Making the adjacency/TOM matrices and dendrogram with Dynamic Tree Cutting (DTC)
#####################################################################################

softPower = 12;
### See example

Ad = adjacency(datExpr, power = 12)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(Ad, TOMType="signed")
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "complete");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_new, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate Eigengenes
MEList1 = moduleEigengenes(datExpr, colors = dynamicColors)
MEs1 = MEList1$eigengenes
# Calculate module eigengene dissimilarity
MEDiss1 = 1-cor(MEs1, use = "pairwise.complete.obs");
# Cluster MEs
METree1 = flashClust(as.dist(MEDiss1), method = "complete");
# Plotting results
sizeGrWindow(7,6)
plot(METree1, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres1 = 0.25
# Plot the threshold
abline(h=MEDissThres1, col = "red")
# Automatic merging
merge1 = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres1, verbose = 3)
# Get merged module colors
mergedColors1 = merge$colors
# New MEs from merged modules
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors1),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors1 = mergedColors1
# Construct numerical labels corresponding to the colors
colorOrder = c('grey', standardColors(50));
moduleLabels = match(moduleColors1, colorOrder)-1;
MEs1 = mergedMEs;
# Save module colors and labels
save(MEs, moduleLabels, moduleColors, geneTree, file = "SKCM-networkConstruction-stepbystep.RData")

#####################################################################################
### WGCNA ANALYSIS - Blockwise module detection for large networks - separate blocks
#####################################################################################

#Can try changing maxBlockSize to get fewer blocks
bwnet = blockwiseModules(datExpr,corType="pearson",
                         maxBlockSize=5500,networkType="signed",power=12,minModuleSize=30,
                         mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                         pamRespectsDendro=FALSE,saveTOMFileBase="SKCMTOM-blockwise",verbose=3)

# Relabel blockwise modules so that their labels
# match those from our previous analysis
moduleLabelsBlockwise = bwnet$colors
#moduleLabelsBlockwise=matchLabels(bwnet$colors,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsBlockwise = labels2colors(moduleLabelsBlockwise)
#moduleColorsBlockwise=labels2colors(moduleLabelsBlockwise)
# measure agreement with single block automatic procedure
#mean(moduleLabelsBlockwise==moduleLabels)

blockNumber=1
# Plot the dendrogram for the chosen block
plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
                    main=paste("Dendrogram and module colors in block",blockNumber),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
blockNumber=2
# Plot the dendrogram for the chosen block
plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
                    main=paste("Dendrogram and module colors in block",blockNumber),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
blockNumber=3
# Plot the dendrogram for the chosen block
plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
                    main=paste("Dendrogram and module colors in block",blockNumber),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

blockNumber=4
# Plot the dendrogram for the chosen block
plotDendroAndColors(bwnet$dendrograms[[blockNumber]],
                    moduleColorsBlockwise[bwnet$blockGenes[[blockNumber]]],"Module colors",
                    main=paste("Dendrogram and module colors in block",blockNumber),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

# Comparing the dentrogram of the one-block and multiblock approaches
sizeGrWindow(12,9)
plotDendroAndColors(geneTree, cbind(moduleColors1, moduleColorsBlockwise),
                    c("Single Block", "4 Blocks"), main = "Gene Dendrogram and Module Colors Compared",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(datExpr, moduleColorsBlockwise)$eigengenes;
single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

#####################################################################################
### WGCNA ANALYSIS - Manual, stepwise module detection
#####################################################################################

# Calculate the weighted adjacency matrix, using the power 12:
A_man = adjacency(datExpr, type = "signed", power = 12)

# Define a dissimilarity based on the topological overlap
dissTOM_man =TOMdist(A_man)
# Hierarchical clustering
geneTree_man = flashClust(as.dist(dissTOM_man),method="complete")
# Define the modules by DTC
moduleLabelsManual1=cutreeDynamic(dendro=geneTree_man,distM=dissTOM_man,
                                  method="hybrid",deepSplit=2,pamRespectsDendro=F,minClusterSize=30)

# Match those from our previous analysis - not used for this analysis because not done above
# moduleLabelsManual2=
#   matchLabels(moduleLabelsManual1,moduleLabelsAutomatic)

# Convert labels to colors for plotting
moduleColorsManual2=labels2colors(moduleLabelsManual1)

# Calculate eigengenes
MEList_man=moduleEigengenes(datExpr,colors=moduleColorsManual2)
MEs_man = MEList_man$eigengenes
# Add the weight to existing module eigengenes
#MET=orderMEs(cbind(MEs,weight))
MET_man=orderMEs(MEs_man)
# Plot the relationships among the eigengenes and the trait
plotEigengeneNetworks(MET_man,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

table(moduleLabelsManual1)
# Convert numeric lables into colors
moduleColorsManual2=labels2colors(moduleLabelsManual1)
table(moduleColorsManual2)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_man, moduleColorsManual2, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge highly correlated modules - not necessary but good next step
merge_man=mergeCloseModules(datExpr,moduleColorsManual2,
                        cutHeight=MEDissThres1)
# Resulting merged module colors
moduleColorsManual3 = merge$colors
# Eigengenes of the newly merged modules:
MEsManual = merge_man$newMEs

# Show the effect of module merging by plotting the
# original and merged module colors below the tree
datColors=data.frame(moduleColorsManual3,moduleColorsManual2)
#datColors=data.frame(moduleColorsManual3,moduleColorsAutomatic,
#                     moduleColorsBlockwise)
View(datColors)
plotDendroAndColors(geneTree_man,colors=datColors,
                    groupLabels=c("Manual Hybrid","DTC"),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

datColors_all = data.frame(moduleColorsManual3, mergedColors1, moduleColorsBlockwise)
plotDendroAndColors(geneTree_man,colors=datColors_all,
                    groupLabels=c("manual hybrid","merged Dynamic","Blockwise"),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

# Observation - Merged Dynamic and Manual Hybrid are the same

datColors_2 = data.frame(moduleColorsManual3, moduleColorsBlockwise)
plotDendroAndColors(geneTree_man,colors=datColors_2,
                    groupLabels=c("manual hybrid", "Blockwise"),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

# check the agreement between manual and automatic module labels
mean(moduleColorsManual3==moduleColorsBlockwise) 

# Rename to moduleColors
moduleColorsManual3 = mergedColors_man
# Construct numerical labels corresponding to the colors
colorOrder_man = c('grey', standardColors(50));
moduleLabels_man = match(mergedColors_man, colorOrder_man)-1;
MEs_man_merge = mergedMEs_man;
# Save module colors and labels
save(MEs_man_merge, moduleLabels_man, mergedColors_man, geneTree_man, file = "SKCM-networkConstruction-MHmerged.RData")


#####################################################################################
### WGCNA ANALYSIS - Network Visualization
#####################################################################################

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Calculate the Topological Overlap again
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 12);
# Transform dissTOM to make connections visible
plotTOM = dissTOM^7;
# Set diag to NA for plot
diag(plotTOM) = NA;
# Plot it
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors1, main = "Network heatmap plot, all genes")

### Visualize the eigengene network ###

MEs2 = moduleEigengenes(datExpr, moduleColors1)$eigengenes
# Calculate eigengenes
MEList2=moduleEigengenes(datExpr,colors=moduleColors1)
MEs2 = MEList$eigengenes
MET2=orderMEs(MEs)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.2)
plotEigengeneNetworks(MET2, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

### Principal Component Analysis (PCA) of eigengenes
View(MET2)
MET3 = MET2[complete.cases(MET2),]
dim(MET2)
dim(MET3)
fit = princomp(MET3, cor = TRUE)
summary(fit)
loadings(fit)
fit$scores
plot(fit, type="lines")
par(mar = c(1.5, 1.5, 1.5, 1.5))
biplot(fit)

### Ideas ###

# Cluster or PCA eigengenes into groups and look deeper into those groups
# Consensus clustering of patients in each (important) module
# Seeing if clusters correspond to oncogenic mutations
