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
library(car) # For regression analysis

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
moduleColorsManual3 = merge_man$colors
# Eigengenes of the newly merged modules:
MEsManual = merge_man$newMEs

MEsManualList=moduleEigengenes(datExpr, align == "along average", colors=moduleColorsManual3)
signif(cor(MEsManualList$eigengenes, use="p"),2)
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
save(A_man, dissTOM_man, moduleLabelsManual1, MEList_man, MEs_man, MET_man, 
     moduleColorsManual2, moduleColorsManual3, merge_man, MEsManual, MEs_man_merge, 
     moduleLabels_man, mergedColors_man, geneTree_man, 
     file = "SKCM-networkConstruction-MHmerged.RData")



########## NEW STUFF #################
### From Tutorial 3.6

#####################################################################################
### WGCNA ANALYSIS - Comparison of Eigengenes across patient samples within modules
#####################################################################################

# Get a variable with just the module eigengene expression for all patients across all modules
# (without NaNs), and find the correlation (matrix) between modules
datME = MEsManualList$eigengenes
datME_filtered = datME[complete.cases(datME),]
signif(cor(datME, use="p"),2)

# Making heatmaps of genes (row) and patient samples (columns)
# Green = underexpression; red = overexpression
sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="maroon";
plotMat(t(scale(datExpr[,moduleColorsManual3==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

sizeGrWindow(8,7);
which.module="lightcyan1"
ME=datME_filtered[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColorsManual3==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)

# Looking into the variance within each module (explained by eigengene)
# Idea - larger the variance - more interesting?
MEsManual_Variance = propVarExplained(datExpr, colors = moduleColorsManual3, MEs = MEsManual, corFnc = "cor", corOptions = "use = 'p'")

write.table(MEsManual, file = "ModuleEigengenes.txt", sep = '\t')
write.table(MEsManual_Variance, file = "ModuleVariance.txt", sep="\t", col.names = NA)

# Idea - mean eigengene value + variance could mean something? box plots?
MEs_mean = colMeans(MEsManual, na.rm = T)

# Took spreadsheet into Excel and sorted MEs by variance (descending) - could have done here
MEs_Manual_Variance_sorted = read.table(file = "ModuleVariance_sorted.txt", sep = '\t', header = T)

# Pulling out "PVE" in "PVE(module_color)" string for plotting variance bar chart
colors = substring(MEs_Manual_Variance_sorted$MEs, 4)
# Making barchart of sorting descending variance amoung MEs
pdf("MEvariance.pdf", width = 12, height = 5)
par(mar=c(10,5,2,1))
barplot(MEs_Manual_Variance_sorted$Variance, main = "ME Variance", xlab = "", ylab = "PropVarExplained", 
        ylim = c(0,0.6), col = colors, names.arg = colors, las = 2)
mtext("Modules", side = 1, line = 7, las = 1, cex = 1)
dev.off()

# IMPORTANT PLOTS FOR F31 - Eigengene variability within each module for all patient samples
# Making a single PDF for all modules - eigengene expression (y) vs. patient samples (x)
pdf("MEplots_EEvPS.pdf", width = 12, height = 5)
plot_list = list()
for(i in colors){
  ME=datME_filtered[, paste("ME",i, sep="")]
  par(mar=c(5, 5, 2.5, 0.2))
  p = barplot(ME, col=i, main= paste("ME", i, sep = ""), cex.main=2,
              ylab="Eigengene Expression",xlab="Patient Sample")
  plot_list[[i]] = p
  print(plot_list[[i]])
}
dev.off()


#####################################################################################
### WGCNA ANALYSIS - Intramodular Connectivity (IC) and Module Membership (MM) Compared
#####################################################################################

Alldegrees1 = intramodularConnectivity(A_man, moduleColorsManual3)
head(Alldegrees1)

# Relationship between gene significance and intramodular connectivity - if have patient clinical data
#colorlevels = unique(moduleColorsManual3)
#sizeGrWindow(9,6)
#par(mfrow = c(2, as.integer(0.5+length(colorlevels)/2)))
#par(mar = c(4,5,3,1))
#for (i in c(1:length(colorlevels))){
#  whichmodule = colorlevels[[i]];
#  restrict1 = (moduleColorsManual3 == whichmodule)
#  verboseScatterplot(Alldegrees1$kWithin[restrict1],
#                     GeneSignificance[restrict1], col=colorh1[restrict1],
#                     main=whichmodule,
#                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
#}

# Finding intramodular connectivity among all genes measured
datKME = signedKME(datExpr, datME, outputColumnName = "MM.")
## Got "Some genes are constant" error message - introduced NaNs into data matrix
head(datKME)
View(datKME)

# Finding genes with high intramodular connectivity in interesting modules (most variance)
# Insert any module here - MM.module_color - or loop over all?
FilterGenes = abs(datKME$MM.maroon)>0.8 #& abs(GS1)>0.2 # for patient data
table(FilterGenes)
FilterGenes = FilterGenes[complete.cases(FilterGenes)] # Removing NaNs
dimnames(data.frame(datExpr))[[2]][FilterGenes] # Getting gene names that pass above threshold

# Finding the relationship between module membership (MM) and intramodular connectivity in interesting
# modules
pdf("MMvIC_top4modules.pdf")
#sizeGrWindow(8,6)
par(mfrow=c(2,2))
# Choose certain number of modules and set up so fits in matrix mfrow (2x2=4)
# Need to change this into a loop
which.color1="maroon";
restrictGenes=moduleColorsManual3==which.color1
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color1, sep="")])^6,
                   col=which.color1,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color1="darkorange2";
restrictGenes=moduleColorsManual3==which.color1
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color1, sep="")])^6,
                   col=which.color1,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color1="thistle2";
restrictGenes=moduleColorsManual3==which.color1
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color1, sep="")])^6,
                   col=which.color1,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color1="skyblue3";
restrictGenes=moduleColorsManual3==which.color1
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color1, sep="")])^6,
                   col=which.color1,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
dev.off()

save(MEsManual, MEsManualList, datME, datME_filtered, ME, MEsManual_Variance, MEs_Manual_Variance_sorted,
    MEs_mean, colors, Alldegrees1, datKME, FilterGenes, file = "July262016.RData")


#####################################################################################
### WGCNA ANALYSIS - Output Gene Lists Per Module and Enrichment
#####################################################################################

### GET GENE LISTS ###

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv"); # REPLACE WITH TCGA GENE ANNOTATION FILE
# Match probes in the data set to the probe IDs in the annotation file 
annot = read.table(file = "TCGA.hg19.June2011.gaf", sep = '\t')
probes = names(datExpr)
probes_to_annot = match(probes, annot$V2)
allLLIDs = annot$V2[probes_to_annot]
#probes2annot = match(probes, annot$substanceBXH) # Probably not right - need to see what IDs to look for
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot]; # Probably not right - need to see what IDs to look for
# Choose interesting/all modules and get gene list - make sure to specify from what method and if merged/not
# OUTPUT should be list of modules (columns) and associated genes (rows)
#intModules = c("maroon", "darkorange2", "thistle2") 

modList = list()
for (module in moduleColorsManual3)
{
  # Select module probes
  modGenes = (moduleColorsManual3==module)
  # Get their entrez ID codes
  #modLLIDs = allLLIDs[modGenes];
  modCodes = probes[modGenes]
  # Write them into a list
  modList[[module]] = modCodes
#  fileName = paste("/Users/Corey/Documents/QuarantaLab/WGCNA/GE_SKCM/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/moduleFiles/modCodes_", module, ".txt", sep="");
#  write.table(as.data.frame(modCodes), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
}
getGeneName = function (geneChar, split="|", pos=1)
  unlist(sapply(strsplit(geneChar, split = split, fixed = TRUE), function(x) (x[pos]), simplify = FALSE)) 


indx = sapply(modList, length)
res = as.data.frame(do.call(cbind, lapply(modList, 'length<-', max(indx))))
#colnames(res) = names(modList[[which.max(indx)]])
#library(stringr)
#str_split_fixed()
for (i in 1:nrow(res_test2)){
  for (i in 1:length(res_test2)) {
    res_codes = sapply(strsplit(res_test2, '|', fixed = TRUE), function(x) (x[2]))
  }
  }

write.table(res, file = "modCodes_all.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

# As background in the enrichment analysis, we will use all probes in the analysis.
# fileName = paste("LocusLinkIDs-all.txt", sep="");
# write.table(as.data.frame(allLLIDs), file = fileName,
#             row.names = FALSE, col.names = FALSE)

# Splitting the strings up so that can use - need to remove NAs

for (i in length(res_test)) {
  a = sapply(strsplit(as.character(res_test), '|', fixed = TRUE), function(x) (x[2]))}

for (i in 1:nrow(res_test2)) {
  a = sapply(strsplit(as.character(res_test2), '|', fixed = TRUE), function(x) (x[2]))}

res_test2_mat = as.matrix(res_test2)
res_test2_mat_sub = res_test2_mat[1:5,]

for (i in 1:dim(res_test2_mat_sub)[1])
{
  for (j in 1:dim(res_test2_mat_sub)[2])
  {
    a_sub = sapply(strsplit(res_test2_mat_sub, split = '|', fixed = TRUE), function(x) (x[2]))
    a_sub_nl = cat(paste0("\t", as.matrix(a_sub[1:]), "\n"))
    
  }
}

within(res_test2, foo <- data.frame(do.call('rbind', strsplit(as.character(res_test2), "|", fixed = TRUE))))
  #a = unlist(strsplit(res_test, "|", fixed = TRUE))[2]}
head(a)


### MODULE ENRICHMENT ###

# Enrichment analysis - runs for a while - outputs top 10 enriched GO terms with each module (colors)
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10);
# Assign enrichment terms for each module to variable
tab = GOenr$bestPTerms[[4]]$enrichment
# Get _______
names(tab)
# Write table to csv file - search for genes of interest in modules
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

# Figure out what columns (components) within table are most important to you
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

#####################################################################################
### WGCNA ANALYSIS - Network Visualization (MDS, PCA)
#####################################################################################

library(cluster)
options(stringsAsFactors = FALSE);

### Multidimensional Scaling (MDS) plot - version of PCA - DON'T DO IN RSTUDIO, WILL CRASH - takes 3+ hours
#pdf("MDS_SKCM.pdf")
cmd1 = cmdscale(as.dist(dissTOM_man),2)
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(cmd1, col = as.character(moduleColorsManual3), main = "MDS Plot", 
     xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
#dev.off()

### Principal Component Analysis (PCA) - outside of WGCNA - rough
View(MET2)
MET3 = MET2[complete.cases(MET2),]
dim(MET2)
dim(MET3)

fit = princomp(MET3, cor = TRUE, scores = TRUE) # cor = TRUE)
summary(fit)
loadings(fit)
fit$scores
plot(fit, type="lines")

par(mar = c(1.5, 1.5, 1.5, 1.5))
biplot(fit)

library(devtools)
install_github("fawda123/ggord")
library(ggord)
ord = prcomp(fit[,1:25])

install.packages("rgl")
library(rgl)
plot3d(fit$scores[,1:3])
text3d(fit$scores[, 1:3], texts = rownames(MET3))
text3d(fit$loadings[,1:3], texts = rownames(fit$loadings), col = "red")
cords = NULL
for (i in 1:nrow(fit$loadings)) {
  cords = rbind(cords, rbind(c(0,0,0), fit$loadings[i, 1:3]))
}
lines3d(cords, col = "red", lwd = 4)

### TOM Plots - Takes lots of time and not very informative for large datasets
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


#####################################################################################
### WGCNA ANALYSIS - Other Tutorial Stuff
#####################################################################################

meanExpressionByArray = apply(datExpr,1, mean, ra.rm = T)
head(meanExpressionByArray)
NumberMissingByArray = apply(is.na(data.frame(datExpr)), 1, sum)
head(NumberMissingByArray)

# Looking at overall gene expression per sample - consistency
sizeGrWindow(9,5)
barplot(meanExpressionByArray, xlab = "Sample", ylab = "Mean Expression",
        main = "Mean Expression across Samples", names.arg = c(1:474), cex.names = 0.7)

keepArray = NumberMissingByArray<500
table(keepArray)
datExpr = datExpr[keepArray,]
#y = y[keepArray]
#ArrayName[keepArray]

# Other potential subsetting techniques
NumberMissingByGene = apply(is.na(data.frame(datExpr)),2,sum)
summary(NumberMissingByGene)
variancedatExpr = as.vector(apply(as.matrix(datExpr), 2, var, na.rm = T))
no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)), 2, sum))
table(no.presentdatExpr)
keepGenes = variancedatExpr>0 & no.presentdatExpr>=4
table(keepGenes)
datExpr = datExpr[, keepGenes]
#GeneName = GeneName[keepGenes]

#sizeGrWindow(9,5)
#plotClusterTreeSamples(datExpr = datExpr, y = keepArray)

#datME = moduleEigengenes(datExpr, colorh1)$eigengenes
#signif(cor(datME, use = "p"), 2)
#signif(cor(MEs_man_merge, use = "p"), 2)

## Already generated a similar plot earlier (above) - just ME labels instead of module labels
#dissimME = (1-t(cor(datME, use = "pairwise.complete.obs", method = "p")))/2
#flashClustMEs_man = flashClust(as.dist(dissimME), method = "complete")
#par(mfrow = c(1,1))
#plot(flashClustMEs_man, main = "Clustering Tree Based on Module Eigengenes")

#sizeGrWindow(8,9)
#par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
#which.module="turquoise"; 
#plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
#        clabels=T,rcols=which.module,
#        title=which.module )
# for the second (blue) module we use
#which.module="blue";  
#plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
#        clabels=T,rcols=which.module,
#        title=which.module )
#which.module="brown"; 
#plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
#        clabels=T,rcols=which.module,
#        title=which.module )

#####################################################################################
### Ideas
#####################################################################################

# Cluster or PCA eigengenes into groups and look deeper into those groups
# Consensus clustering of patients in each (important) module
# Seeing if clusters correspond to oncogenic mutations

### Visualize the eigengene network ###

#MEs2 = moduleEigengenes(datExpr, moduleColors1)$eigengenes
## Calculate eigengenes
#MEList2=moduleEigengenes(datExpr,colors=moduleColors1)
#MEs2 = MEList$eigengenes
#MET2=orderMEs(MEs2)
## Plot the dendrogram
#sizeGrWindow(6,6);
#par(cex = 1.0)
#plotEigengeneNetworks(MET2, "Eigengene dendrogram", marDendro = c(0,4,2,0),
#                      plotHeatmaps = FALSE)
#
## Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
#par(cex = 1.2)
#plotEigengeneNetworks(MET2, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
#                      plotDendrograms = FALSE, xLabelsAngle = 90)