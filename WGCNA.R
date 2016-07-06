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
setwd('/Users/Corey/Documents/WeaverLab/RNASeq_Jan/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/')

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
enableWGCNAThreads()

#####################################################################################
### GETTING RNASeq DATA SET IN CORRECT FORMAT
#####################################################################################

# Make data file with level 3 analysis (TCGA metric, has high-level RNASeq data)
# only using the normalizes files (all the fileIDs with genes.normalized.results file
# type - they are tab-delimited. The first column (seqIDcol) has the gene names, and 
# we want to keep the the column named "normalized_count (column 2) - ignore idCols
# - and minMatchToMerge is the portion of gene names per file that need to match to 
# merge multiple datasets (there are 500+ genes_normalized_results files)
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
barcode = as.data.frame(read.csv("/Users/Corey/Documents/WeaverLab/RNASeq_Jan/METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_HNSC.IlluminaHiSeq_RNASeqV2.1.10.0.sdrf.txt", header = TRUE, sep = '\t'))
View(barcode)
dim(barcode)

# Subset the barcode matrix for the row values associated with normalized gene
# expression values and assign it to a variable - 3396 --> 566 patients
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

### Oscar's way - makes rows as barcodes and columns as genes

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


#####################################################################################
### CLEANING THE DATA SET
#####################################################################################

# NOT ESSENTIAL FOR CODE - used to look for unique cancer types - find a way to use 
# unique function to look at tumor types (inside the barcode)
#cancer.phenotypes <- unique()

# Remove rows from rearranged count dataset that match a certain argument pattern - 
# not actually sure what the biological relevance is here, but it changes the number
# of patients from 566 -- 520
myCleanDataSet=myCountDataRearranged[-grep('(TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-[1]|TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}-06)', rownames(myCountDataRearranged)), ]
length(rownames(myCleanDataSet))
View(myCleanDataSet)

# Simplify the row names of the new clean dataset by just keeping the first 4 digits
# of the TCGA barcode so can match later -- does not change number of patients
rownames(myCleanDataSet)=sub("-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}", "", rownames(myCleanDataSet))
View(myCleanDataSet)


#####################################################################################
### PREPARING THE CLINICAL DATASET
#####################################################################################

# Change working directory - data will save this way - and read in the clinical data
getwd()
setwd("/Users/Corey/Documents/WeaverLab/Clinical_Jan/Clinical/Biotab")
clinical.data = as.data.frame(read.csv("nationwidechildrens.org_clinical_patient_hnsc.txt", header = TRUE, sep = '\t'))
View(clinical.data)

# Because multiple header rows before data, removed additional headers and id numbers;
# IMPORTANT - clinical data has many important information columns - I chose certain
# columns that fit this dataset, but should be tailors to cancer type and what the
# user wants to look at!
# Should make more efficient way of concatenating columns rather than by column 
# names - maybe by column headers as strings
# KEEP patientID, barcode, number lymph nodes examined, vital status,
# days to death, tumor stage, tumor grade, and clinical stage
# Not all data has values - introduces many NaNs (not a number/value)
clinical.data.cleaned = clinical.data[-c(1,2), c(1,2,21,26,28,34,36,56)] 
View(clinical.data.cleaned)

# Attempt to do concatenate colukmns by column names
#column.data <- within(column.data, id <- paste(bcr_patient_uuid, bcr_patient_barcode, lymph_nodes_examined_he_count, vital_status, death_days_to, ajcc_pathologic_tumor_stage, tumor_grade, clinical_stage, sep = "\t")
#clinical.data.columns = clinical.data.rowscleaned[ , column.data]

# Determine the number of rows in each data set, so they can be matched and 
# clustered later - they are different (520 - 528)
length(rownames(myCleanDataSet))
length(rownames(clinical.data.cleaned))

# Match patient barcodes from clinical and RNA seq data (on x axis of table)
# and replace the barcode in patient dataset with the predetermined, simplified
# row names
ordered.set = match(clinical.data.cleaned$bcr_patient_barcode, rownames(myCleanDataSet))
matched.set = clinical.data.cleaned[ordered.set,]
View(clinical.data.cleaned)
View(ordered.set)
View(matched.set)

# Remove those patients (rows) from the patient data set who do not match - should
# get same row length at end of step (number of patients in each data set = 520)
matched.set.removena = matched.set[-which(is.na(ordered.set)), ]
View(matched.set.removena)
length(rownames(matched.set.removena))
length(rownames(myCleanDataSet))

# Determine patient data set that are HPV- by importing auxilary data set, and 
# select the patients with a "Negative" string value for the hpv_status column
# (i.e. remove HPV+ patients). Take matched patient dataset and match the columns
# of the patient and auxilary datasets corresponding to hpv status - remove HPV+
# patients and NaNs (the number of HPV- patients from auxiliary dataset is 427; for
# matched patient dataset = 420)

HPV.negative = as.data.frame(read.csv('/Users/Corey/Documents/WeaverLab/Clinical_Jan/Clinical/Biotab/nationwidechildrens.org_auxiliary_hnsc.txt', header = TRUE, sep = '\t'))
HPV.negative.true = HPV.negative[HPV.negative$hpv_status == "Negative", ]
View(HPV.negative.true)

matched.set.removeHPV = match(matched.set.removena$bcr_patient_barcode, HPV.negative.true$bcr_patient_barcode)
matched.set.HPV = matched.set.removena[-which(is.na(matched.set.removeHPV)), ]
View(matched.set.HPV)

length(rownames(HPV.negative.true))
length(rownames(matched.set.HPV))

#For future reference, may want to determine number of HPV positive patients 
# and how many RNASeq matches exist to clinical data -- DIFFERENT ANALYSIS


# Import recurrence information (clinical follow up) 
# KEEP uuid, barcode, treatment outcome (first), last follow up, treatment outcome
# (follow up), new tumor (y/n), and days to get tumor -- not all these have info for
# all patients ; delete headers ; match the clinical follow up dataset barcode to 
# clinical dataset barcode; delete rows/patients with NaNs and delete non HPV- 
# patients with no recurrence data, within the clinical dataset
recurrence.data = as.data.frame(read.csv('/Users/Corey/Documents/WeaverLab/Clinical_Jan/Clinical/Biotab/nationwidechildrens.org_clinical_follow_up_v1.0_hnsc.txt', header = TRUE, sep = '\t'))
recurrence.data.revised = recurrence.data[-c(1,2), c(1,2,9,11,14,15,16)]
recurrence.data.matched = match(recurrence.data.revised$bcr_patient_barcode, matched.set.HPV$bcr_patient_barcode)
recurrence.data.deleted = which(is.na(recurrence.data.matched))
recurrence.data.new = recurrence.data.revised[-recurrence.data.deleted]
View(recurrence.data.new)

# For loop to make list of multiple recurrent tumors in patients - delete these
# patients from the dataset - keep only first recurrence of tumor in HPV- patients;
# recurrence data now only consists of 188 patients
# This loop needs to be optimized
duplicate.info = duplicated(recurrence.data.new$bcr_patient_barcode)
samples.to.delete = NULL
for (i in which(duplicate.info)) {
  if (recurrence.data.new$new_tumor_event_dx_days_to[i-1] > 0 & recurrence.data.new$new_tumor_event_dx_days_to[i] =='Not Applicable') {
    samples.to.delete = c(samples.to.delete, i)
} else if (recurrence.data.new$new_tumor_event_dx_days_to[i] > 0 & recurrence.data.new$new_tumor_event_dx_days_to[i-1] =='[Not Applicable]') {
    samples.to.delete = c(samples.to.delete,i-1)
} else if (recurrence.data.new$new_tumor_event_dx_days_to[i-1] > recurrence.data.new$new_tumor_event_dx_days_to[i] ) {
    samples.to.delete = c(samples.to.delete,i)
} else {samples.to.delete = c(samples.to.delete,i)}
  }
recurrence.data.new1 = recurrence.data.new[-samples.to.delete,]
View(recurrence.data.new1)
dim(recurrence.data.new1)

# Adding the recurrence data (days to tumor) to matched patient dataset by matching 
# barcodes, ordering the matrix, and adding a new column to the patient dataset 
# matrix - many NaNs introduced in new recurrence column (43/420 patients, only 
# 43/188 from recurrence dataset introduced to patients)
data.ordered = match(matched.set.HPV$bcr_patient_barcode, recurrence.data.new1$bcr_patient_barcode)
recurrence.data.new2 = recurrence.data.new1[data.ordered, ]
View(recurrence.data.new2)
dim(recurrence.data.new2)
matched.set.HPV$Recurrence = as.numeric(recurrence.data.new2$new_tumor_event_dx_days_to)

# Deleting patients that do not have RNAseq data and/or clinical data
# NOTE: not actually sure the point of matching the matrices in reverse, 
# but essentially removing NaNs from matched patient dataset (420--414 patients)

View(myCleanDataSet)
View(matched.set.HPV)
sum(!is.na(matched.set.HPV$Recurrence)) #43 patients with recurrence data
clinical.rows = match(rownames(myCleanDataSet), matched.set.HPV$bcr_patient_barcode)
clinical.rows.reverse = match(matched.set.HPV$bcr_patient_barcode, rownames(myCleanDataSet))
View(clinical.rows)
View(clinical.rows.reverse)
NAindex1 = which(is.na(clinical.rows))
NAindex2 = which(is.na(clinical.rows.reverse))

if (length(NAindex2)==0){
  matched.set.HPV.1 = matched.set.HPV
} else {
  matched.set.HPV.1 = matched.set.HPV[-NAindex2,]
}
View(matched.set.HPV)
dim(matched.set.HPV) # 420 rows
View(matched.set.HPV.1)
dim(matched.set.HPV.1) # 414 rows

# Match RNASeq dataset (520-414 patients after remove NaNs) to patient dataset by
# row names (aka barcodes)
data.counts.matched = myCleanDataSet[-NAindex1,]
clinical.rows.matched = match(rownames(data.counts.matched), matched.set.HPV.1$bcr_patient_barcode)
matched.set.HPV.rearranged = matched.set.HPV.1[clinical.rows.matched,]
View(matched.set.HPV.rearranged)
dim(matched.set.HPV.rearranged) # 414 rows

#Change categorical data to numerical in patient dataset
matched.set.HPV.rearranged$vital_status = as.numeric(recode(matched.set.HPV.rearranged$vital_status,"'Dead'=1;'Alive'=0"))
matched.set.HPV.rearranged$clinical_stage = as.numeric(recode(matched.set.HPV.rearranged$clinical_stage,"'Stage I'=1;'Stage
                                                       II'=2;'Stage III'=3;'Stage IV'=4;'Stage IVA'=4;'Stage IVB'=4;'Stage IVC'=4 "))
matched.set.HPV.rearranged$ajcc_pathologic_tumor_stage =
  as.numeric(recode(matched.set.HPV.rearranged$ajcc_pathologic_tumor_stage,"'Stage I'=1;'Stage II'=2;'Stage III'=3;'Stage
                    IV'=4;'Stage IVA'=4;'Stage IVB'=4;'Stage IVC'=4 "))
matched.set.HPV.rearranged$death_days_to = as.numeric(matched.set.HPV.rearranged$death_days_to)
matched.set.HPV.rearranged$lymph_nodes_examined_he_count = as.numeric(matched.set.HPV.rearranged$lymph_nodes_examined_he_count)
matched.set.HPV.rearranged$tumor_grade = as.numeric(recode(matched.set.HPV.rearranged$tumor_grade,"'G1'=1;'G2'=2;'G3'=3;'G4'=4"))
rownames(matched.set.HPV.rearranged) = matched.set.HPV.rearranged$bcr_patient_barcode

# Final cleaning steps for patient dataset (414 patients, 8 patient outcomes)
matched.set.HPV.rearranged.1 = matched.set.HPV.rearranged[,-1] # remove uuid column
View(matched.set.HPV.rearranged.1)
dim(matched.set.HPV.rearranged.1) # 414 rows, 8 columns

# Deleting genes that have less than 10 counts in more than 90% of the samples --
# Get final RNA Seq gene expression dataset
genes.to.delete = NULL
for (i in 1:ncol(data.counts.matched)) {
  if(length(data.counts.matched[,i][data.counts.matched[,i] < 10]) > 0.9*length(data.counts.matched[,i]))
  {genes.to.delete = c(genes.to.delete, i)
  }
}
final.data = data.counts.matched[,-genes.to.delete]
View(final.data)
dim(final.data) # 414 patients and 16350 genes

# Order the rows of matched.set.HPV.rearranged.1 so that
# they match those of final.data:
Patients=rownames(final.data)
traitRows = match(Patients, matched.set.HPV.rearranged.1$bcr_patient_barcode) # match barcodes in patient and expression dataset
datTraits = matched.set.HPV.rearranged.1[traitRows, -1] # remove additional barcode column and store in new variable
dim(datTraits)
View(datTraits)
sum(!is.na(datTraits$Recurrence)) #43 patients with recurrence data
rownames(datTraits) = matched.set.HPV.rearranged.1[traitRows, 1]
table(rownames(datTraits)==rownames(final.data)) # show that row names agree

#####################################################################################
### END CLEANING ; START WGCNA ANALYSIS - Dendrogram and Trait Heatmap
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
sampleTree = flashClust(dist(final.data), method = "complete");
sizeGrWindow(25,20) # set size of graphical window (inches)
#par sets graphical parameters
par(cex = 1); # 1.5 x magnification of graphical parameters
par(mar = c(0,6,2,0)) # numerical vector for 4 sides of plot 
plot(sampleTree, main = "Clustering to Determine Patient Relatedness", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
traitColors_total = numbers2colors(datTraits, signed = TRUE);
dim(traitColors_total)
#View(traitColors_total)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors_total,
                    groupLabels = colnames(datTraits),
                    main = "Dendrogram and Trait Heatmap")

"""
#Clustering Steps - Use depending on subjective knowledge of similarity between patients
#according to gene expression you want to achieve

# Plot a line to show the cut
abline(h = 2000000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2000000, minSize = 10)
table(clust) # 0=unassigned, 1=largest cluster, 2=next largest cluster
#?cutreeStatic
# clust 1 contains the samples above above line in sample tree
keepSamples = (clust==1)
datExpr = final.data[keepSamples, ]
datTraits.new = datTraits[keepSamples,]
View(datTraits.new)
dim(datTraits.new)
dim(datExpr)
# Re-cluster samples
#library(flashClust)
sampleTree2 = flashClust(dist(datExpr), method = "average")
dim(sampleTree2)
View(sampleTree2)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits.new, signed = FALSE);
### SOMETHING WRONG HERE - used code from own implementation
dim(traitColors)
View(traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits.new),
                    main = "Dendrogram and trait heatmap")
"""

#####################################################################################
### WGCNA ANALYSIS - Choosing a soft threshold to match scale independence and mean connectivity
#####################################################################################


# Choose a set of soft-thresholding powers (1:10 by 1, 12:20 by 2) and use the data
# to pick the soft threshold where SFT > 0.90 and connectivty seems to flatten
# Default Beta for Signed Network (Soft Threshold Power = 12)

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
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of different powers
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")

"""
# sample network based on squared Euclidean distance
# note that we transpose the data
A=adjacency(t(datExpr),type="distance")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold
thresholdZ.k=-2.5 # often -2.5
# the color vector indicates outlyingness (red)
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flahsClust or hclust
#library(flashClust)
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation:
# where red indicates high values
traitColors=data.frame(numbers2colors(datTraits.new,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits.new),"C",sep="")
datColors=data.frame(outlierC=outlierColor,traitColors)
View(outlierColor) #401
View(traitColors) #414
View(datTraits.new) #414
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")

# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
# the following 2 lines differ from what is written in the book
datExpr.new=datExpr[!remove.samples,]
datTraits.new1=datTraits.new[!remove.samples,]
View(datExpr.new)
View(datTraits.new1)
dim(datTraits.new1)
# Recompute the sample network among the remaining samples
A1=adjacency(t(datExpr.new),type="distance")
# Let's recompute the Z.k values of outlyingness
k=as.numeric(apply(A1,2,sum))-1
Z.k=scale(k)

###Additional Clustering Step from Oscar's Code - divided number of genes in half
# Pick soft threshold requires you to disableWGCNAtheads() because it does not work with parallel processing enabled

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
?pickSoftThreshold()
sft = pickSoftThreshold(datExpr.new, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower = 6
#Soft connectivity graph
k.dataOne= softConnectivity(datExpr.new,power=softPower)-1
?softConnectivity
# creates a scale free topology plot for data set
sizeGrWindow(9, 5)
par(mfrow = c(1,1));
cex1 = 0.9;
scaleFreePlot(k.dataOne, main=paste("data HNC, power=",softPower), truncated=F);
#Only uses the 8000 more variable genes. I compare between the 8000 and the ~16000 genes. The results are the same,
#but the computation is faster.
kCut = 8000 # number of most connected genes that will be considered
kRank = rank(-k.dataOne)
vardataOne=apply(datExpr.new,2,var)
# the most connected probesets with non-zero variance
restk = kRank <= kCut & vardataOne>0
datExpr.new.1 = datExpr.new[,restk]

# 16000 genes --> 8000 genes

"""
### Choose a Beta power - not specific to this analysis
# If using R through Graphical User Interface (i.e. RStudio), you may need
# to disable parallel processing with disableWGCNAThreads(). But, make sure it
# is enabled for all further processes
#disableWGCNAThreads()

# Choose a set of soft thresholding powers
#powers=c(1:20) # in practice this should include powers up to 20.
#View(powers)
# choose power based on SFT criterion
#?pickSoftThreshold
#sft=pickSoftThreshold(final.data,powerVector=powers)
#View(sft)

#enableWGCNAThreads()

#####################################################################################
### WGCNA ANALYSIS - Network Construction and Module Detection
#####################################################################################

# Automatic network construction and module detection by Dynamic Tree Cutting

dataExpr = final.data

# Set the threshold for merging closely related modules
mergingThresh = 0.90

# Automatic network construction and module detection on large expression datasets 
# in a block-wise manner (using blockwisemodules function). Input the cleaned
# expression dataset and group by pearson correlation. Create the correlation matrix
# using the signed network equation - beta power = 12. Minimum and maximum module sizes
# (correlations) set, and merge closely associated modules. Module colors output 
# in numerical form.

net = blockwiseModules(datExpr,corType="pearson",
                       maxBlockSize=16500,networkType="signed",power=12,minModuleSize=30,
                       mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                       pamRespectsDendro=FALSE,saveTOMFileBase="HNSCC_TOM")
moduleLabelsAutomatic=net$colors
# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# Obtain module eigengenes in a dataframe
MEsAutomatic=net$MEs
# Convert clinical dataset to dataframe
clinical.traits = as.data.frame(datTraits)

#recurrence
#names(recurrence)="recurrence"
# Next use this trait to define a gene significance variable
#GS.recurrence=as.numeric(cor(final.data,recurrence,use="p"))
# This translates the numeric values into colors
#GS.recurrenceColor=numbers2colors(GS.recurrence,signed=T)

blocknumber=1
# If looking at specific GS, add ,GS.recurrenceColor after moduleColorsAutomatic
datColors=data.frame(moduleColorsAutomatic)[net$blockGenes[[blocknumber]],]
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[blocknumber]],colors=datColors,
                    groupLabels=c("Module colors"),dendroLabels=FALSE,
                    hang=0.03,addGuide=TRUE,guideHang=0.05)

### Blockwise module detection for large networks - separate blocks

#Can try changing maxBlockSize to get fewer blocks
bwnet = blockwiseModules(datExpr,corType="pearson",
                         maxBlockSize=5500,networkType="signed",power=12,minModuleSize=30,
                         mergeCutHeight=mergingThresh,numericLabels=TRUE,saveTOMs=TRUE,
                         pamRespectsDendro=FALSE,saveTOMFileBase="HNSCCTOM-blockwise",verbose=3)

# Relabel blockwise modules so that their labels
# match those from our previous analysis
moduleLabelsBlockwise=matchLabels(bwnet$colors,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsBlockwise=labels2colors(moduleLabelsBlockwise)
# measure agreement with single block automatic procedure
mean(moduleLabelsBlockwise==moduleLabelsAutomatic)

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

### Manual, stepwise module detection

# We now calculate the weighted adjacency matrix, using the power 6:
A = adjacency(datExpr, type = "signed", power = 12)
# Digression: to define a signed network choose
#A = adjacency(datExprFemale, power = 12, type="signed")
#define a dissimilarity based on the topological overlap
dissTOM =TOMdist(A)
#hierarchical clustering
geneTree = flashClust(as.dist(dissTOM),method="average")
# here we define the modules by cutting branches
moduleLabelsManual1=cutreeDynamic(dendro=geneTree,distM=dissTOM,
                                  method="hybrid",deepSplit=2,pamRespectsDendro=F,minClusterSize=30)
# Relabel the manual modules so that their labels
# match those from our previous analysis
moduleLabelsManual2=
  matchLabels(moduleLabelsManual1,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsManual2=labels2colors(moduleLabelsManual2)

# Calculate eigengenes
MEList=moduleEigengenes(datExpr,colors=moduleColorsManual2)
MEs = MEList$eigengenes
# Add the weight to existing module eigengenes
#MET=orderMEs(cbind(MEs,weight))
MET=orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)

### Trying something from Oscar's code - using hclust to see if any difference
### in dynamic tree cutting

Ad = adjacency(final.data, power = 12)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(Ad, TOMType="signed")
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100;
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


### NEED TO DETERMINE TRAITS OF INTEREST!!! ###

### Comparing module detection methods

# automatically merge highly correlated modules
merge=mergeCloseModules(datExpr,moduleColorsManual2,
                        cutHeight=0.90)
# resulting merged module colors
moduleColorsManual3 = merge$colors
# eigengenes of the newly merged modules:
MEsManual = merge$newMEs
# Show the effect of module merging by plotting the
# original and merged module colors below the tree
datColors=data.frame(moduleColorsManual3,moduleColorsAutomatic,
                     moduleColorsBlockwise)
View(datColors)
plotDendroAndColors(geneTree,colors=datColors,
                    groupLabels=c("manual hybrid","single block","3 block"),
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

# Only manual hybrid method

datColors.MH=data.frame(moduleColorsManual3)
plotDendroAndColors(geneTree,colors=datColors.MH,
                    groupLabels="manual hybrid",
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

# check the agreement between manual and automatic module labels
mean(moduleColorsManual3==moduleColorsAutomatic) 

###Relating Modules to Physiological Traits/Outcomes

# Rename to moduleColors
moduleColorsReal.1 = moduleColorsManual3
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColorsReal.1, colorOrder)-1;
MEs.1 = MEsManual;
# Save module colors and labels for use in subsequent parts
save(MEs.1, moduleLabels, moduleColorsReal.1, geneTree, file = "HNSCC_022516.RData")
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

dim(datExpr)
dim(moduleColorsReal.1)
dim(datTraits)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColorsReal.1)$eigengenes
MEs = orderMEs(MEs0)
dim(MEs)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
View(moduleTraitPvalue)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 12, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               cex.lab = 0.7,
               main = paste("Module-trait relationships"))

### Looking a specific module-trait relationship -- Recurrence with certain color
### module - need to determine from heatmap

recurrence = as.data.frame(datTraits$Recurrence)
View(recurrence)
dim(recurrence)
dim(datExpr)
dim(MEs)
names(recurrence) = "recurrence"
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, recurrence, use="p"))
dim(geneTraitSignificance)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
dim(GSPvalue)
names(geneTraitSignificance) = paste("GS.", names(recurrence), sep="")
names(GSPvalue) = paste("p.GS.", names(recurrence), sep="")
#Bar plot that shows the mean gene significance of all modules with a defined trait
verboseBarplot(geneTraitSignificance$GS.recurrence, moduleColors, main='Module significance-recurrence ', xlab='',
               ylab='Mean GS',
               cex = 0.8, cex.axis = 1.2, cex.lab = 1., cex.main = 1.2,
               col=levels(factor(moduleColors)),KruskalTest = F, las=2,
               ylim=c(-0.3, 0.2))

module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(17, 17);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Recurrence",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

############# ADDITIONAL MODULE VISUALIZATION #######################

# Choose a module assignment
moduleColorsReal=moduleColorsManual3
dim(moduleColorsReal)
# Define numbers of genes and samples
nGenes = ncol(datExpr.new)
dim(datExpr.new)
dim(nGenes)
View(nGenes)
nSamples = nrow(datExpr.new)
dim(nSamples)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr.new,moduleColorsReal)$eigengenes
dim(MEs0)
MEsReal = orderMEs(MEs0)
dim(MEsReal)
dim(datTraits.new1)
dim(datExpr.new)
modTraitCor = cor(MEsReal, datExpr.new, use = "p")
View(modTraitCor)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
View(modTraitP)
#Since we have a moderately large number of modules and traits,
#a suitable graphical representation will help in reading
#the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",
                   signif(modTraitP, 1), ")", sep = "")
View(textMatrix)
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datExpr.new),
               yLabels = names(MEsReal), ySymbols = names(MEsReal),
               colorLabels =FALSE,colors=greenWhiteRed(50),textMatrix=textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Choose a module assignment
moduleColorsHNSCC=moduleColorsAutomatic
# Define numbers of genes and samples
View(final.data)
nGenes = ncol(final.data)
nSamples = nrow(final.data)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(final.data,moduleColorsHNSCC)$eigengenes
View(MEs0)
MEsHNSCC = orderMEs(MEs0)
View(MEsHNSCC)
View(datTraits)
modTraitCor = cor(MEsHNSCC, datTraits, use = "p")
View(modTraitCor)
modTraitP = corPvalueStudent(modTraitCor, nSamples)
View(modTraitP)
#Since we have a moderately large number of modules and traits,
#a suitable graphical representation will help in reading
#the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",
                   signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits),
               yLabels = names(MEsHNSCC), ySymbols = names(MEsHNSCC),
               colorLabels =FALSE,colors=greenWhiteRed(50),textMatrix=textMatrix,
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))


