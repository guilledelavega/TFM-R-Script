######################################################################################################################
################### TFM GUILLERMO DE LA VEGA BARRANCO - MÁSTER EN BIOINFORMÁTICA Y BIOESTADÍSTICA ####################
######################################################################################################################

#Loading required packages
library(biomaRt)
library(edgeR)
library(readxl)
library(ggbiplot)
library(ggplot2)
library(ggrepel)
library(EDASeq)
library(DescTools)
library(corrplot)
library(car)
library(visdat)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(caret)
library(ranger)
library(Boruta)
library(leaps)
library(DESeq2)
library(sva)
library(pheatmap)
library(future)
library(fgsea)
library(msigdbr)
library(gridExtra)



#Read all raw RNAseq data
HUGO_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_HUGO_28/Raw_CountMatrix_STAR_HUGO.tsv", header = TRUE, row.names = 1)

GIDE_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_GIDE_91/counts.tsv", header = TRUE, row.names = 1)

RIAZ_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_RIAZ_109/Raw_CountMatrix_STAR_RIAZ.tsv", header = TRUE, row.names = 1)

LIU_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_LIU&VANALLEN_164/Raw_CountMatrix_STAR_LIU.tsv", header = TRUE, row.names = 1)

VANALLEN_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_LIU&VANALLEN_164/Raw_CountMatrix_STAR_VANALLEN.tsv", header = TRUE, row.names = 1)

ABRIL_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_RIBAS_41/Raw_CountMatrix_STAR_RIBAS.tsv", header = TRUE, row.names = 1)

FREEMAN_RAW = read.table("/home/guvega/Documents/COHORTS PROJECT/RawCountMatrices/rawData_FREEMAN_48/counts.tsv", header = TRUE, row.names = 1)

#Load conversion table from ensemble IDs to gene names
Conversion_table_BioMart = read_csv("/home/guvega/Documents/Conversion_table_BioMart.csv")




#########################################################
# INDIVIDUAL COHORT ANALYSIS - 1) FILTERING COUNT MATRIX
#########################################################

####Filtering by: PAR_Y genes, ensembl IDs without gene names, duplicated genes

##Creating a list with all data together
RAW_COUNTS_list = list(
  HUGO = HUGO_RAW,
  GIDE = GIDE_RAW,
  RIAZ = RIAZ_RAW,
  LIU = LIU_RAW,
  VANALLEN = VANALLEN_RAW,
  ABRIL = ABRIL_RAW,
  FREEMAN = FREEMAN_RAW
)

#Empty list
filtered_counts_list = list()

#For loop for all the datasets
for (dataset_name in names(RAW_COUNTS_list)) {
  
  #Selecting each dataset
  CountMatrix = RAW_COUNTS_list[[dataset_name]]
  
  #Step 1: Total expression per gene
  rows = rowSums(CountMatrix)
  rows = sort(rows, decreasing = TRUE)
  rows_dataframe = as.data.frame(rows)
  
  #Step 2:Asign rownames to a new column 'Full_ensemble'
  rows_dataframe$Full_ensemble = rownames(rows_dataframe)
  
  #Step 3: Remove "_PAR_Y" suffix for some IDs
  rows_dataframe$Full_ensemble = sub("_PAR_Y", "", rows_dataframe$Full_ensemble)
  
  #Step 4: Change column names
  colnames(rows_dataframe) = c("Counts", "ensembl_gene_id_version")
  
  #Step 5: Join with Conversion_table_BioMart obtained from BioMart
  merged_data = merge(rows_dataframe, Conversion_table_BioMart, by = "ensembl_gene_id_version", all.x = TRUE)
  
  #Paso 6: Filter rows with invalid gene name
  merged_data = subset(merged_data, !is.na(hgnc_symbol) & hgnc_symbol != "")
  
  #Paso 7: Remove innecesary columns
  merged_data$entrezgene_id = NULL
  merged_data$ensembl_gene_id = NULL
  
  #Paso 8: order by gene expression and remove duplicates
  merged_data = merged_data %>% arrange(desc(Counts))
  merged_data_filtrated = merged_data[!duplicated(merged_data$hgnc_symbol), ]
  
  #Paso 9: Filter with selected genes
  genes_vector = merged_data_filtrated$ensembl_gene_id_version
  CountMatrix_filtrated = CountMatrix[genes_vector, ]
  
  #Paso 10: Add columns with correct id and gene names to the final matrix
  CountMatrix_filtrated$ensembl_gene_id_version = rownames(CountMatrix_filtrated)
  CountMatrix_filtrated = merge(CountMatrix_filtrated, Conversion_table_BioMart, by = "ensembl_gene_id_version", all.x = TRUE)
  
  #Paso 11: Remove innecesary columns
  CountMatrix_filtrated$entrezgene_id = NULL
  CountMatrix_filtrated$ensembl_gene_id = NULL
  
  #Paso 12: Assign hgnc_symbol as row names
  rownames(CountMatrix_filtrated) = CountMatrix_filtrated$hgnc_symbol
  
  #Paso 13: remove redundant columns
  CountMatrix_filtrated$ensembl_gene_id_version = NULL
  CountMatrix_filtrated$hgnc_symbol = NULL
  
  #Verify duplicated genes
  num_duplicados = sum(duplicated(rownames(CountMatrix_filtrated)))
  print(paste("Duplicated genes in", dataset_name, ":", num_duplicados))
  
  #Save results
  filtered_counts_list[[dataset_name]] = CountMatrix_filtrated
}






#################################################################
# INDIVIDUAL COHORT ANALYSIS - 2) EXPLORATORY DATA ANALYSIS EDA
#################################################################

#Assign colors to each dataset
cohort_colors = c("HUGO" = "#458B74", #aquamarine4
                   "GIDE" = "#8B2323", #brown4
                   "RIAZ" = "#548B54", #palegreen4 
                   "LIU" = "#3A5FCD", #royalblue3
                   "VANALLEN" = "#CD6600", #darkorange3 
                   "ABRIL" = "#8B3A62", #hotpink4
                   "FREEMAN" = "#8B5F65") #lightpink4

#Create an empty list to store DGE objects
dge_filtered_list = list()

#Iterate through each filtered dataset in 'filtered_counts_list' to perform PCA analysis
for (dataset_name in names(filtered_counts_list)) {
  
  #Extract the filtered count matrix
  CountMatrix_filtrated = filtered_counts_list[[dataset_name]]
  
  #Check for NA values and remove them if they exist
  if (any(is.na(CountMatrix_filtrated))) {
    CountMatrix_filtrated = na.omit(CountMatrix_filtrated)
  }
  
  #Create the DGEList object with edgeR
  dge = DGEList(counts = CountMatrix_filtrated)
  
  #Filter low expressed genes with filterByExpr
  keep = filterByExpr(dge)
  dge = dge[keep, , keep.lib.sizes = FALSE]
  
  #Save every object filtered in a list for possible further analysis
  dge_filtered_list[[paste0("dge_", dataset_name)]] = dge
  
  
  ## ---- PLOT 1: LIB SIZES PLOT ##
  png(filename = paste0("/home/guvega/Documents/COHORTS PROJECT/TFM/EDA_PLOTS/Library_sizes_", dataset_name, ".png"), width = 800, height = 600)
  par(mar = c(7, 4, 5, 2))  # c(bottom, left, top, right)
  
  # Plot library sizes
  barplot(dge$samples$lib.size / 1e06, las = 2, ann = FALSE, cex.names = 0.75)
  mtext(side = 2, text = "Library size (millions)", line = 3, cex = 1.5)
  title(paste("Barplot of library sizes - ", dataset_name), cex.main=2)
  dev.off()  
  
  
  
  ## ---- PLOT 2: DISTRIBUTION BOXPLOTS ##
  # Get log2 counts per million
  png(filename = paste0("/home/guvega/Documents/COHORTS PROJECT/TFM/EDA_PLOTS/Distribution_BoxPlot_", dataset_name, ".png"), width = 800, height = 600)
  logcounts = edgeR::cpm(dge,log=T)
  #Check distributions of samples using boxplots
  boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, xaxt = "n", cex.axis = 1.5)
  #Add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(logcounts),col="blue")
  title(paste("Boxplots of logCPMs - ", dataset_name), cex.main=2)
  dev.off() 
  
  
  
  
  ## ---- PLOT 3: PCA with labels (NO POINTS) ---- ##
  #Perform PCA analysis using the prcomp function on the transformed CPMs
  pca_analysis = prcomp(t(edgeR::cpm(dge, log = TRUE)))
  #Convert the PCA results to a data frame for ggplot2
  pca_data = data.frame(pca_analysis$x)
  pca_data$Sample = rownames(pca_data)
  
  #Add the Cohort column to identify the cohort
  pca_data$Cohort = dataset_name
  
  #Plot it without points
  pca_plot_no_points = ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
    geom_text_repel(size = 3) +  # Label size
    scale_color_manual(values = cohort_colors) +  # Assign colors
    ggtitle(paste("PCA Labels for", dataset_name)) +
    theme(
      plot.title = element_text(size = 18),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
  
  #Plot it
  print(pca_plot_no_points)
  ggsave(plot=pca_plot_no_points ,path="/home/guvega/Documents/COHORTS PROJECT/TFM/EDA_PLOTS/", filename = paste0("PCA_Label_", dataset_name, ".png"), width =12 , height =12 )
  
  
  
  
  ## ---- PLOT 4: PCA with points (POINTS) ---- ##
  pca_plot_points = ggplot(pca_data, aes(x = PC1, y = PC2, color = Cohort)) +
    geom_point(shape = 19, size = 8) +
    scale_color_manual(values = cohort_colors) +  #Assign colors
    ggtitle(paste("PCA Points for", dataset_name)) + 
    theme(
      plot.title = element_text(size = 25),  
      axis.title.x = element_text(size = 22), 
      axis.title.y = element_text(size = 22)
    )
  
  #Plot
  print(pca_plot_points)
  ggsave(plot=pca_plot_points ,path="/home/guvega/Documents/COHORTS PROJECT/TFM/EDA_PLOTS/", filename = paste0("PCA_Points_", dataset_name, ".png"), width =12 , height =12 )
  
  
  
  
  ## ---- PLOT 5: CORRELATION MATRIX ##
  cpm_matrix = as.data.frame(edgeR::cpm(dge, log=TRUE))
  
  cor_matrix = cor(cpm_matrix)
  
  pheatmap =pheatmap(cor_matrix,
           show_rownames = F,
           show_colnames = F,
           cluster_rows = T,
           cluster_cols = T,
           main = paste("Correlation Matrix Heatmap - ", dataset_name))
  ggsave(plot=pheatmap ,path="/home/guvega/Documents/COHORTS PROJECT/TFM/EDA_PLOTS/", filename = paste0("CorMatrix_", dataset_name, ".png"), width =12 , height =12 )
  
}







#################################################################
########         METADATA PREPARATION                    ########
#################################################################

#Load the clinical data
SuperMetaData = read_excel("/home/guvega/Documents/COHORTS PROJECT/MetaData/SuperMetaData.xlsx")

#Filtering selected outliers
SuperMetaData_filtered = SuperMetaData %>% 
  filter(!Run %in% c("SRR5088872", "SRR5088840", "SRR5088841", "SRR3083584", "SRR10842358", "SRR10900579", "SRR10900547", "SRR10842068", "SRR10900593", "ERR2208971"))

#Replace "unknown" for NA, function
replace_unknown = function(x) {
  if (is.character(x)) {
    x[x == "Unknown"] = NA
  }
  return(x)
}

#Apply NA function to the data
SuperMetaData_filtered = as.data.frame(lapply(SuperMetaData_filtered, replace_unknown))


###### Prepare response column - minor changes ############

#First, the patient with X in VanAllenn et al. The overall survival is so high that it is considered a responder
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST == "X", "Responder", SuperMetaData_filtered$Response.GVB)

#The ones with UNK = NA
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST == "UNK", "NA", SuperMetaData_filtered$Response.GVB)

#CR (Complete Response), PR (Partial Response), PR/CR and PRCR = Responders
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST %in% c("CR", "PR", "PR/CR", "PRCR"), "Responder", SuperMetaData_filtered$Response.GVB)

#PD (Progressive Disease), SD (Stable Disease), PD/SD = Non Responders
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST %in% c("PD", "SD", "PD/SD"), "Non Responder", SuperMetaData_filtered$Response.GVB)

#Mixed Response = Mixed response
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST == "MR", "Mixed Response", SuperMetaData_filtered$Response.GVB)

#Change UNK response with low survival = NonResponder
SuperMetaData_filtered$Response.GVB = ifelse(SuperMetaData_filtered$RECIST == "UNK" & SuperMetaData_filtered$Response.Paper == "Non-Responder", "Non Responder",SuperMetaData_filtered$Response.GVB)

#Filter to remove mixed response and NA patients
SuperMetaData_filtered = SuperMetaData_filtered %>%
  filter(SuperMetaData_filtered$Response.GVB %in% c("Responder", "Non Responder"))


# Remove columns that are not interesting or are largely incomplete, and change class
SuperMetaData_filtered$BioProject = NULL
SuperMetaData_filtered$BioSample = NULL
SuperMetaData_filtered$Experiment = NULL
SuperMetaData_filtered$SRA.Study =NULL
SuperMetaData_filtered$Removed.Patients.from.Papers = NULL
SuperMetaData_filtered$Removed.Patients.by.me = NULL
SuperMetaData_filtered$Specific_treatment = NULL
SuperMetaData_filtered$Progression.Free.Survival..Days. = NULL
SuperMetaData_filtered$site.biopsy = NULL
SuperMetaData_filtered$IntraCohort.Batch = NULL
SuperMetaData_filtered$NF1.mutated = NULL
SuperMetaData_filtered$NRAS.mutated = NULL
SuperMetaData_filtered$KIT.mutated = NULL
SuperMetaData_filtered$KRAS.mutated = NULL
SuperMetaData_filtered$PTEN = NULL
SuperMetaData_filtered$TP53 = NULL
SuperMetaData_filtered$Response.Paper = NULL
SuperMetaData_filtered$RECIST = NULL
SuperMetaData_filtered$Patient.no. = NULL

#Prepare defined class of the variables
SuperMetaData_filtered[sapply(SuperMetaData_filtered, is.character)] = 
  lapply(SuperMetaData_filtered[sapply(SuperMetaData_filtered, is.character)], as.factor)

SuperMetaData_filtered$Run = as.character(SuperMetaData_filtered$Run)
SuperMetaData_filtered$Age = as.numeric(SuperMetaData_filtered$Age)
SuperMetaData_filtered$Overall.Survival..Days. = as.numeric(SuperMetaData_filtered$Overall.Survival..Days.)

summary(SuperMetaData_filtered)


#Check missing data
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Metadata.png", width = 1600, height = 1000)
vis_miss(SuperMetaData_filtered) +
  theme(
    axis.text.x = element_text(size = 30, angle = 80),  # Increase x axis
    axis.text.y = element_text(size = 22), #Increase y axis
    axis.title.y = element_text(size = 22),  # Increase x axis title
    legend.text = element_text(size=22))
dev.off() 

write_csv(SuperMetaData_filtered, "/home/guvega/Documents/COHORTS PROJECT/MetaData/SuperMetaData_filtered.csv")


#Remove non cutaneous melanoma
SuperMetaData_filtered_2.0 = subset(SuperMetaData_filtered, !(Primary.site %in% c("Acral", "Mucosal", "Uveal")))


#Remove non-immune complete treatments
SuperMetaData_filtered_2.0 = subset(SuperMetaData_filtered_2.0, !(Treatment %in% "PDL1+MEKi"))
SuperMetaData_filtered_2.0$Treatment = droplevels(SuperMetaData_filtered_2.0$Treatment)
summary(SuperMetaData_filtered_2.0)


#Excluding huge missing variables 
SuperMetaData_filtered_2.0 = SuperMetaData_filtered_2.0 %>%
  dplyr::select(-Primary.site, -previous_targeted.therapies, -M_stage, -Sex, - Age)

#Check data
vis_miss(SuperMetaData_filtered_2.0) 
dim(SuperMetaData_filtered_2.0)

#Removing remaining samples with any NA
SuperMetaData_filtered_2.0 = na.omit(SuperMetaData_filtered_2.0)
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Metadata_prepared.png", width = 1600, height = 1000)
vis_miss(SuperMetaData_filtered_2.0) +
  theme(
    axis.text.x = element_text(size = 30, angle = 80),  # Increase x axis
    axis.text.y = element_text(size = 22), #Increase y axis
    axis.title.y = element_text(size = 22),  # Increase y axis title
    legend.text = element_text(size=22))
dev.off()  

write_csv(SuperMetaData_filtered_2.0, "/home/guvega/Documents/COHORTS PROJECT/MetaData/SuperMetaData_filtered_further_noNAs.csv")


#Establish metadata rownames
Samples = SuperMetaData_filtered_2.0$Run
SuperMetaData_filtered_2.0$Run = NULL
rownames(SuperMetaData_filtered_2.0) = Samples








#################################################################
########         DATA COMBINATION
#################################################################

##### PCA #######

#Combine all raw data
CountMatrixRaw_Combined = cbind(filtered_counts_list$HUGO, filtered_counts_list$GIDE, filtered_counts_list$RIAZ, filtered_counts_list$LIU,
                                filtered_counts_list$VANALLEN, filtered_counts_list$ABRIL, filtered_counts_list$FREEMAN)

#Remove samples based on metadata
CountMatrixRaw_Combined_filtered = CountMatrixRaw_Combined[, colnames(CountMatrixRaw_Combined) %in% rownames(SuperMetaData_filtered_2.0)]

#Check
all(rownames(SuperMetaData_filtered_2.0) %in% colnames(CountMatrixRaw_Combined_filtered))


#Create dge object
dge = DGEList(counts = CountMatrixRaw_Combined_filtered, samples = SuperMetaData_filtered_2.0)



############# FILTERING ############
keep = filterByExpr(dge)
dge = dge[keep, , keep.lib.sizes=F]




#####  PCA using ggbiplot ###########
## POINTS - Cohort ##
pca_analysis = prcomp(t(edgeR::cpm(dge, log=TRUE)))

pca_plot = ggbiplot::ggbiplot(pca_analysis, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Cohort, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Combined_PCA_Cohort.png", width = 1300, height = 1000)
pca_plot
dev.off()



## POINTS - RNAseq method ##
pca_plot = ggbiplot::ggbiplot(pca_analysis, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$RNAseq.method, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Combined_PCA_RNAseqMethod.png", width = 1300, height = 1000)
pca_plot
dev.off()



## POINTS - Seq Ins ##
pca_plot = ggbiplot::ggbiplot(pca_analysis, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Sequencing.Instrument, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Combined_PCA_SeqInstru.png", width = 1300, height = 1000)
pca_plot
dev.off()



##### COR MATRIX #######
cpm_matrix = as.data.frame(edgeR::cpm(dge, log=TRUE))

cor_matrix = cor(cpm_matrix)

#Make annotations
annotation = data.frame(
  RNAseq_Method = SuperMetaData_filtered_2.0$RNAseq.method,
  Cohort = SuperMetaData_filtered_2.0$Cohort,
  Treatment =SuperMetaData_filtered_2.0$Treatment,
  Response = SuperMetaData_filtered_2.0$Response.GVB,
  Sequencing.Instrument = SuperMetaData_filtered_2.0$Sequencing.Instrument,
  vital_status = SuperMetaData_filtered_2.0$vital_status,
  biopsy_time = SuperMetaData_filtered_2.0$biopsy_time,
  Previous_CTLA4 = SuperMetaData_filtered_2.0$Previous_CTLA4,
  BRAF.mutated = SuperMetaData_filtered_2.0$BRAF.mutated)

row.names(annotation) = rownames(SuperMetaData_filtered_2.0)


#Plot and save
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Combined_CorMatrix.png", width = 1300, height = 1000)
pheatmap::pheatmap(cor_matrix,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation,
         fontsize_annotation = 14,
         fontsize = 12)
dev.off()










#################################################################
########         REMOVING BATCH EFFECT  #########################
#################################################################
#https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

#There are two functions, Combat and Combat-seq. The last one is specifically created for RNA-seq
#Remember, this is made for visualization or further steps as clustering, machine learning methods, etc
#However, for RNAseq and DE, we would take raw counts! And we will include the batch variable in the design formula



############## COMBAT-SEQ - Cohort ######################
plan("multisession", workers = parallel::detectCores() - 2)
options(future.globals.maxSize= 8 * 1024^3)

corrected_COMBATseq_cohort_counts = as.data.frame(ComBat_seq(dge$counts, batch = SuperMetaData_filtered_2.0$Cohort))
write.csv(corrected_COMBATseq_cohort_counts, "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/corrected_COMBATseq_cohort_counts.csv", row.names = TRUE)

# Return to normal
future::plan("sequential")

##CORRELATION PLOT
cor_matrix_Combatseq_Cohort = cor(corrected_COMBATseq_cohort_counts)

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_CorMatrix_COMBATseq_Cohort.png", width = 1300, height = 1000)
pheatmap::pheatmap(cor_matrix_Combatseq_Cohort,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation,
         fontsize_annotation = 14,
         fontsize = 12)
dev.off()


##PCA plots
#Labeling Cohort 
dge_Combat_Cohort = DGEList(counts = corrected_COMBATseq_cohort_counts)
pca_analysis_Combat_Cohort = prcomp(t(edgeR::cpm(dge_Combat_Cohort, log=TRUE)))

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_Cohort.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_Cohort, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Cohort, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22), 
  ) 
dev.off()

#Labeling RNAseq Method
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_Cohort_labelingRNA.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_Cohort, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$RNAseq.method, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()


#Labeling Sequencing instrument
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_Cohort_labelingSeqInstru.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_Cohort, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Sequencing.Instrument, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()




############## COMBAT-SEQ - RNAseq Method ######################
plan("multisession", workers = parallel::detectCores() - 2)
options(future.globals.maxSize= 8 * 1024^3)

corrected_COMBATseq_RNAseqMethod_counts = as.data.frame(ComBat_seq(dge$counts, batch = SuperMetaData_filtered_2.0$RNAseq.method))
write.csv(corrected_COMBATseq_RNAseqMethod_counts, "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/corrected_COMBATseq_RNAseqMethod_counts.csv", row.names = TRUE)

# Return to normal
future::plan("sequential")

##CORRELATION PLOT
cor_matrix_Combatseq_RNAseqMethod = cor(corrected_COMBATseq_RNAseqMethod_counts)

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_CorMatrix_COMBATseq_RNAseqMethod.png", width = 1300, height = 1000)
pheatmap::pheatmap(cor_matrix_Combatseq_RNAseqMethod,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation,
         fontsize_annotation = 14,
         fontsize = 12)
dev.off()


##PCA PLOTS
#Labeling Cohort 
dge_RNAseqMethod = DGEList(counts = corrected_COMBATseq_RNAseqMethod_counts)
pca_analysis_Combat_RNAseqMethod = prcomp(t(edgeR::cpm(dge_RNAseqMethod, log=TRUE)))

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_RNAseqMethod.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_RNAseqMethod, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Sequencing.Instrument, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()

#Labeling RNASeq method 
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_RNAseqMethod_labelingCohort.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_RNAseqMethod, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Cohort, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()


#Labeling Sequencing instrument 
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_RNAseqMethod_labelingSeqInstru.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_RNAseqMethod, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Sequencing.Instrument, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22), 
  ) 
dev.off()




############## COMBAT-SEQ - SEQUENCING INSTRUMENT ######################
plan("multisession", workers = parallel::detectCores() - 2)
options(future.globals.maxSize= 8 * 1024^3)

corrected_COMBATseq_SeqInstrument_counts = as.data.frame(ComBat_seq(dge$counts, batch = SuperMetaData_filtered_2.0$Sequencing.Instrument))
write.csv(corrected_COMBATseq_SeqInstrument_counts, "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/corrected_COMBATseq_SeqInstrument_counts.csv", row.names = TRUE)

# Return to normal
future::plan("sequential")

##CORRELATION PLOT
cor_matrix_Combatseq_SeqInstrument = cor(corrected_COMBATseq_SeqInstrument_counts)

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_CorMatrix_COMBATseq_SeqInstrument.png", width = 1300, height = 1000)
pheatmap::pheatmap(cor_matrix_Combatseq_SeqInstrument,
         show_rownames = F,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = annotation,
         fontsize_annotation = 14,
         fontsize = 12)
dev.off()


##PCA PLOTS
#Labeling Cohort 
dge_SeqInstrument = DGEList(counts = corrected_COMBATseq_SeqInstrument_counts)
pca_analysis_Combat_SeqInstrument = prcomp(t(edgeR::cpm(dge_SeqInstrument, log=TRUE)))

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_SeqInstrument.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_SeqInstrument, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$RNAseq.method, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()

#Labeling RNASeq method 
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_SeqInstrument_labelingCohort.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_SeqInstrument, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$Cohort, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22), 
  ) 
dev.off()


#Labeling Sequencing instrument 
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/COMBATseq/Combined_PCA_COMBATseq_SeqInstrument_labelingRNAseqMethod.png", width = 1300, height = 1000)
ggbiplot::ggbiplot(pca_analysis_Combat_SeqInstrument, var.axes = FALSE,  groups = SuperMetaData_filtered_2.0$RNAseq.method, point.size = 4) +
  ggtitle("PCA") +  
  theme(
    plot.title = element_text(size = 30),  
    axis.title.x = element_text(size = 25), 
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size=32),
    legend.text = element_text(size = 22),  
  ) 
dev.off()







#####################################################
## METADATA VALIDATION - COLLINEARITY ################
######################################################

#### CATEGORICAL VARIABLES #####
# Filter only factor variables
categoricas = SuperMetaData_filtered_2.0[, sapply(SuperMetaData_filtered_2.0, is.factor)]

#Initialize an empty matrix to store Cramér's V values
n = ncol(categoricas)
matriz_cramer = matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(categoricas), colnames(categoricas)))

#Loop to calculate Cramér's V for each pair of categorical variables
for (i in 1:n) {
  for (j in i:n) {
    #Remove NAs before
    tabla = table(na.omit(categoricas[, c(i, j)]))
    matriz_cramer[i, j] = CramerV(tabla)
    matriz_cramer[j, i] = matriz_cramer[i, j]  
  }
}
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/CorPlot_CRAMER.png", width = 1600, height = 1000)
colores = colorRampPalette(c("white", "firebrick2"))(200)
corrplot(matriz_cramer, method = "number", type = "upper", tl.cex = 1.8, col = colores, tl.col = "black", is.corr = FALSE, cl.lim = c(0, 1), cl.cex = 2, number.cex = 3)
dev.off()  

#Samples before removing NAs
total_observaciones_iniciales = nrow(categoricas)

#Remove NAs
categoricas_sin_na = na.omit(categoricas)

#Samples after removing NA
total_observaciones_finales = nrow(categoricas_sin_na)

#Sample removed
observaciones_eliminadas = total_observaciones_iniciales - total_observaciones_finales

#Results
cat("Observaciones iniciales:", total_observaciones_iniciales, "\n")
cat("Observaciones finales después de eliminar NAs:", total_observaciones_finales, "\n")
cat("Observaciones eliminadas:", observaciones_eliminadas, "\n")                                     




#### NUMERICAL VARIABLES - Spearman #####
variables_numericas = SuperMetaData_filtered_2.0[, sapply(SuperMetaData_filtered_2.0, is.numeric)]

#Test normality
shapiro_bases = shapiro.test(variables_numericas$Bases)
shapiro_survival = shapiro.test(variables_numericas$Overall.Survival..Days.)

#Cat normality test results
cat("Prueba de Shapiro-Wilk para Bases: p-value =", shapiro_bases$p.value, "\n")
cat("Prueba de Shapiro-Wilk para Overall.Survival..Days.: p-value =", shapiro_survival$p.value, "\n")
#No normally distributed because <0,05 ---> spearman 


# Select only numerical variables
variables_numericas = SuperMetaData_filtered_2.0[, sapply(SuperMetaData_filtered_2.0, is.numeric)]

# Make correlation matrix between numerical samples with spearman
correlacion_numericas = cor(variables_numericas, use = "pairwise.complete.obs", method = "spearman")

#Plot it
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/CorPlot_SPEARMAN.png", width = 1600, height = 1000)
colores = colorRampPalette(c("dodgerblue3", "white", "firebrick2"))(200)
corrplot(correlacion_numericas, method = "number", type = "upper", tl.cex = 2.5, col = colores, tl.col = "black", cl.cex = 2.5, number.cex = 3)
dev.off()  # Cerrar el dispositivo de gráficos

#Check missing data for each pair
n = ncol(variables_numericas)
observaciones_validas = matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(variables_numericas), colnames(variables_numericas)))

#Calculate the number of valid observations for each pair of variables
for (i in 1:n) {
  for (j in i:n) {
    # Extract the pair of variables
    datos_pareja = variables_numericas[, c(i, j)]
    
    #Remove rows with NA only for that pair
    datos_pareja_sin_na = na.omit(datos_pareja)
    
    #Store the number of valid observations for this pair
    observaciones_validas[i, j] = nrow(datos_pareja_sin_na)
    observaciones_validas[j, i] = observaciones_validas[i, j] 
  }
}

#Print the matrix
print("Número de observaciones usadas para cada par de variables:")
print(observaciones_validas)




#### NUMERICAL VARIABLES - VIF #####
variables_numericas_2 = SuperMetaData_filtered_2.0[, c("Response.GVB", "Bases", "Overall.Survival..Days.")]

logistic_model = glm(Response.GVB ~ Bases + Overall.Survival..Days., 
                        data = variables_numericas_2, 
                        family = binomial)

vif_values = vif(logistic_model)
print(vif_values)









###############################################
### VARIABLE SELECTION METHODS #################
################################################

#### CRITERION-BASED METHODS ####
#https://www.statology.org/regsubsets-in-r/
#Linear Models with R - Julian J. Faraway

#Remove technical variables
SuperMetaData_filtered_for_decisionTree = SuperMetaData_filtered_2.0 %>%
  dplyr::select(-Sequencing.Instrument, -Bases, -Cohort, - RNAseq.method, -biopsy_time, -vital_status, -Overall.Survival..Days.)

#Exhaustive search to find the best regression model
b = regsubsets(Response.GVB ~ ., data=SuperMetaData_filtered_for_decisionTree)
rs = summary(b)

#The stars (*) at the bottom of the output indicate which predictor variables belong in the best regression model for each possible model with a different number of predictor variables.
rs$outmat
rs$which

## AIC ##
n = nrow(SuperMetaData_filtered_for_decisionTree)
k = length(rs$rss)  # Number of predictor variables
p = k + 1           # Number of parameters
(AIC = n*log(rs$rss/n) + (2:p)*2)

plot(1:k, AIC, ylab="AIC", xlab="Number of Predictors", axes=F)
box(); axis(1,at=1:k); axis(2)
rs$outmat[4,] #Check which variables predict the best model



## R2 ##
rs$adjr2
plot(1:k,rs$adjr2, xlab="Number of Predictors", ylab="Adjusted R-square", axes = F)
box(); axis(1,at=1:k); axis(2)
which.max(rs$adjr2)

rs$outmat[4,] #Check which variables predict the best model



## Mallow’s Cp ##
rs$cp
plot(2:p,rs$cp, xlab="No. of parameters", ylab="Cp Statistic", axes = F)
box(); axis(1,at=2:p); axis(2)
abline(a=0,b=1)

rs$outmat[5,]



## PLOT ALL TOGETHER ##
# Make first plot (AIC)
plot1 = qplot(1:k, AIC, geom="point") +
  labs(x="No. of Parameters", y="AIC") +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),        
    axis.title = element_text(size = 14),       
    plot.title = element_text(size = 16, face = "bold"),  
    panel.grid = element_blank(),               
    axis.line = element_line(color = "black", linewidth = 0.8)  
  )

# Make second plot (Adjusted R-square)
plot2 = qplot(1:k, rs$adjr2, geom="point") +
  geom_point(size = 4) +
  labs(x="No. of Parameters", y="Adjusted R-square") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),        
    axis.title = element_text(size = 14),       
    plot.title = element_text(size = 16, face = "bold"), 
    panel.grid = element_blank(),               
    axis.line = element_line(color = "black", linewidth = 0.8)  
  )

# Make third plot (Cp Statistic)
plot3 = qplot(2:p, rs$cp, geom="point") +
  geom_point(size = 4) +
  labs(x="No. of Parameters", y="Cp Statistic") +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),        
    axis.title = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold"),  
    panel.grid = element_blank(),               
    axis.line = element_line(color = "black", linewidth = 0.8)  
  )
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/VariableSelection/CriterionBased.png", width = 700, height = 1000)
grid.arrange(plot1, plot2, plot3, ncol=1)
dev.off()






#### MACHINE LEARNING METHODS

###### FIRST APPROACH: RPART - CART: Classification And Regression Trees
#https://rpubs.com/jboscomendoza/arboles_decision_clasificacion


#Training and test set
set.seed(1998)  # Set seed
n_total = nrow(SuperMetaData_filtered_for_decisionTree)
training_indices = sample(1:n_total, size = 0.7 * n_total)
training_set = SuperMetaData_filtered_for_decisionTree[training_indices, ]  
test_set = SuperMetaData_filtered_for_decisionTree[-training_indices, ]   

#Model
first_tree = rpart(formula = Response.GVB ~ ., data = training_set)

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/VariableSelection/decisiontree.png", width = 400, height = 400)
arbol_decision = rpart.plot(first_tree)
dev.off()

# Show variable importance
print(first_tree$variable.importance)

#Make predictions with test_set
prediccion_1 = predict(first_tree, newdata = test_set, type = "class")
accuracy = mean(prediccion_1 == test_set$Response.GVB)
print(accuracy)

#Check confusion matrix
confusionMatrix(prediccion_1, test_set[["Response.GVB"]])




###### SECOND APPROACH: BORUTA - RANDOM FOREST
#Make the model
boruta_output = Boruta(Response.GVB ~ ., data = SuperMetaData_filtered_for_decisionTree, doTrace = 2)
print(boruta_output)

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/VariableSelection/Boruta.png", width = 1200, height = 1000)
plot(boruta_output, 
     cex.axis=1,  # Increase the axis text size
     las=2,       # Rotate axis labels for better readability
     xlab="", 
     main="Variable Importance",
     cex.lab=1.5, # Increase the label size (if labels are added)
     cex.main=1.5 # Increase the title size
     
)
dev.off()











###############################################
########      MODEL SELECTION           #######
###############################################

######### DESEQ2 ###########
#Check everything is correct
#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order.

#As this is already filtered, we will use this dataframe for further analysis
CountMatrixRaw_Combined_filtered_corrected = as.data.frame(dge$counts)

#Check
all(rownames(SuperMetaData_filtered_2.0) == colnames(CountMatrixRaw_Combined_filtered_corrected))
CountMatrixRaw_Combined_filtered_corrected = CountMatrixRaw_Combined_filtered_corrected[, rownames(SuperMetaData_filtered_2.0)]
all(rownames(SuperMetaData_filtered_2.0) == colnames(CountMatrixRaw_Combined_filtered_corrected))

#Create the object
dds = DESeqDataSetFromMatrix(countData = CountMatrixRaw_Combined_filtered_corrected,
                              colData = SuperMetaData_filtered_2.0,
                              design = ~ Response.GVB + RNAseq.method + Sequencing.Instrument + Treatment)

#Make the factor levels
dds$Response.GVB = factor(dds$Response.GVB, levels = c("Responder","Non Responder"))
dds$Response.GVB = relevel(dds$Response.GVB, ref = "Responder")

#DE analysis
dds = DESeq(dds)
res = results(dds, contrast = c("Response.GVB", "Non Responder", "Responder")) #The last one would be the reference!
res

#Save the data and apply FDR correction
RES_table = as.data.frame(res)
RES_table$Genes = rownames(RES_table)
RES_table$FDR = p.adjust(RES_table$pvalue, method="BH")

write.csv(RES_table, "/home/guvega/Documents/COHORTS PROJECT/TFM/DE_Results/RES_table.csv")

#Count Matrix normalized by Deseq2
DESeq2_count_matrix =as.data.frame(counts(dds, normalized=TRUE))




######### edgeR - QLM ###########
y = DGEList(counts=CountMatrixRaw_Combined_filtered_corrected, samples = SuperMetaData_filtered_2.0)

# Create TMM normalization factors
y = calcNormFactors(y, method = "TMM")

#Levels
y$samples$Response.GVB = as.factor(y$samples$Response.GVB)
y$samples$Response.GVB = relevel(y$samples$Response.GVB, ref = "Responder")

#Design matrix
designMat = model.matrix(~Response.GVB + RNAseq.method + Sequencing.Instrument + Treatment, data = y$samples) 

#DE analysis
y = estimateDisp(y, design)

fit = glmQLFit(y, designMat, coef= "Response.GVBNon Responder")

qlf = glmQLFTest(fit)

#Save data
edgeR_results= data.frame(qlf$table)
edgeR_results$FDR = p.adjust(edgeR_results$PValue, method="BH")
edgeR_results$Genes = rownames(edgeR_results)

write.csv(edgeR_results, "/home/guvega/Documents/COHORTS PROJECT/TFM/DE_Results/edgeR_results.csv")




######### edgeR - exactTest ###########
# Create TMM normalization factors
y = DGEList(counts = CountMatrixRaw_Combined_filtered_corrected, samples = SuperMetaData_filtered_2.0, group = SuperMetaData_filtered_2.0$Response.GVB)
y = calcNormFactors(y, method = "TMM")

#Desing matrix
SuperMetaData_filtered_2.0$Response.GVB = factor(SuperMetaData_filtered_2.0$Response.GVB)

#Levels
SuperMetaData_filtered_2.0$Response.GVB = relevel(SuperMetaData_filtered_2.0$Response.GVB, ref = "Responder") 
designMat = model.matrix(~ Response.GVB + RNAseq.method + Sequencing.Instrument + Treatment, data = SuperMetaData_filtered_2.0) 

#DE analysis
y = estimateDisp(y, designMat)

etest = exactTest(y, pair = c("Responder", "Non Responder"))

#Save data
exacttest_results = etest$table
exacttest_results$FDR = p.adjust(exacttest_results$PValue, method="BH")
exacttest_results$Genes = rownames(exacttest_results)

write.csv(exacttest_results, "/home/guvega/Documents/COHORTS PROJECT/TFM/DE_Results/exacttest_results.csv")




######### limma-voom ###########
y = DGEList(counts=CountMatrixRaw_Combined_filtered_corrected, samples = SuperMetaData_filtered_2.0)

# Create TMM normalization factors
y = calcNormFactors(y, method = "TMM")
designMat = model.matrix(~Response.GVB + RNAseq.method + Sequencing.Instrument + Treatment, data = y$samples) 

#DE analysis
v = voom(y, designMat, plot = T)
vfit = lmFit(v, designMat)
efit = eBayes(vfit)

#Save data
limma_results = topTable(efit, adjust.method = "BH", coef = "Response.GVBNon Responder", number = Inf)
limma_results$Genes = rownames(limma_results)
write.csv(limma_results, "/home/guvega/Documents/COHORTS PROJECT/TFM/DE_Results/limma_results.csv")




##########  Validation of Statistical Models and Dispersion Estimation  ######### 
#Dispersion plots

#DeSeq2
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Model_Dispersions/DeSeq2.png", width = 800, height = 800)
par(cex.axis = 1.5,   
    cex.lab = 1.8,     
    cex.main = 2,      
    mar = c(5, 5, 4, 2) 
)
plotDispEsts(dds,
             main = "Dispersion Estimates (DESeq2)",  
             xlab = "Mean of Normalized Counts",      
             ylab = "Dispersion")                    
dev.off()




# edgeR
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Model_Dispersions/edgeR.png", width = 800, height = 800)
par(cex.axis = 1.5,    
    cex.lab = 1.8,     
    cex.main = 2,      
    mar = c(5, 5, 4, 2) 
)
plotBCV(y, xlab = "Mean of Normalized Counts", ylab = "Dispersion Estimates (edgeR)", main = "Dispersion Estimates (edgeR)")
dev.off()




#limma
png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Model_Dispersions/limma.png", width = 800, height = 800)
par(cex.axis = 1.5,    
    cex.lab = 1.8,     
    cex.main = 2,      
    mar = c(5, 5, 4, 2) 
)
w = voom(y, designMat, plot = TRUE)
dev.off()




#Cross-method reproducibility:
#Examine the overlap in DEGs between methods. A high overlap suggests consistency, while large discrepancies might indicate model-specific biases.
#Use correlation analysis of log fold changes (logFC) from the three methods to assess concordance.

# Subset DEGs based on pvalue < 0.05
DESeq2_DEGs = subset(RES_table, pvalue < 0.05)
edgeR_DEGs = subset(exacttest_results, Pvalue < 0.05)
limma_DEGs = subset(limma_results, P.Value < 0.05)

#Merge data
merged_results = merge(DESeq2_DEGs[, c("Genes", "log2FoldChange")],
                        edgeR_DEGs[, c("Genes", "logFC")],
                        by = "Genes", all = TRUE)
merged_results = merge(merged_results,
                        limma_DEGs[, c("Genes", "logFC")],
                        by = "Genes", all = TRUE)
colnames(merged_results) = c("Genes", "DESeq2", "edgeR", "limma")


# Pairwise correlation
cor_matrix=cor(merged_results[, 2:4], use = "complete.obs", method = "pearson")

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Model_Dispersions/Corplot_methods.png", width = 800, height = 800)
corrplot(cor_matrix, 
         method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),  # Blue to red gradient
         addCoef.col = "black",  # Add correlation values in black
         tl.col = "black",  # Label color
         tl.srt = 45,  # Rotate labels
         tl.cex = 2.5,  # Increase text size for labels
         number.cex = 3,
         cex.main=3,
         title = "Correlation of LogFC across methods", 
         mar = c(0, 0, 1, 0))  # Adjust margin
dev.off()

png(filename = "/home/guvega/Documents/COHORTS PROJECT/TFM/Model_Dispersions/pairs_methods.png", width = 800, height = 800)
pairs(merged_results[, 2:4], 
      main = "LogFC Correlation Across Methods",  # Title of the plot
      cex.main = 2,  # Increase the size of the main title
      cex.labels = 5,  # Increase the size of axis labels
      cex.axis = 1.2,  # Increase the size of axis ticks
      oma = c(4, 4, 6, 4))  # Add extra space for the title
dev.off()




#COMPARISONS FOR VENN DIAGRAM

#UP genes
edgeR_genes_DOWN = rownames(edgeR_results)[edgeR_results$logFC < 0 & edgeR_results$PValue < 0.05]
DESeq2_genes_DOWN = rownames(RES_table)[RES_table$log2FoldChange < 0 & RES_table$pvalue < 0.05]
limma_genes_DOWN = rownames(limma_results)[limma_results$logFC < 0 & limma_results$P.Value < 0.05]
etest_genes_DOWN = rownames(exacttest_results)[exacttest_results$logFC < 0 & exacttest_results$PValue < 0.05]


#DOWN genes
edgeR_genes_UP = rownames(edgeR_results)[edgeR_results$logFC > 0 & edgeR_results$PValue < 0.05]
DESeq2_genes_UP = rownames(RES_table)[RES_table$log2FoldChange > 0 & RES_table$pvalue < 0.05]
limma_genes_UP = rownames(limma_results)[limma_results$logFC > 0 & limma_results$P.Value < 0.05]
etest_genes_UP = rownames(exacttest_results)[exacttest_results$logFC > 0 & exacttest_results$PValue < 0.05]



#Combined in a table (DOWN)

max_length = max(length(DESeq2_genes_DOWN), length(limma_genes_DOWN), length(etest_genes_DOWN))

# Adjust the lists to have the same length by filling with NA
DESeq2_genes_DOWN = c(DESeq2_genes_DOWN, rep(NA, max_length - length(DESeq2_genes_DOWN)))
limma_genes_DOWN = c(limma_genes_DOWN, rep(NA, max_length - length(limma_genes_DOWN)))
etest_genes_DOWN = c(etest_genes_DOWN, rep(NA, max_length - length(etest_genes_DOWN)))

# Make final dataframe
genes_down_df = data.frame(
  DESeq2_genes_DOWN = DESeq2_genes_DOWN,
  limma_genes_DOWN = limma_genes_DOWN,
  etest_genes_DOWN = etest_genes_DOWN
)


#Combined in a table (UP)
max_length = max(length(DESeq2_genes_UP), length(limma_genes_UP), length(etest_genes_UP))

#Adjust the lists to have the same length by filling with NA
DESeq2_genes_UP = c(DESeq2_genes_UP, rep(NA, max_length - length(DESeq2_genes_UP)))
limma_genes_UP = c(limma_genes_UP, rep(NA, max_length - length(limma_genes_UP)))
etest_genes_UP = c(etest_genes_UP, rep(NA, max_length - length(etest_genes_UP)))

# Make final dataframe
genes_UP_df = data.frame(
  DESeq2_genes_UP = DESeq2_genes_UP,
  limma_genes_UP = limma_genes_UP,
  etest_genes_UP = etest_genes_UP
)

#Save data
write_csv(genes_down_df,"/home/guvega/Documents/COHORTS PROJECT/TFM/down_together.csv" )
write_csv(genes_UP_df,"/home/guvega/Documents/COHORTS PROJECT/TFM/up_together.csv" )









###############################################
### FURTHER BIOLOGICAL ANALYSIS #################
################################################

#Check number of DEGs using a strict filtering (FDR)
table(limma_results$FDR < 0.05)



################## PATHWAYS ##############
#Extract pathways described in GSEA

#Test pathways available
human_gene_sets = msigdbr(species = "Homo sapiens")
print(human_gene_sets)
collections_summary = unique(human_gene_sets$gs_cat)
table(human_gene_sets$gs_subcat, human_gene_sets$gs_cat)

#VERY IMPORTANT: Verify it you need the human or mouse genes!! "Homo sapiens" or "Mus musculus"
#Hallmarks 
h_df_h = msigdbr(species = "Homo sapiens", category = "H")
hallmarks = h_df_h %>% split(x = .$gene_symbol, f = .$gs_name)

#C2: REACTOME
h_df_c2_REACTOME = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
C2_REACTOME = h_df_c2_REACTOME %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: Biocarta,KEGG, REACTOME, PID etc.
h_df_c2_BIOCARTA = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
C2_BIOCARTA = h_df_c2_BIOCARTA %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: Biocarta,KEGG, REACTOME, PID etc.
h_df_c2_KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
C2_KEGG = h_df_c2_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: Biocarta,KEGG, REACTOME, PID etc.
h_df_c2_WIKI = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
C2_WIKI = h_df_c2_WIKI %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: Biocarta,KEGG, REACTOME, PID etc.
h_df_c2_PID = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
C2_PID = h_df_c2_PID %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: Biocarta,KEGG, REACTOME, PID etc.
h_df_c2_CP = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
C2_CP = h_df_c2_CP %>% split(x = .$gene_symbol, f = .$gs_name)
#C2: FULL
h_df_c2 = msigdbr(species = "Homo sapiens", category = "C2")
C2 = h_df_c2 %>% split(x = .$gene_symbol, f = .$gs_name)


#C5: Gene Ontology BP
h_df_c5_BP = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
C5_BP = h_df_c5_BP %>% split(x = .$gene_symbol, f = .$gs_name)
#C5: Gene Ontology CC
h_df_c5_CC = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
C5_CC = h_df_c5_CC %>% split(x = .$gene_symbol, f = .$gs_name)
#C5: Gene Ontology MF
h_df_c5_MF = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
C5_MF = h_df_c5_MF %>% split(x = .$gene_symbol, f = .$gs_name)
#C5: Gene Ontology FULL
h_df_c5 = msigdbr(species = "Homo sapiens", category = "C5")
C5 = h_df_c5 %>% split(x = .$gene_symbol, f = .$gs_name)


#7: Immunologic signature
h_df_c7 = msigdbr(species = "Homo sapiens", category = "C7")
C7 = h_df_c7 %>% split(x = .$gene_symbol, f = .$gs_name)

# create a list for each set of pathways
list_fgesea_pathways = list(hallmarks, C2_REACTOME, C2_BIOCARTA, C2_KEGG, C2_WIKI, C2_PID, C2_CP,C2, C5_BP, C5_CC,C5_MF,C5, C7)
names(list_fgesea_pathways) = c("hallmarks", "C2_REACTOME", "C2_BIOCARTA", "C2_KEGG", 
                                "C2_WIKI", "C2_PID", "C2_CP", "C2", "C5_BP", "C5_CC","C5_MF", "C5","C7")



#######################   GSEA   ##################
library(fgsea)
library(msigdbr)

limma_pathways_list = list("hallmarks" = list() , "C2_REACTOME" = list() , "C2_BIOCARTA" = list() , 
                              "C2_KEGG" = list() , "C2_WIKI" = list() , "C2_PID" = list() , 
                              "C2_CP" = list() , "C2" = list(),"C5_BP" = list() , "C5_CC" = list() ,"C5_MF" = list(),"C5" = list(), "C7" = list())


start_time = Sys.time()

#Perform Fast Gene Set Enrichment Analysis (GSEA) for Multiple Pathways
for (n in 1:length(list_fgesea_pathways)){
  gene_sets = list_fgesea_pathways[[n]]
  nPermSimple = 10000
  
  # iterate through the list of pathways
  # log2FC vector
  original_gene_list = limma_results$logFC
  # name the vector
  names(original_gene_list) = rownames(limma_results)
  # omit any NA values 
  gene_list = na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  # fastGSEA
  fgseaRes = fgsea(pathways=gene_sets, stats=gene_list,  minSize  = 10, maxSize  = 500, nPermSimple = nPermSimple, eps = 0)
  # tidy results
  fgseaResTidy = fgseaRes %>% as_tibble() %>% dplyr::arrange(desc(NES))
  # save results in a list and change the name
  limma_pathways_list[[n]]= fgseaResTidy
}

end_time = Sys.time()
end_time - start_time


#Save data
write_xlsx(limma_pathways_list, "/home/guvega/Documents/COHORTS PROJECT/TFM/Pathways_limma.xlsx")





### PLOT - HALLMARKS

# Filter top 5 pos and neg 
top_pos = limma_pathways_list$hallmarks %>% 
  arrange(desc(NES)) %>% 
  head(5)

top_neg = limma_pathways_list$hallmarks %>% 
  arrange(NES) %>% 
  head(5)

#Combine
top_combined = rbind(top_pos, top_neg) %>% 
  arrange(NES)

#Plot
ggplot(top_combined, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Negative NES", "Positive NES"),
                    name = "") + 
  coord_flip() +
  labs(title = "Top 5 Hallmarks Pathways",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme(
    panel.background = element_blank(), # No grey background
    panel.grid = element_blank(),       # Remove squares
    axis.text.y = element_text(face = "bold", size = 22, color = "black"), 
    axis.text.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(face = "bold", size = 28, color = "black"), 
    axis.title.x = element_text(size = 28, face = "bold", color = "black"), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"), 
    legend.position = "top",
    legend.title = element_blank(), # Remove legend title
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.key = element_blank()
  )








### PLOT - C5

# Filter top 5 pos and neg 
top_pos = limma_pathways_list$C5_MF %>% 
  arrange(desc(NES)) %>% 
  head(5)

top_neg = limma_pathways_list$C5_MF %>% 
  arrange(NES) %>% 
  head(5)

# Combine
top_combined = rbind(top_pos, top_neg) %>% 
  arrange(NES)

#Plot
ggplot(top_combined, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("Negative NES", "Positive NES"),
                    name = "") + 
  coord_flip() +
  labs(title = "Top 5 C5 - Molecular function Pathways",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(),       
    axis.text.y = element_text(face = "bold", size = 22, color = "black"), 
    axis.text.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(face = "bold", size = 28, color = "black"), 
    axis.title.x = element_text(size = 28, face = "bold", color = "black"), 
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5, color = "black"), 
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 25), 
    legend.background = element_blank(),
    legend.key = element_blank()
  )
