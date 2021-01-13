###############################################################################################
##        
##                      DECONVOLUTION OF RNA SEQ WITH CIBERSORT
##
##
###############################################################################################

library(factoextra)
library(ComplexHeatmap)
library(dplyr)
library(circlize)


load("/Users/jameskelly/Downloads/ucec.raw_counts.RData") ## raw counts df and an associated clinical df
UCECRNA <- ucec_raw_counts #make a DF with raw counts only
UCECRNA <- UCECRNA[,-1]
write.table(UCECRNA,"/Users/jameskelly/Documents/CIBER_DIR/UCECRNA.txt", sep = "\t") ###make a tab delimitted text file gene names only.


##make a working directory containing:input; immune_sig_matrix; and source code from CIBERSORT.
setwd("/Users/jameskelly/Documents/CIBER_DIR/")
Mixture <- "UCECRNA.txt" # Make object with mixture filename
ImmuneSignature <- "LM22.txt" # and immune signature filename


##This function automates the code for executing
## CIBERSORT deconvolution for any Immune signature matrix
##  and any cell mixture matrix with the same type of data frame
RUN_CIBERSORT <- function(ImmuneSignature, Mixture){
  library("e1071") 
  library(preprocessCore) 
  source("CIBERSORT.R") 
  ImmuneCells <- CIBERSORT(ImmuneSignature, Mixture, perm = 100, QN = TRUE)
  ImmuneCells
}
Immune_proportion <- RUN_CIBERSORT(ImmuneSignature = ImmuneSignature, Mixture = Mixture)



Immune_proportion <- as.data.frame(Immune_proportion) ## Change from matrix to df
names(Immune_proportion)[names(Immune_proportion) == "P-value"] <- "PVALUES" ## R did'nt like when I tried to filter p values 
## until I changed its name


################################################################################################
##        
##            UNSUPERVISED CLUSTERING AND VISUALISATION WITH COMPLEX HEATMAP
##
##
###############################################################################################


##Subsetting clinical info df to only include molecular subtype,
## immune subtype, and histological subtype
CLINICAL_INFO <- ucec_cdr_immune_clin
CLINICAL_INFO <- CLINICAL_INFO[,-1] ## remove ensembl IDs
CLINICAL_INFO<- CLINICAL_INFO[,c(1,6,88,89)] ## interested in annotating molec subtype, immune subtype, histological subtype


##Making a data frame with immune proportions and its associated molecular subtype,
## immune subtype, histological subtype.
CLINICAL_INFO$sample == rownames(Immune_proportion) ## Checking that both dfs in the same order
joined <- cbind(CLINICAL_INFO,Immune_proportion) ## so they can be joined into a single data frame
names(joined)[names(joined) == "Immune Subtype"] <- "Immune.Subtype" ##R complains about these 
names(joined)[names(joined) == "TCGA Subtype"] <- "TCGA.Subtype"  ##  spaces later so replace with .

##ONLY NOW SHOULD YOU REMOVE SAMPLES WITH LOW P VALUES FROM IMMUNE PROPORTION AND JOINED DATA FRAMES.

fltr_pval <- function(joined){
  joined <- dplyr::filter(joined, PVALUES < .05)
  joined
}

joined <- fltr_pval(joined)
Immune_proportion <- fltr_pval(Immune_proportion) 

Immune_proportion <- t(as.matrix(Immune_proportion[,1:22])) #include immune cell estimations only
##and transpose the matrix
fviz_nbclust(Immune_proportion, kmeans, method = "silhouette") ## k means value should be 2


##Subset the clinical info you want to annotate from the 'joined' df
Annotation_names <- joined[, c("histological_type",
                               "Immune.Subtype",
                               "TCGA.Subtype")]
##Each type is given an associated colour
Annotation_Colours <- list("histological_type" = c("Endometrioid endometrial adenocarcinoma" = "red", "Serous endometrial adenocarcinoma" = "black", "Mixed serous and endometrioid" = "orange"),
                           "Immune.Subtype" = c("C1" ="red", "C2" = "yellow","C3"="green", "C4" = "blue", "C5" = "pink", "C6" = "purple" ),
                           "TCGA.Subtype" = c("UCEC." = "white","UCEC.CN_LOW" = "yellow", "UCEC.CN_HIGH" = "red", "UCEC.POLE" = "purple", "NA" =  "white", "UCEC.MSI" = "orange" ))
top_anno <- ComplexHeatmap::HeatmapAnnotation(df = Annotation_names, which = "col", col = Annotation_Colours)

## Immune cell proportion colour scale is determined using
##minimum as white, median as light blue, maximum as navy
maximum <- as.numeric(max(t(Immune_proportion)))
median <- as.numeric(median(t(Immune_proportion)))
minimum <- as.numeric(min(t(Immune_proportion)))

colour_Immune_Proportion <- circlize::colorRamp2(c(minimum, median, maximum), c("white", "light blue", "navy"))
Heatmap(Immune_proportion, name = "Cell Proportion", col = colour_Immune_Proportion, column_km = 2, column_gap = unit(5, "mm"), border = T, top_annotation = top_anno)

