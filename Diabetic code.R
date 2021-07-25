getwd()
setwd("C:/Users/40256911/Documents/OneDrive_1_6-10-2020/Matrices")
library(dplyr)
library(Seurat)
library(cowplot)
library(Matrix)
library(ggplot2)
library(sctransform)



##############
library1.data<-Read10X("filtered_feature_bc_matrix1")
library3.data<-Read10X("filtered_feature_bc_matrix3")
library5.data<-Read10X("filtered_feature_bc_matrix5")
library7.data<-Read10X("filtered_feature_bc_matrix7")
library9.data<-Read10X("filtered_feature_bc_matrix9")
library11.data<-Read10X("filtered_feature_bc_matrix11")
library13.data<-Read10X("filtered_feature_bc_matrix13")
library14.data<-Read10X("filtered_feature_bc_matrix14")
library15.data<-Read10X("filtered_feature_bc_matrix15")

################
library1 <- CreateSeuratObject(counts = library1.data, min.cells = 10, min.features = 800, project = "retina1")
library1$library <- "1"
library1$DR<- "Diabetic"
library1$timepoint<-"Early"
library3 <- CreateSeuratObject(counts = library3.data, min.cells = 10, min.features = 800, project = "retina3")
library3$library <- "3"
library3$DR<- "Diabetic"
library3$timepoint<-"Mid"
library5 <- CreateSeuratObject(counts = library5.data, min.cells = 10, min.features = 800, project = "retina5")
library5$library <- "5"
library5$DR<- "Diabetic"
library5$timepoint<-"Mid"
library7 <- CreateSeuratObject(counts = library7.data, min.cells = 10, min.features = 800, project = "retina7")
library7$library <- "7"
library7$DR<- "Diabetic"
library7$timepoint<-"Mid"
library9 <- CreateSeuratObject(counts = library9.data, min.cells = 10, min.features = 800, project = "retina9")
library9$library <- "9" 
library9$DR<- "Diabetic"
library9$timepoint<-"Early"
library11 <- CreateSeuratObject(counts = library11.data, min.cells = 10, min.features = 800, project = "retina11")
library11$library <- "11" 
library11$DR<- "Diabetic"
library11$timepoint<-"Early"
library13 <- CreateSeuratObject(counts = library13.data, min.cells = 10, min.features = 800, project = "retina13")
library13$library <- "13"
library13$DR<- "Diabetic"
library13$timepoint<-"Late"
library14 <- CreateSeuratObject(counts = library14.data, min.cells = 10, min.features = 800, project = "retina14")
library14$library <- "14" 
library14$DR<- "Diabetic"
library14$timepoint<-"Late"
library15 <- CreateSeuratObject(counts = library15.data, min.cells = 10, min.features = 800, project = "retina15")
library15$library <- "15" 
library15$DR<- "Diabetic"
library15$timepoint<-"Late"

mouse_dbt <- merge(library1, y = list(library3, library5,  library7,
                                      library9, library11, library13, library14, 
                                      library15), add.cell.ids = c("retina1","retina3", "retina5", "retina7",
                                                                   "retina9", "retina11","retina13", "retina14",
                                                                   "retina15"), merge.data = TRUE, project = "SeuratProject")
mouse_dbt <- PercentageFeatureSet(mouse_dbt, pattern = "^MT-", col.name = "percent.mt")
dim(mouse_dbt)
#############
#######################

length(mouse_dbt@meta.data$percent.mt)
#Visualize QC metrics as violin plots 
VlnPlot(mouse_dbt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

mean("nFeature_RNA") 

median(mouse_dbt@meta.data$nFeature_RNA) 
plot1 <- FeatureScatter(mouse_dbt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_dbt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cowplot::plot_grid(plot1,plot2, labels = c("A","B")) 
ggsave("qualitycontrolcorrplots.tiff", plot = last_plot(), dpi = 800, width = 50, height = 15, units = "cm", device = "tiff") 
########################
summary(mouse_dbt@meta.data$percent.mt)
############
#Subset retina object according to quality control metrics to remove unwanted cells from the dataset based on the number of features and percentage of mitochondrial counts 
mouse_dbt <- subset(mouse_dbt, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5) 
dim(mouse_dbt)

####Normalize##########################
mouse_dbt <- NormalizeData(mouse_dbt, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_dbt <- NormalizeData(mouse_dbt)
#Access normalised values  
mouse_dbt[["RNA"]]@data  

mouse_dbt[["RNA"]]@data[c("mt-Nd1", "mt-Nd4"),] 
VlnPlot(mouse_dbt, features = "mt-Nd1") 
VlnPlot(mouse_dbt, features = "mt-Nd4") 
######Identification of highly variable features####################
mouse_dbt <- FindVariableFeatures(mouse_dbt, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse_dbt), 30)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse_dbt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
class(plot1)
VlnPlot(mouse_dbt, features = "mt-Nd1")
######Scaling the data###################
all.genes <- row.names(mouse_dbt)
all.genes
mouse_dbt <- ScaleData(mouse_dbt, vars.to.regress = "percent.mt")
#Access scaled data 
mouse_dbt[["RNA"]]@scale.data[1:5,1:5] 
##############Linear dimentionality reduction#######################
mouse_dbt <- RunPCA(mouse_dbt, features = VariableFeatures(object = mouse_dbt))

print(mouse_dbt[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouse_dbt, dims = 1:3, reduction = "pca")
DimPlot(mouse_dbt, reduction = "pca")
DimHeatmap(mouse_dbt, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mouse_dbt, dims = 1:15, cells = 500, balanced = TRUE)

mouse_dbt <- RunPCA(mouse_dbt, npcs = 30, verbose = FALSE) #####samples
mouse_dbt <- RunUMAP(mouse_dbt, reduction = "pca", dims = 1:30)
p1 <- DimPlot(mouse_dbt, reduction = "umap")
p2 <- DimPlot(mouse_dbt, reduction = "umap", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2 ############

#######dimensionality test############
mouse_dbt <- JackStraw(mouse_dbt, num.replicate = 100, dims = 50)
mouse_dbt <- ScoreJackStraw(mouse_dbt, dims = 1:20)
JackStrawPlot(mouse_dbt, dims = 1:15)

ElbowPlot(mouse_dbt)

##########Clustering the cells########
mouse_dbt <- FindNeighbors(mouse_dbt, dims = 1:30)
mouse_dbt <- FindClusters(mouse_dbt, resolution = 1.2)

# Look at cluster IDs of the first 5 cells
head(Idents(mouse_dbt), 5)

######UMAP/TSNE#######
mouse_dbt <- RunUMAP(mouse_dbt, dims = 1:30)
mouse_dbt <- RunTSNE(mouse_dbt, check_duplicates = FALSE, dims = 1:30)
DimPlot(mouse_dbt, reduction = "tsne", label = TRUE)
DimPlot(mouse_dbt, reduction = "umap", label = TRUE)
saveRDS(mouse_dbt, file = "###New result for all diabetic filtered.rds")

#######finding differentially expressed features#######
# find all markers of cluster 1
cluster1.markers <- FindMarkers(mouse_dbt, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(mouse_dbt, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.markers <- FindAllMarkers(mouse_dbt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
retina.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
saveRDS(retina.markers, "diabetic_markers.rds")





mouse_dbt <- RenameIdents(mouse_dbt,  '0'= "Rod", '1' = "Bipolar", '2' = "Rod", '3' = "Rod", '4' = "Muller", 
                                  '5' = "Rod", '6' = "Muller", '7' = "Cone", '8' = "Rod", '9' = "Bipolar",'10' = "Amacrine", '11' = "Amacrine", 
                                  '12' = "Amacrine", '13' = "Bipolar", '14' = "Cone", '15' = "Bipolar", '16' = "Bipolar" ,'17' = "Muller", '18' = "Amacrine", 
                                  '19'= "Bipolar", '20' = "Amacrine" ,'21' = "Bipolar" , '22' = "Bipolar", '23' = "Bipolar",
                                  '24' = "Amacrine" , '25' = "Amacrine" , '26' = "Bipolar", '27' = "Amacrine" , '28' = "Ganglion", '29' = "Bipolar" , 
                                  '30' = "Amacrine", '31' = "Amacrine", '32' = "Amacrine" , '33' = "Amacrine", '34' = "Vascular", '35' = "Microglia/Astrocytes" ,
                                  '36'= "Bipolar", '37' = "Horizontal")
TSNEPlot(object = mouse_dbt, pt.size = 0.5)
DimPlot(mouse_dbt, reduction = "umap", label = TRUE, pt.size = 0.5) 

DimPlot(mouse_dbt, reduction = "tsne", label = TRUE, pt.size = 0.5) 
saveRDS(mouse_dbt, file = "New####diabetic_result_celltype.rds")

########################################
install.packages("BisqueRNA")
library(BisqueRNA)
library(Biobase)
library(Seurat)

##See Bisque vignette at https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html

setwd("C:/Users/40256911/Documents/Bisque/Bisque")

#Load single cell data
retina_subset_late_normal <- readRDS("New####diabetic_result_celltype.rds")

#Read in bulk data
#NB: Make sure bulk and single cell data aligned to same genome build with same annotation
#The count matrix had several rows named Mar01 or Mar02, probably becasue of automatic date naming in Escel. 
#As sort term fix remove duplicate rows, but better to get original matrix
bulk_RNA_all <- read.delim(file = "PN0159_Merged_raw_counts_with_samples_Mar01_02_deleted.txt", header = TRUE, sep = "\t", dec = ".")
#Set column 1 to row names
bulk_RNA_all <- data.frame(bulk_RNA_all[,-1], row.names = bulk_RNA_all[,1])
#Select samples to analyse
bulk_RNA <- bulk_RNA_all[,c("D1.1","D1.2","D1.3","D1.4","D1.5","D5.1","D5.2","D5.3","D5.4","D5.5")]
#Convert bulk data to expressionset
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk_RNA))
#Check the bulk ExpressionSet:
sampleNames(bulk.eset)

#Bisque  function to automatically Seurat  object to an ExpressionSet
#NB by looking at the cell barcodes, eg AAACCCATCCATTCAT_16 you can see that the delimiter is underscore, followed by the library number (ie individual mouse) 
sc.eset <- BisqueRNA::SeuratToExpressionSet(retina_subset_late_normal, delimiter="_", position=2, version="v3")

#Check the single cell ExpressionSet:
sc.eset$SubjectName
sc.eset$cellType

#If we had the same sample with both bulk and single cell data (named the same) we could perform a reference-based decomposition
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers = NULL, use.overlap=FALSE)
res
#Call the marker-based decomposition method
res_weighted <- BisqueRNA::MarkerBasedDecomposition(bulk.eset,retina.markers, min_gene = 4, weighted=T)
res_weighted
#Look at results
#These estimates are relative within each cell type, so you cannot immediately compare abundance estimates between cell types.
marker.based.estimates <- res_weighted$bulk.props
knitr::kable(marker.based.estimates, digits = 2)
#output results (***remember to move sample names over 1 column to right***)
write.table(marker.based.estimates, file = "fresh for bulk diabetix relative.proportions_weighted.csv", sep = ",", quote = FALSE)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers= retina.markers, use.overlap = FALSE)
reference.based.estimates <- res$bulk.props
knitr::kable(reference.based.estimates, digits = 2)
write.table(reference.based.estimates, file = "reference.based.proportions.csv", sep = ",")

#############################################MUSIC##################

Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'cellType',
                      samples = 'SubjectName', select.ct = c('Rod', 'Bipolar', 'Cone', 'Amacrine','Muller',
                                                             'Vasculer', 'Horizontal', 'Ganglion'))
rownames(sc.eset)
Estimates <- Est.prop$Est.prop.weighted
knitr::kable(Estimates, digits = 2)
Est.prop
write.table(Estimates, file = "##newmusic cell type for diabetes.csv", sep = ",", quote = FALSE)

Jitter_Est(Estimates, method.name = NULL, title = NULL)
jitter.fig = Jitter_Est(list(data.matrix(Est.prop$Est.prop.weighted),
                             data.matrix(Est.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jitter.fig
