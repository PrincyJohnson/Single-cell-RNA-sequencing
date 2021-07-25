getwd()
setwd("C:/Users/40256911/Documents/OneDrive_1_6-10-2020/Matrices")
library(dplyr)
library(Seurat)
library(cowplot)
library(Matrix)
library(ggplot2)
library(sctransform)

options(future.globals.maxSize = 12000 * 1024^2) #sctransform needs this
install.packages("rlang")



library(rlang)
sessionInfo()

##############

library2.data<-Read10X("filtered_feature_bc_matrix2")

library4.data<-Read10X("filtered_feature_bc_matrix4")

library6.data<-Read10X("filtered_feature_bc_matrix6")

library8.data<-Read10X("filtered_feature_bc_matrix8")

library10.data<-Read10X("filtered_feature_bc_matrix10")

library12.data<-Read10X("filtered_feature_bc_matrix12")



library16.data<-Read10X("filtered_feature_bc_matrix16")
library17.data<-Read10X("filtered_feature_bc_matrix17")
library18.data<-Read10X("filtered_feature_bc_matrix18")
###


library2 <- CreateSeuratObject(counts = library2.data, min.cells = 10, min.features = 800, project = "retina2")
library2$library <- "2" 
library2$DR<- "Control"
library2$timepoint<-"Early"

library4 <- CreateSeuratObject(counts = library4.data, min.cells = 10, min.features = 800, project = "retina4")
library4$library <- "4"
library4$DR<- "Control"
library4$timepoint<-"Mid"

library6 <- CreateSeuratObject(counts = library6.data, min.cells = 10, min.features = 800, project = "retina6")
library6$library <- "6"
library6$DR<- "Control"
library6$timepoint<-"Mid"

library8 <- CreateSeuratObject(counts = library8.data, min.cells = 10, min.features = 800, project = "retina8")
library8$library <- "8"
library8$DR<- "Control"
library8$timepoint<-"Mid"

library10 <- CreateSeuratObject(counts = library10.data, min.cells = 10, min.features = 800, project = "retina10")
library10$library <- "10"
library10$DR<- "Control"
library10$timepoint<-"Early"

library12 <- CreateSeuratObject(counts = library12.data, min.cells = 10, min.features = 800, project = "retina12")
library12$library <- "12"
library12$DR<- "Control"
library12$timepoint<-"Early"

library16 <- CreateSeuratObject(counts = library16.data, min.cells = 10, min.features = 800, project = "retina16")
library16$library <- "16"
library16$DR<- "Control"
library16$timepoint<-"Late"
library17 <- CreateSeuratObject(counts = library17.data, min.cells = 10, min.features = 800, project = "retina17")
library17$library <- "17" 
library17$DR<- "Control"
library17$timepoint<-"Late"
library18 <- CreateSeuratObject(counts = library18.data, min.cells = 10, min.features = 800, project = "retina18")
library18$library <- "18" 
library18$DR<- "Control"
library18$timepoint<-"Late" 

##############

mouse_ctrl <- merge(library2, y = list(library4, library6,
                                      library8, library10, library12, library16,
                                      library17, library18), add.cell.ids = c("retina2","retina4", "retina6", "retina8",
                                                                   "retina10", "retina12","retina16", "retina17",
                                                                   "retina18"), merge.data = TRUE, project = "SeuratProject_ctrl")
################
mouse_ctrl <- PercentageFeatureSet(mouse_ctrl, pattern = "^MT-", col.name = "percent.mt")
dim(mouse_ctrl)
#############
#######################

length(mouse_ctrl@meta.data$percent.mt)
#Visualize QC metrics as violin plots 
VlnPlot(mouse_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

mean("nFeature_RNA") 

median(mouse_ctrl@meta.data$nFeature_RNA) 
plot1 <- FeatureScatter(mouse_ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cowplot::plot_grid(plot1,plot2, labels = c("A","B")) 
ggsave("qualitycontrolcorrplots_ctrl.tiff", plot = last_plot(), dpi = 800, width = 50, height = 15, units = "cm", device = "tiff") 
########################
summary(mouse_ctrl@meta.data$percent.mt)
############
#Subset retina object according to quality control metrics to remove unwanted cells from the dataset based on the number of features and percentage of mitochondrial counts 
mouse_ctrl <- subset(mouse_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5) 
dim(mouse_ctrl)

####Normalize##########################
mouse_ctrl <- NormalizeData(mouse_ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_ctrl <- NormalizeData(mouse_ctrl)
#Access normalised values  
mouse_ctrl[["RNA"]]@data  

mouse_ctrl[["RNA"]]@data[c("mt-Nd1", "mt-Nd4"),] 
VlnPlot(mouse_ctrl, features = "mt-Nd1") 
VlnPlot(mouse_ctrl, features = "mt-Nd4") 
######Identification of highly variable features####################
mouse_ctrl <- FindVariableFeatures(mouse_ctrl, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse_ctrl), 30)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouse_ctrl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
class(plot1)
VlnPlot(mouse_ctrl, features = "mt-Nd1")
######Scaling the data###################
all.genes <- row.names(mouse_ctrl)
all.genes
mouse_ctrl <- ScaleData(mouse_ctrl, vars.to.regress = "percent.mt")
#Access scaled data 
mouse_ctrl[["RNA"]]@scale.data[1:5,1:5] 
##############Linear dimentionality reduction#######################
mouse_ctrl <- RunPCA(mouse_ctrl, features = VariableFeatures(object = mouse_ctrl))
mouse_ctrl <- RunPCA(mouse_ctrl, verbose = F)

print(mouse_ctrl[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouse_ctrl, dims = 1:3, reduction = "pca")
DimPlot(mouse_ctrl, reduction = "pca")
DimHeatmap(mouse_ctrl, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mouse_ctrl, dims = 1:15, cells = 500, balanced = TRUE)

mouse_ctrl <- RunPCA(mouse_ctrl, npcs = 30, verbose = FALSE) #####samples
mouse_ctrl <- RunUMAP(mouse_ctrl, reduction = "pca", dims = 1:30)
p1 <- DimPlot(mouse_ctrl, reduction = "umap")
p2 <- DimPlot(mouse_ctrl, reduction = "umap", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2 ############

#######dimensionality test############
mouse_ctrl <- JackStraw(mouse_ctrl, num.replicate = 100, dims = 50)
mouse_ctrl <- ScoreJackStraw(mouse_ctrl, dims = 1:20)
JackStrawPlot(mouse_ctrl, dims = 1:15)

ElbowPlot(mouse_ctrl)

##########Clustering the cells########
mouse_ctrl <- FindNeighbors(mouse_ctrl, dims = 1:30)
mouse_ctrl <- FindClusters(mouse_ctrl, resolution = 1.2)

# Look at cluster IDs of the first 5 cells
head(Idents(mouse_ctrl), 5)

######UMAP/TSNE#######
mouse_ctrl <- RunUMAP(mouse_ctrl, dims = 1:30)
mouse_ctrl <- RunTSNE(mouse_ctrl, check_duplicates = FALSE, dims = 1:30)
DimPlot(mouse_ctrl, reduction = "tsne", label = TRUE)
DimPlot(mouse_ctrl, reduction = "umap", label = TRUE)
saveRDS(mouse_ctrl, file = "###New result for all diabetic filtered.rds")

#######finding differentially expressed features#######
# find all markers of cluster 1
cluster1.markers_ctrl <- FindMarkers(mouse_ctrl, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers_ctrl, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers_ctrl <- FindMarkers(mouse_ctrl, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers_ctrl, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
retina.markers_ctrl <- FindAllMarkers(mouse_ctrl, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
retina.markers_ctrl %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
saveRDS(retina.markers_ctrl, "control_markers.rds")
################################################

mouse_ctrl <- RenameIdents(mouse_ctrl,  '0'= "Rod", '1' = "Bipolar", '2' = "Rod", '3' = "Rod", '4' = "Muller", 
                                  '5' = "Rod", '6' = "Muller", '7' = "Cone", '8' = "Rod", '9' = "Bipolar",'10' = "Amacrine", '11' = "Amacrine", 
                                  '12' = "Amacrine", '13' = "Bipolar", '14' = "Cone", '15' = "Bipolar", '16' = "Bipolar" ,'17' = "Muller", '18' = "Amacrine", 
                                  '19'= "Bipolar", '20' = "Amacrine" ,'21' = "Bipolar" , '22' = "Bipolar", '23' = "Bipolar",
                                  '24' = "Amacrine" , '25' = "Amacrine" , '26' = "Bipolar", '27' = "Amacrine" , '28' = "Ganglion", '29' = "Bipolar" , 
                                  '30' = "Amacrine", '31' = "Amacrine", '32' = "Amacrine" , '33' = "Amacrine", '34' = "Vascular", '35' = "Microglia/Astrocytes" ,
                                  '36'= "Bipolar" , '37' = "Horizontal" , '38' = "Amacrine", '39' = "Bipolar", '40' = "Ganglion", '41' = "Microglia")
TSNEPlot(object = mouse_ctrl, pt.size = 0.5)
DimPlot(mouse_ctrl, reduction = "umap", label = TRUE, pt.size = 0.5) 

DimPlot(mouse_ctrl, reduction = "tsne", label = TRUE, pt.size = 0.5) 
saveRDS(mouse_ctrl, file = "New### control_result_celltype.rds")
levels(mouse_ctrl)

########################################
install.packages("BisqueRNA")
library(BisqueRNA)
library(Biobase)
library(Seurat)

##See Bisque vignette at https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html

setwd("C:/Users/40256911/Documents/Bisque/Bisque")

#Load single cell data
retina_subset_late_normal <- readRDS("New### control_result_celltype.rds")

#Read in bulk data
#NB: Make sure bulk and single cell data aligned to same genome build with same annotation
#The count matrix had several rows named Mar01 or Mar02, probably becasue of automatic date naming in Escel. 
#As sort term fix remove duplicate rows, but better to get original matrix
bulk_RNA_all <- read.delim(file = "PN0159_Merged_raw_counts_with_samples_Mar01_02_deleted.txt", header = TRUE, sep = "\t", dec = ".")
#Set column 1 to row names
bulk_RNA_all <- data.frame(bulk_RNA_all[,-1], row.names = bulk_RNA_all[,1])
#Select samples to analyse
bulk_RNA <- bulk_RNA_all[,c("C5.1","C5.2","C5.3","C5.4","C5.5", "C9.1","C9.2","C9.3","C9.4","C9.5")]
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

reference.based.estimates <- res$bulk.props
reference.based.estimates
knitr::kable(reference.based.estimates, digits = 2)
write.table(reference.based.estimates, file = "controlreference based proportions without markers.csv", sep = ",", quote = FALSE)

#Call the marker-based decomposition method
res_weighted <- BisqueRNA::MarkerBasedDecomposition(bulk.eset,retina.markers_ctrl, min_gene = 5, max_gene = 200 , weighted=F, w_col = "avg_logFC", unique_markers = TRUE,
                                                    verbose = TRUE)
res_weighted
#Look at results
#These estimates are relative within each cell type, so you cannot immediately compare abundance estimates between cell types.
marker.based.estimates <- res_weighted$bulk.props
knitr::kable(marker.based.estimates, digits = 2)
#output results (***remember to move sample names over 1 column to right***)
write.table(marker.based.estimates, file = "marker baed control relative.proportions_weighted.csv", sep = ",", quote = FALSE)

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers= retina.markers_ctrl, use.overlap = FALSE)
reference.based.estimates <- res$bulk.props
knitr::kable(reference.based.estimates, digits = 2)
write.table(reference.based.estimates, file = "controlreference.based.proportions_with markers.csv", sep = ",")

#############################################MUSIC##################

Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'cellType',
                      samples = 'SubjectName', select.ct = c('Rod', 'Bipolar', 'Cone', 'Amacrine','Muller',
                                                             'Vasculer', 'Horizontal', 'Ganglion', 'Microglia/Astrocytes',
                                                             'Microglia'))
rownames(sc.eset)
Estimates <- Est.prop$Est.prop.weighted
knitr::kable(Estimates, digits = 2)
Est.prop
write.table(Estimates, file = "##newmusic control cell type of mice.csv", sep = ",", quote = FALSE)

Jitter_Est(Estimates, method.name = NULL, title = NULL)
jitter.fig = Jitter_Est(list(data.matrix(Est.prop$Est.prop.weighted),
                             data.matrix(Est.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jitter.fig
