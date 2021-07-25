                                       **# Single-cell-RNA-sequencing**
                                                 **ABSTRACT**
The retina is the utmost important part of the eye which enables visual perception. There are many 
studies undertaken to understand the retinal cells and their functions, but still, the proportions of cell 
types in the damaged retina (diseased) are not accurately studied. The complication of diabetes 
mellitus disease known as diabetic retinopathy is one of the growing diseases in the medical field. The 
vast knowledge of the cell types of the retina is essential to know the impact of diabetic retinopathy. 
An understanding of how diabetes affects cell type proportions would help to understand its influence 
in the retina. This project used the current single-cell RNA sequencing pipeline, to provide deeper 
insights about the retinal cell-types and difference between the healthy and diseased retina. Using the 
“10X chromium” technique, two groups of mice were culled, control and streptozotocin-treated 
diabetic mice (control and diabetic retinas). The control dataset contained 14852 genes and 13039 
cells, while the diabetic mouse matrix remained with 11561 cells and 15106 genes. The quality control 
and clustering process were analysed using Seurat package from R. It has been shown that the diabetic 
dataset contains 38 clusters and control has 42 clusters. The clusters were annotated according to cell 
types and detected that both the datasets contain the more or less same number of cell types except 
a difference in one cell-type.
To gain further insights into the effect of diabetes using many existing bulk datasets, a deconvolution 
process is used to estimate the cell type proportions from these bulk data. Information from singlecell data was used to guide deconvolution. Different deconvolution packages were tried to installed 
and use them, but unfortunately due to the recent R updated version, only two packages were able 
to install and used. They are “Bisque” and “MuSiC” which are recently developed packages in R. These 
packages are designed to detect the accurate cell-type proportions but unluckily for this dataset, and
they were unable to retrieve the accurate cell-type estimations. The bulk data used in this step was 
from a collaborator in Queen’s University, Dr Eleni Beli. This bulk dataset was from type 1 diabetic 
mice and control, which contains six-time points. One of them is ZT5 time point which is same as the 
single-cell data that has been used in this project. So, that particular time point was used in this 
deconvolution process. The results from this step were unsuccessful because there was a considerable 
difference between the cell types from bulk and single-cell data. The single-cell data proportions and 
the deconvolution proportions were compared, except for Rod photoreceptor the other cell-types 
were not adequately estimated. There are significant differences in both proportions. This project 
explains the usefulness of cell-type proportions and the importance of filtering out the low-quality 
cells from the dataset, and it also explains the problem in deconvolution algorithm and barriers in 
single-cell RNA sequencing. ScRNA-seq is one of the promising technology in this field. However, the 
deconvolution process has not satisfied as much as expected, with accurate cell proportions in this 
project.

                                                       **Objectives: ** 
 Optimise the scRNA-seq data approach using Seurat by investigating the effects of QC and  background correction for healthy and a diabetic dataset of the mouse retina. 
 Perform differential gene expression.  
 Perform clustering of the different retinal cell types from the scRNA-Seq data by assessing a  range of publicly available algorithms. 
 Identify the cell-types for both control and diabetic retinal dataset. 
 Try to install different deconvolution packages that suits the project. 
 Use them to identify cell-type proportions from bulk data guided by the scRNAseq data  
 Compare the results to identify the best proportions and the most affected cell-types.

                                                       **Control retina**
![image](https://user-images.githubusercontent.com/54199923/126889830-23c216c9-3261-4131-96b9-e8dc05cf40ae.png)

Figure 1 - Cell cluster UMAP to regulate control retina in the mouse. (Figure 1a - left) shows 42 clusters in a variety of clusters. 
(Figure 1b - right) showing the same cell clusters in 8 different cell types.

                                                       **Diabetic retina**
![image](https://user-images.githubusercontent.com/54199923/126889795-683910f7-6488-482e-a530-efac4a267dfe.png)

Figure 2 - Cell cluster UMAP to regulate diabetic retina in the mouse. (Figure 2a - left) shows 38 clusters in a variety of clusters. (Figure 2b - right) showing the same cell clusters in 8 different cell types. 

                                                **Control (left) Vs Diabetic (right)**
![image](https://user-images.githubusercontent.com/54199923/126889875-b22f1265-7612-45dd-a97a-453d613bdf6e.png)

The above figure: 3 represents the difference between control and diseased dataset. Both the datasets are more or less same, except instead of endothelial cell, horizontal cells are available in control datasets.

**Diabetic Results:**
Cell type 	Number of cells 	Proportions 
Rod 	3846 	0.332670 
Bipolar 	2376 	0.205518 
Muller 	2387 	0.206470 
Cone 	911 	0.078799 
Amacrine 	1888 	0.163307 
Ganglion 	60 	0.005189 
Endothelial 	54 	0.004670 
Microglia/Astrocytes 	39 	0.003373 

**Control Results**
Cell type 	Number of cells 	Proportions 
Rod 	3846 	0.332670 
Bipolar 	2376 	0.205518 
Muller 	2387 	0.206470 
Cone 	911 	0.078799 
Amacrine 	1888 	0.163307 
Ganglion 	60 	0.005189 
Endothelial 	54 	0.004670 
Microglia/Astrocytes 	39 	0.003373 

**Bar chart for both the results**
![image](https://user-images.githubusercontent.com/54199923/126889924-c772f4a5-d4c8-4603-9e2a-01b4371f9e56.png)


**Overall Result table**
 	The cell-type proportion from ScRNA-
Seq analysis 	Bisque deconvolution estimations for control and diabetic 	from 	reference-based 	method 	with markers 	MuSiC estimations for control dataset 
 	Control 	C5.1 	C5.2 	C5.3 	C5.4 	C5.5 	C5.
1 	C5.
2 	C5.
3 	C5.
4 	C5.
5 
Rod 	0.334534 	0 	0.0520
68 	0.5687
26 	0.4612
9 	0.5117
71 	0.7
2 	0.7
4 	0.7
4 	0.7
1 	0.7
5 
Muller 	0.183909 	0.0678
64 	0 	0.1771
11 	0 	0 	0 	0 	0 	0 	0 
Bipolar 	0.209678 	0.1061
36 	0.1014
86 	0.1352
48 	0.5387
1 	0 	0 	0 	0 	0 	0 
Cone 	0.067796 	0.1080
84 	0.2756
38 	0 	0 	0.0804
42 	0.0
2 	0.0
1 	0.0
2 	0.0
2 	0.0
2 
Amacrine 	0.164199 	0.3514
13 	0.4357
24 	0.0809
2 	0 	0 	0.0
2 	0 	0.0
2 	0.0
5 	0.0
1 
Ganglion 	0.028376 	0.0308
02 	0.0846
23 	0.0090
68 	0 	0.0549
48 	0.2
3 	0.2
3 	0.2
1 	0.2
1 	0.2
1 
Horizontal 	0.007362 	0.3357
01 	0.0204
93 	0.0289
27 	0 	0 	0 	0 	0 	0 	0 
Microglia/Ast rocytes 	0.004141 	0 	0.0299
68 	0 	0 	0.3528
38 	0.0
2 	0.0
1 	0.0
2 	0.0
2 	0.0
2 
 	 	Bisque reference-based estimations with markers for diabetic 	MuSiC estimations for diabetic 
 	Diabetic 	D5.1 	D5.2 	D5.3 	D5.4 	D5.5 	D5.
1 	D5.
2 	D5.
3 	D5.
4 	D5.
5 
Rod 	0.332670 	0.39 	0.22 	0.38 	0.59 	0 	0.7
1 	0.7
1 	0.7
2 	0.7 	0.7
2 
Bipolar 	0.205518 	0 	0.2 	0.42 	0.38 	0 	0 	0 	0 	0 	0 
Muller 	0.206470 	0.21 	0.07 	0.02 	0 	0.71 	0.0
4 	0.0
5 	0.0
5 	0.0
4 	0.0
3 
Cone 	0.078799 	0.06 	0.07 	0.09 	0 	0.07 	0.0
1 	0.0
1 	0.0
1 	0.0
1 	0.0
2 
Amacrine 	0.163307 	0.26 	0.32 	0 	0 	0.18 	0.0
4 	0.0
4 	0.0
3 	0.0
5 	0.0
5 
Ganglion 	0.005189 	0 	0.07 	0.09 	0 	0.04 	0.2 	0.1
9 	0.1
8 	0.2 	0.1
9 
Endothelial 	0.004670 	0 	0.04 	0 	0 	0 	0 	0 	0 	0 	0 
Microglia/Ast rocytes 	0.003373 	0.08 	0 	0 	0.02 	0 	0 	0 	0 	0 	0 



