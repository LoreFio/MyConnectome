rm(list = ls())
library(fdaPDE)
library(mgcv)

source("createDataMatrix_by_region.R")
source("pca_by_region.R")

myconn_folder = "../data/Myconnectome/"

labels = read.table(
			paste0(myconn_folder,"labels_hammer_myconnectome_new2.csv"),
			quote="\"", comment.char="")
			
nodes <- read.csv(
				paste0(myconn_folder,"node_myconnectome_17287.csv"),
				header=FALSE)
elems <- read.csv(
		paste0(myconn_folder,"elem_myconnectome_all_tissue_17287.csv"),
		header=FALSE)
elems = elems[,1:4]
mesh <- create.MESH.3D(nodes=nodes, tetrahedrons = elems)

list_PCA = NULL

for(idx in c(1:73,75:83))
{
	print(paste0("ROI = ",idx))
	createDataMatrix_by_region(mesh$nnodes,labels,idx)
	list_PCA[[idx]] = pca_by_region(ROI=idx, var_explained=0.5, 
				 				    n_max_pc_plot=3, want_plots=F)	
}

save(list_PCA, file="../result/resultPCA_by_region.RData")
