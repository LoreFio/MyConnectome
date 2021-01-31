library(fdaPDE)
library(mgcv)

source("createDataMatrix.R")
source("computeCorrelation.R")

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
		
max_valid_label = 83
variance_threshold = 0.85
minimum_corr_threshold = 0.4

for (ROI in c(9:14))
{
  print(paste0("ROI: ",ROI))
	##############################################################
	#
	#					VERY IMPORTANT WARNING!!!
	#
	# This will probably take a HUGE amount of time...comment as
	# soon as it is no more necessary!!
	#createDataMatrix(mesh, labels, ROI)
	##############################################################
	
	computeCorrelation(ROI,labels, max_valid_label, 
					   variance_threshold, minimum_corr_threshold,F)
}