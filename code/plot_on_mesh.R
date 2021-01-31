# plot_on_mesh = function (field, onVTU=FALSE, basename="field")
#
# plots the field passed as parameter on the brain mesh
#
# Input:
#		field:  	vector of length 17287 to be plotted on the brain
#					mesh
#		basename:	basename (without extension) of the file .vtu
#
# Output:
# 		none

plot_on_mesh = function (field, basename="field")
{
	source('fdaPDE.write.vtk.R')
	library(fdaPDE)
	# upload data to create the mesh
	nodes_brain = as.matrix(read.csv(
				"../data/Myconnectome/node_myconnectome_17287.csv",
				header=FALSE))
	elems = as.matrix(read.csv(
				"../data/Myconnectome/elem_myconnectome_all_17287.csv",
				header=FALSE))

	# create mesh
	mesh = create.MESH.3D(nodes = nodes_brain, tetrahedrons = elems)
	
	# create finite element basis					
	FEMbasis <- create.FEM.basis(mesh)
	
	# define field on the finite element basis
	FEM_field = FEM(field,FEMbasis)
	
	# export as basename.vtu
	filename = paste0(basename,".vtu")
	write.vtu(FEM_field, file = paste0("../result/vtu/",filename))
}