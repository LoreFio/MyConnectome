# plot_on_mesh_by_region = function (region_field, 
#									 basename="region_field")
#
# plots the field passed as parameter on the brain mesh
#
# Input:
#		region_field:  	vector of length 83 (one attribute per region)
#						to be plotted on the brain mesh				
#		basename:		basename (without extension) of the file .vtu
#
# Output:
# 		none

plot_on_mesh_by_region = function (region_field,
								   basename="region_field")
{
	source("plot_on_mesh.R")
  
	# upload region labels
  	labels = as.matrix(read.table(
			"../data/Myconnectome/labels_hammer_myconnectome_new2.csv",
			quote="\"", comment.char=""))
  
  	# convert region_field to a vector of length 17287, where each 
  	# point belonging to a region has the corresponding value read
  	# from region_field
	point_field = rep(0, dim(labels)[1])
	
	for (j in 1:length(point_field))
	{
		if (labels[j] > 0 && labels[j] < 84)
		{
			point_field[j] = region_field[labels[j]]
		}
	}
	
  	# plot on mesh the point field
	plot_on_mesh(point_field,basename)
}