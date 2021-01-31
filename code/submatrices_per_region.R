# The function generates from datamatrix (passed as parameter)
# and returns a list of 95 submatrices, each
# one containing the nodes corresponding to a specific brain region,
# identified by its label

submatrices_per_region = function (datamatrix, reg_labels)
{
	matrix_list = list()
	# BE careful maybe you will need to access with the commented line.
	#for (j in 1:max(reg_labels$X0))
	for (j in 1:max(reg_labels$V1))
	  {
		indices = c()
		for (k in 1:dim(reg_labels)[1])
		  #if (reg_labels$X0[k] == j)
		  if (reg_labels$V1[k] == j)
		    indices = append(indices,k)
		matrix_list[[j]] = datamatrix[,indices]
	}
	return(matrix_list)
}