#
# The function creates and saves in folder 
# "../result/CorrMatrices/"
# the correlation matrices
# CorrMatrix_k.csv
# generated, for every scan k, from the time means stored in 
# "../data/time_means/" and setting to zero all the values lower,
# in absolute value, than a certain threshold
#
createCorrelationMatrices_thres = function (thres = 0.05)
{
	time_series_folder = "../data/time_means/"

	# NOTE: sessions 12, 14, 15, 17, 18, 20, 21, 23, 24, 27, 33, 36
	# do not have corresponding behavioral data, thus they are not
	# considered
	names_session = sprintf('%0.3d', 
							c(11,13,16,19,22,25,26,28:32,34,35,
							  37:41,43:51,53:73,75:79,81:89,91:104))

	# list of correlation matrices
	corr_matrix_list = NULL
	
	# first loop: it computes a list of correlation matrices
	for(k in 1:length(names_session))
	{
		corr_matrix = NULL
		
		name = paste(time_series_folder,"ses",names_session[k],
				   "_time_means.csv",sep="")
		time_mean = as.matrix(read.csv(name, header=FALSE))
		time_mean = t(time_mean)
		time_mean = time_mean[,c(1:73,75:83)]
  
		corr_matrix = cor(time_mean)
		corr_matrix_list[[k]] = corr_matrix
		
		print(paste0("First loop - Iteration #",k))
	}
	
	print("#################################################")
	
	# second loop: for every corr_matrix in corr_matrix_list, it sets 
	# to zero every element which is <= maxmin (in absolute value)
	# and saves the matrices to files
	for(j in 1:length(names_session))
	{
		corr_matrix = corr_matrix_list[[j]]
		
		for(h in 1:(dim(corr_matrix)[1]-1))
		{
			for(k in (h+1):dim(corr_matrix)[2])
			{
				if(abs(corr_matrix[h,k]) <= thres)	
				{
					corr_matrix[h,k] = 0.0
					corr_matrix[k,h] = 0.0
				}
			}
		}
		
		corr_matrix = as.data.frame(abs(corr_matrix))
		colnames(corr_matrix) = rownames(corr_matrix)
		
		write.csv(corr_matrix, 
				  file=paste0("../result/CorrMatrices/CorrMatrix_",
				  			  j,".csv"))
		
		print(paste0("Second loop - Iteration #",j))
	}
}