#
# The function creates and saves in folder 
# "../data/Myconnectome/"
# the connectivity matrix 
# datamatrix_mean
# generated from the time series stored in 
# "../data/Myconnectome/time_series/"
# and for the region of interest specified by 
# ROI_number
#
# The file is saved as 
# "connectivityMaps_myconnectome_ROI_number.RData"
#
# NOTE: the following libraries are needed:
#	library(mgcv)
#
createDataMatrix_by_region = function (nnodes, labels, ROI_number)
{
	time_series_folder = "../data/time_means/"

	ROI_index = ifelse(labels == ROI_number, 1, 0)
	ROI_index = as.logical(ROI_index)

	names_session = sprintf('%0.3d', c(11:41,43:51,53:79,81:89,91:104))

	datamatrix_mean = NULL

	for(k in 1:length(names_session))
	{
		name = paste(time_series_folder,"ses",names_session[k],
				   "_time_means.csv",sep="")
		time_mean = as.matrix(read.csv(name, header=FALSE))
		time_mean = t(time_mean)

		ROI_avg = time_mean[,ROI_number]
  
		cor_reg = NULL
		for(j in 1:83)
		{
			if (j == ROI_number)
				cor_reg = c(cor_reg,1)
			else
				cor_reg = c(cor_reg, cor(ROI_avg,time_mean[,j]))
		}
  
		z = 0.5 * log((1 + cor_reg) / (1 - cor_reg))
		z[which(is.infinite(z))] = 3.8002
  
		datamatrix_mean = rbind(datamatrix_mean,z)
		print(k)
	}

	datamatrix_mean[which(is.na(datamatrix_mean))] = 0

	save(datamatrix_mean,
		 file=paste0(myconn_folder,
					 "connectivityMaps_myconnectome_mean_",
					 ROI_number,".RData"))
}