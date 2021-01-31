#
# The function creates and saves in folder 
# "../data/Myconnectome/"
# the connectivity matrix 
# datamatrix
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
createDataMatrix = function (mesh, labels, ROI_number)
{
	myconn_folder = "../data/Myconnectome/"
	time_series_folder = paste0(myconn_folder,"time_series/")

	ROI_index = ifelse(labels == ROI_number, 1, 0)
	ROI_index = as.logical(ROI_index)

	names_session = sprintf('%0.3d', c(11:41,43:51,53:79,81:89,91:104))

	datamatrix = NULL

	for(k in 1:length(names_session))
	{
		name = paste(time_series_folder,"ses",names_session[k],
				   "_time_series.csv",sep="")
		time_series = as.matrix(read.csv(name, header=FALSE))
		time_series = t(time_series)
		time_series = time_series[1:518,]
		ROI_ts = time_series[,ROI_index]
		dim(ROI_ts)
  
		cross_sec_avg_ROI_row = as.vector(rowMeans(ROI_ts))
  
		cor_nodes<-NULL
		for(j in 1:mesh$nnodes)
		{
			cor_nodes = c(cor_nodes,
						  cor(cross_sec_avg_ROI_row,time_series[,j]))
		}
  
		#cor_nodes[which(is.na(cor_nodes))] = 0
  
		r = cor_nodes
		z = 0.5 * log((1+r)/(1-r))
  
		datamatrix = rbind(datamatrix,z)
		print(k)
	}

	datamatrix[which(is.na(datamatrix))] = 0

	save(datamatrix,
		 file=paste0(myconn_folder,
					 "connectivityMaps_myconnectome_",
					 ROI_number,".RData"))
}