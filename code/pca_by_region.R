# Compute the pca for the correlation matrix obtained fixing a ROI
#
# Parameters:
#	ROI:	index of the region of interest
#   var_explained: percentage of variance that we want to explain 
#	want_plots: TRUE if we want to plot results
pca_by_region = function(ROI = 59,
					  var_explained = 0.5, n_max_pc_plot = 4, 
					  want_plots = F)
{
	#################################
	####      SOURCE and LOAD     ###
	#################################
 
	# ATTENTION: in the loaded file the matrix must be saved 
	# as datamatrix_mean
	load(paste0(
			"../data/Myconnectome/connectivityMaps_myconnectome_mean_",
  			ROI,".RData"))
  
	source("plot_on_mesh_by_region.R")
  
	#################################
	####            PCA           ###
	#################################
  
	pc_mean = prcomp(datamatrix_mean)
 	summary(pc_mean)
  
	pc_sdev = pc_mean$sdev
  	loadings = pc_mean$rotation
  	n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > 
  							var_explained ))
  	
  	matrix_result = NULL
	for(i in 1: n_pc_variance) 
	{
		l_var = sqrt(var(loadings[,i]))
		m_var = mean(loadings[,i])
		for(j in 1:83)
		{
			if(loadings[j,i] > m_var + 2*l_var)
				matrix_result = rbind(matrix_result, c(j,1))
			if(loadings[j,i] < m_var - 2*l_var)
				matrix_result = rbind(matrix_result, c(j,-1))
		}
	}
	
	if(want_plots)
  	{
    	n_max_pc_plot = min(n_max_pc_plot,n_pc_variance)
    	x11()
	    par(mar = c(3,4,2,2), mfrow = c(n_max_pc_plot,1))
    	for(i in 1:n_max_pc_plot) 
    		barplot(loadings[,i], ylim = c(-1, 1), 
    				names.arg=c(1:(dim(datamatrix_mean)[2])),
    				cex.names=0.8, main=paste0("pca ROI ",ROI))
	}
  
	for(i in 1: n_pc_variance)
  	{
    	plot_on_mesh_by_region(loadings[,i], 
    				 basename = paste0("loadings_mean_",ROI,"_",i))
  	}
  	
	return(matrix_result)
}