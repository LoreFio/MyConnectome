# Compute the pca for the correlation matrix obtained fixing a ROI
#
# Parameters:
#	ROI:	index of the region of interest
#   var_explained: percentage of variance that we want to explain 
#	want_plots: TRUE if we want to plot results
pca_global = function(ROI = 59,
					  var_explained = 0.5, n_max_pc_plot = 4, want_plots = F)
{
  #################################
  ####      SOURCE and LOAD     ###
  #################################
 
  # ATTENTION: in the loaded file the matrix must be saved as datamatrix
  load(paste0("../data/Myconnectome/connectivityMaps_myconnectome_",
  			  ROI,".RData"))
  
  source("plot_on_mesh_by_region.R")
  
  #################################
  ####            PCA           ###
  #################################
  
  pc_full = prcomp(datamatrix)
  summary(pc_full)
  
  pc_sdev = pc_full$sdev
  loadings = pc_full$rotation
  n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > 
  							var_explained ))
  n_pc_variance
  
  if(want_plots)
  {
    x11()
    par(mar = c(1,4,0,2), mfrow = c(n_max_pc_plot,1))
    for(i in 1:n_max_pc_plot) barplot(loadings[,i], ylim = c(-1, 1))
    
    x11()
    layout(matrix(c(2,3,1,3),2,byrow=T))
    barplot(pc_full$sdev^2, las=2, main='Principal Components', 
    		ylim = c(0,100),
            xlab='number of components', ylab='Variances', xaxt = 'n')
    abline(h=1, col='blue')
    barplot(sapply(as.data.frame(datamatrix), sd), las=2, 
    		main='Original Variables', ylim = c(0,1),
            xlab='number of components', ylab='Variances', xaxt = 'n')
    plot(cumsum(pc_full$sdev^2)/sum(pc_full$sde^2), type='b', 
    	 axes=F, xlab='number of components', 
    	 ylab='contribution to the total variance',
         ylim=c(0,1), xaxt = 'n')
    abline(h=var_explained, col='blue')
    box()
  }
  
  for(i in 1: n_pc_variance)
  {
    plot_on_mesh_by_region(loadings[,i], 
    			 basename = paste0("loadings_global_",ROI,i))
  }
}