# Compute the pca for the correlation matrix obtained fixing a ROI
#
# Parameters:
#	ROI:	index of the region of interest
# 	var_explained: percentage of variance that we want to explain 
#	want_plots: TRUE if we want to plot results
pca_correlation = function(ROI = 59,
						   var_explained = 0.5, 
						   load_thres = 0.5, n_max_pc_plot = 4, want_plots = F)
{
  #################################
  ####      SOURCE and LOAD     ###
  #################################
  
  load(paste0("../data/mean_",ROI,"_correlation.RData"))
  myconn_folder = "../data/Myconnectome/"
  # ATTENTION: in the loaded file the matrix must be saved as corr_matrix
  source("plot_on_mesh_by_region.R")
  
  #################################
  ####            PCA           ###
  #################################

  pc.mean_reg = princomp(corr_matrix, scores = T)
  summary(pc.mean_reg)

  pc_sdev = pc.mean_reg$sdev
  loadings = pc.mean_reg$loadings
  n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > var_explained ))
  n_pc_variance
  
  if(want_plots)
  {
    x11()
    par(mar = c(1,4,0,2), mfrow = c(n_max_pc_plot,1))
    for(i in 1:n_max_pc_plot) barplot(loadings[,i], ylim = c(-1, 1))

    x11()
    layout(matrix(c(2,3,1,3),2,byrow=T))
    barplot(pc.mean_reg$sdev^2, las=2, main='Principal Components', ylim = c(0,2),
            xlab='number of components', ylab='Variances', xaxt = 'n')
    abline(h=1, col='blue')
    barplot(sapply(as.data.frame(corr_matrix), sd), las=2, main='Original Variables', ylim = c(0,1),
            xlab='number of components', ylab='Variances', xaxt = 'n')
    plot(cumsum(pc.mean_reg$sdev^2)/sum(pc.mean_reg$sde^2), type='b', axes=F,
        xlab='number of components', ylab='contribution to the total variance',
        ylim=c(0,1), xaxt = 'n')
    abline(h=var_explained, col='blue')
    box()
  }
  
  vector_regions = NULL
  for(i in 1: n_pc_variance)
  {
    plot_on_mesh_by_region(loadings[,i], basename = paste0("loadings_",ROI,i))
  }
  
  for(i in 1: n_pc_variance)
  {
    vector_regions = c(vector_regions, label_name[which(abs(loadings[,i]) >= load_thres)])
  }
  vector_regions = unique(vector_regions)
  save(vector_regions, ROI,
       file = paste0("../data/mean_vec_regions_",ROI,".RData"))
}
