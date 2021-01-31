# Compute the correlation matrix for the datamatrix produced by
# function createDataMatrix
#
# Parameters:
#	ROI:	index of the region of interest
#	reg_labels: labels read from labels_hammer_myconnectome_new2.csv
#	max_valid_label: all regions with label > max_valid_label will 
#					 be discarded
#	variance_threshold:
# 	minimum_corr_threshold: all regions with correlation below this
#							threshold (in absolute value) will be
#							discarded
#	want_plots: TRUE if we want to plot results
computeCorrelation = function (ROI, 
						  	   reg_labels, 
	 						   max_valid_label = 83,
							   variance_threshold = 0.85, 
						  	   minimum_corr_threshold = 0.4,
						  	   want_plots = FALSE)
{
	#################################
	####      SOURCE and LOAD     ###
	#################################
	source("submatrices_per_region.R")
	source("reduction_matrix.R")
	source("corr_region_filter.R")
	
	load(paste0("../data/Myconnectome/connectivityMaps_myconnectome_",
				ROI,".RData"))

	load("../data/mcshapiro.test.RData")

	#################################
	####        MEAN              ###
	#################################

	matrix_list = submatrices_per_region(datamatrix, reg_labels)

	n_row = dim(matrix_list[[1]])[1]
	n_reg = length(matrix_list)
	mean_reg_mat = as.data.frame(matrix(1:n_row*n_reg, 
								 nrow = n_row, ncol = n_reg ))
	meaningless_reg = c(0,ROI)

	# We perform the mean for each region that have at least 
	# 2 points and whose label is at most max_valid_label
	for(j in 1:n_reg)
	{
		if( j > max_valid_label || is.null(dim(matrix_list[[j]])[2]) 
			|| dim(matrix_list[[j]])[2] <= 1 )
		{
		    meaningless_reg = c(meaningless_reg,j)
		} else {
		    mean_reg_mat[,j] = rowMeans(matrix_list[[j]]) 
		}
	}

	# We filter out meaningless regions
	mean_reg_mat = mean_reg_mat[,-meaningless_reg]
	label_name = 1:n_reg
	label_name = label_name[-meaningless_reg]

	# Analysis of evolutions of means in time.
	regions_time = c(4,9,23,51,53,73)	# these are random num ers
	if (want_plots)
	{
		x11()
		matplot(mean_reg_mat[,label_name[regions_time]],
    		    xlab = "time", ylab = "mean of region", 
        		type = "l", lty = 1, lwd = 2, 
        		col = 	rainbow(length(regions_time)))
		legend('topleft', legend = label_name[regions_time], 
				lwd = 1, col = rainbow(length(regions_time)))
	}
	submatrix = mean_reg_mat[,label_name[regions_time]]

	pc_test = princomp(submatrix, scores = T)
	summary(pc_test)
	pc_test$loadings

	#################################
	####  CORRELATION MATRIX      ###
	#################################

	submatrix_cor = cor(submatrix)

	# Study of the shape of the correlation matrix and isolation 
	# of important groups of regions
	# that are correlated among them (intra) => once we have 
	# selected these groups of regions we may could try see if these 
	# groups are effectively correlated with the ROI (inter)
	matrix_cor = cor(mean_reg_mat)

	if (want_plots)
	{
		x11()
		image(matrix_cor, col = rainbow(20), axes = F)
		axis(1, at=seq(0,1, length.out =  length(label_name)), 
			 labels = label_name)
		axis(2, at=seq(0,1, length.out =  length(label_name)), 
			 labels = label_name)

		x11()
		filled.contour(matrix_cor)
	}
	graphics.off()

	correlated_regions_auto=corr_region_filter(matrix_cor, 
								  			   label_name,
											   minimum_corr_threshold)

	sub_cor_auto = matrix_cor[correlated_regions_auto,
							  correlated_regions_auto]
							  
	save(ROI, correlated_regions_auto, sub_cor_auto,
		 minimum_corr_threshold, 
		 file = paste0("../data/MRPCA/MRPCA_",ROI,".RData"))
}