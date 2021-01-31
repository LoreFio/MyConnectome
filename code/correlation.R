#################################
####        LIBRARIES         ###
#################################

rm(list = ls())
library(fdaPDE)

#################################
####        PARAMETERS        ###
#################################

max_valid_label = 83
min_num_points = 5 # in each region
ROI = 60
var_explained = 0.5
load_thres = 0.5
n_max_pc_plot = 4
want_plots = T
call_pc_global = T
call_pc_mean = T
time_window_segmentation=T

#################################
####      SOURCE and LOAD     ###
#################################

source("submatrices_per_region.R")
source("pca_global.R")
source("pca_correlation.R")

load(paste0("../data/Myconnectome/connectivityMaps_myconnectome_",
            ROI,".RData"))	# It's enought to load any datamatrix 
							# in order to get the meaningless regions

### DANGER
# In order to read properly the label file use read.table instead 
# of read.csv
reg_labels = read.table(
			"../data/Myconnectome/labels_hammer_myconnectome_new2.csv",
            quote = "\"", comment.char = "")

matrix_list = submatrices_per_region(datamatrix, reg_labels)
meaningless_reg = 0
for(j in 1:max_valid_label)
{
  if( j > max_valid_label || is.null(dim(matrix_list[[j]])[2]) 
      || dim(matrix_list[[j]])[2] < min_num_points )
  {
    meaningless_reg = c(meaningless_reg,j)
  } 
}

rm(datamatrix, matrix_list)

# Meaningfull regions
label_name = 1:max_valid_label
label_name = label_name[-meaningless_reg] # 78

if(call_pc_global)
{
  pca_global(ROI = ROI, 
             var_explained = var_explained,
             n_max_pc_plot = n_max_pc_plot,
             want_plots = want_plots)
}

##############################################################
####        MEAN Before TRANSFORMATION  and CORRELATION   ####
##############################################################

setwd("../data/time_means")

names_session=sprintf('%0.3d', c(11:41,43:51,53:79,81:89,91:104))
# empty 52 and 90 ma sono state rinominate male [non sono 42 2 80]

matrix_list=list() # correlation matrix (78x78) for each session

for(k in 1:length(names_session))
{
	name=paste("ses",names_session[k],"_time_means.csv",sep="")
	time_series=t(as.matrix(read.csv(name, header=FALSE))) # 518x83
	if (time_window_segmentation)
	{
		# 14 time windows
    	time_series_shorter=matrix(0L, nrow = 14, 
    							   ncol = max_valid_label) 
    	for (j in 1:max_valid_label) # per ogni colonna (regione)
    	{ 
      		s=1
      		for (i in 1:14) 
      		{ 
		        time_series_shorter[i,j] = mean(time_series[s:s+36,j])
        		s = s + 37
      		}
    	}
    	time_series = time_series_shorter
  	}   
	matrix_list[[k]] = cor(time_series[,label_name])
}

setwd("../../code/")

save(matrix_list, label_name, min_num_points, 
     file = paste0("../data/mean_cov_experiments.RData"))

#######################
####      PLOT     ####
#######################

if(want_plots)
{
  # Evolution of some correlations:
  # regions label_name(13) and label_name(56)
  corr_evol=NULL
  for (k in 1:length(names_session)) {
    corr_evol<-c(corr_evol,matrix_list[[k]][13,56]) 
  }
  x11()
  plot(corr_evol,type = "l",col="blue")

  # regions label_name(27) and label_name(68)
  corr_evol=NULL
  for (k in 1:length(names_session)) {
    corr_evol<-c(corr_evol,matrix_list[[k]][27,68])
  }
  lines(corr_evol,col="red")
}

#######################################
####      CORRELATION with ROI     ####
#######################################

# ROI setted at the beginning
ROI_index = which(label_name %in% ROI)

corr_matrix = matrix(0, nrow = length(names_session), 
					 ncol = length(label_name))
for (k in 1:length(names_session))
{
	for (j in 1:length(label_name)) 
	{
    	corr_matrix[k,j] = cor(matrix_list[[k]][,ROI_index],
    					  	   matrix_list[[k]][,j])
  	}  
}

z = 0.5 * log((1+corr_matrix)/(1-corr_matrix))
z[which(z > 0.99)] = 3.8002 # caso limite: rho=1 !
z[which(z < -0.99)] = -3.8002 # caso limite: rho=-1 !

corr_matrix = z

save(corr_matrix, ROI, ROI_index, label_name,
     file = paste0("../data/mean_",ROI,"_correlation.RData"))

if (call_pc_mean)
{
  pca_correlation(ROI = ROI,
                  var_explained = var_explained,
                  load_thres = load_thres,
                  n_max_pc_plot = n_max_pc_plot,
                  want_plots = want_plots)
}

#################################################
####      CORRELATION analysis over TIME     ####
################################################

# mean and variance over the sessions (time)
mean_corr = diag(length(label_name))
var_corr = matrix(0L, nrow = length(label_name), 
				  ncol = length(label_name)) 

corr_label = NULL
mean_vec = NULL
var_vec = NULL
corr_serie = rep(0,length(names_session)) 	# vect di dim 
											# length(names_session)

for (i in 1:(length(label_name)-1)) 
{
	for (j in (i+1):length(label_name)) 
	{
    	for (k in 1:length(names_session)) 
    	{
		    corr_serie[k] = matrix_list[[k]][i,j] # è simmetrica!
    	}
	  corr_label=c(corr_label,paste0("(",i,",",j,")"))
	  
	  mean_corr[i,j] = mean(corr_serie)
	  mean_corr[j,i] = mean_corr[i,j]
	  mean_vec=c(mean_vec,mean_corr[i,j])
	  
	  var_corr[i,j] = var(corr_serie)
	  var_corr[j,i] = var_corr[i,j]
	  var_vec=c(var_vec,var_corr[i,j])
	}
}

x <- corr_label
F <- mean_vec
L <- mean_vec-2*sqrt(var_vec) #75%
U <- mean_vec+2*sqrt(var_vec) #75%

if(want_plots)
  {
  library(ggplot2)
  x11()
  d=data.frame(correlations=x, mean=F, lower=L, upper=U)
  ggplot() + 
  geom_errorbar(data=d, mapping=aes(x=correlations, ymin=upper, ymax=lower), width=0.2, size=1, color="blue") + 
  geom_point(data=d, mapping=aes(x=correlations, y=mean), size=4, shape=21, fill="white")
}

