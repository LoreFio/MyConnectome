#
# The function creates and saves in folder 
# "../result/CorrMatrices/"
# the correlation matrices
# CorrMatrix_k.csv
# generated, for every scan k, from the time means stored in 
# "../data/time_means/" and setting to zero all the values lower,
# in absolute value, than a certain threshold
#
createCorrelationMatrices_taskBased = function (task_code = 1)
{
  
  if ( task_code == 1)
  {
    task_name = 'breathhold'
    names_session = sprintf('%0.3d', 
                            c(62:64,66:74,76:80))
  }
  
  if ( task_code == 2)
  {
    task_name = 'dotstop'
    names_session = sprintf('%0.3d', 
                            c(81,83:86,88,89,91:93))
  }
  
  if ( task_code == 3)
  {
    task_name = 'nback'
    names_session = sprintf('%0.3d', 
                            c(14,17,20,23,29,32,35,38,40,43,48,51,57,59,61))
  }
  
  time_series_folder = paste0("../data/Myconnectome/TaskBased/",task_name,"/time_means/")
  
  # list of correlation matrices
  corr_matrix_list = NULL
  
  # first loop: it computes a list of correlation matrices
  for(k in 1:length(names_session))
  {
    corr_matrix = NULL
    
    name = paste(time_series_folder,task_name,"_ses",names_session[k],
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
    
    corr_matrix = as.data.frame(abs(corr_matrix))
    colnames(corr_matrix) = rownames(corr_matrix)
    
    write.csv(corr_matrix, 
              file=paste0("../result/CorrMatrices/",task_name,"CorrMatrix_",
                          j,".csv"))
    
    print(paste0("Second loop - Iteration #",j))
  }
}