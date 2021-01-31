rm(list = ls())
load("../result/resultPCA_by_region.RData")

N_reg = 83
link_abs = matrix(0, nrow = N_reg, ncol = N_reg)
link_plus = matrix(0, nrow = N_reg, ncol = N_reg)
link_minus = matrix(0, nrow = N_reg, ncol = N_reg)
vec_reg = c(1:73,75:N_reg)

for(k in vec_reg)
{
	for(i in vec_reg)
	{
		for(j in vec_reg)
		{
			link_abs[i,j] = link_abs[i,j] + 
				ifelse(i %in% list_PCA[[k]][,1] &&
				       j %in% list_PCA[[k]][,1], 1,0)

			link_plus[i,j] = link_plus[i,j] + 
				ifelse(sum(c(i,1) %in% list_PCA[[k]])==2 &&
				       sum(c(j,1) %in% list_PCA[[k]])==2, 1,0)
			
			link_minus[i,j] = link_minus[i,j] + 
				ifelse(sum(c(i,-1) %in% list_PCA[[k]])==2 &&
				       sum(c(j,-1) %in% list_PCA[[k]])==2, 1,0)
		}
	}
}

link_abs = as.data.frame(link_abs)
colnames(link_abs) = paste0(1:83) 

write.csv(link_abs, file="../result/link_abs.csv")
write.table(as.data.frame(link_plus), file="../result/link_plus.csv", 
			sep=",")
write.table(as.data.frame(link_minus), file="../result/link_minus.csv", 
			sep=",")
