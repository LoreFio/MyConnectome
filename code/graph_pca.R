## Download and install the package
#install.packages("igraph")
rm(list = ls())
## Load package
library(igraph)
load("../data/mcshapiro.test.RData")
source("createCorrelationMatrices_thres.R")

n_region = 82
n_session = 78
thres = 0
# the next line can be commented if you have the matrices with 0 threshold already created
#createCorrelationMatrices_thres(0)
session_vector = 1:n_session
reg_names = paste0(1:n_region)
e_centrality = matrix(0, nrow = length(session_vector)
, ncol = n_region)
c_centrality = matrix(0, nrow = length(session_vector)
                      , ncol = n_region)
b_centrality = matrix(0, nrow = length(session_vector)
                      , ncol = n_region)

#createCorrelationMatrices_thres(thres)
#uncomment if necessary

for (session in session_vector) {
  corr_matrix = read.csv(file=paste0("../result/CorrMatrices/CorrMatrix_",
                                     session,".csv"))
  corr_matrix = corr_matrix[,-1]
  for (i in 1:n_region) {
    corr_matrix[i,i] = 0
  }
  # we remove correlation with itself since it would createa graph with cycles
  g = graph.adjacency(as.matrix(corr_matrix), mode = "undirected", weighted = TRUE)
  eig_resul = eigen_centrality(g)
  e_centrality[session,] = eig_resul$vector
  c_centrality[session,] =  closeness(g)
  b_centrality[session,] =  estimate_betweenness(g, cutoff = 5)
}

### test of gaussianity # here is the problem!!!
# Since p > n we can't use this fct.
#mcshapiro.test(e_centrality)
#mcshapiro.test(c_centrality)
#mcshapiro.test(b_centrality)
# we can either test on a subset of the features or we can look at qqplots
x11()
par(mfrow=c(1,3))
qqnorm(e_centrality, main = "Eigen Centrality")
qqline(e_centrality)
qqnorm(c_centrality, main = "Closeness Centrality")
qqline(c_centrality)
qqnorm(b_centrality, main = "Betweenness Centrality")
qqline(b_centrality)

### centralities means, cov
e_centrality.mean<-colMeans(e_centrality)
e_centrality.cov = cov(e_centrality)

c_centrality.mean<-colMeans(c_centrality)
c_centrality.cov = cov(c_centrality)

b_centrality.mean<-colMeans(b_centrality)
b_centrality.cov = cov(b_centrality)

# comparison among different centrality measures
x11()
par(mfrow=c(1,3))
plot(e_centrality.mean, ylab = "Eigenvector Centrality")
text(e_centrality.mean, row.names(reg_names), cex = 2*e_centrality.mean/max(e_centrality.mean),
     pos = 4, col = "green")
plot(c_centrality.mean, ylab = "Closeness Centrality")
text(c_centrality.mean, row.names(reg_names), cex = 2*c_centrality.mean/max(c_centrality.mean),
     pos = 4, col = "red")
plot(b_centrality.mean, ylab = "Betweenness Centrality")
text(b_centrality.mean, row.names(reg_names), cex = 2*b_centrality.mean/max(b_centrality.mean),
     pos = 4, col = "blue")

x11()
par(mfrow=c(1,3))
image(e_centrality.cov,col=rainbow(48))
image(c_centrality.cov,col=rainbow(48))
image(b_centrality.cov,col=rainbow(48))

# the most interesting one seems to be the betweennes one
# we try now to plot the CI for all the regions using the bonferroni's correction
k <- n_region  # 2
alpha = 0.1
n = 78
cfr.t <- qt(1-alpha/(2*k),n-1)

IC.BF.C_Centr <- cbind( c_centrality.mean-cfr.t*sqrt(diag(c_centrality.cov)/n),
                        c_centrality.mean,
                        c_centrality.mean+cfr.t*sqrt(diag(c_centrality.cov)/n) )

dimnames(IC.BF.C_Centr)[[2]] <- c('inf','center','sup')
#IC.BF.B_Centr

IC.BF.B_Centr <- cbind( b_centrality.mean-cfr.t*sqrt(diag(b_centrality.cov)/n),
                        b_centrality.mean,
                        b_centrality.mean+cfr.t*sqrt(diag(b_centrality.cov)/n) )

dimnames(IC.BF.B_Centr)[[2]] <- c('inf','center','sup')
#IC.BF.B_Centr

x11()
plot(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2],
     ylim = c(min(IC.BF.B_Centr), max(IC.BF.B_Centr)), col = "red", 
     cex=0.6, pch = 16, xlab = "Regions", ylab = "Betweennes")
axis(1, at = 1:n_region, labels = reg_names[order(IC.BF.B_Centr[,2])],lwd.ticks = 0.6 )
segments(x0=1:n_region,y0=IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),1],
         y1=IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),3],col="orange")
points(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2], col = "red", cex=0.8, pch = 16)
points(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),1], col = "blue", cex=0.6, pch = 4)
points(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),3], col = "blue", cex=0.6, pch = 4)
text(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2], reg_names[order(IC.BF.B_Centr[,2])],
     cex = 1.5*IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2]/max(IC.BF.B_Centr[,2]), pos = 3, col = "red")

x11()
plot(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2],
     ylim = c(min(IC.BF.B_Centr), max(IC.BF.B_Centr)), col = "red",
     type = "l",xlab = "Regions", ylab = "Betweennes")
axis(1, at = 1:n_region, labels = reg_names[order(IC.BF.B_Centr[,2])],lwd.ticks = 0.6 )
lines(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),1], col = "blue")
lines(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),3], col = "blue")
#text(IC.BF.B_Centr[order(IC.BF.B_Centr[,2]),2],
#     reg_names[order(IC.BF.B_Centr[,2])], cex = 1, pos = 4, col = "black")

# So far we should take a too large alpha to conclude => let's try to compare the centrality of the
# 3 most central regions with the mean of the centraliies of all regions and of all regions but 
# those 3

kbis = 5
cfr_t_bis <- qt(1-alpha/(2*kbis),n-1)

# We create a matrix that chooses the 3 most central regions and the 2 means we are interested in
matrix_A = matrix(0, nrow = kbis, ncol = n_region)
matrix_A[1,n_region-2] = 1
matrix_A[2,n_region-1] = 1
matrix_A[3,n_region] = 1
matrix_A[4,1:(n_region-3)] = rep(1/(n_region-3),n_region-3)
matrix_A[5,] = rep(1/n_region,n_region)
#matrix_A # uncomment this to check that A has been well constructed

c_comparison = c_centrality[,order(IC.BF.C_Centr[,2])] %*% t(matrix_A)
c_comparison_cov = matrix_A %*%
  c_centrality.cov[order(IC.BF.C_Centr[,2]),order(IC.BF.C_Centr[,2])]  %*% t(matrix_A)

mcshapiro.test(c_comparison)

IC_c_comparison = cbind(colMeans(c_comparison)-cfr_t_bis*sqrt(diag(c_comparison_cov)/n),
                        colMeans(c_comparison),
                        colMeans(c_comparison)+cfr_t_bis*sqrt(diag(c_comparison_cov)/n))

b_comparison = b_centrality[,order(IC.BF.B_Centr[,2])] %*% t(matrix_A)
b_comparison_cov = matrix_A %*%
  b_centrality.cov[order(IC.BF.B_Centr[,2]),order(IC.BF.B_Centr[,2])]  %*% t(matrix_A)

mcshapiro.test(b_comparison) #doesn't work properly because the system is singular

IC_b_comparison = cbind(colMeans(b_comparison)-cfr_t_bis*sqrt(diag(b_comparison_cov)/n),
                        colMeans(b_comparison),
                        colMeans(b_comparison)+cfr_t_bis*sqrt(diag(b_comparison_cov)/n))

x11()
par(mfrow=c(1,2))
plot(IC_c_comparison[,2], ylim = c(min(IC_c_comparison), max(IC_c_comparison)), col = "red",
     cex=1, pch = 16, xlab = "Regions", ylab = "Closeness", axes = F)
axis(1, at = 1:kbis, labels = c(reg_names[order(IC.BF.C_Centr[,2])][(n_region-2):n_region],
                                "mean without",
                                "total mean"),lwd.ticks = 0.6 )
axis(2, las = 2)
segments(x0=1:5,y0=IC_c_comparison[,1],y1=IC_c_comparison[,3],col="orange")
points(IC_c_comparison[,2], col = "red", cex=1, pch = 16)
points(IC_c_comparison[,1], col = "blue", cex=0.6, pch = 4)
points(IC_c_comparison[,3], col = "blue", cex=0.6, pch = 4)

plot(IC_b_comparison[,2], ylim = c(min(IC_b_comparison), max(IC_b_comparison)), col = "red",
     cex=1, pch = 16, xlab = "Regions", ylab = "Betweennes", axes = F)
axis(1, at = 1:kbis, labels = c(reg_names[order(IC.BF.B_Centr[,2])][(n_region-2):n_region],
                                "mean without",
                                "total mean"),lwd.ticks = 0.6 )
axis(2, las = 2)
segments(x0=1:5,y0=IC_b_comparison[,1],y1=IC_b_comparison[,3],col="orange")
points(IC_b_comparison[,2], col = "red", cex=1, pch = 16)
points(IC_b_comparison[,1], col = "blue", cex=0.6, pch = 4)
points(IC_b_comparison[,3], col = "blue", cex=0.6, pch = 4)

graphics.off()
want_plots = T
n_max_pc_plot = 4
var_explained = 0.8

#### PCA ####
pc_e_centr = prcomp(e_centrality)
summary(pc_e_centr)

pc_sdev = pc_e_centr$sdev
loadings = pc_e_centr$rotation
n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > var_explained ))
n_pc_variance

if(want_plots)
{
  x11()
  par(mar = c(1,4,0,2), mfrow = c(n_max_pc_plot,1))
  for(i in 1:n_max_pc_plot) barplot(loadings[,i], ylim = c(-1, 1))
  
  x11()
  layout(matrix(c(2,3,1,3),2,byrow=T))
  barplot(pc_e_centr$sdev^2, las=2, main='Principal Components', ylim = c(0,2),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  abline(h=1, col='blue')
  barplot(sapply(as.data.frame(e_centrality), sd), las=2, main='Original Variables', ylim = c(0,1),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  plot(cumsum(pc_e_centr$sdev^2)/sum(pc_e_centr$sde^2), type='b', axes=F,
       xlab='number of components', ylab='contribution to the total variance',
       ylim=c(0,1), xaxt = 'n')
  abline(h=var_explained, col='blue')
  box()
}

pc_b_centr = prcomp(b_centrality)
summary(pc_b_centr)

pc_sdev = pc_b_centr$sdev
loadings = pc_b_centr$rotation
n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > var_explained ))
n_pc_variance

if(want_plots)
{
  x11()
  par(mar = c(1,4,0,2), mfrow = c(n_max_pc_plot,1))
  for(i in 1:n_max_pc_plot) barplot(loadings[,i], ylim = c(-1, 1))
  
  x11()
  layout(matrix(c(2,3,1,3),2,byrow=T))
  barplot(pc_b_centr$sdev^2, las=2, main='Principal Components', ylim = c(0,1.5*max(pc_b_centr$sdev^2)),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  abline(h=1, col='blue')
  barplot(sapply(as.data.frame(b_centrality), sd), las=2, main='Original Variables', ylim = c(0,1000),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  plot(cumsum(pc_b_centr$sdev^2)/sum(pc_b_centr$sde^2), type='b', axes=F,
       xlab='number of components', ylab='contribution to the total variance',
       ylim=c(0,1), xaxt = 'n')
  abline(h=var_explained, col='blue')
  box()
}
