graphics.off()
want_plots = T
n_max_pc_plot = 4
var_explained = 0.8

x11()
pairs(full_df)

cor(full_df)
reduced_df = full_df[,-c(1,2)]
cor(reduced_df)

pc.data = princomp(reduced_df, scores = T)
summary(pc.data)

pc_sdev = pc.data$sdev
loadings = pc.data$loadings
n_pc_variance = min(which(cumsum(pc_sdev^2)/sum(pc_sdev^2) > var_explained ))
n_pc_variance

if(want_plots)
{
  x11()
  par(mar = c(1,4,0,2), mfrow = c(n_max_pc_plot,1))
  for(i in 1:n_max_pc_plot) barplot(loadings[,i], ylim = c(-1, 1))
  
  x11()
  layout(matrix(c(2,3,1,3),2,byrow=T))
  barplot(pc.data$sdev^2, las=2, main='Principal Components', ylim = c(0,2),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  abline(h=1, col='blue')
  barplot(sapply(as.data.frame(reduced_df), sd), las=2, main='Original Variables', ylim = c(0,1),
          xlab='number of components', ylab='Variances', xaxt = 'n')
  plot(cumsum(pc.data$sdev^2)/sum(pc.data$sde^2), type='b', axes=F,
       xlab='number of components', ylab='contribution to the total variance',
       ylim=c(0,1), xaxt = 'n')
  abline(h=var_explained, col='blue')
  box()
}

