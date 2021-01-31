library(fdaPDE)
library(mgcv)
#Set your code directory like for instance
#
myconn_folder = "../data/Myconnectome/"
time_series_folder = paste0(myconn_folder,"time_series/")
ROI_number = 60

# Carico i dati per la mesh dei task based
node_myconnectome_17287 <- read.csv(
				paste0(myconn_folder,"node_myconnectome_17287.csv"),
				header=FALSE)
elem_myconnectome_all_tissue_17287 <- read.csv(
		paste0(myconn_folder,"elem_myconnectome_all_tissue_17287.csv"),
		header=FALSE)
# La colonna 5 di elem_BART_all_tissue_small contiene l'indicazione per vedere a quale regione appartiene un tetraedro
# 3 = liquido esterno
# 4 = gray matter
# 5 = white matter


nodes = node_myconnectome_17287
elems = elem_myconnectome_all_tissue_17287[,1:4]

mesh <- create.MESH.3D(nodes=nodes, tetrahedrons = elems)
plot(mesh)

# # Esempio di come estrarre solo una parte della mesh:
# # Conviene farlo nel caso dopo aver definito la ROI
#
# estraggo gli elementi nella gm:
# elem_only_gm = elem_myconnectome_all_tissue_17287[
# 					(elem_myconnectome_all_tissue_17287[,5] == 4),
# 					1:4]    
#
# e i relativi indici:
# elem_only_gm_indexes = sort(unique(c(as.matrix(elem_only_gm))))                         
#
# trovo quali nodi fanno parte della gm:
# nodes_only_gm = nodes[elem_only_gm_indexes,]  
#
# estraggo i nodi nella gm:
# mesh_only_gm <- create.MESH.3D(nodes = nodes_only_gm, 
# 									tetrahedrons = elem_only_gm)

# Carico il vettore con l'indicazione di a che regione del cervello 
# appartengono i nodi
labels_hammer_myconnectome_small <- read.table(
		paste0(myconn_folder,"labels_hammer_myconnectome_new2.csv"),
		quote="\"", comment.char="")

# Per il BART la ROI ? la 24
#
# Per myConnectome, la ROI scelta ? la 59
#
ROI_index = ifelse(labels_hammer_myconnectome_small == ROI_number, 1, 0)
ROI_index = as.logical(ROI_index)

names_session=sprintf('%0.3d', 11:104)
names_session=names_session[-c(42,80)]

datamatrix=NULL

for(k in 1:length(names_session)){
  name=paste(time_series_folder,"ses",names_session[k],"_time_series.csv",sep="")
  time_series=as.matrix(read.csv(name, header=FALSE))
  time_series=t(time_series)
  time_series=time_series[1:518,]
  ROI_ts=time_series[,ROI_index]
  dim(ROI_ts)
  
  cross_sec_avg_ROI_row=as.vector(rowMeans(ROI_ts))
  
  cor_nodes<-NULL
  for(j in 1:mesh$nnodes){
    cor_nodes<-c(cor_nodes,cor(cross_sec_avg_ROI_row,time_series[,j]))
  }
  
  #cor_nodes[which(is.na(cor_nodes))]<-0
  
  r=cor_nodes
  z = 0.5 * log((1+r)/(1-r))
  
  datamatrix<-rbind(datamatrix,z)
  print(k)
}

datamatrix[which(is.na(datamatrix))]=0
save(datamatrix,file=paste0(myconn_folder,"connectivityMaps_myconnectome_",ROI_number,".RData"))

# Tolgo la media dai dati
data_bar = colMeans(datamatrix, na.rm = TRUE)
data_demean = matrix(rep(data_bar, nrow(datamatrix)), nrow = nrow(datamatrix), byrow = TRUE)
datamatrix_demeaned = datamatrix - data_demean

PCA_mv = prcomp(datamatrix_demeaned, center = FALSE)

graphics.off()
