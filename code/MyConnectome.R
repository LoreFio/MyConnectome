library(fdaPDE)
# Il pacchetto fdaPDE essere la versione scaricata da 
# https://github.com/NegriLuca/fdaPDE-manifold
# per funzionare con dati 3D. 
library(mgcv)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#### Carico i dati per la mesh ##########
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

node_brain_vol = as.matrix(read.csv(
					paste0(myconn_folder,"node_myconnectome_17287.csv"),
					header=FALSE))
csf_brain_vol = as.matrix(read.csv(
				paste0(myconn_folder,"elem_myconnectome_csf_17287.csv"),
				header=FALSE))
gm_brain_vol  = as.matrix(read.csv(
				paste0(myconn_folder,"elem_myconnectome_gm_17287.csv"),
				header=FALSE))
wm_brain_vol  = as.matrix(read.csv(
				paste0(myconn_folder,"elem_myconnectome_wm_17287.csv"),
				header=FALSE))
all_brain_vol = as.matrix(read.csv(
				paste0(myconn_folder,"elem_myconnectome_all_17287.csv"),
				header=FALSE))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#### Creo la mesh #### #### #### #### ###
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
mesh = create.MESH.3D(nodes = node_brain_vol,
						tetrahedrons = all_brain_vol)
plot(mesh)

dim(node_brain_vol)
dim(all_brain_vol)

# ROI per il myConnectome e esportazione della regione per 
# visualizzarla in ParaView

load(paste0(myconn_folder,"ROI_myConnectome.RData"))
# Creo la base elementi finiti a partire dalla mesh
FEMbasis = create.FEM.basis(mesh)
# Creo l'oggetto FEM
FEM_ROI  = FEM(ROI_myConnectome, FEMbasis)
# Salvo la ROI in vtu per poterla aprire con ParaView
source('fdaPDE.write.vtk.R')
write.vtu(FEM_ROI, file = paste0(vtu_folder,"ROI_myConnectome.vtu"))


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#### Carico la datamatrix creata myConnectome ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

load(paste0(myconn_folder,"connectivityMaps_myconnectome.RData"))

# Plot della di correlazione media delle funzioni ricavate dalle 
# fMRI per ogni sessione sul paziente

corr.media = colMeans(datamatrix, na.rm = TRUE)
FEMbasis <- create.FEM.basis(mesh)
FEMcorr.media=FEM(corr.media,FEMbasis)
plot(FEMcorr.media)

# Plot di una connectivity map
FEMscan1 = FEM(datamatrix[1,],FEMbasis)
plot(FEMscan1)
