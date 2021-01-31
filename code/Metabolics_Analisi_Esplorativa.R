library(mgcv)
library(rgl)
library(plot3D)
library(plot3Drgl)
library(devtools)
library(car)
library(mvtnorm)
#setwd("C:/Users/Marco/Desktop/Politecnico/4.Quartanno/STATAPP/Progetto/data/Myconnectome")
#load('mcshapiro.test.RData')

load('../data/mcshapiro.test.RData')
####################################################
####################################################
##                                                ##
##  LETTURA DEL FILE DELLE VARIABILI METABOLICHE  ##---------------------------------------------------
##             ANALISI ESPLORATIVA                ##
##                                                ##
####################################################
####################################################

# ho modificato leggermente il file in modo da trovarlo più leggibile "per me".
# Sul Drive dovrei aver messo il file così come l'ho modificato e intitolato io.
# la modifica più grande che ho fatto è stata togliere tutte le features che erano contrassegnate da un
# solo codice numerico, sia perché non siamo in grado di dare un nome a queste sostanze, sia perché
# il rapporto tra p ed n era già di 5 a 1 circa, almeno così è "solo" di 2 a 1
# (che è comunque troppo).

#metabol <- read.table('metabolomics_raw_data_modd.csv', header=T,sep=',')
metabol <- read.table('../data/metabolomics_raw_data_modd.csv', header=T,sep=',')

row.names(metabol) <- metabol$X #imposto le date delle osservazioni come nomi delle righe
row.names(metabol)
metabol <- metabol[,-1]    #sistemo in modo che i nomi di righe e colonne corrispondano

metabol <- as.data.frame(t(metabol))  #traspongo la matrice dati
head(metabol)
colnames(metabol)   #ora le sostanze sono i nomi degli 
metabol.sc <- scale(metabol)
metabol.sc <- data.frame(metabol.sc)
# le analisi che seguono vengono fatte con i dati riscalati, però nulla vieta di farlo anche
# senza questa opzione

##########################
#          PCA           #----------------------------------------------------------------------------
##########################
# chiaramente se facessimo una PCA con questo dataset dovremmo ottenere un errore
pc.metabol <- princomp(metabol.sc, scores=T) # ERRORE! p>n
# proviamo a separare un po' di osservazioni
# ad esempio possiamo isolare gli acidi (un esempio non più valido di altri,
# è tutta un'analisi esplorativa)

head(metabol.sc)
acidi <- c(5,13,17,18,19,27,28,29,30,33,35,36,41,47,50,52,53,54,56,57,63,64,65,
           67,68,72,74,76,80,83,85,88,89,91,92,94,97,98,99,101,102,103,104)
metabol.sc.acid <- metabol.sc[,acidi]
pc.metabol.acid <- princomp(metabol.sc.acid, scores=T)
pc.metabol.acid
summary(pc.metabol.acid) #le componenti sono vicine, ma comunue l'80% della varianza è spiegata
#dalle prime 10.
#esito ugualmente poco significativo dovuto al fatto che ho poca differenza
#tra p ed n (il rapporto tra n e p dovrebbe essere più alto)

load.meta <- pc.metabol.acid$loadings
load.meta                #chiaramente anche i loadings hanno poco da dire

# Essendo p>n, possiamo effettuare invece una prcomp
pc.metabol=prcomp(metabol.sc,retx=T)
summary(pc.metabol)      #80% della varianza spiegata da 14 componenti... non esattamente una
                         #cosa gradevole
load.met <- pc.metabol$rotation

# forse viene più interessante eliminando le componenti con varianza troppo elevata
rapp=sapply(metabol,sd)/sapply(metabol,mean)   #coefficienti di variazione
met.smooth <- metabol.sc[,which(rapp<0.5)]  
pc.met.sm=prcomp(met.smooth,retx=T)            
summary(pc.met.sm)                          #non cambia praticamente nulla
# dall'altro lato si potrebbe provare con i componenti a coefficiente di variazione elevato
met.rough <- metabol.sc[,which(rapp>0.5)]  
pc.met.rg=princomp(met.rough,scores=T) 
summary(pc.met.rg)
load.met.rg <- pc.met.rg$loadings
load.met.rg

##########################################
#  RIDUZIONE DI p MEDIANTE CORRELAZIONI  #---------------------------------------------------------
#    PCA DELLE MAGGIORMENTE CORRELATE    #
##########################################
# Versione "acidi":forse con le correlazioni riesco ad eliminare qualcosa
rhoac=cor(metabol.sc.acid)
rhoac

colla=which(abs(rhoac)>0.9, arr.ind=TRUE) #elementi collineari (tanto correlati da essere ridondanti)
redundac=c(9,24,10,36,39,11,15,35) #si poteva rendere automatico con una procedura meno bruta, 
#se dovesse essere più utile rifarlo in futuro metterò due righe
#più adatte
redundac <- redundac[order(redundac)]

# eliminazione delle variabili "ridondanti" (non molto esatta come cosa: sarebbe più corretto per lo meno
# moltiplicare le variabili ricondotte a una sola per il numero delle stesse, ma finché l'esplorazione
# non è specifica, possiamo tollerare questa modifica)
metabol.sc.acid <- metabol.sc.acid[,-redundac]

# Versione completa (rianalizzo il dataset con tutte le features):
rho=cor(metabol.sc)   # la matrice 
image(rho,col=rainbow(12))

# piccolo giochino con le più negativamente correlate
min(rho)
negc=which(rho< -0.7, arr.ind=TRUE)  # estrapolazione delle componenti
                                     # si potrebbe fare una funzione che evita di farlo manualmente?
x11()
matplot(metabol.sc[,c(36,82)],type='l',lty=1,lwd=2,col=c('blue','red'))
legend('topright',legend=colnames(metabol.sc)[c(36,82)],lty=1,lwd=2,col=c('blue','red'))
dev.off()

# E anche con le più positivamente correlate
coll=which(rho > 0.9, arr.ind=TRUE)
coll=coll[which(coll[,1]<coll[,2]),]
coll  

gr1=c(18,30,67)
gr2=c(29,35,50)
gr3=c(28,33,94,99,100)
gr4=c(89,92)

x11(width=30,height=16)
par(mfrow=c(2,2))
matplot(metabol.sc[,gr1],type='l',lty=1,lwd=2,col=rainbow(length(gr1)))
legend('topright',legend=colnames(metabol.sc)[gr1],lty=1,lwd=2,col=rainbow(length(gr1)),cex=0.8)
matplot(metabol.sc[,gr2],type='l',lty=1,lwd=2,col=rainbow(length(gr2)))
legend('topright',legend=colnames(metabol.sc)[gr2],lty=1,lwd=2,col=rainbow(length(gr2)),cex=0.8)
matplot(metabol.sc[,gr3],type='l',lty=1,lwd=2,col=rainbow(length(gr3)))
legend('topright',legend=colnames(metabol.sc)[gr3],lty=1,lwd=2,col=rainbow(length(gr3)),cex=0.8)
matplot(metabol.sc[,gr4],type='l',lty=1,lwd=2,col=rainbow(length(gr4)))
legend('topright',legend=colnames(metabol.sc)[gr4],lty=1,lwd=2,col=rainbow(length(gr4)),cex=0.8)
graphics.off()

# si può fare una riduzione più intelligente guardando alla pca di questi gruppi
met.gr1 <- metabol.sc[,gr1]  
pc.met.gr1=princomp(met.gr1,scores=T) 
summary(pc.met.gr1)
load.met.gr1 <- pc.met.gr1$loadings      #la prima componente spiega il 95%
load.met.gr1
v1 <- as.matrix(met.gr1) %*% t(t(load.met.gr1[,1]))     # nx1
v1 <- as.data.frame(v1)

met.gr2 <- metabol.sc[,gr2]  
pc.met.gr2=princomp(met.gr2,scores=T) 
summary(pc.met.gr2)
load.met.gr2 <- pc.met.gr2$loadings      #la prima componente spiega il 96%
load.met.gr2
v2 <- as.matrix(met.gr2) %*% t(t(load.met.gr2[,1]))     # nx1
v2 <- as.data.frame(v2)

met.gr3 <- metabol.sc[,gr3]  
pc.met.gr3=princomp(met.gr3,scores=T) 
summary(pc.met.gr3)
load.met.gr3 <- pc.met.gr3$loadings      #la prima componente spiega il 94%
load.met.gr3
v3 <- as.matrix(met.gr3) %*% t(t(load.met.gr3[,1]))     # nx1
v3 <- as.data.frame(v3)

met.gr4 <- metabol.sc[,gr4]  
pc.met.gr4=princomp(met.gr4,scores=T) 
summary(pc.met.gr1)
load.met.gr4 <- pc.met.gr4$loadings      #la prima componente spiega il 95%
load.met.gr4
v4 <- as.matrix(met.gr4) %*% t(t(load.met.gr4[,1]))     # nx1
v4 <- as.data.frame(v4)


p=dim(metabol.sc)[2]
remove=NULL
for(i in 1:p)
{
  if(i %in% gr1 || i %in% gr2 || i %in% gr3 || i %in% gr4){
    remove <- c(remove,i)
  }
}
metabol.red=cbind(metabol.sc[-remove],v1,v2,v3,v4)
colnames(metabol.red)[c(94,95,96,97)]=c('V1','V2','V3','V4')
colnames(metabol.red)
# e ora... rifaccio
rhored=cor(metabol.red)
x11()
image(rhored,col=rainbow(12))
# mi aspetto in verità che non vi siano più features correlate oltre allo 0.9, perché altrimenti
# sarebbero comparse prima. Avremmo forse potuto rimuovere più valori con una pca che ci faceva
#"spendere di più". Abbiamo salvato circa il 95% delle informazioni quando anche l'80% basterebbe.

coll=which(abs(rhored) > 0.7, arr.ind=TRUE)
coll=coll[which(coll[,1]<coll[,2]),]
coll

cl1=c(10,11,18,76,82,3,38,44,48,51,52,35,95,7,25,28,40,58,93)
cl2=c(12,19,83,84)
cl3=c(26,53,54,80)
cl4=c(30,74,81,96)
cl5=c(86,87,90)

x11(width=30,height=20)
par(mfrow=c(2,3))
matplot(metabol.red[,cl1],type='l',lty=1,lwd=2,col=rainbow(length(cl1)))
legend('topright',legend=colnames(metabol.red)[cl1],lty=1,lwd=2,col=rainbow(length(cl1)),cex=0.6)
matplot(metabol.red[,cl2],type='l',lty=1,lwd=2,col=rainbow(length(cl2)))
legend('topright',legend=colnames(metabol.red)[cl2],lty=1,lwd=2,col=rainbow(length(cl2)),cex=0.6)
matplot(metabol.red[,cl3],type='l',lty=1,lwd=2,col=rainbow(length(cl3)))
legend('topright',legend=colnames(metabol.red)[cl3],lty=1,lwd=2,col=rainbow(length(cl3)),cex=0.6)
matplot(metabol.red[,cl4],type='l',lty=1,lwd=2,col=rainbow(length(cl4)))
legend('topright',legend=colnames(metabol.red)[cl4],lty=1,lwd=2,col=rainbow(length(cl4)),cex=0.6)
matplot(metabol.red[,cl5],type='l',lty=1,lwd=2,col=rainbow(length(cl5)))
legend('topright',legend=colnames(metabol.red)[cl5],lty=1,lwd=2,col=rainbow(length(cl5)),cex=0.6)

#seconda PCA
met.cl1 <- metabol.red[,cl1]  
pc.met.cl1=princomp(met.cl1,scores=T) 
summary(pc.met.cl1)
load.met.cl1 <- pc.met.cl1$loadings      #le prime 4 componenti spiegano l'82% (su 19)
load.met.cl1
w1 <- as.matrix(met.cl1) %*% t(t(load.met.cl1[,c(1,2,3,4)]))     # nx4
w1 <- as.data.frame(w1)
colnames(w1)=c('w11','w12','w13','w14')
head(w1)

met.cl2 <- metabol.red[,cl2]  
pc.met.cl2=princomp(met.cl2,scores=T) 
summary(pc.met.cl2)
load.met.cl2 <- pc.met.cl2$loadings      #le prime 2 componenti spiegano quasi il 90% (su 4)
load.met.cl2
w2 <- as.matrix(met.cl2) %*% t(t(load.met.cl2[,c(1,2)]))     # nx2
w2 <- as.data.frame(w2)
colnames(w2)=c('w21','w22')
head(w2)

met.cl3 <- metabol.red[,cl3]  
pc.met.cl3=princomp(met.cl3,scores=T) 
summary(pc.met.cl3)
load.met.cl3 <- pc.met.cl3$loadings      #le prime 2 componenti spiegano l'88% (su 4)
load.met.cl3
w3 <- as.matrix(met.cl3) %*% t(t(load.met.cl3[,c(1,2)]))     # nx3
w3 <- as.data.frame(w3)
colnames(w3)=c('w31','w32')
head(w3)

met.cl4 <- metabol.red[,cl4]  
pc.met.cl4=princomp(met.cl4,scores=T) 
summary(pc.met.cl4)
load.met.cl4 <- pc.met.cl4$loadings      #le prime 2 componenti spiegano il 94% (su 4)
load.met.cl4
w4 <- as.matrix(met.cl4) %*% t(t(load.met.cl4[,c(1,2)]))     # nx2
w4 <- as.data.frame(w4)
colnames(w4)=c('w41','w42')
head(w4)

met.cl5 <- metabol.red[,cl5]  
pc.met.cl5=princomp(met.cl5,scores=T) 
summary(pc.met.cl5)
load.met.cl5 <- pc.met.cl5$loadings      #la prima componente spiega l'83% (su 3)
load.met.cl5
w5 <- as.matrix(met.cl5) %*% t(t(load.met.cl5[,1]))     # nx1
w5 <- as.data.frame(w5)
colnames(w5)=c('w51')
head(w5)

p=dim(metabol.red)[2]
remove=NULL
for(i in 1:p)
{
  if(i %in% cl1 || i %in% cl2 || i %in% cl3 || i %in% cl4 || i %in% cl5){
    remove <- c(remove,i)
  }
}
metabol.red=cbind(metabol.red[-remove],w1,w2,w3,w4,w5)
colnames(metabol.red)   #ora p è sceso a 74, e ciò è molto buono



# Versioni "liscia" e "ruvida" (con l'originale a 106 features: al momento è solo per esplorazione)
rhosm=cor(met.smooth)
rhorg=cor(met.rough)
x11(width=20,height=10)
par(mfrow=c(1,2))
image(rhosm,col=rainbow(12),xlab='Smooth features')
image(rhorg,col=rainbow(48),xlab='Rough features')



#TEST DI GAUSSIANITA'-----------------------------------------------------------------------------------
mcshapiro.test(metabol.sc.acid) # chiaramente le variabili sono troppe e il p-value è 0.
# proviamo a vedere di variabile in variabile
n <- dim(metabol.sc.acid)[1]
p <- dim(metabol.sc.acid)[2]
pval= rep(0,p)

for(i in 1:p)
  pval[i] <- shapiro.test(metabol.sc.acid[,i])$p.value
pval  #9 variabili potrebbero essere prese come gaussiane senza effettuare trasformazioni
nontransf=which(pval>0.05)
#prova di boxcox
lambda <- powerTransform(metabol.sc.acid, family="bcnPower")  #occhio che non è immediata
mac.bc = metabol.sc.acid

for(i in 1:p)
{
  if(is.element(i,nontrasf)==FALSE)
    mac.bc[,i] <- bcnPower(metabol.sc.acid[,i], lambda$lambda[i], gamma = lambda$gamma[i])
}

#riproviamo lo shapiro (univariato per ora)
pval.bc= rep(0,p)

for(i in 1:p)
  pval.bc[i] <- shapiro.test(mac.bc[,i])$p.value
pval.bc   #la situazione pare migliorata, ora le variabili prendibili singolarmente come gaussiane
#sono 19

#vediamo che comportamento hanno le variabili assumibili come gaussiane
normindex=which(pval.bc>0.05)
normmac=metabol.sc.acid[,normindex]
head(normmac)

mcshapiro.test(normmac)  #pvalue sotto 0.05, ma (appena) sopra 0,01
# è già qualcosa: ad essere molto restrittivi, non c'è tutta l'evidenza per rifiutare l'ipotesi di
# gaussianità, anche se la cosa rimane tutt'altro che soddisfacente. Si potrebbe provare anche
# alzando la soglia del pvalue per le colonne prese singolarmente (a 0.1 o 0.2 ad esempio)
# in questo modo presumibilmente ci saranno meno colonne e i dati avanzati saranno globalmente
# più prendibili come gaussiani.

# PROVA DI CLUSTERING-----------------------------------------------------------------------------------
# Riprendo il dataset precedente a 35 colonne
n_clusters = 10
n_random_starts = 48
max_iter = 15
kac = kmeans(metabol.sc.acid, n_clusters, iter.max = max_iter, nstart = n_random_starts)
kac #dato che i cluster sono in un certo senso dei "levels" per le osservazioni, potrebbe essere
# interessante vedere le interazioni con altre variabili di altri dataset che contengono
# dati categorici.

# ANDAMENTO DELLE SOSTANZE NEL TEMPO--------------------------------------------------------------------
# vediamo l'andamento del tempo delle sostanze (prova con le prime 10 o viene un carnaio)
x11()
matplot(metabol[,1:10],type='l',lty=rep(c(1,2,3,4,5),2),col=rainbow(10))
legend('topright',legend=colnames(metabol)[1:10],lty=rep(c(1,2,3,4,5),2),col=rainbow(10),cex=0.6)
dev.off() 
#bello, ma si potrebbe prima fare un test (per ogni variabile) per verificare se c'è evidenza per
#una variazione di certi componenti.
n=dim(metabol.sc)[1]
p=dim(metabol.sc)[2]
diffrates <- as.data.frame(matrix(0, n-1, p)) 
for(i in 1:(n-1)){
  diffrates[i,] <- (metabol.sc[i+1,] - metabol.sc[i,])
}

colnames(diffrates) <- colnames(metabol.sc)
diffrates[1:8,1:4]

mcshapiro.test(diffrates)     #ERRORE! A quanto pare la matrice è singolare
# si può comunque provare a vedere il p_value per le singole componenti (più interpretabile)
pvald= rep(0,p)
for(i in 1:p)
  pvald[i] <- shapiro.test(diffrates[,i])$p.value
pvald    #la maggior parte delle componenti è gaussiana

#analizziamo quelle assumibili come gaussiane in particolare
alpha=0.01
diffrates <- diffrates[-which(pvald<alpha)] #rimangono 81 componenti su 106 per alpha=0.01

# facciamo ora il test (H0: dc=0, H1: dc!=0) UNIVARIATO (per ogni componente) dove dc è la media
# visto che ho voglia di una forte evidenza, valuterei sempre alpha=0.01
# (ma si può interpretare con altri valori ovviamente)
g=dim(diffrates)[2]
pvalt=rep(0,g)
for(i in 1:g){
  dc <- mean(diffrates[,i])
  sdd <- sd(diffrates[,i])
  tstat <- abs(dc-0)/sdd
  pvalt[i] <- 2*(1-pt(tstat,n-1))
}

pvalt  
min(pvalt) #poiché tutti i p_value sono tutti sopra 0.897, si può concludere con abbastanza
           #chiarezza che nessuna componente varia nel tempo. Questa conclusione è comunque
           #utile per dire che i valori delle componenti non dipendono dal tempo almeno singolarmente

#È allora possibile che certe features siano più costanti (oscillino meno) di altre?
sapply(metabol.sc,mean)
sapply(metabol.sc,sd)
# i dati scalati hanno tutti media 0 e varianza 1 giustamente, proviamo a vedere il rapporto tra
# deviazione standard e media (con i dati non scalati), ovvero il coefficiente di variazione
rapp=sapply(metabol,sd)/sapply(metabol,mean)
rapp
im=which(rapp==min(rapp))
iM=which(rapp==max(rapp))

x11()
matplot(metabol.sc[,c(im,iM)],type='l',lty=c(1,1),lwd=c(2,2),col=c('blue','red'))
legend('topright',legend=colnames(metabol.sc)[c(im,iM)],lty=c(1,1),lwd=c(2,2),col=c('blue','red'))
dev.off()
# i grafici sono molto diversi... in particolare nella sessione 20
# (si può vedere anche la versione scalata dato che i grafici possono altrimenti avere odg differenti)
# vediamo in generale le variabili con maggiore variabilità

iG1=which(rapp>1.5)
iG2=which(rapp>0.75 & rapp<1)
x11(width=20,height=10)
par(mfrow=c(1,2))
matplot(metabol.sc[iG1],type='l',lty=1,lwd=1,col=rainbow(length(iG1)))
legend('topright',legend=colnames(metabol.sc)[iG1],lwd=1,col=rainbow(length(iG1)))
matplot(metabol.sc[iG2],type='l',lty=1,lwd=1,col=rainbow(length(iG2)))
legend('topright',legend=colnames(metabol.sc)[iG2],lwd=1,col=rainbow(length(iG2)),cex=0.6)
dev.off()

# ci sono svariati picchi in certe osservazioni, anche condivisi con altri acidi

############################################
# ANALISI DELLE COMPONENTI "INTERESSANTI"  #
############################################
# questa analisi è sommaria perché non dovrebbe risultare molto utile

#urea, tryptophan, taurine, sucrose, oleic acid, lactic acid, nicotinic acid, glucose, fructose,
#cholesterol, glycerol, succinic acid, alanine, lysine.
interesting <- c(6,8,15,16,17,35,36,49,53,71,75,77,84,95)
metabol.irt <- metabol.sc[,interesting]

#PCA?
pc.irt=princomp(metabol.irt,scores=T)   #le prime 4 componenti su 14 spiegano il 60% dei dati
summary(pc.irt)
load.irt <- pc.irt$loadings      
load.irt

x11()    #loadings delle prime 4 componenti
par(mfrow = c(4,1))
for(i in 1:4) barplot(load.irt[,i], ylim = c(-1, 1), col=rainbow(length(interesting)))
# not very interesting actually...


## CONCLUSIONI------------------------------------------------------------------------------------------

## COSE CHE SI POTREBBERO FARE
#1.ricerca degli outlier
#2.analizzare i rapporti con l'altro dataset, quello delle variabili ambientali
#3.capire quali sostanze si possono effettivamente confrontare con le connectivity maps, ovvero
#  le sostanze che rimangono stabili per lo meno nell'arco di una giornata (non abbiamo fonti)
#4.vedere se possono esserci rapporti con le connectivity maps (prima o poi bisogna farlo)
