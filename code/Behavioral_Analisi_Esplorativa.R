library(mgcv)
library(rgl)
library(plot3D)
library(plot3Drgl)
library(devtools)
library(car)
library(mvtnorm)
library(faraway)
library(Matrix)
source("EM.R")
load('../data/mcshapiro.test.RData')
load('mcshapiro.test.RData')
setwd("C:/Users/Marco/Desktop/Politecnico/4.Quartanno/STATAPP/Progetto/data/Myconnectome")

#####################################################
#####################################################
##                                                 ##
##      BEHAVIORAL DATA: ANALISI ESPLORATIVA       ##---------------------------------------------------
##                                                 ##
#####################################################
#####################################################

#Dataset che contiene le sole features comportamentali
beh <- read.table('behavioral_data_clean.csv', header=T,sep=',')
beh <- read.table('../data/behavioral_data_clean.csv', header=T,sep=',')
colnames(beh)
row.names(beh) <- beh[,1]
row.names(beh)
beh <- beh[,-c(1,2,3,4)]

n=dim(beh)[1]
p=dim(beh)[2]

###################
#       PCA       #-------------------------------------------------------------------------------------
###################
pc.beh <- princomp(beh,scores=T)
summary(pc.beh)
load.beh <- pc.beh$loadings      #la prime 3 componenti spiegano già l'85%. Probabile influenza delle
load.beh                         #componenti fatigue, positive e negative, che hanno valori più alti

#Proviamo a rimuoverle
colnames(beh)
behr <- beh[,-c(31,45,47)]
n=dim(behr)[1]
p=dim(behr)[2]
#rifacciamo la pca
pc.behr <- princomp(behr,scores=T)
pc.summ <- summary(pc.behr)
load.behr <- pc.behr$loadings     #ora le prime 3 spiegano il 60%, per arrivare all'80 ce ne vogliono 11

# loadings di ogni stato d'animo in ordine crescente, per interpretare le componenti
load.behr[order(load.behr[,1]),1]  #positività/energia vs negatività/stanchezza        (40,5%)
load.behr[order(load.behr[,2]),2]  #energia/nervosismo vs stanchezza "non negativa"    (12,5%)
load.behr[order(load.behr[,3]),3]  #sorpresa                                           (6,8%)
load.behr[order(load.behr[,4]),4]  #tensione vs distensione                            (5,1%)
load.behr[order(load.behr[,5]),5]  #pare non spiegare nulla
# in conclusione, le prime quattro componenti della pca sono interpretabili, e questo dà l'idea che
# analizzando assieme ai dati delle fMRI, possa uscire qualcosa di interessante, anche in merito
# ad alcuni degli articoli letti.

# plot dei loadings: CREO MATRICE COLORI
colore=rep(1,p)%*%t(rep(1,4))
for(i in 1:p){
  if(load.behr[i,1]>0.17)    colore[i,1]=3
  if(load.behr[i,1]< -0.08)  colore[i,1]=2
  if(load.behr[i,2]>0.14)    colore[i,2]="orange"
  if(load.behr[i,2]< -0.181) colore[i,2]=4
  if(load.behr[i,3]>0.1)     colore[i,3]=5
  if(load.behr[i,3]< -0.2)   colore[i,3]=6
  if(load.behr[i,4]>0.154)   colore[i,4]="salmon"
  if(load.behr[i,4]< -0.1)   colore[i,4]="forest green"
}

coloreb=rep(8,p)%*%t(rep(1,4))
for(i in 1:p){
  if(load.behr[i,1]>0.17)    coloreb[i,1]=3
  if(load.behr[i,1]< -0.08)  coloreb[i,1]=2
  if(load.behr[i,2]>0.14)    coloreb[i,2]="orange"
  if(load.behr[i,2]< -0.181) coloreb[i,2]=4
  if(load.behr[i,3]>0.1)     coloreb[i,3]=5
  if(load.behr[i,3]< -0.2)   coloreb[i,3]=6
  if(load.behr[i,4]>0.154)   coloreb[i,4]="salmon"
  if(load.behr[i,4]< -0.1)   coloreb[i,4]="forest green"
}
vp=seq(0.68,1.194*p,len=p)

#GRAFICO ORDINATO
x11(height=30,width=60)
par(mar=c(3,3,3,2), mfrow=c(1,4))
for(i in 1:4)
{
  barplot(load.behr[order(load.behr[,i]),i],xlim=c(-1,1),names.arg='',
          col=coloreb[order(load.behr[,i]),i],main=paste0("Component #",i),horiz=T)
  text(y=vp,x=ifelse(load.behr[order(load.behr[,i]),i]>0,-0.07,0.07),
       pos=ifelse(load.behr[order(load.behr[,i]),i]>0,2,4),srt=0,
      labels=names(load.behr[order(load.behr[,i]),i]),col=colore[order(load.behr[,i]),i])
}

#GRAFICO NON ORDINATO
x11(height=30,width=60)
par(mar=c(3,3,3,2), mfrow=c(1,4))
for(i in 1:4)
{
  barplot(load.behr[,i],xlim=c(-1,1),names.arg='',col=coloreb[,i],
          main=paste0("Component #",i),horiz=T)
  text(y=vp,x=ifelse(load.behr[,i]>0,-0.07,0.07),
       pos=ifelse(load.behr[,i]>0,2,4),srt=0,
       labels=names(load.behr[,i]),col=colore[,i])
}


# Il fatto che dopo queste prime 4 componenti solo il 65% della varianza sia spiegato, può essere
# dovuto al fatto che parecchie features rappresentano stati .
# Inoltre, tutte le features assumono un range di valori decisamente limitato e il discernimento
# tra l'uno e l'altro di questi valori ha delle elevate componenti soggettive, per cui questo risultato
# può essere considerato buono.


############################
#   REGRESSIONI LINEARI    #----------------------------------------------------------------------
############################

#analizzo i termini "positive, negative, fatigue" che potrebbero essere semplicemente
#somme di altre features

#regressione positività
positivity <- lm(positive ~ .,data=beh)
summary(positivity)

#regressione negatività
negativity <- lm(negative ~ .,data=beh)
summary(negativity)

#regressione fatica
fatica <- lm(fatigue ~ .,data=beh)
summary(fatica)
#come pensavamo, c'è un fit perfetto (positive, negative e fatigue sono somme di altre features)

##########################
#      CORRELAZIONE      #----------------------------------------------------------------------------
##########################

#Matrice di correlazione
rhob=cor(behr)
image(rhob,col=gray.colors(50,start=0,end=1)) # è una bella matrice a blocchi
coll=which(rhob > 0.9, arr.ind=TRUE)          #per segnare gli indici più positivamente correlati
coll=coll[which(coll[,1]!=coll[,2]),]
coll  
coln=which(rhob < -0.65, arr.ind=TRUE)        #per segnare gli indici più negativamente correlati
coln=coln[which(coln[,1]!=coln[,2]),]
coln
#si può ipotizzare che le informazioni ottenibili da tale tipo di analisi
#siano già spiegate almeno in parte dalla pca: più utile per ora cominciare a unire questo
#dataset a quelli delle connectivity maps.

######################
#     CLUSTERING     #----------------------------------------------------------------------------
######################

#CREO VETTORI PRIME COMPONENTI
behcomp <- as.data.frame(as.matrix(behr)%*%as.matrix(load.behr[,c(1,2,3,4)]))   # 81x4
behcomp[,3] <- -1*behcomp[,3]    #inverto la componente della sorpresa per interpretabilità
colnames(behcomp)=c('positivity','activity','surprise','impulse-response')
head(behcomp)
#riscalo
behcomp.sc=as.data.frame(scale(behcomp))
head(behcomp.sc)

#ALTERNATIVA 1: pesare le colonne delle componenti per la loro proporzione di varianza
compsd <- pc.behr$sdev[1:4]
compvar <- c(0.4053207, 0.1251837, 0.0681795, 0.05141865)
for(i in 1:4) behcomp.sc[,i] <- behcomp.sc[,i]*compvar[i]
head(behcomp.sc)

#ALTERNATIVA 2: ignorare quarta componente
behcomp.sc <- behcomp.sc[,-4]

#vediamo se si possono ottenere dei cluster di osservazioni
library(MASS)
behcomp.e <- dist(behcomp.sc, method='euclidean')

n=dim(behcomp.sc)[1]
#matrice distanze
x11()
image(1:n,1:n,as.matrix(behcomp.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )

#cluster e dendrogrammi
behcomp.es <- hclust(behcomp.e, method='single')
behcomp.ea <- hclust(behcomp.e, method='average')
behcomp.ec <- hclust(behcomp.e, method='complete')

x11()
par(mfrow=c(1,3))
plot(behcomp.es, main='euclidean-single',hang=-0.1,xlab='',labels=F,cex=0.6, sub='') #1 cluster solo
plot(behcomp.ec, main='euclidean-complete',hang=-0.1,xlab='',labels=F,cex=0.6, sub='') #forse 3
plot(behcomp.ea, main='euclidean-average',hang=-0.1,xlab='',labels=F,cex=0.6, sub='') #ancora meglio, 3

#per rappresentare i dati, dovremmo escludere la quarta componente
#se non considerassimo la quarta componente ma solo le prime 3, si individuerebbero in verità 5 cluster
#anziché 3

#calcolo di affidabilità
coph.es <- cophenetic(behcomp.es)
coph.ec <- cophenetic(behcomp.ec)
coph.ea <- cophenetic(behcomp.ea)

es <- cor(behcomp.e, coph.es)
ec <- cor(behcomp.e, coph.ec)
ea <- cor(behcomp.e, coph.ea)

c("Eucl-Single"=es,"Eucl-Compl."=ec,"Eucl-Ave."=ea)  #0.78, 0.85, 0.80 rispettivamente
#Complessivamente sono dei valori molto buoni

#PLOT
plot3d(behcomp.sc, size=6, col='blue', aspect = F)

open3d()
clusters <- cutree(behcomp.ec, 5)
plot3d(behcomp.sc, size=6, col=clusters+1, aspect = F) 
#cluster1 (rosso):   giorni comuni (relativamente positivi e senza elementi di sorpresa)
#cluster2 (verde):   giorni negativi, ma ad alto carico di energia
#cluster3 (blu):     giorni di sorpresa rilassati
#cluster4 (ciano):   giorni di netta sorpresa e positività
#cluster5 (magenta): giorni di generale stanchezza

open3d()
clusters <- cutree(behcomp.ea, 5)
plot3d(behcomp.sc, size=6, col=clusters+1, aspect = F)
#cluster1 (rosso):giorni generici
#cluster2 (verde): come prima
#cluster3 (blu): come il ciano
#cluster4 (ciano): come il magenta
#cluster5 (magenta): un'unica osservazione, un giorno di alta energia caratterizzato da sorpresa

