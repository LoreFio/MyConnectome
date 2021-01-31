library(mgcv)
library(rgl)
library(plot3D)
library(plot3Drgl)
library(devtools)
library(car)
library(mvtnorm)
library(faraway)
library(Matrix)

load('../data/mcshapiro.test.RData')

##############################################
##############################################
##                                          ##
##     RELAZIONI TRA I BEHAVIORAL DATA      ##-----------------------------------------------------------
##     E LE CORRELAZIONI TRA LA REGIONI     ##
##                                          ##
##############################################
##############################################
nses=78
nreg=82

#carico dati correlazioni
corr_list=list()
for(j in 1:nses){
  corr_list[[j]] <- read.csv(paste0("../result/CorrMatrices/CorrMatrix_",j,".csv"))
  corr_list[[j]] <- corr_list[[j]][,-1]
}

#carico i dati comportamentali
beh <- read.table('../data/Behavioral/behavioral_data_clean.csv', header=T,sep=',')

#carico gli scores delle comportamentali
beh.scores <- read.csv("../data/behavioral_scores.csv")
head(beh.scores)
beh.scores <- beh.scores[-which(beh[,3]==0),]
colnames(beh.scores) <- c("ses","positivity","energy","surprise")
# la colonna sorpresa è girata giusta?

#########################
#      REGRESSIONI      #-------------------------------------------------------------------------
#########################


#costruzione di 6 matrici, rispettivamente per positività, energia e sorpresa (anche al quadrato)
#che contengano i p_value dei beta delle regressioni lineari (modello completo)
#una settima matrice (sempre 82x82) con l'R2 come altro indicatore di significatività

#nel modello non includiamo le interazioni, poiché le nostre tre componenti sono già ortogonali
pval_pos <- matrix (0,nreg,nreg)
for(j in 1:nreg){
  pval_pos[j,j] <- 1
}   #poiché non ci interessano i dati sulla diagonale, setto tali valori a 1

pval_ene <- pval_pos
pval_sur <- pval_pos
pval_pos2 <- pval_pos
pval_ene2 <- pval_pos
pval_sur2 <- pval_pos
R2 <- matrix (0,nreg,nreg)
v <- rep(0,nses)

for(i in 1:(nreg-1)){
  for(j in (i+1):nreg){
    for(k in 1:nses){
      v[k] <- corr_list[[k]][i,j]
    }
    dset <- as.data.frame(cbind(v,beh.scores))
    regr = lm(v ~ positivity + energy + surprise + I(positivity^2) + I(energy^2) + I(surprise^2),
              data = dset)
    pval_pos[i,j] = summary(regr)$coefficients[2,4]
    pval_ene[i,j] = summary(regr)$coefficients[3,4]
    pval_sur[i,j] = summary(regr)$coefficients[4,4]
    pval_pos2[i,j] = summary(regr)$coefficients[5,4]
    pval_ene2[i,j] = summary(regr)$coefficients[6,4]
    pval_sur2[i,j] = summary(regr)$coefficients[7,4]
    pval_pos[j,i] = pval_pos[i,j]
    pval_ene[j,i] = pval_ene[i,j]
    pval_sur[j,i] = pval_sur[i,j]
    pval_pos2[j,i] = pval_pos2[i,j]
    pval_ene2[j,i] = pval_ene2[i,j]
    pval_sur2[j,i] = pval_sur2[i,j]
    
    R2[i,j] <- summary(regr)$r.squared
    R2[j,i] = R2[i,j]
  }
  print(paste0("Region #",i))
  print(paste0(nreg-i," regression(s) done in the iteration"))
}   

#esplorazione matrici
x11()
par(mfrow=c(2,3))
image(pval_pos, col = rainbow(180)[1:61])
image(pval_ene, col = rainbow(180)[1:61])
image(pval_sur, col = rainbow(180)[1:61])
image(pval_pos2, col = rainbow(180)[1:61])
image(pval_ene2, col = rainbow(180)[1:61])
image(pval_sur2, col = rainbow(180)[1:61])
#le matrici dei pvalue sono interessanti per quanto riguarda le regioni: si vede che ci sono
#delle regioni la cui correlazione con le altre è in parte spiegata (linearmente o quadraticamente)
#da uno stato d'animo.

image(R2, col = rainbow(180)[1:61])
#la matrice di R2 invece dà dei valori piuttosto bassi. C'era da aspettarselo? Sì, sia perché
#è chiaro che la complessità delle interazioni tra le regioni del cervello va molto oltre questa
#analisi, sia perché le nostre 3 componenti già spiegano solo una parte della variabilità
#degli stati d'animo


