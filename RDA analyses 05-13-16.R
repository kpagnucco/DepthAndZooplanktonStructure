rm(list=ls(all=TRUE))
#set working directory

getwd()
setwd("/Users/ckluk/Desktop/Post doc/PROJECTS/Invert Proxy Paper/BythoRDA/Compare")

# Load required packages
library(ade4)
library(vegan)
library(recluster)
library(ape)

# Script for running stepwise RDA 
# Without fish P/A data

# Import Hellinger transformed size class data (S1 to S4 only)
Biotic <- read.csv("HellBiomSizeClassesLong.csv", row.name = 1)
dim(Biotic)
#include only size classes 1 to 4
ZSS<-Biotic[c(1:4)]
#change names of columns to get rid of Hell
colnames(ZSS)=c("S1","S2","S3","S4")

# Import ENV data
Envdata <- read.csv("Envdatashortlong.csv", row.name = 1)
dim(Envdata)
# includes only CAL, TPL, SIO3, DOCL, pH, INVERT, DEPTHL, AREAL, BYTHO, BYTHO_P/A
# drop BYTHO_P/A for RDA
Envdatshort<-Envdata[c(1:9)]

# For Variation partitioning, set up a Bytho matrix. 
# Best to use Bytho abundances
Bytho<-Envdata[c(10)]

#set up Env.nobytho
Env.no.bytho<-Envdata[c(1:8)]


# Redundancy analysis 
# ------------------------------------
## RDA
ZSS.all <- rda(ZSS ~., Envdatshort[c(1:5,7:9)], scale=T)
#test significance of the canonical relationship
anova(ZSS.all, by = "axis")
#test the significance of each environmental variable
anova(ZSS.all, by = "margin")
ZSS.all
RsquareAdj(ZSS.all)


#with variable selection
ZSS.step.wise = ordistep(rda(ZSS~ 1, Envdatshort[c(1:5,7:9)], scale = T), scope = formula(ZSS.all), pstep = 1000, perm.max = 1000)

#get R.squared
RsquareAdj(ZSS.step.wise)

#RDA of community composition by environmental variables chosen
#Change the columns below based on results of first step (sio3, area, ph, depth)
ZSS.rda.result = rda(ZSS ~., Envdatshort[,c(3,5,7,8)], scale=T)
ZSS.rda.result

#test significance of the canonical relationship
anova(ZSS.rda.result, by = "axis")

#test the significance of each environmental variable
anova(ZSS.rda.result, by = "margin")

#RDA triplot
plot(ZSS.rda.result, scaling = 3, main="Triplot no-fish-RDA")
spe.sc = scores(ZSS.rda.result, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

RsquareAdj(ZSS.rda.result)

# Script for running RDA and Variation Partitioning
# With fish P/A data
# --------------------------------------------------#

# Import Hellinger transformed size class data (S1 to S4 only)
fish.Biotic <- read.csv("ZooSizeClasses_NONHell.csv", row.name = 1)
dim(fish.Biotic)
#include only size classes 1 to 4
fish.SzBiotic<-fish.Biotic[c(2:5)]
#transform size biomass data (Hellinger)
ZSS.fish<-decostand(fish.SzBiotic,"hellinger")

# Import ENV data
fish.EnvMaxTL <- read.csv("EnvFishMaxTL.csv", row.name = 1)
dim(fish.EnvMaxTL)
# includes only INVERT + SIG ENV
fish.Env.Invert<-fish.EnvMaxTL[c(3,5,7,8,6)]
# includes only MaxTL + SIG ENV
fish.Env.MaxTL<-fish.EnvMaxTL[c(3,5,7,8,10)]
# includes only INVERT and MaxTL + SIG ENV
fish.Env.Invert.MaxTL<-fish.EnvMaxTL[c(3,5,7,8,6,10)]
# includes only SIG ENV 
fish.Env<-fish.EnvMaxTL[c(3,5,7,8)]
# includes only MaxTL
fish.MaxTL<-fish.EnvMaxTL[c(10)]
# includes only Invert
fish.Invert<-fish.EnvMaxTL[c(6)]
# includes only Invert and MaxTL
fish.Invert.MaxTL<-fish.EnvMaxTL[c(6,10)]

###CAN ONE INDEX REPLACE ANOTHER###
#RDA1 with INVERT
RDA1 <- rda(ZSS.fish ~., fish.Invert, scale=T)
#test significance of the canonical relationship
anova(RDA1, by = "axis")
#test the significance of each environmental variable
anova(RDA1, by = "margin")
RsquareAdj(RDA1)

#RDA1 triplot
plot(RDA1, scaling = 3, main="RDA1")
spe.sc = scores(RDA1, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

#RDA2 with MaxTL
RDA2 <- rda(ZSS.fish ~., fish.MaxTL, scale=T)
#test significance of the canonical relationship
anova(RDA2, by = "axis")
#test the significance of each environmental variable
anova(RDA2, by = "margin")
RsquareAdj(RDA2)

#RDA2 triplot
plot(RDA2, scaling = 3, main="RDA2")
spe.sc = scores(RDA2, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

#RDA3 with INVERT and MaxTL
RDA3 <- rda(ZSS.fish ~., fish.Invert.MaxTL, scale=T)
#test significance of the canonical relationship
anova(RDA3, by = "axis")
#test the significance of each environmental variable
anova(RDA3, by = "margin")
RsquareAdj(RDA3)

#RDA3 triplot
plot(RDA3, scaling = 3, main="RDA3")
spe.sc = scores(RDA3, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

#RDA4 with Env
RDA4 <- rda(ZSS.fish ~., fish.Env, scale=T)
#test significance of the canonical relationship
anova(RDA4, by = "axis")
#test the significance of each environmental variable
anova(RDA4, by = "margin")
RsquareAdj(RDA4)

#RDA4 triplot
plot(RDA4, scaling = 3, main="RDA4")
spe.sc = scores(RDA4, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

# --------------------------------------------------#

##Var Partitioning to compare all four RDAS
#RDA1 vs RDA2
RDA1.RDA2.vp<-varpart(ZSS.fish, fish.Invert, fish.MaxTL, scale=T)
summary(RDA1.RDA2.vp)
plot(RDA1.RDA2.vp)
showvarparts(2)
RDA1.RDA2.vp
#test significance of MaxTL without Invert: X2barX1
MaxTL.bar.Invert =  rda(ZSS.fish ~ ., data=fish.MaxTL)
anova(MaxTL.bar.Invert)
#test significances of Invert without MaxTL: X1barX2
Invert.bar.MaxTL =  rda(ZSS.fish ~ ., data=fish.Invert)
anova(Invert.bar.MaxTL)

#RDA1 vs RDA3
RDA1.RDA3.vp<-varpart(ZSS.fish, fish.Invert, fish.Invert.MaxTL, scale=T)
summary(RDA1.RDA3.vp)
plot(RDA1.RDA3.vp)
showvarparts(2)
RDA1.RDA3.vp
#test significance of MaxTL without Invert: X2barX1
InvertMaxTL.bar.Invert =  rda(ZSS.fish ~ ., data=fish.Invert.MaxTL)
anova(InvertMaxTL.bar.Invert)
#test significances of Invert without MaxTL: X1barX2
Invert.bar.InvertMaxTL =  rda(ZSS.fish ~ ., data=fish.Invert)
anova(Invert.bar.InvertMaxTL)

#RDA1 vs RDA4
RDA1.RDA4.vp<-varpart(ZSS.fish, fish.Invert, fish.Env, scale=T)
summary(RDA1.RDA4.vp)
plot(RDA1.RDA4.vp)
showvarparts(2)
RDA1.RDA4.vp
#test significance of MaxTL without Invert: X2barX1
Env.bar.Invert =  rda(ZSS.fish ~ ., data=fish.Env)
anova(Env.bar.Invert)
#test significances of Invert without MaxTL: X1barX2
Invert.bar.Env =  rda(ZSS.fish ~ ., data=fish.Invert)
anova(Invert.bar.Env)

#RDA2 vs RDA3
RDA2.RDA3.vp<-varpart(ZSS.fish, fish.MaxTL, fish.Invert.MaxTL, scale=T)
summary(RDA2.RDA3.vp)
plot(RDA2.RDA3.vp)
showvarparts(2)
RDA2.RDA3.vp
#test significance of MaxTL without Invert: X2barX1
Invert.MaxTL.bar.MaxTL =  rda(ZSS.fish ~ ., data=fish.Invert.MaxTL)
anova(Invert.MaxTL.bar.MaxTL)
#test significances of Invert without MaxTL: X1barX2
MaxTL.bar.Invert.MaxTL =  rda(ZSS.fish ~ ., data=fish.MaxTL)
anova(MaxTL.bar.Invert.MaxTL)

#RDA2 vs RDA4
RDA2.RDA4.vp<-varpart(ZSS.fish, fish.MaxTL, fish.Env, scale=T)
summary(RDA2.RDA4.vp)
plot(RDA2.RDA4.vp)
showvarparts(2)
RDA2.RDA4.vp
#test significance of MaxTL without Invert: X2barX1
Env.bar.MaxTL =  rda(ZSS.fish ~ ., data=fish.Env)
anova(Env.bar.MaxTL)
#test significances of Invert without MaxTL: X1barX2
MaxTL.bar.Env =  rda(ZSS.fish ~ ., data=fish.MaxTL)
anova(MaxTL.bar.Env)

#RDA3 vs RDA4
RDA3.RDA4.vp<-varpart(ZSS.fish, fish.Invert.MaxTL, fish.Env, scale=T)
summary(RDA3.RDA4.vp)
plot(RDA3.RDA4.vp)
showvarparts(2)
RDA3.RDA4.vp
#test significance of Invert.MaxTL without Env: X2barX1
Invert.MaxTL.bar.Env =  rda(ZSS.fish ~ ., data=fish.Env)
anova(Invert.MaxTL.bar.Env)
#test significances of Env without Invert.MaxTL: X1barX2
Env.bar.Invert.MaxTL =  rda(ZSS.fish ~ ., data=fish.Invert.MaxTL)
anova(Env.bar.Invert.MaxTL)


# --------------------------------------------------#

###HOW MUCH VARIATION CAN WE EXPLAIN###
#RDA5 with Env, Invert
# includes only Depth, INVERT
fish.Env.Invert<-fish.EnvMaxTL[c(3,5:8)]
RDA5 <- rda(ZSS.fish ~., fish.Env.Invert, scale=T)
#test significance of the canonical relationship
anova(RDA5, by = "axis")
#test the significance of each environmental variable
anova(RDA5, by = "margin")
RsquareAdj(RDA5)

#RDA5 triplot
plot(RDA2, scaling = 3, main="RDA5")
spe.sc = scores(RDA5, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

#RDA6 with Env, MaxTL
fish.Env.MaxTL<-fish.EnvMaxTL[c(3,5,7,8,10)]
RDA6 <- rda(ZSS.fish ~., fish.Env.MaxTL, scale=T)
#test significance of the canonical relationship
anova(RDA6, by = "axis")
#test the significance of each environmental variable
anova(RDA6, by = "margin")
RsquareAdj(RDA6)

#RDA6 triplot
plot(RDA6, scaling = 3, main="RDA6")
spe.sc = scores(RDA6, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

#RDA7 with Env, Invert + MaxTL
fish.Env.Invert.MaxTL<-fish.EnvMaxTL[c(3,5:8,10)]
RDA7 <- rda(ZSS.fish ~., fish.Env.Invert.MaxTL, scale=T)
#test significance of the canonical relationship
anova(RDA7, by = "axis")
#test the significance of each environmental variable
anova(RDA7, by = "margin")
RsquareAdj(RDA7)

#RDA7 triplot
plot(RDA7, scaling = 3, main="RDA7")
spe.sc = scores(RDA7, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(ZSS), col = "red", cex=0.8)

# --------------------------------------------------#

##Var Partitioning to compare all 3 RDAS
#RDA5 vs RDA6
RDA5.RDA6.vp<-varpart(ZSS.fish, fish.Env.Invert, fish.Env.MaxTL, scale=T)
summary(RDA5.RDA6.vp)
plot(RDA5.RDA6.vp)
showvarparts(2)
RDA5.RDA6.vp

#RDA5 vs RDA7
RDA5.RDA7.vp<-varpart(ZSS.fish, fish.Env.Invert, fish.Env.Invert.MaxTL, scale=T)
summary(RDA5.RDA7.vp)
plot(RDA5.RDA7.vp)
showvarparts(2)
RDA5.RDA7.vp

#RDA6 vs RDA7
RDA6.RDA7.vp<-varpart(ZSS.fish, fish.Env.MaxTL, fish.Env.Invert.MaxTL, scale=T)
summary(RDA6.RDA7.vp)
plot(RDA6.RDA7.vp)
showvarparts(2)
RDA6.RDA7.vp

# --------------------------------------------------#

DM = read.csv("DepthMaxTL.csv", header=TRUE)
attach(DM)
DM
stem(DEPTHL)
stem(MaxTL)
plot(DEPTHL,MaxTL)
model = lm(MaxTL ~ DEPTHL)
model
summary(model)
abline(lm(MaxTL ~ DEPTHL))
