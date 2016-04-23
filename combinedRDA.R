# Script for running RDA and Variation Partitioning
# --------------------------------------------------

rm(list=ls(all=TRUE))
#set working directory

getwd()
setwd("/Users/ckluk/Desktop/Post doc/BythoRDA/Compare")

# Load required packages
library(ade4)
library(vegan)
library(recluster)
library(ape)

#Import info on number piscivores and planktivores vs invert#
Correlation <- read.csv("FishSdata.csv", row.name = 1)

# Import Hellinger transformed size class data (S1 to S4 only)
Biotic <- read.csv("HellBiomSizeClassesLong.csv", row.name = 1)
dim(Biotic)
#include only size classes 1 to 4
SizeBiotic<-Biotic[c(1:4)]
#change names of columns to get rid of Hell
colnames(SizeBiotic)=c("S1","S2","S3","S4")

# Import ENV data
Envdata <- read.csv("Envdatashortlong.csv", row.name = 1)
dim(Envdata)
# includes only CAL, TPL, SIO3, DOCL, pH, INVERT, DEPTHL, AREAL, BYTHO, BYTHO_P/A
# drop BYTHO_P/A for RDA
Envdatshort<-Envdata[c(1:9)]

# check correlations between environmental variables
# plot pairwise combinatins of environmental data
# pairs(Envdata, main = "Mean")

# For Variation partitioning, set up a Bytho matrix. 
# Best to use Bytho abundances
Bytho<-Envdata[c(10)]

#set up Env.nobytho
Env.no.bytho<-Envdata[c(1:8)]


# Redundancy analysis 
# ------------------------------------
## RDA
mod <- rda(SizeBiotic ~., Envdatshort, scale=T)
#test significance of the canonical relationship
anova(mod, by = "axis")
#test the significance of each environmental variable
anova(mod, by = "margin")
mod
RsquareAdj(mod)


#with variable selection
step.wise = ordistep(rda(SizeBiotic~ 1, Envdatshort, scale = T), scope = formula(mod), pstep = 1000, perm.max = 1000)

#get R.squared
RsquareAdj(step.wise)

#RDA of community composition by environmental variables chosen
#Change the columns below based on results of first step (sio3, area, invert, ph, depth)
rda.result = rda(SizeBiotic ~., Envdatshort[,c(3,5:8)], scale=T)
rda.result

#test significance of the canonical relationship
anova(rda.result, by = "axis")

#test the significance of each environmental variable
anova(rda.result, by = "margin")

#RDA triplot
plot(rda.result, scaling = 3, main="Triplot no-fish-RDA")
spe.sc = scores(rda.result, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col="red")
text(spe.sc[,1], spe.sc[,2], colnames(SizeBiotic), col = "red", cex=0.8)

RsquareAdj(rda.result)

# variation partitioning
# ------------------------------------
out.vp<-varpart(SizeBiotic, Bytho, Env.no.bytho, scale=T)
summary(out.vp)
plot(out.vp)
showvarparts(2)
out.vp

#test significance of Env without Bytho: X2barX1
rda.env =  rda(SizeBiotic ~ ., data=Env.no.bytho)
anova(rda.env)

#test significances of Bytho without Env: X1 bar X2
rda.bytho =  rda(SizeBiotic ~ ., data=Bytho)
anova(rda.bytho)

# Script for running RDA and Variation Partitioning
# With fish P/A data
# --------------------------------------------------

# Import Hellinger transformed size class data (S1 to S4 only)
fish.Biotic <- read.csv("ZooSizeClasses_NONHell.csv", row.name = 1)
dim(fish.Biotic)
#include only size classes 1 to 4
fish.SzBiotic<-fish.Biotic[c(2:5)]
#transform size biomass data (Hellinger)
fish.SizeBiotic<-decostand(fish.SzBiotic,"hellinger")

# Import ENV data
fish.Envdata <- read.csv("EnvFishlakes.csv", row.name = 1)
dim(fish.Envdata)
# includes only INVERT, CAL, TPL, SIO3, DOCL, pH, DEPTHL, AREAL, BYTHO, BYTHO_P/A
# drop BYTHO_P/A for RDA
fish.Envdatshort<-fish.Envdata[c(7:15)]

# Import Fish data
Fish <- read.csv("FishPAdata.csv", row.name = 1)
Fish <- Fish[,-1] #remove 1st col of CAISNID
dim(Fish)

# Make a big matrix of fish and Env for RDA
fish.env<-cbind(fish.Envdatshort,Fish)

# Check some correlations
cor.p = cor(fish.env$TPL, fish.env$DEPTHL)
plot(fish.env$TPL, fish.env$DEPTHL)
abline(lm(fish.env$TPL~fish.env$DEPTHL))

# ------------------------------------# Redundancy analysis # ------------------------------------

# Stepwise RDA without fish (just environment) on smaller fish-length zoo-env datasets
modnofish <- rda(fish.SizeBiotic ~., fish.Envdatshort, scale=T)
#test significance of the canonical relationship
anova(modnofish, by = "axis")
#test the significance of each environmental variable
anova(modnofish, by = "margin")


## RDA with variable selection
## Biotic matrix is: SizeBiotic
## Explan matrix is: fish.env (includes environment + fish P/A)

fish.mod <- rda(fish.SizeBiotic ~., fish.env, scale=T)
fish.step.wise = ordistep(rda(fish.SizeBiotic~ 1, fish.env, scale = T), scope = formula(fish.mod), pstep = 1000, perm.max = 1000)
#test significance of the canonical relationship
anova(fish.mod, by = "axis")
#test the significance of each environmental variable
anova(fish.mod, by = "margin")
RsquareAdj(fish.mod)

# Scaling = 3 (both species and site scores are scaled by eigenvalues)
plot(fish.mod, scaling = 3, main="Triplot fish-RDA: scaling = 3") #includes all variables
fish.spe.sc3 = scores(fish.mod, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,fish.spe.sc3[,1], fish.spe.sc3[,2], length=0, lty=1, col="red")
text(fish.spe.sc3[,1], fish.spe.sc3[,2], colnames(fish.SizeBiotic), col = "red", cex=0.8)

#-----------------------------------
# using only significant variates from no fish (sio3, area, invert, ph, depth) and fish
# to get pretty graphs
# variates: SizeBiotic ~  DEPTHL + Invert + Area + PH + sio3 + all fish

fish.mod2 = rda(fish.SizeBiotic ~., fish.env[,c(1,5,6,8:20)], scale=T)
fish.mod2
#test significance of the canonical relationship
anova(fish.mod2, by = "axis")
#test the significance of each environmental variable
anova(fish.mod2, by = "margin")
RsquareAdj(fish.mod2)

# Scaling = 3 (both species and site scores are scaled by eigenvalues)
plot(fish.mod2, scaling = 3, main="Triplot fish-RDA signif vars only: scaling = 3")
fish.spe.sc3 = scores(fish.mod2, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,fish.spe.sc3[,1], fish.spe.sc3[,2], length=0, lty=1, col="red")
text(fish.spe.sc3[,1], fish.spe.sc3[,2], colnames(fish.SizeBiotic), col = "red", cex=0.8)

#-----------------------------------
# RDA STEPWISE
fish.mod4 <- rda(fish.SizeBiotic ~., fish.env, scale=T)
#test significance of the canonical relationship
anova(fish.mod4, by = "axis")
#test the significance of each environmental variable
anova(fish.mod4, by = "margin")
fish.mod4
RsquareAdj(fish.mod4)


#with variable selection
step.wise2 = ordistep(rda(fish.SizeBiotic~ 1,fish.env, scale = T), scope = formula(fish.mod4), pstep = 1000, perm.max = 1000)

#get R.squared
RsquareAdj(step.wise2)

# using only significant variates (depth, cisco, rainbow smelt)
# to get pretty graphs
# variates: SizeBiotic ~ DEPTHL + cisco + rainbowsmelt

fish.mod3 = rda(fish.SizeBiotic ~., fish.env[,c(8,13,14)], scale=T)
fish.mod3
#test significance of the canonical relationship
anova(fish.mod3, by = "axis")
#test the significance of each environmental variable
anova(fish.mod3, by = "margin")
RsquareAdj(fish.mod3)

# Scaling = 3 (both species and site scores are scaled by eigenvalues)
plot(fish.mod3, scaling = 3, main="Triplot fish-RDA signif vars only: scaling = 3") #includes all variables
fish.spe.sc3 = scores(fish.mod3, choices=1:2, scaling = 3, display = "sp")
arrows(0, 0,fish.spe.sc3[,1], fish.spe.sc3[,2], length=0, lty=1, col="red")
text(fish.spe.sc3[,1], fish.spe.sc3[,2], colnames(fish.SizeBiotic), col = "red", cex=0.8)

# variation partitioning
# ------------------------------------
fish.out.vp<-varpart(fish.SizeBiotic, Fish, fish.Envdatshort, scale=T)
summary(fish.out.vp)
plot(fish.out.vp)
showvarparts(2)
fish.out.vp

#test significance of Env without Fish: X2barX1
fish.rda.env =  rda(fish.SizeBiotic ~ ., data=Envdatshort)
anova(fish.rda.env)

#test significances of Fish without Env: X1 bar X2
fish.rda.bytho =  rda(fish.SizeBiotic ~ ., data=Fish)
anova(fish.rda.bytho)

#PROCRUSTES ROTATION#
#Perform the procrustes and plot the matrices
procr1<-recluster.procrustes(rda.result,modnofish,num=156)
procr1

rda.proc <- procrustes (modnofish, fish.mod)
rda.proc
plot (rda.proc)
protest (modnofish, fish.mod)

rda.proc2 <- procrustes (fish.mod, fish.mod2)
rda.proc2
plot (rda.proc2)
protest (fish.mod, fish.mod2)

rda.proc3 <- procrustes (fish.mod2, fish.mod3)
rda.proc3
plot (rda.proc3)
protest (fish.mod2, fish.mod3)

#CORRELATION BETWEEN INVERT AND PLANKTOVIRES/PISC#

plank <- cor.test(Correlation$INVERT, Correlation$PlanktivS)
plank

pisc <- cor.test(Correlation$INVERT, Correlation$PiscivS)
pisc


