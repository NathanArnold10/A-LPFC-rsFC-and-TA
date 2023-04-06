#Hierarchical regression in R
#Amygdala-LPFC rsFC does not predict TASC above MR abilities.

#Importing the data
Hierarchical <- read.csv("NA.csv")
print(Hierarchical)
Hierarchical <- Hierarchical[!apply(is.na(Hierarchical), 1, any),]



#Loading packages
install.packages("rmarkdown")
library(rmarkdown)
library("psych")
exists("psych")
library("lsr")
exists("lsr")
library("car")
exists("car")
library("MASS")
exists("car")
library("kader")
exists("kader")
library("ggplot2")

#Descriptive statistics
summary(Hierarchical)
describe(Hierarchical)
iqr_values <- sapply(Hierarchical, IQR, na.rm = T)
iqr_values

#Remove outliers
#Means and SD for z-scores
meanTASC <- mean(Hierarchical$TASC, na.rm = T)
sdTASC <- sd(Hierarchical$TASC, na.rm = T)
meanMR <- mean(Hierarchical$MatrixReasoning, na.rm = T)
sdMR <- sd(Hierarchical$MatrixReasoning, na.rm = T)
meanLPFC <- mean(Hierarchical$Amygdala_LPFC, na.rm = T)
sdLPFC <- sd(Hierarchical$Amygdala_LPFC, na.rm = T)
meanGP <- mean(Hierarchical$Pallidum.LPFC, na.rm = T)
sdGP <- sd(Hierarchical$Pallidum.LPFC, na.rm = T)
meanJLC <- mean(Hierarchical$JLC.Amygdala, na.rm = T)
sdJLC <- sd(Hierarchical$JLC.Amygdala, na.rm = T)
#Calculating Z-scores
z_scoresTASC <- ((Hierarchical$TASC - meanTASC) / sdTASC)
z_scoresTASC
z_scoresMR <- ((Hierarchical$MatrixReasoning - meanMR) / sdMR)
z_scoresMR
z_scoresLPFC <- ((Hierarchical$Amygdala_LPFC - meanLPFC) / sdLPFC)
z_scoresLPFC
Z_scoresGP <- ((Hierarchical$Pallidum.LPFC - meanGP) / sdGP)
Z_scoresGP
z_scoresJLC <- ((Hierarchical$JLC.Amygdala - meanJLC) / sdJLC)
z_scoresJLC
#MR Rows 1 and 74 have z_scores -3.1066456 and -3.6960761
Hierarchical <- Hierarchical[-c(1, 74),]
Hierarchical

#Testing Normality of Distribution
#Normality of Distribution
#TASC is Not Normally Distributed
shapiro.test(Hierarchical$TASC)
#MatrixReasoning is Not Normally Distributed
shapiro.test(Hierarchical$MatrixReasoning)
#Amygdala_LPFC is Normally distributed
shapiro.test(Hierarchical$Amygdala_LPFC)
#JPL.Amygdala is Not Normally Distributed
shapiro.test(Hierarchical$JLC.Amygdala)
#Pallidum.LPFC is Normally Distributed
shapiro.test(Hierarchical$Pallidum.LPFC)
#Test +ive
#MatrixReasoning Transformation of MR - BoxCox
MR1 <- boxcox(lm(Hierarchical$MatrixReasoning ~ 1))
MRnd <- (Hierarchical$MatrixReasoning^3)
shapiro.test(MRnd) #Normally Distributed
TASC1 <- boxcox(lm(Hierarchical$TASC ~ 1))
TASCnd <- sqrt(Hierarchical$TASC)
shapiro.test(TASCnd) #Not Normally Distributed
#Test -ive/+ive
#Transformation of TASC and JLC - Yeo-Johnson
JLCtrans <- yjPower(Hierarchical$JLC.Amygdala, lambda = -1)
shapiro.test(JLCtrans)
#All variables are normally distributed, or transformed
#to normality - with the exception of TASC

#Must Use Spearmans Rank and Kruskal-Wallis 
#Instead to Complete Non-Parametric hierarchical
#Regression

#Spearman's Rank Regressions
#Restricted Model
SR1 <- cor.test(Hierarchical$TASC, Hierarchical$MatrixReasoning, method = "spearman", exact = F)
#Unrestricted Amy-LPF model
SR2 <- cor.test(Hierarchical$TASC, Hierarchical$MatrixReasoning + Hierarchical$Amygdala_LPFC, method = "spearman", exact = F)
#Unrestricted GP-LPFC model
SR3 <- cor.test(Hierarchical$TASC, Hierarchical$MatrixReasoning + Hierarchical$Pallidum.LPFC, method = "spearman", exact = F)
#Unrestricted Amy-JPL model
SR4 <- cor.test(Hierarchical$TASC, Hierarchical$MatrixReasoning + Hierarchical$JLC.Amygdala, method = "spearman", exact = F)
PH5 <- cor.test(Hierarchical$TASC, Hierarchical$Amygdala_LPFC, method = "spearman", exact = F) #Post-Hoc
SR1 #Sig
SR1df <- data.frame(MatrixReasoning = Hierarchical$MatrixReasoning,
                 TASC = Hierarchical$TASC)
ggplot(SR1df, aes(x = MatrixReasoning, y = TASC)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
SR2 #Sig
SR3 #Sig
SR4 #Sig
PH5 #NS
PEARSON1 <- cor.test(Hierarchical$Amygdala_LPFC, MRnd, corr.method = "pearson")
PEARSON1 #MRnd and Amygdala_LPFC do not correlate

# Compare the Correlations Using Kruskal-Wallis
#Test comparison
kruskal.test(list(SR1 = SR1$estimate, SR2 = SR2$estimate)) #NS
#GP control comparison
kruskal.test(list(SR1 = SR1$estimate, SR3 = SR3$estimate)) #NS
#JPL control comparison
kruskal.test(list(SR1 = SR1$estimate, SR4 = SR4$estimate)) #NS

#All Models Were Statistically Significant, Unrestricted Models Were Not Significantly Better Predictors Than
#The Restricted Model, Suggesting Amygdala-LPFC does not explain variance beyond MR, and MR only predicts TASC.

