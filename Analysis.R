
##########################################################################
## Impact of life cycle variation on feeding system musculature in Caudata
## Morgane Taillades
## Statistical analysis
##########################################################################


### Libraries 
library(readxl)      # read_excel
library(ape)         # read.nexus
library(mvMORPH)     # mvgls, mvgls.pca, manova.gls, GIC
library(nlme)        # gls, correlation structures
library(car)         # Anova (Type II)
library(dplyr)       # light wrangling

### Working directory 
setwd("C:/Users/mtaillad/OneDrive - UGent/UGent-PC/Mtaillad/Desktop/Salamandres/Redaction/Taillades et al. 2026/MT/V2/Datas and analysis/")

### Data input 
RawData <- read_excel("Mean_by_species.xlsx", col_names = TRUE)

###  Build the table with functional groups

df <- data.frame(
  RawData[, c(2:11)],
  # Sum of muscle masses by functional group, then convert mass to volume
  Vol_HyoLev   = ((RawData$pds_MQP + RawData$pds_MIM + RawData$pds_MIH + RawData$pds_MIM_MIH + RawData$pds_MIHp + RawData$pds_MGH)* 1e-3) / 1.06,
  Vol_HyoPro   = ((RawData$pds_MGH + RawData$pds_MSbR + RawData$pds_MSbH + RawData$pds_MBh + RawData$pds_Mbhint + RawData$pds_MGG)* 1e-3) / 1.06,
  Vol_HyoRet   = ((RawData$pds_MRC + RawData$pds_MRCpro + RawData$pds_MRCsup)* 1e-3) / 1.06,
  Vol_MouthOp  = ((RawData$pds_MDM + RawData$pds_MDMa + RawData$pds_MDMp)* 1e-3) / 1.06,
  Vol_MouthClo = ((RawData$pds_MAME + RawData$pds_MAMEsup + RawData$pds_MAMEpro + RawData$pds_MAMIproI + RawData$pds_MAMIproII + RawData$pds_MAMIsup + RawData$pds_MAMP + RawData$pds_MPt)* 1e-3) / 1.06,
  
  # Sum of PCSAs by functional group
  pcsa_HyoLev   = (RawData$pcsa_MQP + RawData$pcsa_MIM + RawData$pcsa_MIH + RawData$pcsa_MIM_MIH + RawData$pcsa_MIHp + RawData$pcsa_MGH),
  pcsa_HyoPro   = (RawData$pcsa_MGH + RawData$pcsa_MSbR + RawData$pcsa_MSbH + RawData$pcsa_MBh + RawData$pcsa_Mbhint + RawData$pcsa_MGG),
  pcsa_HyoRet   = (RawData$pcsa_MRC + RawData$pcsa_MRCpro + RawData$pcsa_MRCsup),
  pcsa_MouthOp  = (RawData$pcsa_MDM + RawData$pcsa_MDMa + RawData$pcsa_MDMp), 
  pcsa_MouthClo = (RawData$pcsa_MAME + RawData$pcsa_MAMEsup + RawData$pcsa_MAMEpro + RawData$pcsa_MAMIproI + RawData$pcsa_MAMIproII + RawData$pcsa_MAMIsup + RawData$pcsa_MAMP + RawData$pcsa_MPt))

### Import phylogeny 
phy <- read.nexus(file = "PhySalamanders.nex")

# Reorder df to match the phylogeny tip order
df <- df[unlist(sapply(unique(phy$tip.label), function(x) which(df$Espece == x))), ]


# Factors
df$T_eco    <- as.factor(df$T_eco)
df$T_morpho     <- as.factor(df$T_morpho)
df$Adult_habitat<- as.factor(df$Adult_habitat)

###Log transformation

df <- data.frame(
  df[, 1:5],
  log10(df[, 6:20]))

# Row names as species
rownames(df) <- df$Espece



# Variable groups
variablesVol  <- c("Vol_HyoLev","Vol_HyoPro","Vol_HyoRet","Vol_MouthOp","Vol_MouthClo")
variablespcsa <- c("pcsa_HyoLev","pcsa_HyoPro","pcsa_HyoRet","pcsa_MouthOp","pcsa_MouthClo")


# IMPORTANT: Uncomment the line corresponding to the variable you want to use for the analyses.

variablesEco <- "T_eco"
#variablesEco <- c("T_morpho")
#variablesEco <- c("Adult_habitat")

####------------------------------------------------------------------------------####
####                                 VOLUME                                       ####
####------------------------------------------------------------------------------####

#Build formula 
formVol <- as.formula(paste0("cbind(",paste(variablesVol, collapse = ","),") ~ ",variablesEco," * HL"))

### Multivariate phylogenetic models (mvgls) 
# Model set: BM, OU, EB, lambda
fitBMVol     <- mvgls(formVol, phy, data = df, model = "BM",     method = "LL")
fitOUVol     <- mvgls(formVol, phy, data = df, model = "OU",     method = "LL")
fitEBVOLVol  <- mvgls(formVol, phy, data = df, model = "EB",     method = "LL")
fitLAMBDAVol <- mvgls(formVol, phy, data = df, model = "lambda", method = "LL")

# Compare models using GIC (smaller is better)
GIC(fitBMVol) 
GIC(fitOUVol) 
GIC(fitEBVOLVol) 
GIC(fitLAMBDAVol) 

#IMPORTANT: Choose the evolutionary model with the lowest GIC, 
#           and modify "model = 'X'" in the formula below accordingly.

### Phylogenetic PCA
fitACPVol <- mvgls(cbind(Vol_HyoLev,Vol_HyoPro,Vol_HyoRet,Vol_MouthOp,Vol_MouthClo) ~ HL, phy, data = df, model = "OU", method = "LL")
ACP_Vol   <- mvgls.pca(fitACPVol, plot = TRUE)

# Variance explained per axis (%)
print((ACP_Vol$values / sum(ACP_Vol$values)) * 100)

# Variable loadings per axis
PCA_loadings_Vol <- ACP_Vol$vectors
rownames(PCA_loadings_Vol) <- c("Vol_HyoLev","Vol_HyoPro","Vol_HyoRet","Vol_MouthOp","Vol_MouthClo")
print(PCA_loadings_Vol)

# Components for plotting
PCA_components_Vol <- data.frame(ACP_Vol$scores[, 1:4], df[[variablesEco]])
colnames(PCA_components_Vol)[5] <- "Eco"
PCA_components_Vol$Eco<- as.factor(PCA_components_Vol$Eco)

#IMPORTANT: Choose the evolutionary model with the lowest GIC, 
#           and modify fit'X'Vol in the formula below accordingly.

### MANCOVA (Wilks)
man_Vol <- manova.gls(fitOUVol, test = "Wilks", nperm = 999, verbose = TRUE, type = "II")
print(man_Vol)
print(effectsize(man_Vol))

#IMPORTANT: In the gls function, choose the correlation structure that most closely matches 
#           the evolutionary model previously determined.

### ANCOVA-like (univariate PGLS per variable) 
DFAncVol <- data.frame(df, species = row.names(df))
ancova_Vol <- list()

formVol <- as.formula(paste0("cbind(",paste(variablesVol, collapse = ","),") ~ ",variablesEco," * HL"))


for (var in variablesVol) {
  # Univariate PGLS with Martins correlation (fixed)
  form <- as.formula(paste(var, "~", variablesEco, "+ HL"))
  model <- gls(form,
               correlation = corMartins(value = 1, phy = phy, fixed = TRUE, form = ~ species),
               data = DFAncVol, na.action = na.exclude)
  
  A <- Anova(model) # Type II
  B <- as.data.frame(A)
  ancova_Vol[[var]] <- B
}

print(ancova_Vol)

####------------------------------------------------------------------------------####
####                                   PCSA                                       ####
####------------------------------------------------------------------------------####

formpcsa <- as.formula(paste0("cbind(",paste(variablespcsa, collapse = ","),") ~ ",variablesEco," * HL"))

### Multivariate phylogenetic models (mvgls) 
# Model set: BM, OU, EB, lambda
fitBMpcsa     <- mvgls(formpcsa, phy, data = df, model = "BM",     method = "LL")
fitOUpcsa     <- mvgls(formpcsa, phy, data = df, model = "OU",     method = "LL")
fitEBpcsapcsa  <- mvgls(formpcsa, phy, data = df, model = "EB",     method = "LL")
fitLAMBDApcsa <- mvgls(formpcsa, phy, data = df, model = "lambda", method = "LL")

# Compare models using GIC (smaller is better)
GIC(fitBMpcsa) 
GIC(fitOUpcsa) 
GIC(fitEBpcsapcsa) 
GIC(fitLAMBDApcsa) 

#IMPORTANT: Choose the evolutionary model with the lowest GIC, 
#           and modify "model = 'X'" in the formula below accordingly.

### Phylogenetic PCA
fitACPpcsa <- mvgls(cbind(pcsa_HyoLev,pcsa_HyoPro,pcsa_HyoRet,pcsa_MouthOp,pcsa_MouthClo) ~ HL, phy, data = df, model = "OU", method = "LL")
ACP_pcsa   <- mvgls.pca(fitACPpcsa, plot = TRUE)

# Variance explained per axis (%)
print((ACP_pcsa$values / sum(ACP_pcsa$values)) * 100)

# Variable loadings per axis
PCA_loadings_pcsa <- ACP_pcsa$vectors
rownames(PCA_loadings_pcsa) <- c("pcsa_HyoLev","pcsa_HyoPro","pcsa_HyoRet","pcsa_MouthOp","pcsa_MouthClo")
print(PCA_loadings_pcsa)

# Components for plotting
PCA_components_pcsa <- data.frame(ACP_pcsa$scores[, 1:4], df[[variablesEco]])
colnames(PCA_components_pcsa)[5] <- "Eco"
PCA_components_pcsa$Eco<- as.factor(PCA_components_pcsa$Eco)

#IMPORTANT: Choose the evolutionary model with the lowest GIC, 
#           and modify fit'X'pcsa in the formula below accordingly.

### MANCOVA (Wilks)
man_pcsa <- manova.gls(fitOUpcsa, test = "Wilks", nperm = 999, verbose = TRUE, type = "II")
print(man_pcsa)
print(effectsize(man_pcsa))

#IMPORTANT: In the gls function, choose the correlation structure that most closely matches 
#           the evolutionary model previously determined.

### ANCOVA-like (univariate PGLS per variable) 
DFAncpcsa <- data.frame(df, species = row.names(df))
ancova_pcsa <- list()

for (var in variablespcsa) {
  # Univariate PGLS with Martins correlation (fixed)
  form <- as.formula(paste(var, "~", variablesEco, "+ HL"))
  model <- gls(form,
               correlation = corMartins(value = 1, phy = phy, fixed = TRUE, form = ~ species),
               data = DFAncpcsa, na.action = na.exclude)
  
  A <- Anova(model) # Type II
  B <- as.data.frame(A)
  ancova_pcsa[[var]] <- B
}

print(ancova_pcsa)

