# 26/05/2020
# patricia.sebastian@ivirma.com


# Settings ----------------------------------------------------------------

setwd("/home/patricia/IVI/POF/")

library(Rmisc)
library(readxl)  # read excel sheets


# Read data ---------------------------------------------------------------

datos = read.delim("150 Muestras FOP FOO y CONTROLES_anonimizado.csv", header = T, sep = "\t", as.is = T)
datos = datos[, -c(14, 15)]
str(datos)

table(datos$TIPO.FO, useNA = "ifany")
# CONTROL     FOO     FOP 
#      31      83      35 

# Group case and control
datos$GRUPO = ifelse(datos$TIPO.FO == "CONTROL", "CONTROL", "CASE")
table(datos$GRUPO)
# CASE CONTROL 
#  118      31 

# Analyze variables -------------------------------------------------------

#### EDAD ####
summarySE(data = datos, measurevar = "EDAD", groupvars = "GRUPO", na.rm = T)
# GRUPO   N     EDAD       sd        se        ci
# 1    CASE 118 35.68644 3.398590 0.3128656 0.6196138
# 2 CONTROL  31 35.06452 2.803991 0.5036117 1.0285122

# Test normalidad
shapiro.test(datos$EDAD) 
# p-value = 1.648e-05 -> No normal data

# Evaluar las diferencias
wilcox.test(EDAD ~ GRUPO, datos)
# p-value = 0.2417

####### Si incluimos los tres grupos
# summarySE(data = datos, measurevar = "EDAD", groupvars = "TIPO.FO", na.rm = T)
# kruskal.test(TIPO.FO ~ EDAD, datos) # si sale significativo --> pairwise
# pairwise.wilcox.test(datos$EDAD, datos$TIPO.FO)

#### IMC ####
summarySE(data = datos, measurevar = "IMC", groupvars = "GRUPO", na.rm = T)
shapiro.test(datos$IMC) 
# p-value = 3.375e-08  -> No normal data
wilcox.test(IMC ~ GRUPO, datos)
# p-value = 0.1859
####### Si incluimos los tres grupos --> Como en edad

#### Origen ####
table(datos$origen, datos$GRUPO)
# CASE CONTROL
# bilbao          29       0
# la fe            8       0
# pamplona        59      31
# san sebastian    2       0
# vitoria         20       0

# Al ser variable discreta la evaluamos con test de Fisher
fisher.test(table(datos$origen, datos$GRUPO))
# p-value = 3.62e-06

##### Menarquia ####
summarySE(datos, measurevar = "MENARQUIA", groupvars = "GRUPO", na.rm = T)
shapiro.test(datos$MENARQUIA)
# p-value < 2.2e-16 -> No normal data
wilcox.test(MENARQUIA ~ GRUPO, datos)
# p-value = 0.6173
####### Si incluimos los tres grupos --> Como en edad

#### RNV ####
summarySE(datos, measurevar = "RNV", groupvars = "GRUPO", na.rm = T)
shapiro.test(datos$RNV)
# p-value = 1.171e-13 -> No normal data
wilcox.test(RNV ~ GRUPO, datos)
# p-value = 0.001855
####### Si incluimos los tres grupos --> Como en edad

#### RFA ####
summarySE(data = datos, measurevar = "RFA", groupvars = "GRUPO", na.rm = T)
# GRUPO   N       RFA       sd        se        ci
# 1    CASE 117  3.649573 3.190277 0.2949412 0.5841682
# 2 CONTROL  30 12.766667 3.410767 0.6227181 1.2736014
shapiro.test(datos$RFA)
# p-value = 2.265e-08 -> No normal data
wilcox.test(RFA ~ GRUPO, datos)
# p-value = 8.542e-16 
####### Si incluimos los tres grupos --> Como en edad

#### FSH ####
datos$FSH = as.numeric(gsub("<", "", gsub(",", ".", datos$FSH))) # Ñapa
summarySE(datos, measurevar = "FSH", groupvars = "GRUPO", na.rm = T)
# GRUPO   N       FSH       sd        se       ci
# 1    CASE 112 28.911429 32.92015 3.1106617 6.163984
# 2 CONTROL  29  5.723103  3.13316 0.5818133 1.191790
shapiro.test(datos$FSH)
# p-value = 7.75e-16 -> No normal data
wilcox.test(FSH ~ GRUPO, datos)
# p-value = 1.658e-08
####### Si incluimos los tres grupos --> Como en edad

#### AMH ####
datos$AMH = as.numeric(gsub("<", "", gsub(",", ".", datos$AMH))) # Ñapa
summarySE(datos, measurevar = "AMH", groupvars = "GRUPO", na.rm = T)
# GRUPO   N       AMH        sd         se         ci
# 1    CASE 115 0.2830696 0.2553781 0.02381413 0.04717561
# 2 CONTROL  30 3.6220000 2.7693662 0.50561477 1.03409832
shapiro.test(datos$AMH)
# p-value < 2.2e-16 -> No normal data
wilcox.test(AMH ~ GRUPO, datos)
# p-value < 2.2e-16   
####### Si incluimos los tres grupos --> Como en edad

#### AMH grouped ####
datos$AMH_grouped = ifelse(datos$AMH<0.68, "Baja", 
                           ifelse(datos$AMH>2.7, "normal", "suficiente"))
table(datos$AMH_grouped, datos$GRUPO)
fisher.test(table(datos$AMH_grouped, datos$GRUPO))
# p-value < 2.2e-16

#### E2 ####
datos$E2 = as.numeric(gsub("<", "", gsub(",", ".", datos$E2))) # Ñapa
summarySE(datos, measurevar = "E2", groupvars = "GRUPO", na.rm = T)
# GRUPO  N       E2       sd        se       ci
# 1    CASE 56 59.41518 51.71251  6.910374 13.84870
# 2 CONTROL  7 62.14286 31.42944 11.879211 29.06738
shapiro.test(datos$E2)
# p-value = 3.762e-06 -> No normal data
wilcox.test(E2 ~ GRUPO, datos)
# p-value = 0.4772


# for A vs B --------------------------------------------------------------

setwd("workspace/projects/covid19/")
load("../../matrix.Rda")

datos = read_xlsx("../pof/pof demographic data.xlsx", sheet = 2 )
datos = datos[, -c(14)]
str(datos)
datos$FSH <- as.numeric(datos$FSH)
datos$FSH
datos$RFA <- as.numeric(datos$RFA)
datos$AMH <- as.numeric(datos$AMH)
datos <-
  as.data.frame(sapply(mat_df, as.numeric))

table(datos$`TIPO FO`, useNA = "ifany")
# FOO FOP 
# 83  35 
# CONTROL     FOO     FOP 
#      31      83      35 

# Group case and control
datos$GRUPO = ifelse(datos$Subtipo == "A", "B")
datos$GRUPO = datos$Subtipo
table(datos$GRUPO)
# A   B 
# 17 101 

#### RFA ####
summarySE(data = datos, measurevar = "RFA", groupvars = "GRUPO", na.rm = T)
# GRUPO   N      RFA       sd        se        ci
# 1     A  17 4.117647 3.314407 0.8038618 1.7041109
# 2     B 100 3.570000 3.179019 0.3179019 0.6307863
shapiro.test(datos$RFA)
# p-value = 6.156e-07 -> No normal data  # (low pvalue)
wilcox.test(RFA ~ GRUPO, datos)
# W = 946.5, p-value = 0.4544


#### FSH ####
datos$FSH = as.numeric(gsub("<", "", gsub(",", ".", datos$FSH))) # Ñapa
summarySE(datos, measurevar = "FSH", groupvars = "GRUPO", na.rm = T)
# GRUPO  N      FSH       sd        se        ci
# 1     A 15 42.49400 43.46888 11.223616 24.072262
# 2     B 97 26.81103 30.72829  3.119985  6.193122
shapiro.test(datos$FSH)
# p-value = 9.576e-13 -> No normal data
wilcox.test(FSH ~ GRUPO, datos)
# W = 891.5, p-value = 0.1625


#### AMH ####
datos$AMH = as.numeric(gsub("<", "", gsub(",", ".", datos$AMH))) # Ñapa
summarySE(datos, measurevar = "AMH", groupvars = "GRUPO", na.rm = T)
# GRUPO  N       AMH        sd         se         ci
# 1     A 17 0.2694118 0.2484822 0.06026578 0.12775774
# 2     B 94 0.2941809 0.2590850 0.02672258 0.05306574
shapiro.test(datos$AMH)
# p-value < p-value = 2.985e-09 -> No normal data
wilcox.test(AMH ~ GRUPO, datos)
# p-value = 0.6371
####### Si incluimos los tres grupos --> Como en edad

df <- datos[ , c("ID","AMH","FSH","RFA")]
df[all, ]

rownames(df) <- NULL

df

fit$labels
correct <- fit$labels[all]
correct

rownames(datos) <- datos$ID
datos[correct,]
datos
df <- datos[correct , c("ID","AMH","FSH","RFA")]
t(df)

write.table(x = t(df), "../pof/clinical_data_dendrogram_order.tsv", sep = "\t", quote = F,
            col.names = T, row.names = F)
