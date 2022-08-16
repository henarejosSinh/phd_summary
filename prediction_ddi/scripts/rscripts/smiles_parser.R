
# Parser test -------------------------------------------------------------

## Script para parsear SMILES y extraer matrices de similitud (por Tanimoto) y de distancias entre fármacos
rm(list = ls())
setwd("/home/pgarcia/Escritorio/")
install.packages("rcdk")
library(fingerprint)
library(gplots)
library(ggplot2)
library(rcdk)
cdk.version()
#bpdata #Example of data.frame to parser
# https://cran.r-project.org/web/packages/rcdk/vignettes/using-rcdk.html (scripts de ejemplo)

## Data.frame con rownames = nombre del fármaco y SMILES
drugbank = read.delim("../../remote/drugbank/2020/2020_07_02/drugbank_info_2020_07.txt", header = T, sep = "\t", as.is = T)
drugbank = drugbank[grep("\\bapproved\\b", drugbank$groups),]
query_smiles = data.frame(ID = drugbank$X.ID, SMILES = drugbank$smiles, stringsAsFactors = F)
query_smiles = query_smiles[query_smiles$SMILES != "",] ## Eliminamos los que no tienen SMILE
#query_smiles$SMILES[2022] ## NULL no lo lee bien al hacer el parseo con la función parse.smiles. Las demás sí
query_smiles = query_smiles[query_smiles$SMILES != "OC1=CC=CC(=C1)C-1=C2\\CCC(=N2)\\C(=C2/N\\C(\\C=C2)=C(/C2=N/C(/C=C2)=C(\\C2=CC=C\\-1N2)C1=CC(O)=CC=C1)C1=CC(O)=CC=C1)\\C1=CC(O)=CC=C1",] ## Eliminamos este fármaco porque no reconoce el SMILE (fila 2022) y peta al parsear
query_smiles = data.frame(SMILES = query_smiles$SMILES, stringsAsFactors = F, row.names = query_smiles$ID)

head(bpdata)
head(query_smiles)

mols <- parse.smiles(query_smiles$SMILES, omit.nulls = T) ## Otra opción para no tener que cargarnos la fila 2022 en el paso anterior
# mols <- parse.smiles(bpdata$SMILES) ## Ejemplo
mols <- parse.smiles(query_smiles$SMILES) ## Con nuestros datos, función para convertirlo a objeto que nos pide el paquete
mols
fps <- lapply(mols, get.fingerprint, type='circular') ## Cálculo de bits a partir de SMILES
str(fps)
head(fps)

## Calculamos las matrices de similitud y la de distancias (la inversa a la de similitud)
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto') ## Matriz de similitud
fp.sim <- fingerprint::fp.sim.matrix(t_m2, method='tanimoto') ## Matriz de similitud
rownames(fp.sim) = rownames(query_smiles)
colnames(fp.sim) = rownames(query_smiles)
head(fp.sim)
View(fp.sim) ## Aunque parezca que no está haciendo la matriz por el head(), al hacer el View, vemos que sí la hace
feature(t_m2)
fp.dist <- 1 - fp.sim ## Matriz de distancias
rownames(fp.dist) = rownames(query_smiles)
colnames(fp.dist) = rownames(query_smiles)
head(fp.dist)
View(fp.dist) ## Lo mismo ocurre con la matriz de distancias

## Podemos filtrar por una similitud mayor que X para quedarnos con las que más se parezcan entre sí
