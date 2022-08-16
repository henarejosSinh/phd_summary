########### Description #######################################################
# 30 july 2020
# ihc.europa@gmail.com
# limma-voom for rna-seq COVID GEO 
# In vivo antiviral host response to SARS-CoV-2 by viral load, sex, and age [dataset I]
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075
# GSE152075

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

# specific use
library(NOISeq)
library(tidyverse)
library(edgeR)
library(limma)
library(gsrmUtils)
library(data.table)

# general use

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)  # for color palettes
library(clusterProfiler)  # for enrichment

display.brewer.all()

colours_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "Gray", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",  # light green
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


########### Enviroment #########################################################

setwd("workspace/projects/")
getwd()

# save.image("scripts/rdata/limma-voom_03_08_20.Rda")
# load("scripts/rdata/limma-voom_03_08_20.Rda")

# Create variable for working directory
dir <- "workspace/projects/"


# Functions to be used ----------------------------------------------------

# ggplot theme settings
bar.zyx.theme <-  theme_light( base_size = 18 ) +
  theme(
    axis.text.x = element_text(hjust = 0.65),
    axis.title.x = element_text(vjust = -1.0),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", 
                                    colour = NA, inherit.blank = T),
    plot.background = element_rect( colour = "white" , inherit.blank = T),
    panel.grid = element_line( colour = "gray88" , 
                               arrow = F, inherit.blank = T),
    panel.grid.minor = element_line( size = rel(0.25),
                                     arrow = F, inherit.blank = T),
    panel.border = element_rect( colour = "gray70" , size = rel(1), 
                                 inherit.blank = T),
    axis.line = element_blank() ) 


# load information from geo GSE152075 -------------------------------------------

load("covid19/scripts/rdata/annot_table_GSE150275.Rda", verbose = T)
head(ed)

# annot should be gene annotations
# experimental desing must have the same rownames as the rest of my data, 
# in same order, if there's a col named samples
rownames(ed) <- ed$title
# ed$sample <- NULL

# load count matrix
# rownames= genes, colnames = samples
counts <- read.delim("covid19/data/geo_series/GSE152075_raw_counts_GEO.txt", 
                     row.names = 1, sep = " ")
head(counts)

all(rownames(ed) == colnames(counts))
# if needed to reorder ed[colnames(counts),]

# create an object of DGEList 
dge <- DGEList(counts = counts)
dge$samples
dim(dge)

# we need to do an exp design mat
design <- model.matrix(~0+sars.cov.2.positivity, data = ed)
rownames(design) <- rownames(ed)
head(design)

 # filter low counts or count 0 
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = F]
dim(dge)

# alternative way to remove low counts
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
d <- dge[-drop,] 
dim(d)

# TMM Normalization
# Trimmed mean of M values (TMM) normalization estimates sequencing
# depth after excluding genes for which the ratio of counts between a pair
# of experiments is too extreme or for which the average expression is too
# extreme. The edgeR software implements a TMM normalization.

dge <- calcNormFactors(dge)
head(dge@.Data)

# # When the library sizes are quite variable between samples, then the voom approach 
# is theoreticallymore  powerful  than  limma-trend.   In  this  approach,  
# the  voom  transformation  is  applied to  thenormalized and filteredDGEListobject:
v <- voom(dge, design, plot = T)
v <- voom(dge, design, plot = T, normalize.method = "quantile")

# alternative
v <- voom(d, design, plot = T)
v <- voom(d, design, plot = T, normalize.method = "quantile")
# v <- voom(dge, design, plot = T, normalize.method = "quantile")

names(v) # E normalized values
head(v$E)

# repersent with boxplot
dat <- v$E
dat <- normalizeQuantiles( dat, ties = T)
plot_Boxplot(dat = dat, ed, condition = "sars.cov.2.positivity")

# the usual limma pipelines for differential expression can be applied, for example
fit <- lmFit(v, design = design)
fit <- eBayes(fit)
topTable(fit, coef = ncol(design))

# to give more weight to fold-changes in the ranking, one could use say:
fit <- treat(fit, lfc = log2(1.2))
topTreat(fit, coef = ncol(design))


# alternative using NOISEQ ------------------------------------------------
# https://rpubs.com/jberghout/525119
# https://bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
# load ed and counts
load("projects/covid19/scripts/rdata/annot_table_GSE150275.Rda", verbose = T)
head(ed)
rownames(ed) <- ed$title  # for downstream analysis, ed needs to have rownames
# = colnames of count matrix (normally, the samples from the study)

counts <- 
  read.delim("projects/covid19/data/geo_series/GSE152075_raw_counts_GEO.txt", 
                     row.names = 1, sep = " ")
head(counts)
# exclude bad quality samples (check preprocessing_RNAseq...R)
counts <- 
  counts[,!names(counts) %in% c("POS_074", "POS_231", "POS_233", "POS_320")]
ed <- ed[!ed$title %in% c("POS_074", "POS_231", "POS_233", "POS_320"), ]

design <- model.matrix(~0+sars.cov.2.positivity, data = ed)

# NOISEQ biomart fetch ----------------------------------------------------

# check data(marioni) to load objects use as example for the package
data(Marioni)

str(mybiotypes)
mychroms
str(mychroms)
mylength
str(mylength)  # mychroms start-end
mycounts  # ok 
myfactors  # ok from ed (matrix desing(~0+...))
mygc

# we need to obtain biological information (gene length, biotype...) for NOISEQ
# that means using biomart
# GSE used v96 of biomart for mapping
# extract gene names and start from there
# if it does not work properly, convert genenames to ENG using corres and them
# retrieve information from biomart
write.table("gse150275_genenames.txt", x = rownames(counts), quote = F,
            col.names = F, row.names = F)
length(rownames(counts))  # 35784 genes
mart <- fread("projects/covid19/data/geo_series/GSE152075_mart_export.txt")
length(mart$`Gene name`)  # 31747 genes, 131 anyDuplicated by gene name
head(mart)
mart2 <- mart[!duplicated(mart$`Gene name`) ,]
mart2[mart2$`Gene name` %in% rownames(counts), ]

# subset counts by mart info
# first, we create a new col with ensg ID
counts$ensg <- unlist(lapply(rownames(counts), function(name){
  print(name)
  if (!name %in% mart2$`Gene name`) {
    return(NA)
  } else {
    translated <- mart2[mart2$`Gene name` == name, ]$`Gene stable ID`
    print(class(translated))
    return(translated)
  }}))
# we use that col to subset count matrix, which means going from
# ~35,000 genes to 27,000
counts2 <- counts[!is.na(counts$ensg), ]

# now we interchange rownames with colnames
rownames(counts2) <- counts2$ensg
counts2$ensg <- NULL

## METHOD USING BIO INFO FROM BIOMART
# obtain rest of bioinfo (skip this step if you obtained everything)
write.table(row.names(counts2), file = "ensg_toget.txt", quote = F,
            row.names = F, col.names = F)

mart_plus <- 
  read.delim("projects/covid19/data/geo_series/GSE152075_ensg_more_info.txt",
             row.names = 1, stringsAsFactors = F)
mart_plus
# create biotypes named vector
biotypes <-  setNames(mart_plus$Gene.type, rownames(mart_plus))
head(biotypes)
# create gc named vector
gcseq <- setNames(mart_plus$Gene...GC.content, rownames(mart_plus) )
head(gcseq)
# create factors
factors <- data.frame(row.names =  rownames(ed), 
                      infectivity = ed$sars.cov.2.positivity,
                      batch = ed$sequencing_batch,
                      age = ed$age,
                      gender = ed$gender
                      )
factors
rownames(factors)
# subset mart2 info considering counts ensg:
mart_noiseq <- mart2[mart2$`Gene stable ID` %in% rownames(counts2), ]
rownames(mart_noiseq) <- mart_noiseq$`Gene stable ID`  ## needed for NOISEQ
mart_noiseq$`Gene stable ID` <- NULL
mart_noiseq$`Gene name` <- NULL
mart_noiseq
# create length
ensg_len <- setNames(
(as.numeric(mart_noiseq$`Gene end (bp)`) - as.numeric(mart_noiseq$`Gene start (bp)`)),
 rownames(mart_noiseq))

myData_all <- readData(data = counts2, factors = factors,
                      chromosome = mart_noiseq,
                      length = ensg_len,
                      biotype = biotypes,
                      gc = gcseq)  # ALSO HAVE OPTIONS FOR GC 

myData_all
# in this function, parameter k refers to the minimum number of counts 
# required for the gene to be counted as detected. When `factor = NULL` 
# these are calculated for each sample. If instead factor = a string that 
# matches the name of one of your columns in your `myfactors` object,
# the samples and data within that condition are aggregated. Try `factor = "Tissue"`.
mybiodetection <- dat(myData_all, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2), cex = 1)
explo.plot(mybiodetection, samples = c(1, 480), plottype = "persample")

#Sensitivity plot allows you to see the average expression and expression 
# distribution across all samples. 
mycountsbio = dat(myData_all, factor = "infectivity", type = "countsbio")
mycountsbio = dat(myData_all, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = "protein_coding", samples = c(1:5,475:480), 
           plottype = "boxplot")

#RNA-seq can preferentially detect longer genes (versus shorter genes) that
# are expressed in the same molecular quantities. NOISeq includes a function 
# that is able to normalize for gene length by sorting genes into bins using 
# gene length annotation data (`mylength`) and fitting a cubic spline regression.
#In the function below, `factor="Tissue"` generates a separate 
# normalization/spline fitting for each Tissue condition. This is optional.
mylengthbias = dat(myData_all, factor = "infectivity", 
                   type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")
# # If the numberof  samples  or  conditions  to  appear  in  the  plot  is  2 
# or  less  and  no  biotype  is  specified  (toplot  = “global”),  a10
# # diagnostic test is provided.  A cubic spline regression model is fitted
# to explain the relationship between lengthand expression.  Both the model
# p-value and the coefficient of determination (R2) are shown in the plot
# as wellas  the  fitted  regression  curve.   If  the  model  p-value
# is  significant  and  R2  value  is  high  (more  than  70%),  
# theexpression depends on the feature length and the curve shows the type of dependence


# The “GCbias” plot describes the relationship between the feature GC content
# and the expression values.
myGCbias = dat(myData_all, factor = "infectivity", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
# If the model p-value is significant and R2 value is high (more than 70%), 
# the expression will depend onthe feature GC content and the curve will show the type of dependence.

# PCA is one of the most useful visualizations you can do in order to determine in an 
# unsupervised way whether your samples are clustering according with your experimental
# design. It is also valuable to determine whether technical noise due to batch effects # or other artifacts are contributing meaningful variation that you’ll need to remove.
myPCA = dat(myData_all, type = "PCA")
explo.plot(myPCA, factor = "infectivity")
explo.plot(myPCA, factor = "age")
explo.plot(myPCA, factor = "gender")
explo.plot(myPCA, factor = "batch")

# apply batch effect correction and normalization
mydata2corr1 = ARSyNseq(myData_all, factor = "batch", batch = TRUE, norm = "rpkm", logtransf = F)
mydata2corr1 = ARSyNseq(myData_all, factor = "batch", batch = TRUE, norm = "tmm", logtransf = F)
myPCA = dat(mydata2corr1, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA, factor = "infectivity")
explo.plot(myPCA, factor = "batch")

# save.image(file = "noiseq.Rda")
# load("")
# export
#to see in console
QCreport(myData_all, samples = NULL, factor = "infectivity", norm = FALSE)
#to export as a pdf with a custom name
QCreport(myData_all, file = "QCreportForRNASeqWithNOISeq.pdf", samples = NULL, 
         factor = "infectivity", norm=FALSE)

mynoiseqBio <- noiseqbio(myData_all,  norm = "tmm", plot = T, 
                         factor = "infectivity",
                         k = 0.5, lc = 1, r = 20, adj = 1.5,
                         a0per = 0.9, random.seed = 12345, filter = 1)

mynoiseqBio.deg = degenes(mynoiseqBio, q = 0.95, M = NULL)
"11167 differentially expressed features"
# plot
par(mfrow = c(1, 2))
DE.plot(mynoiseqBio, q = 0.95, graphic = "expr", log.scale = TRUE)
# 4662 DE with q = 0.95
# save
write.csv(mynoiseqBio.deg, file = "DEGs_NOISeqBio_allbio.csv")

## ALTERNATIVE JUST USING COUNT AND FACTOR DATA
# NoISEQ readData function. needs count matrix and experimental design as 
# factors
myData <- readData(data = counts[,-481], factors = factors)
myData

# filtering
myfilt <- filtered.data(counts[,-481], factor = factors$infectivity, norm = F, 
                        method = 1, 
                        p.adj = "fdr", cpm = 1, cv.cutoff = 100)
dim(myfilt)
# [1] 2770  480

# normalization (can be applied separately)
myTMM = tmm(assayData(myData)$exprs, long = 1000, lc = 0)
head(myTMM[, 1:4])
head(counts2[, 1:4])

# DE
# technical replicates (from the same sample)
mynoiseqTech = noiseq(myData, k = 0.5, norm = "tmm", 
                      # factor = "sars.cov.2.positivity",
                      # conditions = ed$sars.cov.2.positivity,
                      pnr = 0.2, nss = 5, v = 0.02, lc = 1, 
                      replicates = "technical")
# biological replicates (different samples)
mynoiseqBio <- noiseqbio(myData,  norm = "tmm", plot = F, 
                          factor = "infectivity",
                          k = 0.5, lc = 1, r = 20, adj = 1.5,
                          a0per = 0.9, random.seed = 12345, filter = 1)

# higher probability = more chances of changes due to different treatment of 
# biological samples:
head(mynoiseqBio@results[[1]])

# DE bio 
mynoiseqBio.deg = degenes(mynoiseqBio, q = 0.95, M = NULL)
# plot
par(mfrow = c(1, 2))
DE.plot(mynoiseqBio, q = 0.95, graphic = "expr", log.scale = TRUE)
# 4662 DE with q = 0.95
# save
write.csv(mynoiseqBio.deg, file = "DEGs_NOISeqBio.csv")

# load information from geo GSE149601 -------------------------------------------
# all cases

load("covid19/scripts/rdata/annot_table_GSE149601.Rda", verbose = T)
head(ed)
# all cases

# annot should be gene annotations
# experimental desing must have the same rownames as the rest of my data, 
# in same order, if there's a col named samples
rownames(ed) <- ed$title
# ed$sample <- NULL

# load count matrix
# rownames= genes, colnames = samples
counts <- read.delim("covid19/data/geo_series/GSE149601_Boson_raw_read_counts.tsv", 
                     row.names = 1, sep = "\t")
head(counts)

all(rownames(ed) == colnames(counts))
# if needed to reorder ed[colnames(counts),]

# need to convert rownames of count matrix from ENGS to genename
# rownames ENSG to genename -----------------------------------------------

# Identifiers in ENSG --> Change to GeneName --> Use mean
load("general_scripts/conversion_to_HGNC/conversion_to_HGNC_2020.01.30.RData")
corres <- conversion_to_HGNC$ensembl_gene_id

head(corres) #To check which column is gene and wich is probe
dim(corres) # 42499   2

# Hay filas duplicadas?
nrow(corres) == nrow(corres[!duplicated(corres),]) # TRUE = no hay filas duplicadas
corres = corres[!duplicated(corres),] # Si las hay las quitamos
dim(corres) # 42499   2

# Hay genes sin ENSEMBL (de normal no) y ENSEMBL sin gene (de normal sí)? Los quitamos
sum(corres$ensembl_gene_id=="") #0
sum(corres$hgnc_symbol=="")  #0
(nrow(corres) - nrow(corres[corres$hgnc_symbol!="" & corres$ensembl_gene_id!="", ])) == sum(corres$ensembl_gene_id=="")
corres_static = corres[corres$hgnc_symbol!="" & corres$ensembl_gene_id!="", ]
dim(corres_static) # [1] 42499     2

# Hay ENSEMBL IDs (ENSG son genes, no transcritos) que mapean a más de un gen? 
# Las quitamos, porque no nos fiamos.
sum(duplicated(corres_static$ensembl_gene_id)) # 7
dup_names = as.character(corres_static[duplicated(corres_static$ensembl_gene_id),"ensembl_gene_id"])
length(dup_names) == sum(duplicated(corres_static$ensembl_gene_id)) # Debe ser TRUE
corres_static = corres_static[!corres_static$ensembl_gene_id %in% dup_names,] # Quitamos los ENSGs que mapean a más de 1 gen
sum(duplicated(corres_static$ensembl_gene_id)) # Debe ser 0
dim(corres_static) # 42485     2
# Hay que comprobar que los genes que quito porque el mismo ENSG mapea a varios
# ya no están en corres_static. EJEMPLO:
# GEN   ENSG
#  A    E1  --> Lo quiero (es único) 
#  B    E2  --> NO lo quiero (Mismo ENSG E2  asociado a diferente gen) 
#  C    E2  --> NO lo quiero (Mismo ENSG E2 asociado a diferente gen) 
#  D    E3  --> Lo quiero (Mismo gen codificado por ENSG distintos)
#  D    E4  --> No lo quiero (A pesar de mismo gen codificado por ENSG distintos, E4 también codifica otro gen E)
#  E    E4  --> No lo quiero (Mismo ENSG E4  asociado a diferente gen) 
corres_static[corres_static$hgnc_symbol=="PLEKHG7",]
corres_static[corres_static$ensembl_gene_id=="ENSG00000230417",]

# Filtrar corres_static por los ENSGs que están en la matriz de datos
# matriz de datos = rownames de matriz de conteos 
rownames(counts)
corres = corres_static[corres_static$ensembl_gene_id %in% rownames(counts),]
dim(corres)
# [1] 31931     2

# Pasamos a Gen
length(unique(corres$hgnc_symbol)) # 31924
sum(duplicated(corres$hgnc_symbol)) # 7
# next function last for like 5-10 min
mat_annot = do.call("rbind",lapply(unique(corres$hgnc_symbol),function(x){
  ensgs <- as.character(corres$ensembl_gene_id[corres$hgnc_symbol == x])
  summ_expr <- apply(as.matrix(dat[ensgs,]),2,mean) # No la mediana porque cada gen sólo tiene 2 ENSG máximo
}))
head(mat_annot)
dim(mat_annot) #[1] 31924    480
rownames(mat_annot) = unique(corres$hgnc_symbol)

# now I have an object with the the genenames in the rownames and must
# transform rownames of original count matrix
length(rownames(counts)) == length(rownames(mat_annot))

rownames(mat_annot) # genes translated, but not all
rownames(counts) 

aux <- corres[corres$ensembl_gene_id %in% rownames(counts), ]
aux <- aux[!duplicated(aux$hgnc_symbol), ]
counts$genename <- unlist(lapply(rownames(counts), function(name){
  print(name)
  if (!name %in% aux$ensembl_gene_id) {
    return(NA)
  } else {
  translated <- aux[aux$ensembl_gene_id == name, ]$hgnc_symbol
  return(translated)
}}))
counts$genename
rownames(counts)

# subset counts to match translated genenames
counts2 <- counts[!(is.na(counts$genename)), ]
rownames(counts2) <- counts2$genename
head(counts2)
counts2$genename <- NULL

# bear in mind: we went from 48,445 rows to 31,924 

# create DGELIST object ---------------------------------------------------


# create an object of DGEList 
dge_gse149 <- DGEList(counts = counts2)
dge_gse149$samples
dim(dge_gse149)

# make sure ed belongs to corresponding GSE
# we need to do an exp design mat
design_gse149 <- model.matrix(~0+sars.cov.2.positivity, data = ed)
rownames(design) <- rownames(ed)
head(design)

# filter low counts or count 0 
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = F]
dim(dge)

# alternative way to remove low counts
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
d <- dge[-drop,] 
dim(d)

# TMM Normalization
# Trimmed mean of M values (TMM) normalization estimates sequencing
# depth after excluding genes for which the ratio of counts between a pair
# of experiments is too extreme or for which the average expression is too
# extreme. The edgeR software implements a TMM normalization.

dge <- calcNormFactors(dge)
head(dge@.Data)

# # When the library sizes are quite variable between samples, then the voom approach 
# is theoreticallymore  powerful  than  limma-trend.   In  this  approach,  
# the  voom  transformation  is  applied to  thenormalized and filteredDGEListobject:
v <- voom(dge, design, plot = T)
v <- voom(dge, design, plot = T, normalize.method = "quantile")

# alternative
v <- voom(d, design, plot = T)
v <- voom(d, design, plot = T, normalize.method = "quantile")
# v <- voom(dge, design, plot = T, normalize.method = "quantile")

names(v) # E normalized values
head(v$E)

# repersent with boxplot
dat <- v$E
dat <- normalizeQuantiles( dat, ties = T)
plot_Boxplot(dat = dat, ed, condition = "sars.cov.2.positivity")

# the usual limma pipelines for differential expression can be applied, for example
fit <- lmFit(v, design = design)
fit <- eBayes(fit)
topTable(fit, coef = ncol(design))

# to give more weight to fold-changes in the ranking, one could use say:
fit <- treat(fit, lfc = log2(1.2))
topTreat(fit, coef = ncol(design))


