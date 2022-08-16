########### Description #######################################################
# 23 oct 2020

# 1000genomes extraction of specific variants and subsequent analysis
# ihc.europa@gmail.com

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

# specific

library(gsrmUtils)
library(gsrmNGSAnalysis)
library(RWeka)

# general use

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(gplots)  # for heatmap2 function
library(RColorBrewer)  # for color palettes
library(clusterProfiler)  # for enrichment
# for visualization of dendrograms
library(factoextra)
library(dendextend)
# for getting clusters n
library(NbClust)
library(ggdendro)
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

# save.image("workspace/projects/pof/scripts/rdatas/# ")
load("scripts/rdata/de_gse150720_gsrmutils.Rda")

# Create variable for working directory
dir.data <- "workspace/projects/"
getwd()
setwd("")


# tabix command order -----------------------------------------------------

# first, load 66 variants 

v66 <- read.delim(
  "results/14_filtered_variants/variants_66_with_info.tsv", header = T,
           sep = "\t", stringsAsFactors = F)
v66$id
pos_targets <- do.call("rbind", lapply(1:nrow(v66), function(i){
  cat(i, "\n") # iterate every column (subpath)
  id <- v66[["id"]][[i]]
  pos <- unlist(strsplit(x = id, split = "_")[[1]][2])
  chr <- gsub(unlist(strsplit(x = id, split = "_")[[1]][1]), pattern = "chr",
              replacement = "")
  print(pos)
  print(chr)
  target <- paste0(chr, ":", pos, "-", pos)
  print(target)
  data.frame(id = target, chr = chr, stringsAsFactors = F)
}))
head(pos_targets)


# load tabix (in rstudio terminal, ml HTSlib)
as.vector(pos_targets[ pos_targets$chr == "1", ]$id )

tabix(regions = as.vector(pos_targets[ pos_targets$chr == "1", ]$id ), 
      file.name = chr1)

# will return a dataframe with corresponding genotypes and samples from
# 1000g vcf


# clinical_values_prediction ----------------------------------------------

# abc model to order rows
abc <- read.delim("results/20_for_paper/nbclusters_66_3v_AvsBvsC.csv", 
                  header = T, sep = ",", stringsAsFactors = F, row.names = 1)
abc

# clinical values for prediction dataframe
amh_afc <- read.delim2("results/21_predictions_test_models/pof_amh_afc_CACO.csv",
           sep = ",", header = T, stringsAsFactors = F, row.names = 1)
amh_afc

# change order to assing correct class
amh_afc_abc <- amh_afc[rownames(abc),]
amh_afc_abc
amh_afc_abc$class <- abc$class
# remove NA (either AMH/AFC missing values) samples
amh_afc_abc <- amh_afc_abc[!is.na(amh_afc_abc$AMH), ]
amh_afc_abc

# ADD a col for samples or remove them...if not it will trigger an error in weka!
write.table(amh_afc_abc, 
            file = "results/21_predictions_test_models/pof_amh_afc_CACO_abc.csv",
            sep = ",", quote = F, col.names = T)

# fsh

fsh <- read.delim2("results/21_predictions_test_models/FSH.csv", 
                  header = T, sep = ",", stringsAsFactors = F, row.names = 1)

# same process (order by abc, assing corresponding class...)
fsh_abc <- fsh[rownames(abc),]
fsh_abc
fsh_abc$class <- abc$class
fsh_abc <- fsh_abc[!is.na(fsh_abc$FSH), ]
fsh_abc
fsh_abc <- fsh_abc[!fsh_abc$class == "C", ]
fsh_abc

write.table(fsh_abc, 
            file = "results/21_predictions_test_models/FSH_ab.csv",
            sep = ",", quote = F, col.names = T, row.names = F)

# amh + afc + 66v genotypes

abc <- read.delim("results/20_for_paper/nbclusters_66_3v_AvsBvsC.csv", 
                  header = T, sep = ",", stringsAsFactors = F, row.names = 1)
abc

# clinical values for prediction dataframe
amh_afc_missings <- 
  read.delim2("results/21_predictions_test_models/pof_amh_afc_CACO_with_missings.csv",
                       sep = ",", header = T, stringsAsFactors = F, row.names = 1)
amh_afc_missings

# replace missing values with mean of the numerical distribution 
# AMH
amh_afc_missings$AMH
mean(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AMH), ]$AMH))  # mean is more affected by outliers 0.97
median(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AMH), ]$AMH)) # median is less affected by outliers 0.3

amh_afc_missings[is.na(amh_afc_missings$AMH),]$AMH <- mean(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AMH), ]$AMH))  
amh_afc_missings$AMH
amh_afc_missings
# AFC
amh_afc_missings$AFC
mean(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AFC), ]$AFC))  # mean is more affected by outliers 5.45
median(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AFC), ]$AFC)) # median is less affected by outliers 4

amh_afc_missings[is.na(amh_afc_missings$AFC),]$AFC <- 
  median(as.numeric(amh_afc_missings[!is.na(amh_afc_missings$AFC), ]$AFC))  
amh_afc_missings$AFC
amh_afc_missings

# change order to assing correct class
amh_afc_abc <- amh_afc_missings[rownames(abc),]
amh_afc_abc
amh_afc_abc$class <- abc$class

aux <- cbind(abc[, -ncol(abc)], amh_afc_abc)

# remove NA (either AMH/AFC missing values) samples
# aux <- aux[!is.na(aux$AMH), ]
# aux <- aux[!is.na(aux$AFC), ]

write.table(aux, 
            file = "results/21_predictions_test_models/66v_amh_afc_150.csv",
            sep = ",", quote = F, col.names = T, row.names = F)

# write.arff(aux, 
#            file = "results/21_predictions_test_models/66v_amh_afc.arff")


# IGSR_TestModel_results --------------------------------------------------

res = fread(file = "results/21_predictions_test_models/06_66_variants/weka_results_model_ABC_1000g.csv", sep = ",")
res
classPred = lapply(res$predicted, function(x){
  res = unlist(strsplit(x, split = ":", fixed = T)[[1]])[2]
})
res$classPred = unlist(classPred)
res
str(res)
score = 0.8
nA = res[res$classPred == "A" & res$prediction >= score,]
nA
nB = res[res$classPred == "B" & res$prediction >= score,]
nB

nrow(nA)/1271
nrow(nB)/1271
