###########SCRIPT DESCRIPTION#######################
# jan 15,2020
# obtain network genes related to targeted hypothesis

#working directory
# dir=('/home/ihenarejos/workspace/')
# setwd(dir)
# setwd("/home/ihenarejos/workspace/")

#clean global enviroment if needed
# rm(list=ls())

########## LOAD LIBRARIES####################################

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(gplots) 
library(gsrmUtils)

########### save_load enviroment ----------------------------------------------------

save.image("workspace/projects/pof/scripts/rdatas/poi_variants_network.Rda")
load("workspace/projects/pof/scripts/rdatas/poi_variants_network.Rda")

poigenes <- as.character(read.table(file = "workspace/projects/pof/lists_genes/genespoi611.tsv")$V1)

bioinfo_possible_variants <- fread(file = "workspace/projects/pof/results/14_toconfirm/bioinfo_possible_variants.txt")
# find if there's any gene from poi list in the possible variants to confirm

# get genes in possible variants
genes_in_variants <- as.data.frame(unique(unlist(strsplit(bioinfo_possible_variants$gene,split = ","))))

# are poi genes in these genes?
intersect(genes_in_variants$`unique(unlist(strsplit(bioinfo_possible_variants$gene, split = ",")))`,poigenes$V1) # 3 poi genes

# clinvar variants --------------------------------------------------------

clinvar_info <- fread("workspace/projects/pof/results/13_study_candidates/clinvar.txt" , sep = "\t" , quote = F)

clinvar_info$id <- unlist(lapply(clinvar_info$id, function(element){chr5_79950724
  value <- paste0("chr",element)
  return(value)}))

variants_tocheck <- read.table(file = "workspace/projects/pof/data/txt_tsv/variants.to.check.txt", col.names = "id")

clinical <- clinvar_info[clinvar_info$id %in% variants_tocheck$id, ]

# or, to keep variants that do not find data in clinvar
clinical <- left_join(variants_tocheck, clinvar_info, by = "id")

# multi allelic check in clinvar and rownames
multi_check <- clinvar_info[grep("chr5_79950724", x = clinvar_info$id, fixed = F)]

clinical <- rbind(clinical,multi_check)