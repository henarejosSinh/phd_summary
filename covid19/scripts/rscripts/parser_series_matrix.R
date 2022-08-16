########### Description #######################################################
# 30 july 2020

# Parse series from GEO

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

library(data.table)
library(foreach)

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
# load("workspace/projects/pof/scripts/rdatas/")

# Create variable for working directory
dir <- "workspace/projects/"

getwd()
setwd("../covid19/")

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


# Proper start ------------------------------------------------------------

args <- commandArgs(TRUE)

#### Input files ####

# Series Matrix File
SERIES_FILE <- "covid19/data/geo_series/GSE152075_series_matrix_484_samples.txt"
SERIES_FILE <- "covid19/data/geo_series/GSE149601_series_matrix.txt"
PLATFORM_FILE <- "covid19/data/geo_series/"

#### Parameters ####
# Column containing gene IDs in platform file #
# Generally "ENTREZ_GENE_ID" or "GENE"
GENE_COLUMN <- if(!is.na(args[3])) args[3] else "ENTREZ_GENE_ID"

#### Output files ####
EXPRESSION_OUTPUT_FILE <- if(!is.na(args[4])) args[4] else "expression.csv"
ANNOTATION_OUTPUT_FILE <- if(!is.na(args[5])) args[5] else "annotation.csv"

# Read characteristics
con <- file(SERIES_FILE, "r")
characteristics <- c()
while(TRUE) {
  line <- readLines(con, n = 1)
  if(length(line) == 0) {
    break
  } else if(startsWith(line, "!Sample_title")) {
    titles <- unlist(strsplit(line, "\t"))[-1]
    titles <- gsub("\\\"", "", titles)
  } else if(startsWith(line, "!Sample_characteristics")) {
    characteristics <- c(characteristics, line)
  } else if(startsWith(line, "!Sample_geo_accession")) {
    accession <- unlist(strsplit(line, "\t"))[-1]
    accession <- gsub("\\\"", "", accession)
  }
}
close(con)

# Parse characteristics
anno <- data.frame(lapply(characteristics, function(x) {
  values <- unlist(strsplit(x, "\t"))[-1]
  values <- gsub("\\\"", "", values)
  parts <- strsplit(values, ": ")
  
  name <- parts[[1]][[1]]
  values <- sapply(parts, function(x) x[2])
  
  out <- list()
  out[[name]] <- values
  return(out)
}))

ed <- data.table(sample = accession, title = titles, anno)

# save datatable if needed
save(ed, file = "covid19/scripts/rdata/annot_table_GSE149601.Rda")

# Read probe-level expression data
D <- fread(SERIES_FILE, header=TRUE, skip="\"ID_REF\"", fill=TRUE, na.strings=c("","NA","null"))
D <- D[1:(nrow(D)-1),] # remove table end marker


# process platform file -------------------------------------------------------

ref <- read.table(PLATFORM_FILE, header = TRUE, 
                  sep = "\t", quote = "", comment.char ="#", fill = TRUE)[,c("ID",GENE_COLUMN)]
colnames(ref) <- c("ID_REF","entrez")
ref$entrez <- as.character(ref$entrez)
ref <- subset(ref, !is.na(entrez) & entrez != "")

entrez_split <- strsplit(ref$entrez, " /// ")
ref <- data.frame(
  ID_REF=rep(ref$ID_REF, sapply(entrez_split, length)),
  entrez=unlist(entrez_split)
)

# Merge tables to map Entrez genes ids
m <- data.table(merge(ref, D, all=FALSE)[,-1])

# Aggregate duplicate genes by median
m <- m[, lapply(.SD, median(na.rm=TRUE)), by=entrez]

# Extract gene and sample names
genes <- m$entrez
samples <- colnames(m)[-1]

# Transpose matrix, add sample column and set gene names as column names
m <- data.table(samples, transpose(m[,-1]))
colnames(m) <- c("sample", as.character(genes))

# Write results to separate expression and annotation files
fwrite(m, file=EXPRESSION_OUTPUT_FILE, sep=",")
fwrite(anno, file=ANNOTATION_OUTPUT_FILE, sep=",")
