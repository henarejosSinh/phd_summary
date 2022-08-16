########### Description #######################################################
# 11 aug 2020

# DE covid
# ihc.europa@gmail.com

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

# specific

library(gsrmUtils)

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


# functions ---------------------------------------------------------------

diffExprAnalysis = function(dat, ed, condition, paired = FALSE, pair = NULL){
  # Create experimental design
  group = ed[colnames(dat), condition]
  if (paired){
    pair = ed[colnames(dat), pair]
    design <- model.matrix(~ 0+group+pair)
    rownames(design) <- colnames(dat)
    fit <- lmFit (dat, design= design)
    res.limma <- eBayes (fit)
    res = list(topTable(res.limma, coef=1, number = nrow((dat))))
    names(res) = paste(gsub("group", "", colnames(design)[1:2]), collapse="-")
  }else{
    design <- model.matrix(~ 0+group)
    rownames(design) <- colnames(dat)
    # Create contrasts
    combs = combn(x=colnames(design), m=2, simplify = TRUE)
    cons2 = apply(combs,2,function(y){
      paste(y,collapse="-")
    })
    cons = gsub("group","",cons2)
    contrasts <- makeContrasts(contrasts=cons2, levels = design)
    
    # Fit linear model
    fit <- lmFit (dat, design= design)
    fit.cont <- contrasts.fit (fit, contrasts)
    res.limma <- eBayes (fit.cont)
    # Create results
    res = list()
    for(i in 1:length(cons2)){
      res[[cons2[i]]] = topTable(res.limma, coef = i, number = nrow(dat))
    }
    names(res) = gsub("group","",names(res))
  }
  
  # Include FDRs column
  res2 = lapply(names(res), function(x){
    # cat(x, "\n")
    aux = unlist(strsplit(x, split = "-", fixed = T))
    condition1 = aux[1]
    condition2 = aux[2]
    res[[x]]$FC = unlist(lapply(rownames(res[[x]]),function(y){
      # cat(y, "\n")
      A = mean(as.numeric(dat[y, group==condition1]))
      B = mean(as.numeric(dat[y, group==condition2]))
      FC = 2^(A-B)
      if (FC<1){
        FC = -1/FC
      }
      FC
    }))
    res[[x]]
  })
  names(res2) = names(res)
  
  return(res2)
}

plot_volcanoPlot = function(DEG_results, threshold = 0.05, up_color = "#D55E00", down_color = "#0072B2",
                            no_sig_color = "dimgray", labels = T, title = ""){
  toplot = DEG_results[,c("logFC","adj.P.Val")]
  toplot$p.value = -10*log10(toplot$adj.P.Val)
  toplot$sig = ifelse(toplot$adj.P.Val<=threshold,"SIG","NO_SIG")
  toplot$sense = ifelse(toplot$logFC>0,"UP","DOWN")
  toplot$sig2 = toplot$sig
  toplot$sig2[toplot$sig2=="SIG"] = paste(toplot$sig2[toplot$sig2=="SIG"],toplot$sense[toplot$sig2=="SIG"],sep="_")
  toplot$names = rownames(DEG_results)
  toplot$names[toplot$sig=="NO_SIG"] = ""
  
  if (labels){
    p = ggplot(toplot,aes(x=logFC,y=p.value,color=sig2,label=names))+
      geom_point(size=2)+
      geom_text(hjust=0, vjust=0)
  }else{
    p = ggplot(toplot,aes(x=logFC,y=p.value,color=sig2,label=names))+
      geom_point(size=2)
  }
  p = p+
    geom_hline(yintercept=-10*log10(threshold),colour="cadetblue4",linetype="longdash",size=0.2)+
    theme_bw()+
    xlab("log2 fold change")+
    ylab("-10*log10(p-value)")+
    ggtitle(title)+
    scale_color_manual(values=c(SIG_DOWN=down_color,SIG_UP=up_color,NO_SIG=no_sig_color))+
    theme(legend.position="none",
          axis.title = element_text(size=18),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22, hjust=0.5))
  
  return(p)
}



# Differential expression start -------------------------------------------

gse <- "GSE150275" # 484 samples from covid CA/CO (4 extracted because BQ)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075
file <- paste0("dat_normalized_",gse,".Rda")
load(paste0("scripts/rdata/", file), verbose = T)                               
load("scripts/rdata/dat_normalized_GSE150275.Rda")

nrow(dat_norm$targets) # 480 samples (4 excluded)
head(dat_norm$targets) 
head(dat_norm$E)
head(dat_norm$weights)
head(dat_norm$design)
dim(dat_norm$targets)

# important: rownames from ed must be colnames from count matrix
rownames(ed) <- ed$title
ed


# pca analysis ---------------------------------------------------------------------

# do an overview of Sex, Age, n1_ct, sequencing batch, infection
groupAge
summary(as.numeric(ed$age))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    27.0    43.0    42.6    59.0    77.0 
# Does not work well...there's people with 80~ years, and some
# have 90+ years...
ed$age <- as.character(ed$age)


groupAge_2 <- function(age){
  num <- as.numeric(age)
  res <- ifelse(age == "90+", yes = "> 90", 
                no = ifelse(age == "Unknown", 
                            yes = "Unknown", 
                            no = ifelse(num < 59, ifelse(num < 43, ifelse(num < 27, "< 27", "[27, 43]"), "[43, 59]"), "[59, 90]")))
}
res <- groupAge_2(ed$age)
ed$grouped_age <- res
plot_PCAscores

summary(as.integer(as.character(levels(ed$n1_ct))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   12.00   18.00   21.00   20.69   23.00   30.00       2 

# ct explanation:
# ct stands for "cycle threshold" which is the tempeture that was necessary
# to amplify the DNA after RNA conversion. Low ct means that a lot 
# of RNA existed from the start.
# N1 SARS is a fragment of the virus that is target for PCR
ed$n1_ct
ed$n1_ct[ed$n1_ct == "N/A"] <- NA
ed$n1_ct[ed$n1_ct == "Unknown"] <- NA
ed$n1_ct

nrow(ed[ed$sars.cov.2.positivity == "pos" & is.na(ed$n1_ct), ])
# 17 cases don't have CT (neither age, gender... sequencing batch T,R)

groupCT <- function(ct){
  num <- as.numeric(as.character(ct))
  res <- ifelse( is.na(ct), yes = "NA", 
                 no = ifelse(num < 23, ifelse(num < 20.69, ifelse(num < 18, "< 18", "[18, 20.69]"), "[20.69, 23]"), "> 23"))
}
res <- groupCT(ed$n1_ct)
ed$n1_ct_grouped <- res

# repeat for age, gender, sars.cov.2.positivity, sequencing_batch
plot_PCAscores(dat = dat_norm$E, ed = ed, condition1 = "n1_ct_grouped", 
               components = c(1,2))  
# CACO
plot_PCAscores(dat = dat_norm$E, ed = ed, condition1 = "sars.cov.2.positivity", 
               components = c(1,2), colors = c("lightsteelblue4", "firebrick"))  


# DE ----------------------------------------------------------------------

# make a subset to exclude non Female samples from study

ed_2 <- ed[ed$gender == "F", ]
ed_2
ed_2 = as.data.frame(ed_2)
ed_2
dim(ed_2) #231 women
ed_2
rownames(ed_2) <- ed_2$title
rownames(ed_2) 

# how many women?

table(ed_2[ed_2$sars.cov.2.positivity == "neg", ]$gender)  #30
table(ed_2[ed_2$sars.cov.2.positivity == "pos", ]$gender)  #201
norm_counts_female <- dat_norm$E[, colnames(dat_norm$E) %in% ed_2$title ]
ncol(norm_counts_female) # 231

# what is the mean of n1_ct, age?
mean(as.numeric(ed_2$age))
median(as.numeric(ed_2$age))
mean(as.numeric(ed_2[ed_2$sars.cov.2.positivity == "pos", ]$n1_ct))
median(as.numeric(ed_2[ed_2$sars.cov.2.positivity == "pos", ]$n1_ct))

# CHANGE POS/NEG levels of factor (29/03/2021)
levels(ed_2$sars.cov.2.positivity) = c("pos", "neg")
levels(ed_2$sars.cov.2.positivity) 
# save for correlation matrix
# write.table(x = norm_counts_female, "coun_matrix_normalised_only_women.tsv",
# sep = "\t",  quote = F,  row.names = T,  col.names = T)

# REPEAT PCA after excluding samples
plot_PCAscores(dat = norm_counts_female, ed = ed_2, condition1 = "grouped_age", 
               components = c(1,2))  
# CACO
plot_PCAscores(dat = norm_counts_female, ed = ed_2, 
               condition1 = "sars.cov.2.positivity", 
               components = c(1,2), colors = c("lightsteelblue4", "firebrick"))  
# pos = firebrick, neg = lightsteelblue4
# follow instructions grsmutils wiki
# In the case of having two groups, the diffExprAnalysis() function should be 
# used. It is the case, is a CA/CO

# The diffExprAnalysis() function performs a tipical limma R-package 
# differential expression analysis between each pair of groups defined by the
#interest variable to identify genes that are signficantly differentially 
#expressed (DEGs) between the groups of each comparison. 
# detach("package:data.table", unload = T) # wont work otherwise?
DE_results = diffExprAnalysis(dat = norm_counts_female, ed = as.data.frame(ed_2), 
                              condition = "sars.cov.2.positivity", pair = NULL,
                              paired = F)
# male and unknown included samples ~4000 DE genes > 2800~ FDR

# only female 5056 > 2580 padj by FDR
dim(DE_results$`neg-pos`[DE_results$`neg-pos`$adj.P.Val <= 0.05, ])
dim(DE_results$`pos-neg`[DE_results$`pos-neg`$adj.P.Val <= 0.05, ])


# The diffExprAnalysis() function returns a list with all the possible comparisons.
# Each element of the list is named as the comparison
# (early_secretory vs late_secretory, etc.) and includes a data.frame with the
# p-value of the comparison, the adjusted p-value
# (typically calculated with the Benjamini & Hochberg method (8)
#   known as False Discovery Rate (FDR))
# and the fold-change (FC) for each gene in the dataset. Results for
# each comparison can be plotted using Volcano plots with the function
plot_volcanoPlot(DEG_results = DE_results$`neg-pos`) 
plot_volcanoPlot(DEG_results = DE_results$`pos-neg`) 
dim(DE_results$`neg-pos`)

write.table(DE_results$`pos-neg`, file = "de_results_caco.tsv", quote = F, row.names = T,
            col.names = T, sep = "\t")
# neg-pos comparison: blue genes are expressed lower in controls than cases,
# while red genes are genes expressed higher in control than cases

# check expression of known covid genes
check <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", "CTSL")
DE_results$`neg-pos`["TMPRSS2", ]
DE_results$`neg-pos`["GPX3", ]
View(DE_results$`neg-pos`[check,])
# TMPRSS4 is not significant
check_s <- c("TMPRSS2", "BSG", "MX1", "CTSB", "CTSL")
# FURIN, ACE2 not in DE
plot_volcanoPlot(DEG_results = DE_results$`neg-pos`[check, ]) 

# res
significant_genes <- lapply(DE_results, function(res){
  out = res[res$adj.P.Val <= 0.05,]
  return(out)
}) # 2895
lapply(significant_genes, nrow)

# heatmap
plot_HeatmapSelectedGenes(dat = dat_norm$E[, c(1:20,450:470)], 
                          sel_genes = rownames(significant_genes$`neg-pos`[check_s,]),  #[1:10], 
                          ed = ed, condition = "sars.cov.2.positivity")

plot_genesProfile(genes = check, dat = dat_norm$E, 
                  ed = ed, condition = "sars.cov.2.positivity")

# COMPARISON WAS CO/CA
# We should change FC sign to make it more understandable:
DE_results$`neg-pos`$FC_CACO <- -(DE_results$`neg-pos`$FC)
DE_results$`neg-pos`$logFC_CACO <- -(DE_results$`neg-pos`$logFC)
head(DE_results$`neg-pos`)
sum(DE_results$`neg-pos`$adj.P.Val <= 0.05) # 2580

write.table(
  "results/1-covid/de_gse150720_gsrmutils_voom_cmp_quartile_females_fc_corrected.txt", 
            sep = "\t",
            x = DE_results$`neg-pos`, quote = F, row.names = T, col.names = T )
# save enviroment
save.image("scripts/rscripts/de_gse150720_gsrmutils_27oct20.Rda")


# plot expression ---------------------------------------------------------
# Remove experiment effect
# Evaluate expression USING NORMALIZED DATASET!!
genes_of_interest <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", 
                       "CTSL") # "ACE2", "FURIN" LOST IN PREPROCESSING
DE_results$`neg-pos`[genes_of_interest, ]

dat = melt(norm_counts_female[genes_of_interest,])  # norm data
head(dat)
# new col of genes before melt

dat$covid = ed_2[dat$Var2, "sars.cov.2.positivity"]  #ed must have rownames
head(dat)
colnames(dat) = c("gene", "sample", "expression_value", "covid")
table(is.na(dat))
# View(dat2)
limits = quantile(dat$expression_value, c(0, 0.1, 0.5, 1), na.rm = F)
toplot = summarySE(dat, measurevar = "expression_value", 
                   groupvars = c("gene", "covid"))
toplot

toplot$gene = factor(toplot$gene, levels = genes_of_interest)
lapply(genes_of_interest, function(gene){
  fc <- round(DE_results$`neg-pos`[gene, ]$FC, digits = 2)
  ymax_ref = (max(toplot[toplot$gene == gene, ]$expression_value) +
                max(toplot[toplot$gene == gene, ]$ci) ) + 1 # depends on case
  # print(ymax_ref)
  # adjpval <- round(DE_results$`neg-pos`[gene, ]$adj.P.Val, digits = 2)
  p = ggplot(toplot[toplot$gene == gene,], aes(x = covid, y = expression_value, color = covid)) +
    ylim(0, ymax_ref) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = expression_value - ci, ymax = expression_value + ci),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c(neg = "lightsteelblue4", pos = "firebrick")) +
    # geom_hline(yintercept = limits, linetype = 2, color = "grey60")+
    xlab("") + ylab("") +
    ggtitle(paste0(gene , " expression change in COVID-19, FC: ", fc)) +
    # " adj-pval: ", adjpval)) +
    theme_light() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "none")
  p
})


# diff expr hands out -----------------------------------------------------
diffExprAnalysis
# first, subset the dataframe from ed considering samples from count
# matrix and condition
ed_2$sars.cov.2.positivity = factor(ed_2$sars.cov.2.positivity, levels = c("pos", "neg"))
group <- ed_2[colnames(norm_counts_female), "sars.cov.2.positivity"]
# create design matrix
design <- model.matrix(~0 + group)
# label rownames with colnames
rownames(design) <- colnames(norm_counts_female)
head(design)
# generate model matrix of combinations between designs
combs <- combn(x = colnames(design), m = 2, simplify = T)
head(combs)
# contrasts
combs2 <- apply(combs, 2, function(y){
  paste(y, collapse = "-")
})
# correct name of contrasts and design
contrast <- gsub("group","", combs2)
colnames(design) <- gsub("group","", colnames(design))
head(design)
contrasts <- makeContrasts(contrasts = contrast, levels = design)
head(contrasts)
# fit linear model and bayes
fit <- lmFit(norm_counts_female, design = design)
fit_contrast <- contrasts.fit(fit, contrasts)
res_limma <- eBayes(fit_contrast)

# create list to put results
res <- list()
for (i in 1:length(combs2)){
  res[[combs2[i]]] <- topTable(res_limma, coef = i, 
                               number = nrow(norm_counts_female))
}
names(res) <- gsub("group","", names(res))
head(res)
# fix gene table and calculate FC
res2 <- lapply(names(res), function(x){
  aux <- unlist(strsplit(x, split = "-", fixed = T))
  condition1 <- aux[1]
  condition2 <- aux[2]
  res[[x]]$FC <- unlist(lapply(rownames(res[[x]]), function(y) {
    A = mean(as.numeric(norm_counts_female[y, group == condition1]))
    B = mean(as.numeric(norm_counts_female[y, group == condition2]))
    FC = 2^(A - B)
    if (FC < 1) {
      FC = -1/FC
    }
    FC
  }))
  res[[x]]
})
names(res2) <- names(res)

plot_volcanoPlot(DEG_results = res2$`pos-neg`) 
write.table(res2$`pos-neg`, file = "de_results_caco.tsv", quote = F, row.names = T,
            col.names = T, sep = "\t")

# pca3D -------------------------------------------------------------------


# install.packages("rgl")
library(rgl)
pc = prcomp(t(norm_counts_female))
colors = c(rep("firebrick", 201), rep("skyblue", 30))  # los que quieras 
Var = round(summary(pc)$importance[2,]*100,1)

plot3d(pc$x[, 1:3],
       col = colors, 
       type = "s",
       size = 2,
       xlab = "PC1: 11.2%",  # Esto son los valores de Var (manualmente)
       ylab = "PC2: 8.5%",
       zlab = "PC3: 7%", box=F, axes=T, cex = 1.5)

# correlations ------------------------------------------------------------

# dame azucar

# check that correlations script worked well
# ZNF24   JUP     -0.26472226080733297
# ZNF26   JUP     -0.23245460607898555
# PRPSAP1 JUP     -0.10353007497504496
# LARP1   JUP     0.2856263608992222
# RDH10   JUP     0.09107915453475797
# CREBBP  JUP     0.16617760548392907
# GSTM2   JUP     -0.11745085318390792
# JUN     JUP     0.12172696441638026
# PRPF38B JUP     -0.3684249910566091
# JAG1    JUP     0.13987287003179888
# PRPF38A JUP     -0.1297781445093985
# ITSN2   NEBL    -0.1540575043189309
# TFRC    NEBL    -0.16816136161908585
# TXNIP   OPTN    0.18118714158211685
# CLOCK   OPTN    -0.012404316723930822
# ITSN2   SF1     0.014652432316846079
# TFRC    SF1     -0.10880814535281999


cor( x =  norm_counts_female[c("ZNF24"), ], y = norm_counts_female["JUP", ] )
cor( x =  norm_counts_female[c("JUP"), ], y = norm_counts_female["ZNF24", ] )
cor( x =  norm_counts_female[c("PRPSAP1"), ], y = norm_counts_female["JUP", ] )
cor( x =  norm_counts_female[c("ITSN2"), ], y = norm_counts_female["NEBL", ] )
cor( x =  norm_counts_female[c("TFRC"), ], y = norm_counts_female["NEBL", ] )
cor( x =  norm_counts_female[c("TXNIP"), ], y = norm_counts_female["OPTN", ] )
cor( x =  norm_counts_female[c("TFRC"), ], y = norm_counts_female["SF1", ] )

correlation_table <- fread("/home/ihenarejos/workspace/projects/covid19/results/1_DE_covid/correlation_results.tsv")

hist(correlation_table$correlation)

# pos
nrow(correlation_table[correlation_table$correlation >= 0.5, ])
# [1] 44250
nrow(correlation_table[correlation_table$correlation >= 0.6, ])
# 11,114
nrow(correlation_table[correlation_table$correlation >= 0.55, ])
# 22,794
nrow(correlation_table[correlation_table$correlation >= 0.7, ])
# 2290

# nega
nrow(correlation_table[correlation_table$correlation <= -0.5, ])
# 3000
nrow(correlation_table[correlation_table$correlation <= -0.6, ])
# 63
nrow(correlation_table[correlation_table$correlation <= -0.55, ])
# 505
nrow(correlation_table[correlation_table$correlation <= -0.7, ])

num <- 0.7

dfaux <- correlation_table[correlation_table$correlation >= num, ]
dfaux2 <- correlation_table[correlation_table$correlation <= -num, ]

dftosave <- rbind(dfaux, dfaux2)

write.table(dftosave, quote = F, file = "/home/ihenarejos/workspace/projects/covid19/results/1_DE_covid/correlation_results_filtered.tsv", sep = "\t",
            row.names = F, col.names = T)

# ATT

detosave <- DE_results$`neg-pos`[,c("adj.P.Val", "FC_corrected")]
detosave$genes <- rownames(detosave)

detosave$adj_true <- ifelse(detosave$adj.P.Val <= 0.05, 1, 0)

write.table(
  detosave, quote = F, 
  file = "results/1-covid/att_covid_for_networks_corrected.tsv", 
  sep = "\t",
  row.names = F, col.names = T)



# network analysis --------------------------------------------------------

# read node table and do mean of degree
# and check which nodes have degree over mean + high betweennes centrality

table <- read.delim("projects/covid19/results/2_DE_Endometrium/node_table_endometrium_corr_filtered.csv", header = T, sep = ","  )

table$mean_degree <- mean(table$Degree)
head(table)
# subset by those node which Degree > mean degree
# hubs
table_degree <- table[table$Degree > table$mean_degree, ]
head(table_degree[
  order(table_degree[,c("Degree")], decreasing = T)
  , c("shared.name", "Degree")])


# influential genes by betweeness centrality and closeness centrality
head(table[
  order( table[,"ClosenessCentrality"], table[,"BetweennessCentrality"], 
         decreasing = T)
  , c("shared.name", "ClosenessCentrality", "BetweennessCentrality")] )
