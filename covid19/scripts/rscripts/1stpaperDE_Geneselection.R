# 3 sep 2021

# DE endometrium controls, selection and visualization of important genes
# ihc.europa@gmail.com

#clean global enviroment if needed
rm(list = ls())

# load("scripts/rdata/endometrial_control_joined.Rda", verbose = T)

library(limma)
library(dbplyr)
library(reshape2)
library(Rmisc)
library(ggplot2)

# pre-processing -----------------------------------------------------------

files = dir("data/expression_data/", full.names = T)

dats2 = lapply(files, function(f){
  load(f)
  return(dat)
})
common_genes = Reduce("intersect", lapply(dats2, rownames))
length(common_genes)
# [1] 13437

# check genes of interest in join dataset

# Endometrium datasets processing
dats2 = do.call("cbind", lapply(dats2, function(x){
  x[common_genes,]
}))
dim(dats2)

# Experimental designs
eds2 = do.call("rbind", lapply(files, function(f){
  load(f)
  ed = data.frame(sample = rownames(ed),
                  time = ed$time,
                  experiment = gsub(
                    "data//", "", unlist(strsplit(
                      f, split = "_", fixed = T))[2], fixed = T),
                  stringsAsFactors = F)
  return(ed)
}))

# add rownames to eds2
rownames(eds2) = eds2$sample

# Fix times
eds2$time[eds2$time=="LH+2"] = "ESE"
eds2$time[eds2$time=="LH+8"] = "MSE"
eds2$time[eds2$time=="proliferative"] = "PF"
eds2$time[eds2$time=="mid-secretory"] = "MSE"
eds2$time[eds2$time=="Proliferative"] = "PF"
eds2$time[eds2$time=="Secretory"] = "MSE"
eds2$time[eds2$time%in%c("-1", "-3", "-12", "-7", "-14", "-8", "-5", "-4")] = "PF"
eds2$time[eds2$time%in%c("1", "2", "3", "4")] = "ESE"
eds2$time[eds2$time%in%c("6", "7")] = "MSE"
eds2$time = factor(eds2$time, levels = c("PF", "ESE", "MSE", "LSE"))

# Remove Altmae and Bradley outliers
dats2 = dats2[, !colnames(dats2) %in% c("GSM2593180", "GSM2593181", "GSM742078")]
eds2 = eds2[!rownames(eds2) %in% c("GSM2593180", "GSM2593181", "GSM742078"), ]

# do it by phase of the cycle
dim(dats2)
# 13,437 genes 109 samples
table(eds2$time)
# PF ESE MSE LSE 
# 29  29  43   8 
dim(dats2[,rownames(eds2[eds2$time == "PF",])])  #29
dim(dats2[,rownames(eds2[eds2$time == "MSE",])])  #43
dats2
eds2

# visualization -----------------------------------------------------------

genes_of_interest = c('COBL', 'GPX3', 'SOCS3', 'DOCK2', 'SLC2A3')
dats2[, eds2[!eds2$time %in% c('MSE') ,]$sample]

dat = melt(dats2[genes_of_interest,])
dat
dat$time = eds2[dat$Var2, "time"]
dat
colnames(dat) = c("gene", "sample", "expression_value", "time")
dat = dat[dat$time != 'MSE',]

limits = quantile(dats2, c(0, 0.1, 0.5, 1))
toplot_p = summarySE(dat, measurevar = "expression_value", groupvars = c("gene", "time"))
toplot_p
toplot_p$gene = factor(toplot_p$gene, levels = genes_of_interest)
toplot_p$group <- "A"
# Plot of expression values for each phase by gene
gene = 'SOCS3'

plot_expression_values = function(toplot, gene){
  ggplot(toplot[toplot$gene==gene,], aes(x = time, y = expression_value, color = time, group=group))+
  geom_line()+
  geom_point(size = 3, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = expression_value - ci, ymax = expression_value + ci),  
                width = 0.2, position = position_dodge(width = 0.5))+
  scale_color_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MSE = "chartreuse4",
                                LSE = "lightpink3"))+
  # geom_hline(yintercept = limits, linetype = 2, color = "grey60")+
  xlab("") + ylab("") +
  ggtitle("")+
  theme_void()+
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none")
}
# 7.3 , 4.51

toplot_p
plot_expression_values(toplot_p, gene)

