# 4 sep 2020

# DE endometrium controls
# ihc.europa@gmail.com

#clean global enviroment if needed
rm(list = ls())

# load("scripts/rdata/endometrial_control_joined.Rda", verbose = T)

library(limma)
library(reshape2)
library(Rmisc)
library(ggplot2)

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

# check covid DE genes in endometrium datasets
load("scripts/rdata/de_gse150720_gsrmutils.Rda")
genes_of_interest <- rownames(DE_results$`neg-pos`)  # 5056 genes in expr matrix
genes_of_interest <- rownames(DE_results$`neg-pos`[DE_results$`neg-pos`$adj.P.Val
                                                   < 0.05, ])  
table(genes_of_interest %in% common_genes)  # 4596 TRUE 461 FALSE
genes_of_interest <- rownames(DE_results$`neg-pos`[DE_results$`neg-pos`$adj.P.Val
                                                   < 0.05, ])  
table(genes_of_interest %in% common_genes)  # 2611 TRUE 284 FALSE

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

# remove batch effect for correlation networks ------------------------------------------------------

# REMOVE BATCH EFFECT (experiment + TIME, for menstrual cycle network)
dats_exp_time = removeBatchEffect(dats2, batch = eds2$experiment, 
                            batch2 = eds2$time )

# save general expression matrix for all phases correlation
write.table("endometrium_controls_normalised_counts_matrix_menstrual_cycle.tsv", 
            sep = "\t",
            x = dats_exp_time, quote = F, row.names = T, col.names = T )

# REMOVE BATCH EFFECT (EXPERIMENT)
dats_exp = removeBatchEffect(dats2, batch = eds2$experiment)
dats3 = dats_exp

# save and upload it to server
save(dats_exp_time, dats_exp, eds2, 
     file = "scripts/rdata/endometrium_joinControl_counts_ed.Rda")

# save tsv for each phase
for (phase in unique(eds2$time)) {
  write.table(paste0("endometrium_controls_normalised_counts_matrix_", phase,
                     ".tsv"),
              sep = "\t",
              x = dats3[,rownames(eds2[eds2$time == phase,])], quote = F, 
              row.names = T, col.names = T )
}

# DE (JOIN DATASET) -------------------------------------------------------
# 
# load("scripts/rdata/endometrium_joinControl_counts_ed.Rda")
# dats2 = read.table("results/2-endometrium/endometrium_controls_normalised_counts_matrix_MSE.tsv",
#            sep = "\t", header = T, row.names = 1)

head(dats_exp)
# eds2_phase = eds2[eds2$time == "MSE",]

# Differential analysis 
# Condition "time" (not corrected previously) to check the differences in expression
# through phases. To check a particular phase against the other, rename ed
# so is the study phase vs the others
ed_exp = eds2
ed_exp$time = factor(ed_exp$time, levels = c("Rest", "MSE", "LSE", "ESE", "PF"))
ed_exp
ed_exp$time[!ed_exp$time == "MSE"] = factor("Rest")
ed_exp$time = factor(ed_exp$time, levels = c("MSE", "Rest"))

DEGs = diffExprAnalysis(dat = dats_exp, ed = ed_exp, condition = "time") # condition experiment for original analysis
sigs = lapply(DEGs, function(deg){
  rownames(deg)[deg$adj.P.Val<=0.05]
})
to_remove = unique(unlist(sigs))
length(to_remove) # 0 --> No genes to be removed (JOIN Dataset)

# sigs of study phase vs Rest
DEGs$`MSE-Rest`[to_remove,]

DEGs$`MSE-Rest`$gene = rownames(DEGs$`MSE-Rest`)
View(DEGs$`MSE-Rest`)

write.table(DEGs$`MSE-Rest`, file = "results/2-endometrium/MSE_vs_rest_expression_matrix.tsv",
            sep = "\t", row.names = F, col.names = T)

# plot expression ---------------------------------------------------------
genes_of_interest = c('ACE2', 'TMPRSS2', 'TMPRSS4', 'CTSL', 'CTSB', 'BSG', 'FURIN', 'MX1')
# Evaluate expression
dat = melt(dats3[genes_of_interest,])
dat$time = eds2[dat$Var2, "time"]
colnames(dat) = c("gene", "sample", "expression_value", "time")
limits = quantile(dats3, c(0, 0.1, 0.5, 1))
toplot = summarySE(dat, measurevar = "expression_value", groupvars = c("gene", "time"))
toplot$gene = factor(toplot$gene, levels = genes_of_interest)

gene = 'ACE2'

checks = eds2[eds2$time == 'MSE',]$sample
mean(dats_exp['ACE2', checks]) 
toplot[toplot$time == 'MSE' & toplot$gene == 'ACE2',]

# Plot of expression values for each phase by gene

p = ggplot(toplot[toplot$gene==gene,], aes(x = time, y = expression_value, color = time))+
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
p
# Evaluate differences
anova = anovaAnalysis(dat = dats3[genes_of_interest,], ed = eds2, condition = "time")
pairwise.t.test(dats3[gene,], eds2$time)



# correlation analysis ----------------------------------------------------

correlation_table <- fread("/data/network/nas01/fivibio-data/projects/covid/results/1-DE/correlation_endometrium_results.tsv")

hist(correlation_table$correlation)

# check that correlations script worked well
# LURAP1L	GSTM5	0.7577578713600622
# FAM19A2	GSTM5	0.7961742853663709
# SHMT2	PKIG	0.8595563833276877
# TARSL2	PKIG	0.8369100003411295

cor( x =  dats2[c("LURAP1L"), ], y = dats2["GSTM5", ] )
cor( x =  dats2[c("FAM19A2"), ], y = dats2["GSTM5", ] )
cor( x =  dats2[c("SHMT2"), ], y = dats2["PKIG", ] )
cor( x =  dats2[c("PKIG"), ], y = dats2["TARSL2", ] )
cor( x =  dats2[c("MRPL9"), ], y = dats2["RPL30", ] )
cor( x =  dats2[c("MX1"), ], y = dats2["CTSL", ] )


cor( x =  dats3[c("MX1"), ], y = dats3["CTSL", ] )
cor( x =  dats3[c("SHMT2"), ], y = dats3["PKIG", ] )

dats2[c("MRPL9", "RPL30"), ]
dats2[c("MX1", "CTSL"), ]

# pos
nrow(correlation_table[correlation_table$correlation >= 0.95, ])
# [1] 994420
nrow(correlation_table[correlation_table$correlation >= 0.98, ])
# 21,538
nrow(correlation_table[correlation_table$correlation >= 0.99, ])
# 473

# nega
nrow(correlation_table[correlation_table$correlation <= -0.99, ])

num <- 0.99

dfaux <- correlation_table[correlation_table$correlation >= num, ]
dfaux2 <- correlation_table[correlation_table$correlation <= -num, ]

dftosave <- rbind(dfaux, dfaux2)

write.table(dftosave, quote = F, file = "/data/network/nas01/fivibio-data/projects/covid/results/1-DE/endometrium_corr_filtered.tsv", sep = "\t",
            row.names = F, col.names = T)

# ATT

detosave <- dats2[,c("adj.P.Val", "FC")]
detosave$genes <- rownames(detosave)

detosave$adj_true <- ifelse(detosave$adj.P.Val <= 0.05, 1, 0)

write.table(detosave, quote = F, file = "/home/ihenarejos/workspace/projects/covid19/results/1_DE_covid/att_correlation.tsv", sep = "\t",
            row.names = F, col.names = T)

# intersect ERA GENES -----------------------------------------------------

era <- read.delim("RIF_signature_4 - DiazGimeno_signature.tsv", header = T,
           sep = "\t" )

genes = unique(c(correlation_table$gene1, correlation_table$gene2))
era_in <- intersect(era$Gene.symbol, correlation_table$gene1)
era_in_2 <- intersect(era$Gene.symbol, correlation_table$gene2)
# same, 187 genes from 238
era_in_3 = intersect(genes, era$Gene.symbol)
length(era_in_3)

# aux <- correlation_table[correlation_table$gene1 %in% era_in, ]
# aux2 <- correlation_table[correlation_table$gene2 %in% era_in, ]
# aux3 <- rbind(aux,aux2)
# aux3 <- aux3[!duplicated(aux3)]

aux4 = correlation_table[(correlation_table$gene1 %in% era_in_3 & correlation_table$gene2 %in% era_in_3), ]

hist(aux4$correlation)

# same tendency in PF,ESE,MSE except LSE