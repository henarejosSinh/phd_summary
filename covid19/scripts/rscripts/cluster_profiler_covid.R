############ Description #######################################################
# 16 feb 2021

# clusterprofiler COVID Modules

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# BiocManager::install("clusterProfiler")
# 
# R.version
library(ggplot2)
library(clusterProfiler)
library(data.table)
library(reshape2)
library(VennDiagram)
library(RColorBrewer)

# tosav = read.csv("results/5 - network_study/Mid_Secretory_phase/mse_cytohubba_res.csv", header = T)
# write.table(tosav, file = "results/5 - network_study/Mid_Secretory_phase/mse_cytohubba_res.tsv", sep = "\t",
            # col.names = T, row.names = F, quote = F)
load("scripts/rdata/go_bp.Rda") # also mf, cc

# functions ---------------------------------------------------------------
enrich_barplot_gsea = function(phase, ontology, n_modules, minGSSize, maxGSSize){
  load(paste0("/home/sinh/workspace/local/covid19/scripts/rdata/", ontology, "_data.Rda"))  # load ontology files
  
  set.seed(1234)
  out = vector("list", n_modules)
  
  for (i in 1:n_modules){
    dir = "/home/sinh/workspace/local/covid19/results/6 - modules/"
    file = paste0(dir, phase, "_network_module", i, ".csv")
    
    print(file)
    df = read.csv(file, header = T)
    df = df[df$selected == "true",]  # select subset
    print("number of nodes")
    print(nrow(df))
    
    # we want the log2 fold change 
    original_gene_list <-df$logFC_CACO
    # name the vector
    names(original_gene_list) <- df$name
    head(original_gene_list)
    # omit any NA values 
    gene_list<-na.omit(original_gene_list)
    
    # sort the list in decreasing order (required for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)
    gene_list
    
    # from string to object
    # After GOeval
    TERM2GENE = eval(parse(text=paste0("aux_", ontology, "_corr")))
    TERM2NAME= eval(parse(text=paste0("aux_", ontology)))
    
    # originals
    # TERM2GENE = eval(parse(text=paste0(ontology, "_corr")))
    # TERM2NAME= eval(parse(text=paste0(ontology, "_ont")))
    
    y = GSEA(gene_list, TERM2GENE= TERM2GENE, TERM2NAME= TERM2NAME, minGSSize = minGSSize, maxGSSize = maxGSSize, pAdjustMethod =  "fdr", pvalueCutoff = 0.05, seed = T, nPerm = 10000)
    write.table(y@result, file = paste0(dir, ontology, "_module" , i, ".tsv"), sep = "\t", quote = F)
    # print(y@result)
    # View(y@result)
    if (nrow(y@result) == 0){
      print("no res")
      next
      }
    
    # PLOTS
    # cnetplot = cnetplot(y, foldChange=gene_list)
    # heatplot = heatplot(y, foldChange=gene_list)
    # emapplot = emapplot(y)
    # 
    # # dotplot
    # require(DOSE)
    # dotplot = dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign)
    
    # barplot enrich results
    
    print("Number of enriched functions")
    print(nrow(y@result))
    
    # NES CHECK
    print("NES RESULTS")
    print(sum(y@result$NES > 0))
    print(sum(y@result$NES < 0))
    
    limit_rows = ifelse(nrow(y@result) > 20, 20, nrow(y@result))
    
    enrich_table = y@result[1:limit_rows,]
    enrich_table$sig = ifelse(enrich_table$p.adjust <= 0.05, "sig", "no_sig")
    enrich_table$updown = ifelse(enrich_table$NES > 0, "up", "down")
    
    enrich_table$NES <- unlist(lapply(enrich_table$NES, function(x){
      return(round(x,2))
    }))
    
    # General plot - Enriched functions barplot (x = -10*log10(p-value))
    toplot = data.frame(term = paste0(enrich_table$Description, " (", enrich_table$ID, ")"),
                        x = -10*log10(enrich_table$p.adjust),
                        sig = enrich_table$sig,
                        NES = enrich_table$NES,
                        stringsAsFactors = F)
    toplot$term = factor(toplot$term, levels = rev(toplot$term))
    toplot$term
    
    
    # only sigs
    type = enrich_table$updown
    print(type)
    type
    toplot = toplot[toplot$sig == "sig", ]
    p1 = ggplot(toplot, aes(x = x, y = term, fill = type, label = NES))+
      geom_bar(stat = "identity") # fill = black
    p1 = p1 + geom_vline(xintercept =-10*log10(0.05), linetype = 2, color = "#D55E00")
    p1 = p1 + xlab("-10*log10(p-value)") + ylab("Functions") + ggtitle(paste0(phase," Module ", i))+
      theme_light()+
      theme(legend.position="none", plot.title = element_text(size = 22, hjust = 0.5),
            axis.text = element_text(size = 14), axis.title = element_text(size = 18))
    
    if (all(type == "up")){
    p1 = p1 + scale_fill_manual(values = c(down = "royalblue"))
    }
    if (all(type == "down")){
    p1 = p1 + scale_fill_manual(values = c(down = "royalblue"))
    }
    else {
    p1 = p1 + scale_fill_manual(values = c(up = "firebrick", down = "royalblue"))
    }
    p1 = p1 +   geom_text(aes(label=NES),position="stack",vjust=.5, hjust = -0.1, size = 5)
    
    p1
    
    ggsave(paste0("results/6 - modules/barplot_module_", i, ".png"), plot = p1, width = 14, height = 10, units = "cm")
    file_name = paste0("results/6 - modules/barplot_module_", i, ".tiff")
    tiff(file_name, width=600, height=500)
    print(p1)
    dev.off()
    # out[[paste0("module_",i)]] = list(barplot = p1, dotplot = dotplot, cnetplot = cnetplot, heatplot = heatplot,
    #                    emapplot = emapplot)
  }
  return(out)
}

enrich_barplot_ora = function(enrichment_results, threshold = 0.05, onlySigs = T, maxN = 20, general_color = "black", DEG_results = NULL, up_color = "#D55E00", down_color = "#0072B2", title = ""){
  # Add a columnn indicating if it is significant or not
  enrich_table = enrichment_results[1:maxN,]
  enrich_table$sig = ifelse(enrich_table$padjust<=threshold, "sig", "no_sig")
  
  # General plot - Enriched functions barplot (x = -10*log10(p-value))
  toplot = data.frame(term = paste0(enrich_table$Name, " (", enrich_table$ID, ")"),
                      x = -10*log10(enrich_table$padjust),
                      sig = enrich_table$sig,
                      stringsAsFactors = F)
  toplot$term = factor(toplot$term, levels = rev(toplot$term))
  
  if (onlySigs){
    toplot = toplot[toplot$sig == "sig", ]
    p1 = ggplot(toplot, aes(x = term, y = x))+
      geom_bar(stat = "identity", fill = general_color)
  }else{
    color2_rgb = col2rgb(general_color)
    sec_color = rgb(color2_rgb[1], color2_rgb[2], color2_rgb[3], alpha = 0.5)
    p1 = ggplot(toplot, aes(x = term, y = x, fill = sig))+
      geom_bar(stat = "identity")+
      scale_fill_manual(values = c(sig = general_color, no_sig = sec_color))
  }
  p1 = p1 +
    geom_hline(yintercept = -10*log10(threshold), linetype = 2, color = "#D55E00")+
    coord_flip()+
    xlab("") + ylab("-10*log10(p-value)") + ggtitle(title)+
    theme_light()+
    theme(legend.position="none", plot.title = element_text(size = 20, hjust = 0.5))
  
  # DEG plot - Enriched functions barplot (x = # of UP or DOWN genes)
  if (!is.null(DEG_results)){
    toplot2 = data.frame(term = paste0(enrich_table$Name, " (", enrich_table$ID, ")"),
                         stringsAsFactors = F)
    toplot2$term = factor(toplot2$term, levels = rev(toplot2$term))
    genes = lapply(enrich_table$geneID, function(x){
      unlist(strsplit(x, split = ", ", fixed = T))
    })
    n_UPs = unlist(lapply(genes, function(x){
      sum(DEG_results[x, "FC"]>=0)
    }))
    n_DOWNs = unlist(lapply(genes, function(x){
      sum(DEG_results[x, "FC"]<0)
    }))
    toplot2 = cbind(rbind(toplot2, toplot2),
                    n_genes = c(n_UPs, -n_DOWNs),
                    type = c(rep("UP", length(n_UPs)), rep("DOWN", length(n_DOWNs))))
    limits = c(min(toplot2$n_genes[toplot2$type=="DOWN"]),
               max(toplot2$n_genes[toplot2$type=="UP"]))
    limits[1] = limits[1] - limits[1]%%10
    limits[2] = limits[2] + (10-limits[2]%%10)
    p2 = ggplot(toplot2, aes(x = term, y = n_genes, fill = type))+
      geom_bar(stat = "identity")+
      coord_flip()+
      scale_fill_manual(values = c(UP = up_color, DOWN = down_color))+
      scale_y_continuous(labels = c(seq(abs(limits[1]), 0, -10), seq(10, limits[2], 10)),
                         breaks = c(seq(abs(limits[1]), 0, -10), seq(10, limits[2], 10)))+
      xlab("") + ylab("Number of genes") + ggtitle(title)+
      theme_light()+
      theme(legend.position = "none")
    
  }else{
    p2 = NULL
  }
  
  results = list(general_plot = p1,
                 DEGs_plot = p2)
  
  return(results)
  
}

ora_net_modules = function(phase, module,ontology){
  
  bggenes = read.table(paste0("results/5 - network_study/", phase,"_network_genes.txt"),
                       stringsAsFactors = F)[[1]]
  
  # MODULE GENES
  df = read.csv(paste0("results/6 - modules/", phase,"_network_module", module,".csv"), header=TRUE, stringsAsFactors = F)
  df = df[df$selected == "true",]
  
  load(paste0("scripts/rdata/go_", ontology, ".Rda"))
  
  
  ORA_BP = enricher(gene=df$name,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    universe = bggenes,
                    minGSSize = 5,
                    maxGSSize = 500,
                    TERM2GENE = GOannot,
                    TERM2NAME = GO_id_name)
  
  # View(ORA_BP@result[ORA_BP@result$p.adjust <= 0.05,])
  
  
  enrich_table = ORA_BP@result[ORA_BP@result$p.adjust <= 0.05,]
  
  if (nrow(enrich_table) == 0) {
    print("no enrichment")
    break}
  
  colnames(enrich_table)[6] = "padjust"
  colnames(enrich_table)[2] = "Name"
  
  enrich_barplot_ora(enrichment_results = enrich_table, maxN = ifelse(nrow(enrich_table) >= 20, 20, nrow(enrich_table)), title = paste0("ORA ", toupper(ontology)," ",toupper(phase)," ", module ," module"))
}

ora_net_modules_tables = function(phase, module,ontology){
  
  bggenes = read.table(paste0("results/5 - network_study/", phase,"_network_genes.txt"),
                       stringsAsFactors = F)[[1]]
  
  # MODULE GENES
  df = read.csv(paste0("results/6 - modules/", phase,"_network_module", module,".csv"), header=TRUE, stringsAsFactors = F)
  df = df[df$selected == "true",]
  
  load(paste0("scripts/rdata/go_", ontology, ".Rda"))
  
  
  ORA_BP = enricher(gene=df$name,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    universe = bggenes,
                    minGSSize = 5,
                    maxGSSize = 500,
                    TERM2GENE = GOannot,
                    TERM2NAME = GO_id_name)
  
  # View(ORA_BP@result[ORA_BP@result$p.adjust <= 0.05,])
  
  
  enrich_table = ORA_BP@result[ORA_BP@result$p.adjust <= 0.05,]
  
  if (nrow(enrich_table) == 0) {
    print("no enrichment")}
  
  colnames(enrich_table)[6] = "padjust"
  colnames(enrich_table)[2] = "Name"
  enrich_table
}

# launch function ---------------------------------------------------------

# res = enrich_barplot(phase = "mse", ontology = "bp", n_modules = 9, 5, 200)
# 
# View(res)
# res$module_1$barplot
# res$module_2$barplot
# res$module_3$barplot
# res$module_4$barplot

res = list()
for (i in 1:2){
res[[i]] = ora_net_modules(phase = "mse", ontology = "cc", module = i)
print(res$general_plot)
}

for (i in 1:2){
  print(res[[i]])
  Sys.sleep(2)
}

ontology = "bp"

res = list()
for (i in 1:9){
  res[[i]] = ora_net_modules_tables(phase = "mse", ontology = ontology, module = i)
  print(res)
}

str(res)
for (i in 1:9){
  if (nrow(res[[i]]) >= 1){
  write.table(res[[i]], file = paste0("results/6 - modules/mse_ora_module_",ontology, "_", i, ".tsv"), sep = "\t"
              ,quote = F, row.names = F, col.names = T)}
}
View(res[[1]])
View(res[[2]])

# prepare ontology --------------------------------------------------------

# If user have GO annotation data (in data.frame format with first column of gene 
# ID and second column of GO ID), they can use enricher() and gseGO() functions 
# to perform over-representation test and gene set enrichment analysis.
# PD: Cols are inversed internally

# # bp_corr = read.csv("Downloads/bingo/bingo_gene_GO_BP.txtgene_corr.tsv", sep = "\t", header = F, colClasses = "character" )
# bp_corr = read.csv("results/5 - network_study/bingo/BP_gene_corr.tsv", sep = "\t", header = F, colClasses = "character" )
# 
# head(bp_corr)
# colnames(bp_corr) = c("GO", "Gene")
# 
# # If genes are annotated by direction annotation, it should also annotated by 
# # its ancestor GO nodes (indirect annation). If user only has direct annotation, 
# # they can pass their annotation to buildGOmap function, which will infer indirection
# # annotation and generate a data.frame that suitable for both enricher() and gseGO().
# 
# bp_ont = read.csv("Downloads/bingo/bp_ont.tsv", stringsAsFactors = F, colClasses = "character" , sep = "\t", header = F)
# bp_ont = read.csv("results/5 - network_study/bingo/BP_ont.tsv", stringsAsFactors = F, colClasses = "character" , sep = "\t", header = F)
# 
# head(bp_ont)
# colnames(bp_ont) = c("GO","Gene")
# bp_ont$Gene <- factor(bp_ont$Gene)
# 
# # buildgomap function (http://yulab-smu.top/clusterProfiler-book/chapter5.html)
# # provided by a data.frame of GO (column 1) and gene (column 2) direct 
# # annotation this function will building gene to GO and GO to gene mapping, 
# # with directly and undirectly (ancestor GO term) annotation.
# 
# aux_bp = buildGOmap(bp_ont)
# head(aux_bp)
# 
# aux_bp_corr = buildGOmap(bp_corr)
# 
# save.image("scripts/rdata/bp_data.Rda")
# load("scripts/rdata/bp_data.Rda")

# prepare data ------------------------------------------------------------

# tutorial followed:
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# reading in data from deseq2
df = read.csv("results/6 - modules/mse_network_module1.csv", header=TRUE)
df = df[df$selected == "true",]
# df = read.csv("Downloads/mse_network_module2.csv", header=TRUE)
head(df)
# we want the log2 fold change 
original_gene_list <-df$logFC_CACO

# name the vector
names(original_gene_list) <- df$name
head(original_gene_list)
# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list

# all genes network
# write.table("results/5 - network_study/bingo/all_network_genes.txt", x = df$name, quote = F,
#             col.names = F, row.names = F)

# gene syms ---------------------------------------------------------------
gene <- names(gene_list)[abs(gene_list) > 2]
gene
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

### GSEA All trees/ontologies
# rm(gse)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# gse <- gseGO(geneList = gene_list, 
#              ont = "BP", 
#              keyType = "SYMBOL", 
#              nPerm = 1000, 
#              minGSSize = 5, 
#              maxGSSize = 200, 
#              pvalueCutoff = 1, 
#              verbose = TRUE, 
#              OrgDb = organism,   # organism
#              pAdjustMethod = "fdr")
# View(gse@result)
# plotGOgraph(gse)

# GSEA Podded Tree --------------------------------------------------------


# Propagated podded tree
# str(aux_bp_corr)
head(bp_ont)
str(bp_ont)

# gse_custom <- GSEA(geneList = gene_list, 
#                    TERM2GENE = bp_corr,
#                    TERM2NAME = bp_ont,
#                    nPerm = 1000, 
#                    minGSSize = 30, 
#                    maxGSSize = 200, 
#                    pvalueCutoff = 0.05, 
#                    verbose = TRUE,
#                    pAdjustMethod = "fdr")
# gse_custom
# View(gse_custom@result)

# # Enricher DEG
# deg <- names(gene_list)[abs(gene_list) > 1.5]
# deg
# x = enricher(deg, TERM2GENE=bp_corr, TERM2NAME=bp_ont, 
#              minGSSize = 30, maxGSSize = 200, pAdjustMethod = "fdr", pvalueCutoff = 0.05)
# x
# head(summary(x))
# View(x@result[x@result$p.adjust < 0.05, ])
# barplot(x)

# Alternative GSEA
y = GSEA(gene_list, TERM2GENE=aux_bp_corr, TERM2NAME= aux_bp, minGSSize = 5, 
         maxGSSize = 50, pAdjustMethod =  "fdr", pvalueCutoff = 0.05, n = 20000)
if (nrow(y@result) == 0){print("no res")}
head(summary(y))
# gseaplot(y, "0006412")
View(y@result)

# write.table(y@result, file = "results/6 - modules/bp_all.tsv", sep = "\t", quote = F)

cnetplot(y, foldChange=gene_list)
heatplot(y, foldChange=gene_list)
emapplot(y)




# bar plot funt. enrich. --------------------------------------------------

# check UP/DOWN 
# GeneRatio in clusterProfiler::dotplot() is calculated as: count / setSize
# https://support.bioconductor.org/p/130447/

# require(DOSE)
# dotplot(y, showCategory=10, split=".sign") + facet_grid(.~.sign)
# 
# res1 = read.table("Downloads/bingo/gsea_bp_module1.tsv", sep = "\t")
# View(res1)
# sum(res1$NES > 0)
# # 36
# sum(res1$NES < 0)
# # 21


# https://guangchuangyu.github.io/2016/12/dotplot-for-gsea-result/
# GeneRatio = core_enrichment/setsize

nrow(y@result)
enrich_table = y@result[1:12,]
enrich_table$sig = ifelse(enrich_table$p.adjust <= 0.05, "sig", "no_sig")
enrich_table$updown = ifelse(enrich_table$NES > 0, "up", "down")

# setsize/total genes
# colnames(enrich_table)
# enrich_table$generatio = apply(enrich_table, 1, function(x){
#   core_genes = x[["core_enrichment"]] # this is the leading genes
#   setsize = x[["setSize"]]
#   print(core_genes)  
#   core_genes = length(unlist(strsplit(core_genes, split = "/")))
#   res = as.numeric()/as.numeric(setsize)
#   res = round(res, 2)
#   return(res)
# })
# enrich_table$generatio

# NES rounded
enrich_table$NES <- unlist(lapply(enrich_table$NES, function(x){
  print(x)
  return(round(x,2))
}))
enrich_table$NES

# General plot - Enriched functions barplot (x = -10*log10(p-value))
toplot = data.frame(term = paste0(enrich_table$Description, " (", enrich_table$ID, ")"),
                    x = -10*log10(enrich_table$p.adjust),
                    sig = enrich_table$sig,
                    NES = enrich_table$NES,
                    stringsAsFactors = F)
toplot$term = factor(toplot$term, levels = rev(toplot$term))
toplot$term

# only sigs
type = enrich_table$updown
type
toplot = toplot[toplot$sig == "sig", ]
p1 = ggplot(toplot, aes(x = x, y = term, fill = type, label = NES))+
  geom_bar(stat = "identity") # fill = black
p1
p1 = p1 + geom_vline(xintercept =-10*log10(0.05), linetype = 2, color = "#D55E00")
# coord_flip()

p1

p1 = p1 + xlab("-10*log10(p-value)") + ylab("Functions") + ggtitle("MSE Module 3")+
  theme_light()+
  theme(legend.position="none", plot.title = element_text(size = 22, hjust = 0.5),
        axis.text = element_text(size = 12), axis.title = element_text(size = 18))

p1

p1 = p1 + scale_fill_manual(values = c(up = "firebrick", down = "royalblue"))

p1

p1 = p1 +   geom_text(aes(label=NES),position="stack",vjust=.5, hjust = -0.1, size = 5)
p1 

# boxplots topological features -------------------------------------------

# read each phase
phase = read.csv("results/5 - network_study/Proliferative_phase/pf_cytohubba_res.csv")
head(phase)
phase_name = "PF"

par(cex.lab = 1.25)
par(cex.axis = 1.25)

# Degree
box_plot = boxplot(phase$Degree, data=phase, xlab = paste0(phase_name," genes"), ylab = "Degree" )
probs = quantile(phase$Degree, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
probs
probs[[7]][1]
# 0%   25%   50%   75%   90%   95%   99%  100% 
#   1.0   5.0  15.0  38.5  65.0  92.0 144.0 209.0 
sum(phase$Degree >= probs[[7]][1]) # mse 21 genes, ese 23
box_plot$stats[[5]] # maximum
sum(phase$Degree >= box_plot$stats[[5]])

# Betweenness
box_plot = boxplot(phase$Betweenness, data=phase, xlab = paste0(phase_name," genes"), ylab = "Betweenness" )
probs = quantile(phase$Betweenness, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
probs
probs[[7]][1]
# 0%        25%        50%        75%        90%        95%        99%       100% 
#   0.0000   315.5832  2556.3121  7324.9274 14446.0043 20120.7184 34824.4027 65973.4272 
sum(phase$Betweenness >= probs[[7]][1]) # mse 21 genes 23 ese
box_plot$stats[[5]]
sum(phase$Betweenness >= box_plot$stats[[5]])

# Bottleneck
box_plot = boxplot(phase$BottleNeck, data=phase, xlab = paste0(phase_name," genes"), ylab = "Bottleneck" )
probs = quantile(phase$BottleNeck, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
probs
probs[[7]][1]
# 0%        25%        50%        75%        90%        95%        99%       100% 
#   0.0000   315.5832  2556.3121  7324.9274 14446.0043 20120.7184 34824.4027 65973.4272 
sum(phase$BottleNeck >= probs[[7]][1]) # mse 28 genes, ese 24
box_plot$stats[[5]]
sum(phase$BottleNeck >= box_plot$stats[[5]])


# boxplot maximum relatives -----------------------------------------------

par(cex.lab = 1.25)
par(cex.axis = 1.25)

box_deg = boxplot(mse$Degree, data= mse$Degree, xlab = "MSE", ylab = "Degree")
quantile(mse$Degree, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
(box_deg)
# stats [[5]] # max relative
# hubs > 88
# bottleneck > 3
# betweenness > 17826.3536

box_bot = boxplot(mse$BottleNeck, data= mse$BottleNeck, xlab = "MSE", ylab = "BottleNeck")
# quantile(mse$Degree, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
(box_bot$stats)

box_betw = boxplot(mse$Betweenness, data= mse$Betweenness, xlab = "MSE", ylab = "BottleNeck")
# quantile(mse$Degree, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
(box_betw$stats)

sum(mse$Degree >= 88)
sum(mse$BottleNeck >= 3)
sum(mse$Betweenness >= 17826.3536)

# venn topological features --------------------------------------------------------------------

mse = read.csv("results/5 - network_study/Mid_Secretory_phase/mse_cytohubba_res.csv", stringsAsFactors = F)
head(mse)


net = mse
feature = "Degree"
feature = "Betweenness"
feature = "BottleNeck"

res = quantile(net[,feature], probs = 0.99)[[1]][1]
res

degree_g = net$node_name[net$Degree >= res]
degree_g
length(degree_g)
# 21

betw_g = net$node_name[net$Betweenness >= res]
betw_g
# 21

bottle_g = net$node_name[net[feature] >= res]
bottle_g
# 28



# max rel
# 88
# bottleneck > 3
# betweenness > 17826.3536

degree_g = net$node_name[net$Degree >= 88]
degree_g

betw_g = net$node_name[net$Betweenness >= 17826.3536]
betw_g

bottle_g = net$node_name[net[feature] >= 3]
bottle_g


### VENN

myCol <- brewer.pal(3, "Pastel2")
name = "quantile"

venn.diagram(x = list(degree_g,betw_g,bottle_g),
             category.names = c("Degree", "Betweenness", "BottleNecks"),
             filename = paste0("results/5 - network_study/Mid_Secretory_phase/venn_", name), output = T ,
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             col = myCol,
             
             # Numbers
             cex = 2,
             # fontface = "sans",
             fontfamily = "sans",
             
             # col=c("#440154ff", '#21908dff', '#fde725ff'),
             # fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
             
             # Names
             cat.cex = 1.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             cat.col = myCol,
             rotation = 1)
