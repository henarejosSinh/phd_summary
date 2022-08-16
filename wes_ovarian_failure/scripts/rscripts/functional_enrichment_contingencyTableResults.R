###########Descripción del script#######################
#18 de octubre 2019

#appy functional enrichment by ORA to results from fisher tests (variants)


#clean global enviroment if needed
rm(list=ls())

########## USO DE LIBRERIAS ####################################

library(clusterProfiler)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(biomaRt)
library(gsrmUtils)
#if needed to load functions from other scripts
# source("/home/ihenarejos/workspace/proyectos/pof/scripts/functions.R")
source("/home/ihenarejos/workspace/scripts/rscripts/functional_annotation_GO.R")
source("/home/ihenarejos/workspace/scripts/rscripts/functional_annotation_KEGG.R")

########### save_load enviroment ----------------------------------------------------

# save.image("workspace/proyectos/pof/scripts/rdatas/enrichmentCaseControl.Rda")
# load("workspace/proyectos/pof/scripts/rdatas/enrichmentCaseControl.Rda")


# load variants tsv ------------------------------------------------------

# load('workspace/scripts/rdatas/KEGG_annotation_January2019.RData') #for annot table and more
# load("workspace/scripts/rdatas/GO_annotation_January2019.RData")

exome_genes<-fread("workspace/proyectos/pof/lists_genes/all_exome.txt", header = F)
exome_genes<-as.vector(exome_genes$V1)
length(exome_genes)

allgg<-fread("workspace/projects/pof/results/09_functionalEnrichment/all_genes_variants_protein_coding.txt")
tsvCC<-fread("workspace/proyectos/pof/results/09_functionalEnrichment/poi_2395_variants_genes.txt")

allgg_genes<-as.vector(allgg$genes)  #bg genes
genes<-as.vector(tsvCC$genes) #selected genes

poigenes<-read.table("workspace/proyectos/pof/lists_genes/genespoi611.tsv")
length(intersect(genes,poigenes$V1))


# checking gene numbers in our data ---------------------------------------

# intersect exome genes (protein coding) with genes in our variants (protein coding)

mart_proteinC_genes<-fread("workspace/proyectos/pof/lists_genes/gene_regions_grch37_no_patches_protein_coding.tsv")

hgnc_mart<-mart_proteinC_genes$`HGNC symbol`

hgnc_mart<-hgnc_mart[hgnc_mart != ""]
length(hgnc_mart)

length(intersect(exome_genes,hgnc_mart))

length(intersect(allgg$genes,hgnc_mart))

length(intersect(exome_genes,allgg$genes))

genes_poi<-read.table("workspace/proyectos/pof/lists_genes/genespoi611.tsv")

#check how many genes in our variants we have from poi
length(intersect(allgg$genes,genes_poi$V1))

#check how many genes from poi list we found in vcf file created with restrictive filters
df_vcf_poi<-fread("workspace/proyectos/pof/results/09_functionalEnrichment/5%cases__poi_genes.txt")

length(intersect(df_vcf_poi$genes,genes_poi$V1))


# using translate function ------------------------------------------------

translation_go<-translate_genes(gene_list = genes,ini_gene_id = "genename",final_gene_id = "goid")
translation_go_all<-translate_genes(gene_list = allgg_genes,ini_gene_id = "genename",final_gene_id = "goid")
# translation_go_all_exome<-translate_genes(gene_list = exome_genes,ini_gene_id = "genename",final_gene_id = "goid")

translation_kegg<-translate_genes(gene_list = genes,ini_gene_id = "genename",final_gene_id = "keggid")
translation_kegg_all<-translate_genes(gene_list = allgg_genes,ini_gene_id = "genename",final_gene_id = "keggid")
# translation_kegg_all_exome<-translate_genes(gene_list = exome_genes,ini_gene_id = "genename",final_gene_id = "keggid")

# exome has genes that are not in snpEff annotations (nomenclature error? add those anyways?)

# functional enrichment go ------------------------------------------------

#go annotation works with genename and hugo nomenclature directly

# electronic 
auxgo<-functional_annotation_GO(gene_id = "goid",bg_genes = translation_go_all$final, 
                                selected_genes = translation_go$final, ontology = "BP", experimental = F, propagated = F)
aux_go_electronic<-functional_enrichment(auxgo)

# propagated and experimetal (min 50 max 200)
auxgo<-functional_annotation_GO(gene_id = "goid",bg_genes = translation_go_all$final, 
                                selected_genes = translation_go$final, ontology = "BP", experimental = T, propagated = T, min_genes =50, max_genes = 200 )
auxgo_propagated<-functional_enrichment(auxgo)

#only experimental
auxgo<-functional_annotation_GO(gene_id = "goid",bg_genes = translation_go_all$final, 
                                selected_genes = translation_go$final, ontology = "BP", experimental = T, propagated = F, min_genes =1, max_genes = 1000 )
auxgo_experimental<-functional_enrichment(auxgo)

(plotgo<-plot_funcEnrichmentResults(auxgo_propagated,threshold = 0.05,onlySigs = T,maxN = 150, title = "Functional Enrichment Gene Ontology BP"))
(plotgo$general_plot + default_theme())
(plotgo$general_plot + theme_light()+ theme(axis.title = element_text(size = 18), 
                                          axis.text = element_text(size = 7.5),  plot.title = element_text(size = 22, 
                                        hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 14)))



# functional enrichment kegg ----------------------------------------------

auxkegg<-functional_annotation_KEGG(bg_genes = translation_kegg_all$final,selected_genes = translation_kegg$final , gene_id = "keggid",diseases = F)

auxkegg2<-functional_enrichment(auxkegg)

(keggPlot<-plot_funcEnrichmentResults(auxkegg2,title = "Functional Enrichment Kegg"))

(plotkegg<-keggPlot$general_plot + default_theme() )

# keggPlot<-plotFunctEnrich(auxkegg2,title = "Functional Enrichment Kegg")
# 
# dfconversion<-fread("workspace/proyectos/pof/results/09_functionalEnrichment/hgcn_entrezSelection.txt", stringsAsFactors = F)
# 
# dfconversion$hsaIDs = paste0("hsa:", dfconversion$`NCBI gene ID`)
# 
# auxkegg2$genename<-as.character(apply(auxkegg2,1,function(row){
#   id_value<-row[["geneID"]]
#   genes<-""
#   id_values<-(strsplit(id_value,split = ",")[[1]])
#   for (hsa in id_values){
#     genes<-paste0(genes,(dfconversion[ dfconversion$hsaIDs %in% hsa]$`HGNC symbol`)," ")
#   } 
#   return(genes)
# }
# ))


# cluster profiler --------------------------------------------------------

#  these are all the annots for the entire database used (GO, BP, propagated and experimental, 50-200)
term2gene <- auxgo$TERM2GENE
term2name <- auxgo$TERM2NAME

# enricher
aux_enrich <- enricher(gene = translation_go$final, universe = translation_go_all$final,
                       TERM2GENE =term2gene , TERM2NAME = term2name, pAdjustMethod = "fdr", 
                       minGSSize = 50, maxGSSize = 200 )

View(aux_enrich@result)  # compare with View of gsrmutils -> go ok
View(auxgo_propagated[auxgo_propagated$padjust <0.05,])

(plot_dot_enricher <-dotplot(aux_enrich))


# term2gene_go<-sort(term2gene$term, decreasing = T)  required by cluster profiler
aux <- enrichGO(universe = translation_go_all$final,gene = translation_go$final, ont = "BP",
                minGSSize = 50, maxGSSize = 200, pAdjustMethod = "fdr", keytype = "SYMBOL",
                OrgDb ="org.Hs.eg.db" ) 


plot_cluster <- plotGOgraph(aux)
plot_network <- enrichMap( vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai, aux)
plot_genes <- cnetplot(aux)
(plot_dot <- dotplot(aux))

results_cluster <- aux@result

# compare results ??
compareCluster(results_cluster$ID, auxgo_propagated$ID)

# intersect
length(setdiff(results_cluster$ID, auxgo_propagated$ID))

# kegg

kegg_all_cp <- gsub(translation_kegg_all$final, pattern = "hsa:", replacement = "")
kegg_cp <- gsub(translation_kegg$final, pattern = "hsa:", replacement = "")


# enricher 

term2gene <- auxkegg$TERM2GENE
term2name <- auxkegg$TERM2NAME

aux_enrich_kegg <-  enricher(gene = translation_kegg$final, universe = translation_kegg_all$final, 
                             TERM2GENE =term2gene , TERM2NAME = term2name, pAdjustMethod = "fdr" )

View(aux_enrich_kegg@result[aux_enrich_kegg@result$p.adjust <= 0.05 ,])
View(auxkegg2[auxkegg2$padjust <= 0.05,])

plot_dot_enricher_kegg <- dotplot(aux_enrich_kegg)
# enrichkegg
aux_2 <- enrichKEGG(universe = kegg_all_cp, gene = kegg_cp , organism = "hsa" , pAdjustMethod = "fdr")

plot_network <- enrichMap( vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai, aux_2)
plot_genes <- cnetplot(aux_2)
(plot_dot <- dotplot(aux_2))

results_cluster_2 <- aux_2@result

