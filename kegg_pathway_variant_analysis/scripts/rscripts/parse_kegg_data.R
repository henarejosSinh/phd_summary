# 27/02/2020
# patricia.sebastian@ivirma.com

# Script to parse KEGG information to be used to pathways variant analysis in POF project


# Settings ----------------------------------------------------------------

rm(list=ls())

# Set working directory
setwd("/data/network/nas01/fivibio-data/projects/pof/pathway_variant_analysis/")

# Libraries
library(XML)
library(gplots)

# Get biomart information
load("/data/network/nas01/fivibio-data/data/biomart/conversion_to_HGNC/conversion_to_HGNC_2020.01.30.RData", 
     verbose = T)
corres = conversion_to_HGNC$entrezgene_id
corres$entrezgene_id = paste0("hsa:", corres$entrezgene_id)


# KEGG data ---------------------------------------------------------------

# Getting KGML files paths 
mypath = "/data/network/nas01/fivibio-data/data/kegg/keggannot_2020/pathway"
myfiles = dir("/data/network/nas01/fivibio-data/data/kegg/keggannot_2020/pathway")
length(myfiles)

# Selection of genes done by Ismael Henarejos
pathways_selection = read.delim("selection_kegg_pathways.csv", header = T, as.is = T, 
                                colClasses = c("character", "character", "character", "character", 
                                               "logical", "character"))
selected_pathways = paste0("hsa", pathways_selection$id[pathways_selection$variant_analysis])

# Check if all selected nodes have an associated KGML
venn(list(selected_pathways = selected_pathways, kegg_xgmls = myfiles))
setdiff(selected_pathways, myfiles)
# [1] "hsa04624" "hsa04320" "hsa04361"
# They are pathways not available in human!! --> They have to be removed

# Files to be parsed
myfiles = intersect(selected_pathways, myfiles)
length(myfiles)

# Parse pathways information
pathwaysInfo = lapply(myfiles, function(xx){
  
  # cat(xx, "\n")

  x = paste(mypath,"/",xx,sep="")
  
  # Read file
  xmlpath = xmlParse(x)
  

  ##### Entries = ATTS #####
  # Get entries information
  en = getNodeSet(xmlpath,"//pathway/entry")
  
  # Parse entries information
  entries = do.call("rbind",lapply(1:length(en),function(i){
    
    # Type info
    type = xmlGetAttr(en[[i]], "type")
    
    # Only genes and groups are included
    if (!(type %in%c("gene", "group"))){
      return(NULL)
    }
    
    # Network id info
    id = xmlGetAttr(en[[i]], "id")
    
    # Genes included in the node info
    name = xmlGetAttr(en[[i]], "name")
    
    # Get genenames
    if (name !=  "undefined"){ # If it is not a group
      genes = unlist(strsplit(name, split = " "))
      if (all(genes%in%corres$entrezgene_id)){
        genes_names = paste(corres$hgnc_symbol[corres$entrezgene_id%in%genes], collapse = " | ")
      }else{ # If a gene has not correspondence hsa identifier is retained
        genes_names = paste(c(corres$hgnc_symbol[corres$entrezgene_id%in%genes],
                        setdiff(genes, corres$entrezgene_id)), collapse = " | ")
      }
    }else{ # If it is a group
      genes_names = "--" # Will be fixed later
    }
    
    # Add OR information
    name = gsub(" ", " | ", name)
    
    # Name to be displayed
    nodename = xmlGetAttr(getNodeSet(en[[i]], "graphics")[[1]],"name")
    if (is.null(nodename)){
      nodename = "--"# Necessary for groups (nodename = NULL). Will be fixed later
    }
    nodename = unlist(strsplit(nodename, split = ", "))[1]
    
    # x & y coordinates
    xx = xmlGetAttr(getNodeSet(en[[i]], "graphics")[[1]],"x")
    yy = xmlGetAttr(getNodeSet(en[[i]], "graphics")[[1]], "y")
    
    # Retain groups components in the name. Will be fixed later
    if (type == "group"){
      comps = getNodeSet(en[[i]],"component")
      name = paste(unlist(lapply(1:length(comps),function(j){
        xmlGetAttr(comps[[j]],"id")
      })),collapse=" & ")
    }
    
    # data.frame with node information
    res = data.frame(id=id, genes=name, type=type, genes_names = genes_names,
               nodename=nodename, x=xx, y=yy, stringsAsFactors = F)
    return(res)
    
  }))
  
  # Fix ... in nodename
  entries$nodename = gsub("...", "", entries$nodename, fixed = T)
  
  # Nodes finally included in the pathway network
  nodes = unique(entries$id)



  ##### Relations = SIF #####
  
  # Get relations information
  rels = getNodeSet(xmlpath,"//pathway/relation")
  
  # Parse relations information
  relations = do.call("rbind",lapply(1:length(rels),function(i){
    
    # Get entry1 & entry2 information
    entry1 = xmlGetAttr(rels[[i]],"entry1")
    entry2 = xmlGetAttr(rels[[i]], "entry2")
    
    # Keep only relations where entry1 and entry2 arre included in nodes in relations
    if (!all(c(entry1, entry2) %in% nodes)){
      return(NULL)
    }
      
    # Get type information
    type = xmlGetAttr(rels[[i]], "type")
    
    # Get subtype information
    subt = getNodeSet(rels[[i]],"subtype")
    if (length(subt) == 0){
      subtype_name = subtype_value = "Unknown" # Relations without subtype (to be checked manually)
    }else{
      subtype_name = paste(unlist(lapply(1:length(subt), function(j){
        xmlGetAttr(subt[[j]],"name")
      })),collapse=", ")
    }
    
    # data.frame with edge information
    res = data.frame(entry1=entry1,entry2=entry2,type=type,subtype=subtype_name, stringsAsFactors = F)
    return(res)

  }))
  
  # There are relations but not associated with any node in entries (genes or groups)
  if (is.null(relations)){
    cat("Pathway", xx, "has not relations associated with nodes defined in entries (genes or groups)\n")
    return(NULL)
  }
  
  ##### Fix obtained information #####
  
  # Fix binding/association as a group
  last_entry = max(as.numeric(nodes))
  ba_indexes = grep("binding/association", relations$subtype)
  entries_toadd = do.call("rbind", lapply(ba_indexes, function(i){
    last_entry <<- last_entry+1
    res = data.frame(id = last_entry, 
                     genes = paste(relations$entry1[i], relations$entry2[i], sep = " & "),
                     type = "group",
                     genes_names = "--",
                     nodename = "--",
                     x = entries$x[entries$id == relations$entry1[i]],
                     y = entries$y[entries$id == relations$entry1[i]],
                     stringsAsFactors = F)
    return(res)
  }))
  entries = rbind(entries, entries_toadd)
  if(length(ba_indexes)>0){
    relations = relations[-ba_indexes, ]
  }

  # Fix groups information (genes, nodename and edges)
  gr_indexes = which(entries$type == "group")
  nodes2remove = NULL
  for (i in gr_indexes){
    
    # Modify entries information
    gr_nodes = unlist(strsplit(entries$genes[i], split = " & "))  
    entries$genes[i] = paste(paste0("(", entries$genes[entries$id%in%gr_nodes], ")"), collapse = " & ")
    entries$genes_names[i] = paste(paste0("(", entries$genes_names[entries$id%in%gr_nodes], ")"), collapse = " & ")
    entries$nodename[i] = paste(entries$nodename[entries$id%in%gr_nodes], collapse = " & ")
    
    # Modify relations information
    relations$entry1[relations$entry1%in%gr_nodes] = entries$id[i]
    relations$entry2[relations$entry2%in%gr_nodes] = entries$id[i]
    
    # Add to the list of nodes to remove (nodes that are included in a group)
    nodes2remove = c(nodes2remove, gr_nodes)
    
  }
  # # Remove entries associated to a group
  # entries = entries[!entries$id%in%nodes2remove, ]
  
  # Remove from entries nodes not in relations
  nodes_rels = unique(c(relations$entry1, relations$entry2))
  entries = entries[entries$id %in% nodes_rels,]
  
  # Remove self-loops
  relations = relations[relations$entry1 != relations$entry2, ]
  
  # Return a list with entries (ATTS) and relations (SIF) information
  res = list(SIF = relations,
             ATTS = entries)
  return(res)
  
})
names(pathwaysInfo) = myfiles

# Remove pathways witout relations
pathwaysInfo = pathwaysInfo[unlist(lapply(pathwaysInfo, function(x){!is.null(x)}))]

# FALTA DECIDIR SI QUITAMOS LOS QUE TIENEN POCAS RELACIONES O NO
# unlist(lapply(pathwaysInfo, function(x){nrow(x$SIF)})) # evaluar número de relaciones
# unlist(lapply(pathwaysInfo, function(x){nrow(x$ATTS)})) # evaluar número de nodos


# Save information --------------------------------------------------------

# RData
save(pathwaysInfo, file = "pathways_information.RData")


lapply(names(pathwaysInfo), function(x){
  print(paste0("SIFs/", x, "_SIF.txt"))
  print(head(pathwaysInfo[[x]]$SIF))
  write.table(pathwaysInfo[[x]]$SIF, file = paste0(getwd(),"/SIFs/", x, "_SIF.txt"), 
              col.names = T, row.names = F, 
              quote = F, sep = "\t")

  write.table(pathwaysInfo[[x]]$ATTS, file = paste0(getwd(),"/ATTs/", x, "_ATT.txt"), 
              col.names = T, row.names = F, 
              quote = F, sep = "\t")
})

