########### Description #######################################################
# 25 feb 2020

# DDI correspondences
# Ismael Henarejos Castillo
# ihc.europa@gmail.com
# ismael.henarejos@ivirma.com

#working directory
# dir=('/home/ihenarejos/workspace/')
# setwd(dir)
# setwd("/home/ihenarejos/workspace/")

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)  # for color palettes  
library(gsrmUtils)
library(readxl)  # read excel sheets
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


########## Enviroment ##########################################################

dir <- "workspace/projects/ddi/data/dbs/"

# save.image("scripts/rdata/ddi_03_08_2020.Rda")
load("scripts/rdata/ddi_25_09_2020.Rda")

########## Functions ####################################

myTile_theme <- (theme_light( base_size = 30) + 
                   theme(legend.title = element_blank(), 
                         axis.text.x = element_text(angle = 0,size = rel(0.8), vjust = 2) ,
                         axis.text.y = element_text(size = rel(0.8)) ,
                         axis.text.x.top = element_text(),
                         axis.title.x = element_text("drug", size = rel(1) ),
                         axis.title.y = element_text("guideline", size = rel(1)),
                         legend.text  = element_blank(),
                         legend.position = "none",
                         strip.background.x = element_blank(),
                         strip.background.y = element_blank(),
                         # axis.ticks.x = element_blank(),
                         # axis.ticks.y = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor  = element_blank(),
                         axis.line.x = element_blank(),
                         axis.line.y = element_blank()
                   ) )

# load database -----------------------------------------------------------

stitch <- fread("workspace/projects/ddi/data/dbs/drugbank.aliases.matches.tsv",
                header = F,
                col.names = c("CIDm","CIDs","synonym","db"),
                stringsAsFactors = F)

sider_atc <- read.delim("workspace/projects/ddi/data/dbs/sider_atc.tsv", 
                        header = F,
                        col.names = c("cid","atc"),
                        stringsAsFactors = F)
sider_names <- read.delim("workspace/projects/ddi/data/dbs/sider_names.tsv", 
                          header = F, col.names = c("cid","alias"),
                          stringsAsFactors = F)
sider_all <-
  read.delim("data/dbs/sider_all_se.tsv", header = F,
             col.names = c("cid","cid2","ef_id","type","ef_id2","description"),
             stringsAsFactors = F)

# sider_all <-
#   read.delim("data/dbs/sider_all_label_with_ref.tsv", header = F,
#              col.names = c("ref","cid","cid2","ef_id","type","ef_id2","description"),
#              stringsAsFactors = F)

# filter by Preferred Term on sider
sider_all <- sider_all[sider_all$type == "PT", ]
sider_all <- sider_all[c("cid","cid2", "description")]

write.table("results/06_vilar_ADE/sider_pt_ddi_ref.tsv", sep = "\t",
            quote = F, row.names = F, col.names = F, x = sider_all)

all_effects <- unique(sider_all$description)

write.table("results/06_vilar_ADE/sider_effects.tsv", sep = "\t",
            quote = F, row.names = F, col.names = F, x = all_effects)

# add atc for both drugs of sider interactions
sider_to_map <- merge(sider_all,sider_atc, by.x = c("cid"))

# add name 
sider_to_map2 <- merge(sider_all, sider_names, by.x = c("cid"))

# transform sider to a list of lists where each element is the side effects
# of an atc id (a drug)
by_effect <- split(sider_to_map[,c("atc","description","cid","cid2")], f = sider_to_map$description)
by_atc <- split(sider_to_map[,c("atc","description","type","cid","cid2")], f = sider_to_map$atc)
by_names <- split(sider_to_map2[,c("alias","description","type","cid","cid2")], f = sider_to_map2$alias)

all(by_atc$B05CA02$description == by_atc$A01AB03$description) # check if side
# effects are the same for different ATC from same drugs


# map sider to drugbank ---------------------------------------------------------
library(data.table)
drugbank <- fread("data/dbs/drugbank_info_2020_07.txt", header = T)

# drugbank <- drugbank[grepl(drugbank$groups, pattern = "approved") , ]
# get a new col in drugbank with child atc:
drugbank$atc <- (apply(drugbank, 1, function(row){
  atcs <- row["atcCodes"]  # get all atcs
  value <- strsplit(x = atcs, split = ";")[[1]]  # first delimiter
  value2 <- strsplit(x = value, split = ",") # second delimiter, will
  # return a list of lists in case there's more than one atc for that drug
  vector <- list()
  for (i in value2) {  # we need to iterate for each possible atc
    child <- i[1]
    vector <- append(vector,child)
  }
  # print(value2)
  # print(list_atc)
  return(vector)
}))
names(drugbank)
head(drugbank$atc)
drugbank = drugbank[grep("\\bapproved\\b", drugbank$groups),]
drugbank_set = drugbank[drugbank$smiles != "" | drugbank$targets != "" | 
                          drugbank$enzymes != "" | drugbank$transporters != ""
                        | drugbank$carriers != "",]

drugbank_set$atc <- vapply(drugbank_set$atc, paste, collapse = " ", character(1L))
drugbank_set <-drugbank_set[c("#ID", "atc")]
drug_and_atc <-drugbank_set[,c("#ID","atc"),]

drug_and_atc
View(drug_and_atc)
write.table(x = drug_and_atc, quote = F, sep = " ",
            row.names = F, col.names = F, file = "results/06_vilar_ADE/drugbank_id_atc_2.tsv")

new_names <- lapply(names(by_atc), function(n) {
  # print(n)  will print code atc
  
  pos_res <- lapply(1:nrow(drugbank), function(pos){
    
    atc_list <- unlist(drugbank$atc[pos])
    if (n %in% atc_list) {
      return(pos)
    }else{
      return(NA)
    }
  })
  
  pos_res <- unlist(pos_res[!is.na(pos_res)])
  return(drugbank[pos_res, ]$name)
})

# modify old names with new names
drugs_sides <- by_atc
names(drugs_sides) <- new_names
drugs_sides$Chlorhexidine
rm(by_atc)

length(unique(plyr::compact(new_names))) # 1146


# map sider in drugbank synonyms (worse)---------------------------------------------

# proof of concept
aliases <- tolower(unlist(strsplit(drugbank$alias, split = ";", perl = T)))

names(by_names) <- tolower(names(by_names))
new_names2 <- lapply(names(by_names), function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(drugbank), function(index){
    
    syns_list <- tolower(unlist(strsplit(drugbank$alias[index], 
                                         split = ";", perl = T)))
    if (n %in% syns_list) {
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  return(drugbank[index_res, ]$name)
  
  # pos_res <- lapply(1:nrow(drugbank), function(pos){
  # atc_list <- unlist(drugbank$atc[pos])
}
)

names(by_names) <- new_names2
library(plyr)
length(unique(plyr::compact(new_names2))) # 1023

# sider>drugbank>twosides (by drug name) ----------------------------

# twosides 2012 with pvalues

twosides <- fread("workspace/projects/ddi/data/dbs/twosides_subset_confidence.tsv",
                  header = T)

vector_names <- tolower(names(drugs_sides))
names(twosides)
twosides$by_name <- unlist(apply(twosides,1,function(row){
  d1 <- row["drug1"]
  d2 <- row["drug2"]
  
  check_drug <- F
  
  if (tolower(d1) %in% vector_names && 
      tolower(d2) %in% vector_names) {
    check_drug <- TRUE
    # print("hola")
  }
  if (check_drug) {
    # relation is ok as both drug exits in sider/drugbank
    return(T)} else {
      return(F)
    }
}
)
)

# which are true rel
twosides_by_name <- twosides[twosides$by_name == T, ]
# 3,678,309 from 4,650,...
length(unique(union(twosides$drug2,twosides$drug1)))
# 645 drugs
length(unique(union(twosides_by_name$drug2,twosides_by_name$drug1)))
# 548 drugs (~100 less mapping with names )

# sider>drugbank>twosides (by cid) ----------------------------------------

# we need to convert fisrt column CID1 to CID0 

list_df <- lapply(drugs_sides, function(drug){
  # # drug is the name of the df
  
  cid_mod <-
    gsub(drug$cid[1], # must change 1 to 0 after CID
         pattern = "([CID]+)([1])", perl = T, replacement = "CID0" )
  cid_mod <-
    gsub(cid_mod, 
         pattern = "([CID]+)([0]+)", perl = T, replacement = "" )
  # cid2
  cid2_mod <-
    gsub(drug$cid2[1], 
         pattern = "([CID]+)([0]+)", perl = T, replacement = "" )
  df <- data.frame(cid_mod = cid_mod, cid2_mod = cid2_mod)
  # add as new col
  return(df)
})

cids <- rbindlist(list_df)
cids <- unique(cids)

twosides$stitch_id1_mod <- unlist(apply(twosides,1,function(row){
  s1 <- row["stitch_id1"]
  
  value <- gsub(s1, 
                pattern = "([CID]+)([0]+)", perl = T, replacement = "" )
  
  return(value)
}
)  
)

twosides$stitch_id2_mod <- unlist(apply(twosides,1,function(row){
  s2 <- row["stitch_id2"]
  
  value <- gsub(s2, 
                pattern = "([CID]+)([0]+)", perl = T, replacement = "" )
  
  return(value)
}
)  
)

twosides$by_cid <- unlist(apply(twosides,1, function(row){
  s1 <- row["stitch_id1_mod"]
  s2 <- row["stitch_id2_mod"]
  
  found_s1 <- F
  found_s2 <- F
  # print(s1)
  if (s1 %in% cids$cid_mod |
      s1 %in% cids$cid_mod2 )  {
    found_s1 <- T
  }
  if (s2 %in% cids$cid_mod |
      s2 %in% cids$cid_mod2 )  {
    found_s2 <- T
  }
  # print(found_s1, found_s2)
  if (found_s1 & found_s2) {
    # relation is ok as both drugs exits in sider/drugbank
    return(T)} else {
      return(F)
    }
}
)
)

twosides_by_cid <- twosides[twosides$by_cid == T, ]
# 3,407,078 from 4,650,...

length(unique(union(twosides_by_cid$drug2,twosides_by_name$drug1)))
# 611 drugs from 645

# map infertily drugs to dbs ----------------------------------------------

infertile <- fread("workspace/projects/ddi/data/infertily_drugs/drugs_28_02_2020.txt", header = F)
infertile <- tolower(infertile$V1)

# sider names 15/18
length(infertile[infertile %in% tolower(unique(names(drugs_sides)))])

# twosides cid 6/18
length(infertile[infertile %in% tolower(unique(twosides_by_cid$drug1)) | 
                   infertile %in%  tolower(unique(twosides_by_cid$drug2)) ])

names(drugs_sides) <- tolower(names(drugs_sides))
list_alias_art_drugbank <- lapply(infertile, function(n) {
  print(n)
  
  
  index_res <- lapply(1:nrow(drugbank), function(index){
    
    name_db <- tolower(drugbank$name[index])
    syns_list <- tolower(unlist(strsplit(drugbank$alias[index], 
                                         split = ";", perl = T)))
    # print(syns_list)
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  name_db <- drugbank[index_res, ]$name
  alias <- drugbank[index_res, ]$alias
  
  df <- data.frame( name = name_db, syns = alias, stringsAsFactors = F)
  return(df)
  
  # pos_res <- lapply(1:nrow(drugbank), function(pos){
  # atc_list <- unlist(drugbank$atc[pos])
}
)
alias_art <- do.call( "rbind", list_alias_art_drugbank)

# check drugs indv:
# aux <- which(drugbank$name == "Follitropin")


# map art sider -----------------------------------------------------------

# using sider subset after mapping in drugbank with sider atcs
length(infertile[infertile %in% tolower(unique(names(drugs_sides)))]) # 15/18

list_art_sider <- lapply(names(drugs_sides), function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)

# 15/18 
## check
list_art_sider <- unique(list_art_sider)
list_art_sider <- unlist(list_art_sider)
length(list_art_sider)
for (i in setdiff(alias_art$name, list_art_sider)){print(paste0(i))}

# from sider map directly using sider names
sider_names$alias <- tolower(sider_names$alias)

list_art_sider_directly <- unlist(lapply(sider_names$alias, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)
)
list_art_sider_directly <- unlist(unique(list_art_sider_directly))
length(unique(list_art_sider_directly))
for (i in setdiff(alias_art$name, list_art_sider_directly)){print(paste0(i))}
# 13/18 lost Ganirelix and Buserelin

# try to recover lost art drugs -------------------------------------------

drugbank[drugbank$name == "Follitropin",]
# at sider, Follitropin == Metrodin?
# Metrodin == Urofollitropin at drugbank (used also for infertility)

sider_names[grep(sider_names$alias, pattern = "metrodin"), ]
sider_atc[sider_atc$atc == "G03GA04", ]  # is in df of drugbank
drugs_sides$urofollitropin ## OK! To add

sex_related_sider <- sider_atc[grep(sider_atc$atc , pattern = "^G03") , ]$atc

sex_related_sider_which <- unlist(lapply(sex_related_sider, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(drugbank_app), function(index){
    
    atcs <- unlist(drugbank_app$atc[index])
    if (n %in% atcs) {
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(drugbank_app[index_res, ]$name)
}
)
)
# print
for (i in unique(sex_related_sider_which)){print(paste0(i))}


# offsides art mapping----------------------------------------------------------------

offsides <- fread("workspace/projects/ddi/data/dbs/OFFSIDES.csv",
                  header = T)

offsides <- fread("data/dbs/offsides_subset.tsv",
                  header = T)

# pval threshold == 0.05
offsides_aliases <- unique(unlist(lapply(offsides$drug, function(n) {
  n <- tolower(strsplit(n, split = ",")[[1]][1])
  return(n)
}
)
)
)

art_off <- unique(unlist(lapply(offsides_aliases, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)
)
)
length(art_off) #96


# missing drugs in offsides:
for (i in setdiff(alias_art$name, art_off)){print(paste0(i))}

# twosides 2012 art mapping (drug names) ---------------------------------------

# twosides is combinations of drug 1 and drug 2
twosides_aliases_1_2012 <- unique(unlist(lapply(twosides$drug1, function(n) {
  n <- tolower(strsplit(n, split = ",")[[1]][1])
  return(n)
}
)
)
)

twosides_aliases_2_2012 <- unique(unlist(lapply(twosides$drug2, function(n) {
  n <- tolower(strsplit(n, split = ",")[[1]][1])
  return(n)
}
)
)
)

twosides_names_2012 <- union(twosides_aliases_1_2012, twosides_aliases_2_2012)

art_two_2012 <- unique(unlist(lapply(twosides_names_2012, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)
)
)
length(art_two_2012)  # 27/70  27/81


# twosides 2014 art mapping -----------------------------------------------

twosides_new <- fread("workspace/projects/ddi/data/dbs/TWOSIDES.csv",
                      header = T)

# twosides is combinations of drug 1 and drug 2
twosides_aliases_1 <- unique(unlist(lapply(twosides_new$drug_1_concept_name, function(n) {
  n <- tolower(strsplit(n, split = ",")[[1]][1])
  return(n)
}
)
)
)

twosides_aliases_2 <- unique(unlist(lapply(twosides_new$drug_2_concept_name, function(n) {
  n <- tolower(strsplit(n, split = ",")[[1]][1])
  return(n)
}
)
)
)

twosides_names <- union(twosides_aliases_1, twosides_aliases_2)

art_two <- unique(unlist(lapply(twosides_names, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  # df <- data.frame(name_db = alias_art[index_res, ]$name,
  #                  offiside_name = n, stringsAsFactors = F)
  return(alias_art[index_res, ]$name)
}
)
)
)
length(art_two)  # from 2014 curated release


# drugbank to kegg --------------------------------------------------------

# check how many unique sider drugs map to kegg from drugbank
# 1146 from 1560 of sider map to drugbank
drugbank_sider <- drugbank[drugbank$name %in% unique(plyr::compact(new_names)), ]
drugbank_sider <- drugbank_sider[grepl(drugbank_sider$groups, pattern = "approved") , ]

head(drugbank_sider$targets)

# need to transform targets to entrez and keep those that are in hsa pathways
targets <- unique(unlist(strsplit(drugbank_sider$targets, split = ";")))

kegg_targets <- 
  translate_genes(targets, ini_gene_id = "genename", final_gene_id = "keggid")
# 1038 translated / 13,93 % lost
ids <- drugbank_kegg[, c("#ID")]
ids$kegg <- T
write.table(ids, file = "drugs_with_targets_in_kegg.tsv" , sep = "\t", quote = F, row.names = F,
            col.names = T)


# tiles visualizations ----------------------------------------------------

guidelines <- read.table("workspace/projects/ddi/data/dbs/drugs_for_infertility_guidelines_NUMBER_OF_GUIDELINES.csv",
                         header = T, sep = ",")
guidelines <- melt(guidelines)

tile_plot <- 
  function(df, varx, vary, colour, filler) {
    ggplot(df, aes(x = varx, y = vary, fill= filler)) +
      geom_tile( colour = colour, size = 0.7) +
      scale_fill_manual(values = c("white", "lightgreen"))
  }

levels_drugs <- as.character(guidelines$NAME_DRUGBANK[!guidelines$NAME_DRUGBANK %in% c("Acetylsalicylic acid",
                                                                                       "Heparin", "Bromocriptine")])
levels_drugs <- append(levels_drugs, c("Acetylsalicylic acid",
                                       "Heparin", "Bromocriptine"))
factor(guidelines$NAME_DRUGBANK, ordered = T, levels = levels_drugs )
(base <- tile_plot(guidelines, varx = guidelines$variable, 
                   vary = factor(guidelines$NAME_DRUGBANK, ordered = T,
                                 levels = levels_drugs ) ,
                   colour = "gray50", filler = as.factor(guidelines$value)))

base + theme_light( base_size = 30) + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 0,size = rel(0.8), vjust = 2) ,
        axis.text.y = element_text(size = rel(0.8)) ,
        axis.text.x.top = element_text(),
        axis.title.x = element_text("drug", size = rel(1) ),
        axis.title.y = element_text("guideline", size = rel(1)),
        legend.text  = element_blank(),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank()
  ) + 
  scale_x_discrete(position = "top")  +
  labs(x = NULL, y = NULL)

atc_infertility <- 
  read.table("workspace/projects/ddi/data/dbs/DRUGS FOR INFERTILITY - All infertility drugs (clinical use)_ATC_CODE.csv",
             header = T, sep = ",", stringsAsFactors = F)

atc_infertility$atcs <- apply(atc_infertility, 1, function(row){
  atcs <- row["ATC_CODE"]  # get all atcs
  value <- strsplit(x = atcs, split = ";")[[1]]  # first delimiter
  value2 <- strsplit(x = value, split = ",")
  
  vector <- unique(unlist(lapply(value2, function(elem){
    return(elem[length(elem)])
  })))
  # vector <- c()
  # for (i in value2) {  # we need to iterate for each possible atc
  # child <- rev(i)[1]
  # vector <- append(vector, child)
  # }
  # print(unique(vector))
  # print(value2)
  # print(atcs)
  # break
  return(paste(vector, collapse = ","))
}
)

atc_infertility <- atc_infertility[c("NAME_DRUGBANK", "atcs")]

index2evaluate = grep(",", atc_infertility$atcs)

toadd = do.call("rbind", lapply(index2evaluate, function(i){
  data.frame(NAME_DRUGBANK = atc_infertility$NAME_DRUGBANK[i],
             atcs = unlist(strsplit(atc_infertility$atcs[i], split = ",")),
             stringsAsFactors = F)
}))
atc_infertility = atc_infertility[-index2evaluate,]
atc_infertility = rbind(atc_infertility, toadd)
atc_infertility$value <- "1"
atc_infertility <-
  reshape::add.all.combinations(atc_infertility,vars = c("NAME_DRUGBANK","atcs"))
atc_infertility$value[is.na(atc_infertility$value) ] <- 0

(base_atc <- tile_plot(atc_infertility, varx = atc_infertility$atcs, vary = atc_infertility$NAME_DRUGBANK ,
                       colour = "gray50", filler = as.factor(atc_infertility$value)))

base_atc + theme_light( base_size = 25) + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 0,size = rel(0.8)) ,
        axis.text.y = element_text(size = rel(0.8)) ,
        axis.text.x.top = element_text(),
        axis.title.x = element_text("drug", size = rel(1), vjust = .9 ),
        axis.title.y = element_text("guideline", size = rel(1)),
        legend.text  = element_blank(),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank()
  ) + 
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL)

#  number of drugs mapped in each database

mapped <- read.table("workspace/projects/ddi/data/infertily_drugs/DRUGS FOR INFERTILITY - ART drugs - databases.csv",
                     header = T, sep = ",")
mapped <- melt(mapped)

mapped

(base <- tile_plot(mapped, varx = factor(mapped$variable, ordered = F,
                                         levels = c("DRUGBANK", "KEGG", "SIDER",
                                                    "INTERACTOME", 
                                                    "OFFSIDES", "TWOSIDES")), 
                   vary = factor(mapped$NAME_DRUGBANK, 
                                 levels = sort(unique(mapped$NAME_DRUGBANK),
                                               decreasing = T)) ,
                   colour = "gray50", filler = as.factor(mapped$value)))

base + theme_light( base_size = 30) + 
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 0,size = rel(0.8), vjust = 2) ,
        axis.text.y = element_text(size = rel(0.8)) ,
        axis.text.x.top = element_text(),
        # axis.title.x = element_text("drug", size = rel(1) ),
        # axis.title.y = element_text("guideline", size = rel(1)),
        legend.text  = element_blank(),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank()
  ) + 
  scale_x_discrete(position = "top")  +
  labs(x = NULL, y = NULL)


# Mapping V2 ---------------------------------------------------------

# redo plots using novel added drugs after meeting decision
# load new data

getwd()

my_data <- list()
for (i in 1:9) {
  my_data[[i]] <- read_xlsx(
    "data/infertily_drugs/DRUGS FOR INFERTILITY_23_09.xlsx", 
    sheet = i, 
    col_names = T
  )
}

# merge each df (list element) in one big df
merged_excel <- as.data.frame(do.call(my_data, what = rbind))
merged_excel_sel <- merged_excel[c("ID_DRUGBANK","NAME_DRUGBANK","GUIDELINE_TYPE")]

# check duplicates/uniques
duplicated(merged_excel$ID_DRUGBANK)
unique(merged_excel$ID_DRUGBANK)

df_drugs <- merged_excel[merged_excel$STATUS == "C", ]
# drop duplicates based on ID
df_drugs <- df_drugs[!duplicated
                         (df_drugs[, c("ID_DRUGBANK")]), ]

# order of filtering matters
# df_drugs <- df_drugs[df_drugs$STATUS == "C", ]
selection <- df_drugs[c("ID_DRUGBANK","NAME_DRUGBANK","GUIDELINE_TYPE")]
selection
length(unique(selection$ID_DRUGBANK))
# 100

# map infertily drugs to dbs V2 ----------------------------------------------

infertile <- tolower(df_drugs$NAME_DRUGBANK)
length(infertile)

# sider names 55/91 46/70  // 48/81  61/100
length(infertile[infertile %in% tolower(unique(names(drugs_sides)))])

# twosides cid 29/91 25/81
length(infertile[infertile %in% tolower(unique(twosides_by_cid$drug1)) | infertile %in%  tolower(unique(twosides_by_cid$drug2)) ])

names(drugs_sides) <- tolower(names(drugs_sides))
list_alias_art_drugbank <- lapply(infertile, function(n) {
  print(n)
  
  
  index_res <- lapply(1:nrow(drugbank), function(index){
    
    name_db <- tolower(drugbank$name[index])
    syns_list <- tolower(unlist(strsplit(drugbank$alias[index], 
                                         split = ";", perl = T)))
    # print(syns_list)
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  name_db <- drugbank[index_res, ]$name
  alias <- drugbank[index_res, ]$alias
  
  df <- data.frame( name = name_db, syns = alias, stringsAsFactors = F)
  return(df)
  
  # pos_res <- lapply(1:nrow(drugbank), function(pos){
  # atc_list <- unlist(drugbank$atc[pos])
}
)
alias_art <- do.call( "rbind", list_alias_art_drugbank)

# check drugs indv:
# aux <- which(drugbank$name == "Follitropin")


# map art sider V2-----------------------------------------------------------

# using sider subset after mapping in drugbank with sider atcs
length(infertile[infertile %in% tolower(unique(names(drugs_sides)))]) # 61/100

list_art_sider <- lapply(names(drugs_sides), function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)
length(list_art_sider)
## check
list_art_sider <- unique(list_art_sider)
list_art_sider <- unlist(list_art_sider)
length(list_art_sider) #47/70 # 49/81  # 62/100

for (i in setdiff(alias_art$name, list_art_sider)) {print(paste0(i))}
setdiff(tolower(list_art_sider), infertile[infertile %in% tolower(unique(names(drugs_sides)))])

# from sider map directly using names from sider database
sider_names$alias <- tolower(sider_names$alias)

list_art_sider_directly <- unlist(lapply(sider_names$alias, function(n) {
  print(n)
  
  index_res <- lapply(1:nrow(alias_art), function(index){
    
    name_db <- tolower(alias_art$name[index])
    syns_list <- tolower(unlist(strsplit(alias_art$syns[index], 
                                         split = ";", perl = T)))
    
    if (n %in% syns_list | n %in% name_db) {
      print("found")
      return(index)
    }else{
      return(NA)
    }
  })
  
  index_res <- unlist(index_res[!is.na(index_res)])
  
  return(alias_art[index_res, ]$name)
}
)
)
list_art_sider_directly <- unlist(unique(list_art_sider_directly))
length(unique(list_art_sider_directly)) 
for (i in setdiff(alias_art$name, list_art_sider_directly)){print(paste0(i))}
# 50/91  o 41/70 o 42/81 or 57/100

# Now, check the sections to map to off/twosides
# offides 47/91
# twosides 18/91


# map to kegg V2------------------------------------------------------

# find how many drugs from selected have either targets, carrier or enzymes

# kegg
# load kegg
load("../pathway_variant_analysis/kegg_info_2020_no_annotation.RData")

# extract all genenames in pathways

kegg_genenames <- do.call(what = "rbind", (lapply(pathwaysInfo, function(path){
  genes <- gsub("^\\s+|\\s+$|[()]", "", unique(unlist(strsplit(path$ATTS$genes_names, split = "&|\\|"))))
  print(genes)
  print(path)
  aux <- append(aux_genenames, values = genes)
  df <- data.frame(values = genes, stringsAsFactors = F)
  return(df)
} )))
unique(kegg_genenames)

# or load annot table from gsrmutils

annot_kegg <- gsrmUtils::annot_table
annot_kegg$entrez <- gsub(annot_kegg$gene, pattern = "hsa:", replacement = "")

load("../pathway_variant_analysis/conversion_to_HGNC_2020.01.30.RData")

head(conversion_to_HGNC$entrezgene_id$hgnc_symbol)
entrez_hgnc <- conversion_to_HGNC$entrezgene_id

correspondence_kegg <- 
  entrez_hgnc[entrez_hgnc$entrezgene_id
              %in% annot_kegg$entrez, ]


to_kegg_targets <- 
  drugbank[(drugbank$`#ID`) %in% (selection$ID_DRUGBANK) , ]
head(to_kegg_targets[, c("name", "targets", "transporters", "carriers", "enzymes")])
targets <- (to_kegg_targets[,  c("name", "targets", "transporters", "carriers", "enzymes")])
# length(targets[targets == "" | transporters == "" | carriers == "" | enzymes == "", ])  # 67/70

# 85 drugs from 91 have targets. How many of them have those targets in kegg?
targets_in_kegg <- unique(unlist(strsplit(targets$targets, split = ";")))
kegg_targets <- 
  translate_genes(targets_in_kegg, ini_gene_id = "genename", 
                  final_gene_id = "keggid")
# 186 translated / 4,12 % lost
ids <- drugbank_kegg[, c("#ID")]

# map to interactome V2 -----------------------------------------------------------------
# load interactome 

int <- read.table("data/dbs/interactome_gene_name_final.txt", 
                  stringsAsFactors = F,header = T
)

int$node1
int$node2 

# is either targets/enzymes... from selected drugs in node 1 or node 2?
# save in the targets dataframe as new cols

targets$in_kegg <- unlist(apply(targets, 1, function(row){
  
  value <- 0
  
  target <- row[["targets"]]
  targets_split <- unique(unlist(strsplit(target, split = ";")))
  # print(targets_split)
  for (i in targets_split) {
    # print(tar)
    if ((i) %in%  correspondence_kegg$hgnc_symbol) {
      value <- 1
    }
  }
  enzymes <- row[["enzymes"]]
  enzymes_split <- unique(unlist(strsplit(enzymes, split = ";")))
  for (i in enzymes_split) {
    if ((i) %in% correspondence_kegg$hgnc_symbol) {
      value <- 1
    }
  }
  
  carriers <- row[["carriers"]]
  carriers_split <- unique(unlist(strsplit(carriers, split = ";")))
  for (i in carriers_split) {
    if ((i) %in% correspondence_kegg$hgnc_symbol) {
      value <- 1
    }
  }
  
  transporters <- row[["transporters"]]
  transporters_split <- unique(unlist(strsplit(transporters, split = ";")))
  for (i in transporters_split) {
    if ((i) %in% correspondence_kegg$hgnc_symbol) {
      value <- 1
    }
  }
  return(value)
}))

targets$in_int <- unlist(apply(targets, 1, function(row){
  
  value <- 0
  
  target <- row[["targets"]]
  targets_split <- unique(unlist(strsplit(target, split = ";")))
  # print(targets_split)
  for (i in targets_split) {
    # print(tar)
    if ((i) %in% int$node1 | (i) %in% int$node2) {
      value <- 1
    }
  }
  enzymes <- row[["enzymes"]]
  enzymes_split <- unique(unlist(strsplit(enzymes, split = ";")))
  for (i in enzymes_split) {
    if ((i) %in% int$node1 | (i) %in% int$node2) {
      value <- 1
    }
  }
  
  carriers <- row[["carriers"]]
  carriers_split <- unique(unlist(strsplit(carriers, split = ";")))
  for (i in carriers_split) {
    if ((i) %in% int$node1 | (i) %in% int$node2) {
      value <- 1
    }
  }
  
  transporters <- row[["transporters"]]
  transporters_split <- unique(unlist(strsplit(transporters, split = ";")))
  for (i in transporters_split) {
    if ((i) %in% int$node1 | (i) %in% int$node2) {
      value <- 1
    }
  }
  return(value)
}))

targets <- as.data.frame(targets)

# number of drugs by database ---------------------------------------------


# build df with presence/absence of selected drugs in each database
# example

mapped
selection$DRUGBANK <- 0
selection$SIDER <- 0
selection$OFFSIDES <- 0
selection$TWOSIDES <- 0
selection$KEGG <- 0
selection$INTERACTOME <- 0
head(selection)
length(selection$DRUGBANK)

head(selection)
head(selection$DRUGBANK, n = 100)
selection$NAME_DRUGBANK[83]

# grep to find index of drugname to be modified
grep("esterified", drugbank$name)  # 8372
drugbank$name[8372] <- "Esterified estrogens"

# order to create list DRUGBANK SIDER OFFSIDES TWOSIDES KEGG INTERACTOME
list_of_drugs <- list(drugbank$`#ID`, list_art_sider, art_off, art_two, targets$name)
str(list_of_drugs)

selection$DRUGBANK <- unlist(apply(selection, 1, function(row){
  drug <- row[["ID_DRUGBANK"]]
  value <- row[["DRUGBANK"]]
  if (drug %in% drugbank$`#ID`) {
    value <- 1
  } 
  return(value)
}
))

selection$SIDER <- unlist(apply(selection, 1, function(row){
  drug <- row[["NAME_DRUGBANK"]]
  value <- row[["SIDER"]]
  if (tolower(drug) %in% tolower(list_art_sider)) {
    value <- 1
  } 
  return(value)
}
))

selection$OFFSIDES <- unlist(apply(selection, 1, function(row){
  drug <- row[["NAME_DRUGBANK"]]
  value <- row[["OFFSIDES"]]
  if (tolower(drug) %in% tolower(art_off)) {
    value <- 1
  } 
  return(value)
}
))

selection$TWOSIDES <- unlist(apply(selection, 1, function(row){
  drug <- row[["NAME_DRUGBANK"]]
  value <- row[["TWOSIDES"]]
  if (tolower(drug) %in% tolower(art_two)) { # in 2014 version map more
    value <- 1
  } 
  return(value)
}
))

selection$KEGG <- 0
for (i in selection$NAME_DRUGBANK) {
  
  if (tolower(i) %in% tolower(targets$name)) {
    value <- targets[tolower(targets$name) == tolower(i), ]$in_kegg
    
    selection[tolower(selection$NAME_DRUGBANK) == tolower(i),]$KEGG <- value
  } else {
    selection[tolower(selection$NAME_DRUGBANK) == tolower(i),]$KEGG <- 0 } 
} 

selection$INTERACTOME <- 0
for (i in selection$NAME_DRUGBANK) {

  if (tolower(i) %in% tolower(targets$name)) {
    value <- targets[tolower(targets$name) == tolower(i), ]$in_int

    selection[tolower(selection$NAME_DRUGBANK) == tolower(i),]$INTERACTOME <- value
    } else {
      selection[tolower(selection$NAME_DRUGBANK) == tolower(i),]$INTERACTOME <- 0 } 
    } 


View(selection)

write.table(selection, 
            file = "results/01_mapping/table_drugs_in_databases_23_09_2020_sel",
            quote = F, row.names = F, col.names = T, sep = "\t")


# bar plot drugs by databases ----------------------------------------------------------------

selection_melted <- melt(selection,id.vars = colnames(selection[,1:3]))

ggplot(selection_melted, (aes(x = as.factor(variable), 
                         y = as.numeric(value), fill = variable))) +
  geom_bar(stat = "identity") + 
  xlab("Databases") + ylab("Number of drugs") +
  theme_light(base_size = 18) + 
  # theme(legend.text = element_text( size = rel(0.8)),
  theme(legend.text = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.title.x = element_text(vjust = -1.0)) +
  scale_fill_brewer(palette = "Dark2") 
  

# recheck if we lost some drugs  ----------------------------------------------------------------

# "acetylsalicylic acid" ("Aspirin") is on sider



# create atcs network -----------------------------------------------------

# load drugbank 2020
library(data.table)
drugbank <-fread("/home/sinh/workspace/local/ddi/data/dbs/drugbank_info_2020_07.txt", sep = "\t", header = T)

# see how is the object of drugbank
head(drugbank)
names(drugbank)
# get id, name, atcCods in a subset
drugAtc <- drugbank[, c("#ID", "name", "atcCodes", "groups")]
head(drugAtc)
drugAtc[7,] # example of a problematic one
test_atc <- drugAtc[7,]$atcCodes
# test_atc[nchar(strsplit(test_atc , split = ";" )) == 1]
# aux <- strsplit(test_atc , split = ";" )[[1]]
# lapply(aux, function(i){
#   cat(i, "\n")
#   # print(i)
# })
c <- vector()
for (i in aux) {
  print(i)
  aux_list <- strsplit(i, split = ",")
  res <- aux_list[[1]][(grepl("(?=^.{1}$)", aux_list[[1]], perl = T))]
  print(res)
  c <- append(c,res)
  final_res <- unique(c)
}
final_res
aux2 <- gsub(paste(final_res, ",", collapse = ""), pattern = " ", replacement = "")
sub(",$", replacement = "", aux2)
test_atc = "L02AE51,L02AE,L02A,L02,L;L02AE02,L02AE,L02A,L02,L;ABA,A;ADWED,B;ADAEAD,A"

drugAtc$atc_compressed <- unlist(apply(drugAtc, 1, function(row){
  atcs = row[["atcCodes"]]
  print(atcs)
  if (atcs == "") {
    return(NA)
  }
  aux = strsplit(atcs , split = ";" )[[1]]
  c = vector()
  for (i in aux) {
    # print(i)
    aux_list = strsplit(i, split = ",")
    res <- aux_list[[1]][(grepl("(?=^.{1}$)", aux_list[[1]], perl = T))]
    # print(res)
    c = append(c,res)
    final_res = unique(c)
  }
  print(final_res) # "V" "H"
  aux2 = gsub(paste(final_res, ",", collapse = ""), pattern = " ", 
               replacement = "")
  end = sub(",$", replacement = "", aux2)
  print(end) # "V,H"
  return(end)
}))
names(drugAtc)

drugAtc_approved <- drugAtc[
  grepl(drugAtc$groups, pattern = "approved", fixed = T ), ]

write.table(drugAtc[,c("#ID","atc_compressed")], 
            "results/03_atc_net/drugAtcs_2020.tsv", quote = F, sep = "\t",
            row.names = F)

check_atcs <- read.table("results/03_atc_net/drugbank_sif_atcs.tsv",
                         sep = "\t", header = T)
head(check_atcs)
length(unique(check_atcs$id))
# 1964 # unique

# sider correspondence ----------------------------------------------------

sider_all <-
  read.delim("data/dbs/sider_all_se.tsv", header = F,
             col.names = c("cid","cid2","ef_id","type","ef_id2","description"),
             stringsAsFactors = F)

sider_atc <- read.delim("data/dbs/sider_atc.tsv", 
                        header = F,
                        col.names = c("cid","atc"),
                        stringsAsFactors = F)
head(sider_all)
head(sider_atc)

