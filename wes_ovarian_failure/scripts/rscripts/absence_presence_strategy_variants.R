########### Description #######################################################
# 25 nov 2019

# Presence/absence in cases/control groups of certain variant

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(gplots)  # for heatmap2 function
library(RColorBrewer)  # for color palettes
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

save.image("workspace/projects/pof/scripts/rdatas/
# absence_presence_variants-")
load("workspace/projects/pof/scripts/rdatas/fisher_test_strategy_variants.Rda")

# Create variable for working directory
dir <- "workspace/projects/pof/"


# Functions to be used ---------------------------------------------------------

# Function that filters variants over a certain % of cases and below a % of 
# controls
filter_proportions <- function(df , control_value, case_value){
  res <- df %>% 
    filter(df$controlpercent <= as.integer(control_value)
           & 
             df$casepercent >= as.integer(case_value)) 
  
  return(res)
}

# ggplot theme settings
bar.zyx.theme <-  theme_light( base_size = 18 ) +
  theme(axis.text.x = element_text(hjust = 0.65),
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

# CHECK
#  plot_effect + default_theme()

# Create matrix heatmap-like with counts over and below cases/controls ---------

# load fisher test results (we could also load the python output) as it has
# the number of variants for each group
fisher_cc_power <-
  read.table(paste0(dir,"results/13_study_candidates/bio_fisher_pow.txt"), 
                                               header = T)

# Create columns with proportions of cases and controls
fisher_cc_power$controlpercent <- 
  ((fisher_cc_power$control/32) * 100)
fisher_cc_power$casepercent <-
  ((fisher_cc_power$cases/118) * 100)


## create matrix 100 * 100 <= x controls, => y cases
size <- 101  # so it actually adds the 0% 
m_variants <- matrix(nrow = size, ncol = size, 
                     data = 0)  # we create a matrix of the size desired
                                # and fill it with zeros
# Rownames and colnames should be the % intervals
rownames(m_variants) <- as.numeric(0:(size - 1)) 
colnames(m_variants) <- as.numeric(0:(size - 1))

# matrix zero substitution by actual values from number of variants for each
# interval
for (i in 1:101) {
  # iterate rows
  for (j in 1:101) {
    # iterate cols
    m_variants[i,j] <- 
      nrow(filter_proportions(fisher_cc_power,case_value = i - 1,
                              control_value = j - 1)) 
    # we would start by adding for example the number of variants below 0% ctr
    # and over 0% cases; for when we have i = 15 and j = 15, we would be filling
    # that index in the matrix with the number of variants below 15% of ctr
    # and over 15% of cases
  } 
  
}
# We subset the matrix creating intervals of five by five. (dim goes down to
# 21 x 21 from 101 x 101)
mv_intervals <- m_variants[as.numeric(rownames(m_variants)) %% 5 == 0, 
                  as.numeric(colnames(m_variants)) %% 5 == 0]

# plot the matrix
(heat_matrix <- heatmap.2(
  mv_intervals, dendrogram = "none", 
                Rowv = FALSE, 
                Colv = FALSE, trace = "none" ,
                cellnote=ifelse(mv_intervals ==0, NA, mv_intervals),
                # main = "Number of variants over % cases and below % controls ",
                # xlab = " Over % Controls",
                # ylab = " Below % Cases",
                # Row/Column Labeling
                RowSideColors,
                notecol = "black",
                colRow = "black",
                colCol = "black",
                ColSideColors,
                margins =  c(5,5),
                sepcolor ="black",
                sepwidth=c(0.01,0.01),
                hline=median(breaks),
                vline=median(breaks),
                offsetRow = 0.001,
                offsetCol = 0.001,
                keysize = 1.5,
                notecex = 1.2,
                cexRow = 2.0,
                cexCol = 2.0,
                # color key + density info
                density.info = "histogram",
                key = FALSE,
                # col = colorRampPalette(c("darkorange4","yellow", "white"))
                col = colorRampPalette(rev
                                       (c("paleturquoise1","paleturquoise",
                                          "paleturquoise2","paleturquoise3", 
                                          "paleturquoise4")))
                (n = 1000)))

# check distributuion of cases in a subset of variants --------------------

# must load base data for contingency tables from python output
# subset of variants not in controls 
filter.v.0c <- 
  fisher_cc_power[fisher_cc_power$control == 0,]

barplot(table(filter.v.0c$cases),
        col = "coral",
        main = "Number of variants not in controls by presence in cases", 
        xlab = "Number of cases",
        ylab = "Number of variants") 

# lets put a ylim an substract only 1 case values
barplot(table(filter.v.0c$cases[filter.v.0c$cases > 1 ]),
        col = "brown1",
        main = "Number of variants not in controls by presence in cases", 
        xlab = "Number of cases affected > 1",
        ylab = "Number of variants",
        ylim = c(0,5000)) 

barplot(table(filter.v.0c$cases[filter.v.0c$cases > 2 ]),
        col = "cadetblue",
        main = "Number of variants not in controls by presence in cases", 
        xlab = "Number of cases affected > 2",
        ylab = "Number of variants",
        ylim = c(0,2000)) 

# check distribution of cases in all variants --------------------------------------

fisher_disease <- 
  read.table(paste0(
    dir,"results/07_test_fisher_rest/contingency_table_targeted.txt"),
    header = T)

# get a boleean value column if variants are in targeted hypothesis or not.
fisher_cc_power$disease <- fisher_cc_power$variant %in% 
  fisher_disease$variant

# transform boleean value to string
fisher_cc_power$disease[fisher_cc_power$disease == FALSE] <- 
  "No-POI"
fisher_cc_power$disease[fisher_cc_power$disease == "TRUE"] <- 
  "POI"

# now, do a subset of variants not in controls 
filter.v.0c <-
  fisher_cc_power[fisher_cc_power$control == 0,]

# alternative : do subset of variants not in cases
filter.v.0c <- 
  fisher_cc_power[fisher_cc_power$cases == 0,]

# plot
# change thresholds at will
ggplot(filter.v.0c, (aes(x = as.factor(cases), 
                             fill = as.factor(disease)))) +
  geom_bar(stat = "count",
           colour = "black") + 
  xlab("Número de casos afectados") + ylab("Número de variantes") +
  theme_light(base_size = 18) + 
  theme(legend.text = element_text( size = rel(0.8)),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(hjust = 0.7),
        axis.title.x = element_text(vjust = -1.0)) +
  scale_fill_brewer(palette = "Set1") 


# visualization by effect at protein level -----------

bioeffect <- 
  fread(
    file = paste0(dir,"results/13_study_candidates/",
                  "info_non_targeted_with_n_genes_aff.txt"), header = T,
    stringsAsFactors = F)

# add all columns from fisher from their respective variant
bioeffect <- left_join(fisher_cc_power,  bioeffect, by = "variant")
# select those that appear in 0 controls
bioeffect <- bioeffect[bioeffect$control   == 0,]

# split and select 1st list effect on protein level annotated
bioeffect$biotype_1 <- unlist(lapply(bioeffect$bioEff, function(row){
  
  split <- unlist(strsplit(row, split = "," ,fixed = T ))[[1]][1]
  if (split == "synonymous_variant") {
    split_syn <- unlist(strsplit(row, split = "," ,fixed = T ))[[2]]
    split2 <- unlist(strsplit(split_syn, split = "&" ,fixed = T ))[[1]][1]
  } else {
    split2 <- unlist(strsplit(split, split = "&" ,fixed = T ))[[1]][1]
  }
  return(split2)
}))

# effect on protein list
target_change_list <- c("protein_protein_contact" ,"frameshift_variant" ,
                       "stop_gained" ,"start_lost"
                       ,"splice_acceptor_variant" ,"splice_donor_variant" ,
                       "structural_interaction_variant"
                       ,"3_prime_UTR_variant" ,"5_prime_UTR_variant" ,
                       "disruptive_inframe_insertion" ,"missense_variant")

target_change_list_del <- 
  c("splice_acceptor_variant","splice_donor_variant",
    "stop_gained","stop_lost","start_lost","protein_protein_contact",
    "frameshift_variant")

# add a column separating most deleterious variants from missense (as factors) 
# && group splice and stop/start variants
bioeffect$deleterious <- unlist(lapply(bioeffect$biotype_1, 
                                       function(type){
  res <- "missense&UTR"
  if (type %in% target_change_list_del) {
    res <- type
  }
  return(as.factor(res))
}))


bioeffect$deleterious <- unlist(lapply(bioeffect$deleterious, 
                                             function(biotype){
  res <- biotype
  if (biotype == "splice_acceptor_variant" || biotype == "splice_donor_variant") {
    res <- "splice_variant"
  } 
  if (biotype == "stop_lost" || biotype == "stop_gained" || biotype == "start_lost") {
    res <- "stop/start_variant"
  }
  return(as.factor(res))
  
}))

# graphical visualization
# dark2 or set1 recommended colours for visualization
plot_effect <- ggplot(bioeffect[bioeffect$deleterious != "missense&UTR", ] , 
       aes(x = as.factor(
         bioeffect[bioeffect$deleterious != "missense&UTR" , ]$cases ) , 
           fill = as.factor(
             bioeffect[
               bioeffect$deleterious != "missense&UTR" , ]$deleterious) ) ) +
  geom_bar( stat = "count" , 
           colour = "black" ,
           position = "stack") +  # If we would want to put the number of v
  # in graph:
  # geom_text(aes(x = as.factor(bioeffect$cases), label = ..count..  ) , 
  #           data =  NULL, 
  #           inherit.aes = T,
  #           stat = "count", 
  #           size = 5, 
  #           position_stack(vjust = 0.5)
  #           ) +
  # facet_zoom(as.numeric(cases)) +
  # facet_grid(~poi) +
  # facet_wrap(~poi,drop = T, ) +
  xlab("Número de casos afectados") + 
  ylab("Número de variantes") + 
  scale_x_discrete(breaks = c(1:31)) +
  theme_light( base_size = 18)  +
  scale_fill_brewer(palette = "Set1") 

# add barplot theme
plot_effect + bar.zyx.theme

# presence absence of mutations in each sample---------------------------------

presence_absence <- 
  fread(file = paste0(dir,"/results/13_study_candidates/",
                           "presence_absence.txt"), header = T)

# how many variants (0>1/2) can we found for each patient?
presence.by.sample <- colSums(presence_absence[,-1], na.rm = F )

# select variants by presence in a % of groups ----------------------------

# get variants in at least 2 controls 0 cases
variants.to.study <- 
  filter_proportions(fisher_cc_power, control_value =  10,
                                        case_value = 0)
# get over 5% cases 0 controls
variants.to.study <- 
  filter_proportions(fisher_cc_power, control_value =  0,
                                        case_value = 5)
# get over 10% cases 0 controls
variants.to.study <- filter_proportions( fisher_cc_power, control_value = 0, 
                                        case_value = 10)
# save to process in python
write.table(
  variants.to.study$variant,
  file = paste0(dir,"results/14_filtered_variants/variants_5p_cases.txt"), 
  quote = F, col.names = F, row.names = F)

# check coverage of genes obtained  ---------------------------------------

cover <- read.delim(
  "pof/results/20_for_paper/coverage_per_gene.tsv", sep = " ", header = T)
cover

genes <- read.delim("pof/results/20_for_paper/variants_66_info_paper.txt", header = T,
                    sep = "\t")$Genes
all <- read.delim("pof/results/20_for_paper/variants_66_info_paper.txt", header = T,
                    sep = "\t", stringsAsFactors = F)

res <- cover[cover$gene %in% genes, ]
min(res$mean)
rownames(res) <- res$gene

setdiff(res$gene, genes)  # genes have more than one variant

all$Coverage <- unlist(lapply(all$Genes, function(x){
   # gene <- x["Genes"] 
   # print(gene)
   coverage <- res[x,]$mean
   return(coverage)
}))
all$Coverage
all$Coverage <- round(all$Coverage)
all
# C1orf101 not found -> is in reallity CATSPERE
# CATSPERE 151.15480987055

write.table(all, "pof/results/20_for_paper/66_info_coverage.tsv", sep = "\t",
            row.names = F, col.names = T, quote = F)
