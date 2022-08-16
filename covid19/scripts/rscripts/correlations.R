# 4 sep 2020

# DE endometrium controls
# ihc.europa@gmail.com

#clean global enviroment if needed
rm(list = ls())

# save/load env
# save.image("scripts/rda/correlation_analysis_V3.Rda")
# load("scripts/rda/correlation_analysis_V2.Rda")

defile <- 
  "results/1-covid/de_gse150720_gsrmutils_voom_cmp_quartile_females.txt"
de_covid <- read.table(defile, header = T, sep = "\t")
head(de_covid)

# packages
library(ggplot2)
library(data.table)
library(igraph)


# data --------------------------------------------------------------------
genes_of_interest <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", 
                       "CTSL", "ACE2", "FURIN") 

# functions ---------------------------------------------------------------

# write function to subset network by correlation threshold and create igraph
# object undirected
subset_network <- function(network, threshold){
  sub <- network[
    abs(network$correlation) >= threshold, ]
  graph <- graph_from_data_frame(d = sub, directed = F)
}

# function to get number of correlations, fit square, number of genes...
# depending on correlations thresholds
power_fit_seq <- function(cor_df, bottom_threshold, upper_threshold, intervals, 
                          df_name, padj_value){
  # res_r2 <- list()
  # res_plots <- list()
  res = do.call("rbind", lapply(
    seq(bottom_threshold, upper_threshold, intervals), function(x){
      # obtain a subset of the respective df with padj and correlation value
      sub <- cor_df[
        cor_df[["p.adj"]] <= padj_value & abs(cor_df[["correlation"]]) >= x,]
      # keep nrows as a separate value
      num <- nrow(sub)
      # keep number of genes as a separate value
      genes <- length(unique(append(sub$gene1, sub$gene2)))
      # create a graph from dataframe
      sub_graph <- graph_from_data_frame(d = sub, directed = F)
      # decovid genes
      de <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                             rownames(de_covid)))
      de_sig <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                                 rownames(de_covid[de_covid$FC <= 0.05, ])))
      # intersect genes from covid infectivity
      infectivity <- paste0(intersect(unique(append(sub$gene1, sub$gene2)), 
                                      genes_of_interest), collapse = ",")
      # plot
      # plot_degree_distribution(sub_graph)
      # res_plots[[x]] <- p
      # obtain alpha and R²
      sub_res <- fit_power_law_b(sub_graph, name = df_name, name2 = x )
      # print(sub_res)
      data.frame( threshold = x,
                  Rsq = sub_res,
                  n_corrs = num, 
                  n_genes = genes,
                  n_de = de,
                  n_de_sig = de_sig,
                  infect = infectivity)
    }))
  return(res) 
}

# plot degree
plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", 
       ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
}

# from graph calculate power law fit
fit_power_law_b = function(graph, name, name2) {
  # dir check
  ifelse(dir.exists(paste0("results/4-correlation_thresholds/",name)), yes = 
           "donothing", dir.create(paste0("results/4-correlation_thresholds/",name)))
  savedir <- paste0("results/4-correlation_thresholds/",name,"/")
  
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  alpha_res <- paste0("Alpha =", round(alpha, 3))
  print(paste("R square =", round(R.square, 3)))
  rsquare_res <- paste0("R²=", round(R.square, 3))
  # res <- rsquare_res
  res <- round(R.square, 3)
  # plot
  png(paste0(savedir,name2,".png"), 1080, 920)
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)",
       ylab = "Probability (log)",
       col = 1, main = paste0("Degree Distribution ", name, "_", name2))
  curve(power.law.fit, col = "red", add = T, n = length(d))
  dev.off()
  return(res)
}

# correlation analysis (general) ----------------------------------------------------

# load endometrium join matrix

load("../../../fivibio-data/data/datasets/endometrium_joinControl_counts_ed.Rda", 
     verbose = T)

correlation_table_all_phases <- fread("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/corr_menstrual_cycle.tsv")
correlation_table_all_pf <- fread("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/corr_PF.tsv")
correlation_table_all_ese <- fread("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/corr_ESE.tsv")
correlation_table_all_mse <- fread("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/corr_MSE.tsv")
correlation_table_all_lse <- fread("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/corr_LSE.tsv")

hist(correlation_table$correlation)
hist(correlation_table_all_phases$pvalue)

# all phases

# LURAP1L  PKIA -0.07585158 0.433090495
# 2: FAM19A2  PKIA  0.25054531 0.008597350
# 3:     BBX  PKIA -0.21231766 0.026661341
# 4: FAM19A5  PKIA  0.10081583 0.296914117
# 5: WDR83OS  PKIA -0.30278890 0.001373149
# 6:   GSTK1  PKIA -0.22310530 0.019705554

# check if it's ok
cor.test( x =  dats3[c("FAM19A2"), ], y = dats3["PKIA", ] )

correlation_table_all_phases$p.adj <- p.adjust(correlation_table_all_phases$pvalue, 
                                               method = "fdr" )
correlation_table_all_pf$p.adj <- p.adjust(correlation_table_all_pf$pvalue ,
                                           method = "fdr")
correlation_table_all_ese$p.adj  <- p.adjust(correlation_table_all_ese$pvalue ,
                                             method = "fdr")
correlation_table_all_mse$p.adj  <- p.adjust(correlation_table_all_mse$pvalue ,
                                             method = "fdr")
correlation_table_all_lse$p.adj <- p.adjust(correlation_table_all_lse$pvalue ,
                                            method = "fdr")

# save in lists
corrs <- list(correlation_table_all_phases, correlation_table_all_pf,
              correlation_table_all_ese, correlation_table_all_mse, 
              correlation_table_all_lse)
names(corrs)
# add names of corrs
names(corrs) <- c("all", "pf", "ese", "mse", "lse")
names(corrs) 

# example of what is inside the list
# all:Classes ‘data.table’ and 'data.frame':	90269766 obs. of  5 variables:
#   ..$ gene1      : chr [1:90269766] "GSTM3" "NEBL" "TMEM254-AS1" "DHFR" ...
# ..$ gene2      : chr [1:90269766] "KCNB1" "KCNB1" "KCNB1" "KCNB1" ...
# ..$ correlation: num [1:90269766] 0.763 0.538 0.859 0.614 0.685 ...
# ..$ pvalue     : num [1:90269766] 5.32e-22 1.62e-09 7.57e-33 1.29e-12 2.11e-16 ...
# ..$ p.adj      : num [1:90269766] 1.59e-21 2.47e-09 5.30e-32 2.26e-12 4.51e-16 ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   $ pf :Classes ‘data.table’ and 'data.frame':	90269766 obs. of  5 variables:
#   ..$ gene1      : chr [1:90269766] "ORAI2" "VOPP1" "SHMT2" "TARSL2" ...
# ..$ gene2      : chr [1:90269766] "MBD4" "MBD4" "MBD4" "MBD4" ..

# to tests scripts
coraux <- list(head(correlation_table_all_phases), head(correlation_table_all_pf),
               head(correlation_table_all_ese), head(correlation_table_all_mse),
               head(correlation_table_all_lse))

# plot distribution of the p values adjusted by fdr
hist(correlation_table_all_phases$p.adj)
hist(correlation_table_all_pf$p.adj)
hist(correlation_table_all_ese$p.adj)
hist(correlation_table_all_mse$p.adj)
hist(correlation_table_all_lse$p.adj)

# summary of the padj data
menstrual_cycle_summary <- summary(correlation_table_all_phases$p.adj)
PF_summary <- summary(correlation_table_all_pf$p.adj)
ESE_summary <- summary(correlation_table_all_ese$p.adj)
MSE_summary <- summary(correlation_table_all_mse$p.adj)
LSE_summary <- summary(correlation_table_all_lse$p.adj)
menstrual_cycle_summary
PF_summary
ESE_summary
MSE_summary
LSE_summary

# number of padj correlations?
sum(correlation_table_all_phases$p.adj <= 0.05)
sum(correlation_table_all_pf$p.adj <= 0.05)
sum(correlation_table_all_ese$p.adj <= 0.05)
sum(correlation_table_all_mse$p.adj <= 0.05)
sum(correlation_table_all_lse$p.adj <= 0.05)


# filter correlations by covid genes --------------------------------------

# first, load DE COVID-19 results:
cov <- read.table("results/1 - covid/de_gse150720_gsrmutils_voom_cmp_quartile_females_fc_corrected.txt",
                  header = T)
sum(cov$adj.P.Val <= 0.05) #2580 sig genes
head(cov)
cov$gene <- rownames(cov)

# subset covid genes 
cov_paper1 <- cov[genes_of_interest, ]  # FURIN AND ACE ARE LOST
cov_paper1 <- cov_paper1[!is.na(cov_paper1$logFC), ]


# subset by padj
cov_sig <- cov[cov$adj.P.Val <= 0.05, ]
cov_sig$gene <- rownames(cov_sig)
head(cov_sig)

# RBIND WITH COVPAPER1
cov_sig <- rbind(cov_sig, cov_paper1)

# tranform to att and save
cov_att <- cov_sig[,c("gene", "FC_CACO", "logFC_CACO", "adj.P.Val")]
head(cov_att)
sum(cov_att$adj.P.Val <= 0.05) #2580
sum(duplicated(cov_att)) # 4
cov_att <- cov_att[!duplicated(cov_att), ]
nrow(cov_att) # 2582 genes after removing duplicates (2 not padj from paper1)
cov_att$adj_true <- 1
# cov_att$adj_true <- cov_att[cov_att$adj.P.Val <= 0.05]

write.table(cov_att, "results/5 - network_study/att_covid.tsv",
            sep = "\t", quote = F, col.names = T, row.names = F)


# how many of the COVID-19 DE genes can be found in the join endometrium dataset
# with ~13,000 genes?
# load whichever expression matrix
countmatrix <- read.table(
  "results/2 - endometrium/endometrium_controls_normalised_counts_matrix_MSE.tsv",
  header = T)
head(countmatrix)
nrow(countmatrix)
# 13437 genes in join dataset
length(intersect(rownames(cov_sig), countmatrix$genes))
# 2322 genes
total_cov <- 2322

# now, subset the big correlation matrix with correlations where
# only DE COVID-19 are in gene1 and gene2
# cov_net <- corrs[["pf"]][corrs[["pf"]]$gene1 |
#                            corrs[["pf"]]$gene2 , ]
# aux_cors <- head(corrs[["pf"]], n = 10000)

# test 1
# cov_net <- corrs[["pf"]][corrs[["pf"]]$gene1 %in% rownames(cov_sig) &
#                            corrs[["pf"]]$gene2 %in% rownames(cov_sig), ]

# loop for all

cov_nets <- list()
names(cov_nets) <- c("all", "pf", "ese", "mse", "lse")
phase_names <- c("all", "pf", "ese", "mse", "lse")

# in the next loop we extract correlations were in both cols we have genes
# from the expression matrix of CACO covid-dataset
# here we get the same number of correlations for the COVID genes,
# so nrows will be the same, but corr value will be different depending of phase
toadd <- lapply(names(corrs), function(x){
  print(x)
  toadd = corrs[[x]][corrs[[x]]$gene1 %in% rownames(cov_sig) &
                       corrs[[x]]$gene2 %in% rownames(cov_sig), ]
})
# same names as original list
names(toadd) 
names(toadd) <- names(corrs)
names(toadd) 
View(toadd)

# now, we filter by padj value
tofilter <- lapply(names(toadd), function(x){
  print(x)
  if (x == "lse") {# since lse does not have padj correlations; 
    threshold = 1
  } else {
    threshold = 0.05
  }
  print(threshold)
  tofilter = toadd[[x]][toadd[[x]]$p.adj <= threshold, ]
})
# again, assign names
names(tofilter)
names(tofilter) <- names(toadd)
names(tofilter)
View(tofilter)
str(tofilter)

# We also want to include correlations of CACO-diff expr genes with the genes 
# from paper 1: "TMPRSS2" "TMPRSS4" "BSG"     
# "MX1"     "CTSB"    "CTSL"    "ACE2"    "FURIN"  
# first, get all the correlations of paper1 genes from the big matrix of 90~MM corrs
toadd_cov <- lapply(names(corrs), function(x){
  toadd_cov = corrs[[x]][corrs[[x]]$gene1 %in% genes_of_interest |
                           corrs[[x]]$gene2 %in% genes_of_interest, ]
})
names(toadd_cov) <- names(tofilter)
str(toadd_cov)
"ACE2" %in% toadd_cov[["mse"]]$gene1 # T
"ACE2" %in% toadd_cov[["mse"]]$gene2 # T

# get a vector of covid genes and concatenate with the genes from the 
# CACO expression matrix 
targets <- lapply(names(tofilter), function(x){
  cov_in = unique(c(tofilter[[x]]$gene1, tofilter[[x]]$gene2))
  targets = c(cov_in,genes_of_interest)
})
names(targets) <- phase_names
str(targets)
"ACE2" %in% targets[["mse"]]

# we generate a subset were we have all target genes desired, but applying a 
# pvalue threshold
toadd_cov_mod <- lapply(names(toadd_cov), function(x){
  toadd_cov_mod = toadd_cov[[x]][toadd_cov[[x]]$gene1 %in% targets[[x]] &
                                   toadd_cov[[x]]$gene2 %in% targets[[x]], ]
  # INSERT PVALUE THRESHOLD
  toadd_cov_mod = toadd_cov_mod[toadd_cov_mod[[x]]$pvalue <= 0.05 ,]
})
names(toadd_cov_mod) <- phase_names
str(toadd_cov_mod)
"ACE2" %in% toadd_cov_mod[["pf"]]  #FALSE, filtered by pvalue

# second, add correlations to the filtered networks, keeping unique rels
toadd_net <- lapply(names(toadd_cov_mod), function(x){
  res = unique(rbind(tofilter[[x]], toadd_cov_mod[[x]]))
})
names(toadd_net) <- phase_names
str(toadd_net)

# now, subset by value of corr (by phase)
tocorr <- lapply(names(toadd_net), function(x){
  print(x)
  if (x == "all" ) {
    value = 0.98
  } else { 
    if (x == "lse") {
      value = 0.9
    } else {
      value = 0.65
    }
  }
  print(value)  # IN ABSOLUTE VALUE SO WE DON'T LOSE -CORRS
  tocorr = toadd_net[[x]][abs(toadd_net[[x]]$correlation) >= value, ]
})
names(tocorr)
names(tocorr) <- phase_names
names(tocorr)
str(tocorr)
# View(tocorr)

# finally, add a new col for each network with correlation value in absolute 
tonet <- lapply(names(tocorr), function(x){
  print(x)
  tocorr[[x]]$abs_corr <- round(abs(tocorr[[x]]$correlation), 3)
  tocorr[[x]]$correlation <- round(tocorr[[x]]$correlation, 3)
  tocorr[[x]]
})
names(tonet) <- phase_names
names(tonet)
str(tonet)

# recheck genes of interest
"BSG" %in% tonet[["mse"]]$gene1
"TMPRSS4" %in% tonet[["all"]]$gene2
"TMPRSS4" %in% tonet[["all"]]$gene1
"TMPRSS2" %in% tonet[["lse"]]$gene1
"ACE2" %in% tonet[["ese"]]$gene1
"ACE2" %in% tonet[["ese"]]$gene2
"ACE2" %in% tonet[["mse"]]$gene1
"ACE2" %in% tonet[["mse"]]$gene2
"ACE2" %in% tonet[["lse"]]$gene2
"ACE2" %in% tonet[["lse"]]$gene1
"ACE2" %in% tonet[["pf"]]$gene1
"ACE2" %in% tonet[["pf"]]$gene2


# calculate power fit law for each subset CALCULATE NUMBER OF GENES REMAINING FROM COVID
tograph <- lapply(names(tonet), function(x){
  print(x)
  graph = graph_from_data_frame(d = tonet[[x]], directed = F)  # creates graph
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  alpha_res <- paste0("Alpha =", round(alpha, 3))
  print(paste("R square =", round(R.square, 3)))
  rsquare_res <- paste0("R²=", round(R.square, 3))
  res = round(R.square, 3)
  rels = nrow(tonet[[x]])
  genes = length(intersect(unique(c(tonet[[x]]$gene1, tonet[[x]]$gene2)), cov_sig$gene))
  data.frame(fit_law = res, num_rels = rels, genes = genes, ratio = genes/total_cov)
}) 
names(tograph) <- names(toadd)
str(tograph)
# SAVE THE STR TO A TXT!!

phase_to_save <- "pf" 
write.table(tonet[[phase_to_save]], 
            file = paste0("results/5 - network_study//",phase_to_save,"_network.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T, append = F)


# Pvalue adjusted + correlation threshold decision ----------------------

# loop to obtain FDR correlations
listdf <- list()
plots <- list()
mens <- c("menstrual_cycle", "PF", "ESE", "MSE", "LSE")

for (i in (1:length(corrs))) {
  
  toplot = do.call("rbind", lapply(seq(0, 1, 0.05), function(x){
    data.frame(threshold = x,
               value = sum(corrs[[i]][["p.adj"]] <= 0.05 
                           & abs(corrs[[i]][["correlation"]]) >= x))
  }))
  
  listdf[[mens[i]]] <- toplot
  
  # plotdf
  p <- ggplot(toplot, aes(y = value/1e6, x = threshold, group = 1)) +
    geom_point() +
    geom_line() +
    theme_light(base_size = 18) + 
    theme(plot.title = element_text(paste0("n: ", mens[i])))
  # + ggtitle(paste0("n: ", mens[i]))
  print(p)
  # save plots in object
  plots[[mens[i]]] <- p
}


for (i in 1:length(plots)) {
  print(plots[[i]])
}

for (i in 1:length(plots)) { 
  # dir check
  ifelse(dir.exists(paste0("results/4-correlation_thresholds/padj")), yes = 
           "donothing", 
         dir.create(paste0("results/4-correlation_thresholds/padj")))
  svg(paste0("results/4-correlation_thresholds/padj/", mens[[i]], ".svg", 
             sep = ""), 1080, 920, onefile = T, pointsize = 18)
  print(plots[[i]])
  dev.off()
}

# igraph power fit law and degree distribution ----------------------------------------------------

# subset by threshold
mc_subset <- correlation_table_all_phases[
  abs(correlation_table_all_phases$correlation) >= 0.9, ]

# load a graph from a dataframe
mc_graph <- graph_from_data_frame(d = mc_subset, directed = F)

# write function to subset network by correlation threshold and create igraph
# object undirected
subset_network <- function(network, threshold){
  sub <- network[
    abs(network$correlation) >= threshold, ]
  graph <- graph_from_data_frame(d = sub, directed = F)
}

mc_graph <- subset_network(correlation_table_all_pf, threshold = 0.98)
pf_graph <- subset_network(correlation_table_all_pf, threshold = 0.75)

# degree and degree distribution

# For degree a numeric vector of the same length as argument v.
# For degree_distribution a numeric vector of the same length as the maximum 
# degree plus one. The first element is the relative frequency zero degree 
# vertices, the second vertices with degree one, etc.

d <- degree(mc_graph, mode = "all")
hist(d)
dd <- degree.distribution(mc_graph, mode = "all", cumulative = FALSE)

# write a function to plot the degree distribution
plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", 
       ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
}

plot_degree_distribution(mc_graph)


# plot and fit the power law distribution
fit_power_law_b = function(graph, name, name2) {
  # dir check
  ifelse(dir.exists(paste0("results/4-correlation_thresholds/",name)), yes = 
           "donothing", dir.create(paste0("results/4-correlation_thresholds/",name)))
  savedir <- paste0("results/4-correlation_thresholds/",name,"/")
  
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  alpha_res <- paste0("Alpha =", round(alpha, 3))
  print(paste("R square =", round(R.square, 3)))
  rsquare_res <- paste0("R²=", round(R.square, 3))
  # res <- rsquare_res
  res <- round(R.square, 3)
  # plot
  png(paste0(savedir,name2,".png"), 1080, 920)
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)",
       ylab = "Probability (log)",
       col = 1, main = paste0("Degree Distribution ", name, "_", name2))
  curve(power.law.fit, col = "red", add = T, n = length(d))
  dev.off()
  return(res)
}

mc_pwfl <- fit_power_law_b(mc_graph, "example", "0.98")
mc_pwfl

# for 0.98  threshold
# [1] "Alpha = 1.181"
# [1] "R square = 0.873"


# obtain all values for criteria in a loop ---------------------------------------------------

genes_of_interest <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", 
                       "CTSL", "ACE2", "FURIN") 

power_fit_seq <- function(cor_df, bottom_threshold, upper_threshold, intervals, 
                          df_name, padj_value){
  # res_r2 <- list()
  # res_plots <- list()
  res = do.call("rbind", lapply(
    seq(bottom_threshold, upper_threshold, intervals), function(x){
      # obtain a subset of the respective df with padj and correlation value
      sub <- cor_df[
        cor_df[["p.adj"]] <= padj_value & abs(cor_df[["correlation"]]) >= x,]
      # keep nrows as a separate value
      num <- nrow(sub)
      # keep number of genes as a separate value
      genes <- length(unique(append(sub$gene1, sub$gene2)))
      # create a graph from dataframe
      sub_graph <- graph_from_data_frame(d = sub, directed = F)
      # decovid genes
      de <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                             rownames(de_covid)))
      de_sig <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                                 rownames(de_covid[de_covid$FC <= 0.05, ])))
      # intersect genes from covid infectivity
      infectivity <- paste0(intersect(unique(append(sub$gene1, sub$gene2)), 
                                      genes_of_interest), collapse = ",")
      # plot
      # plot_degree_distribution(sub_graph)
      # res_plots[[x]] <- p
      # obtain alpha and R²
      sub_res <- fit_power_law_b(sub_graph, name = df_name, name2 = x )
      # print(sub_res)
      data.frame( threshold = x,
                  Rsq = sub_res,
                  n_corrs = num, 
                  n_genes = genes,
                  n_de = de,
                  n_de_sig = de_sig,
                  infect = infectivity)
    }))
  return(res) 
}

res_mc <- power_fit_seq(corrs[[1]], 0.98, 0.995, 0.005, "mc", padj_value = 0.05)
res_pf <- power_fit_seq(corrs[[2]], 0.75, 0.95, 0.05, "pf", padj_value = 0.05)
res_ese <- power_fit_seq(corrs[[3]], 0.75, 0.95, 0.05, "ese", padj_value = 0.05)
res_mse <- power_fit_seq(corrs[[4]], 0.75, 0.95, 0.05, "mse", padj_value = 0.05)
res_lse <- power_fit_seq(corrs[[5]], 0.95, 0.99, 0.01,  "lse", padj_value = 1)  
# plots are saved at fivibio-data/projects/covid/results/4-correlation_thresholds

l_res <- list(res_mc, res_pf, res_ese, res_mse)

# plot threshold and fit square
res_plots <- list()
lapply(l_res, function(x){
  # plotdf
  pr <- ggplot(x, aes(y = x$Rsquare, x = x$threshold, group = 1)) +
    geom_point() +
    geom_line() +
    theme_light(base_size = 18)  +
    theme(axis.title.x = element_text("threshold"), 
          axis.title.y = element_text("R²"))
  # theme(plot.title = element_text(paste0("n: ", mens[i])))
  # + ggtitle(paste0("n: ", mens[i]))
  res_plots[[deparse(substitute(x))]] <- pr
  print(pr)
  # save plots in object
  
})

ggplot(res_lse, aes(y = Rsquare, x = threshold, group = 1)) +
  geom_point() +
  geom_line() +
  theme_light(base_size = 18)  +
  theme(axis.title.x = element_text("threshold"), 
        axis.title.y = element_text("R²"))


# choose criteria for networks -------------------------------------------

genes_of_interest <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", 
                       "CTSL", "ACE2", "FURIN") 

# power fit law must be > 0.8
# preferible if it includes all genes_of_interest
# at least 2000 posible COVID genes mapped from which ~1500 p.adj?

grid_search <- function(cor_df, bottom_threshold, upper_threshold, intervals, 
                        df_name, padj_value){
  
  res = do.call("rbind", lapply(
    seq(bottom_threshold, upper_threshold, intervals), function(x){
      # obtain a subset of the respective df with padj and correlation value
      sub <- cor_df[
        cor_df[["p.adj"]] <= padj_value & abs(cor_df[["correlation"]]) >= x,]
      # keep nrows as a separate value
      num <- nrow(sub)
      # keep number of genes as a separate value
      genes <- length(unique(append(sub$gene1, sub$gene2)))
      # create a graph from dataframe
      sub_graph <- graph_from_data_frame(d = sub, directed = F)
      # decovid genes
      de <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                             rownames(de_covid)))
      de_sig <- length(intersect(unique(append(sub$gene1, sub$gene2)),
                                 rownames(de_covid[de_covid$FC <= 0.05, ])))
      # intersect genes from covid infectivity
      infectivity <- paste0(intersect(unique(append(sub$gene1, sub$gene2)), 
                                      genes_of_interest), collapse = ",")
      # plot
      # plot_degree_distribution(sub_graph)
      # res_plots[[x]] <- p
      # obtain alpha and R²
      sub_res <- fit_power_law_b(sub_graph, name = df_name, name2 = x )
      # print(sub_res)
      
      # ifelse check
      if ( sub_res >= 0.8 & de >= 2000 & de_sig >= 1000 & num <= 100000 ) {
        pass = T
      } else {
        pass = F
      }
      data.frame( threshold = x,
                  Rsq = sub_res,
                  n_corrs = num, 
                  n_genes = genes,
                  n_de = de,
                  n_de_sig = de_sig,
                  infect = infectivity,
                  pass = pass)
    }))
  return(res) 
}

# hyperparameter search based on results of seq approach
res_mc_grid <- grid_search(corrs[[1]], 0.965, 0.97, 0.0001, "mc", padj_value = 0.05)
res_pf_grid <- grid_search(corrs[[2]], 0.7, 0.8, 0.01, "pf", padj_value = 0.05)
res_ese_grid <- grid_search(corrs[[3]], 0.75, 0.8, 0.01, "ese", padj_value = 0.05)
res_mse_grid <- grid_search(corrs[[4]], 0.7, 0.75, 0.01, "mse", padj_value = 0.05)
res_lse_grid <- grid_search(corrs[[5]], 0.90, 0.95, 0.01, "lse", padj_value = 1)

# check if corr values are correct--------------------------------------------------------

## All phases
correlation_table_all_phases[correlation_table_all_phases$gene1 == "ORAI2" & 
                               correlation_table_all_phases$gene2 == "B3GLCT", ]
# ORAI2	B3GLCT	0.59
## PF
correlation_table_all_pf[correlation_table_all_pf$gene1 == "ORAI2" & 
                           correlation_table_all_pf$gene2 == "B3GLCT", ]
# ORAI2	B3GLCT	-0.4279988
## ESE
correlation_table_all_ese[correlation_table_all_ese$gene1 == "ORAI2" & 
                            correlation_table_all_ese$gene2 == "B3GLCT", ]
# -0.4425059
## MSE
correlation_table_all_mse[correlation_table_all_mse$gene1 == "ORAI2" & 
                            correlation_table_all_mse$gene2 == "B3GLCT", ]
# -0.1366914
## LSE
correlation_table_all_lse[correlation_table_all_lse$gene1 == "ORAI2" & 
                            correlation_table_all_lse$gene2 == "B3GLCT", ]
# -0.06564081

# OK?
cor( x =  dats3[c("ORAI2"), ], y = dats3["B3GLCT", ] ) # OK

# threshold decision
summary(correlation_table_all_phases[
  correlation_table_all_phases$correlation >= 0,  ]$correlation)
summary(correlation_table_all_phases$correlation < 0)


# subset networks by criteria and save SIFs ------------------------------------------

# define threshold value
num <- 0.97
padj <- 0.05  # 0.05 for all, except LSE (1)
c <- 1  # 1 == All phases, 2 == PF, 3 == ESE, 4 == MSE, 5 == LSE

subset_net_mc <- corrs[[1]][abs(corrs[[1]]$correlation) >= num &
                              corrs[[1]]$p.adj <= padj ] 
subset_net_pf <- corrs[[2]][abs(corrs[[2]]$correlation) >= num &
                              corrs[[2]]$p.adj <= padj ] 
subset_net_ese <- corrs[[3]][abs(corrs[[3]]$correlation) >= num &
                               corrs[[3]]$p.adj <= padj ] 
subset_net_mse <- corrs[[4]][abs(corrs[[4]]$correlation) >= num &
                               corrs[[4]]$p.adj <= padj ] 
subset_net_lse <- corrs[[5]][abs(corrs[[5]]$correlation) >= num &
                               corrs[[5]]$p.adj <= padj ] 



sub_nets <- list()

###
# ADD genes

# loop genes of interest not in subset network
# find them in corresponding network, (c) and add them as rbind to 
# subset network

# should introduce genes of interest correlations anyways
genes_of_interest <- c("TMPRSS2", "TMPRSS4", "BSG", "MX1", "CTSB", 
                       "CTSL", "ACE2", "FURIN")
# [1] "TMPRSS2" "TMPRSS4" "BSG"     "MX1"     "CTSB"    "CTSL"   
# [7] "ACE2"    "FURIN"  
# find then in the corresponding correlation network (either A,B) and include
# all in the filtered network (considering if it's already in and also
# B>A directions...)

names(corrs) <- c("mc", "pf", "ese", "mse", "lse")

nt <- "mc"

write.table(subset_net, quote = F, 
            file = paste0(
              "/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/", 
              "filtered_corr_test_alternative_", nt, ".tsv"), 
            sep = "\t",
            row.names = F, col.names = T)



# for loop for several networks
# save filtered correlation and export to cytoscape
# TESTING

# filtered_list <- list()
# 
# for (i in (1:length(corrs))) {
# 
# to_add <- corrs[[i]][
#   abs(corrs[[i]]$correlation) >= num, ]
# 
# filtered_list[[mens[i]]] <- to_add
# 
# }
# 
# for (i in mens) {
# write.table(i, quote = F, file = paste0("/data/network/nas01/fivibio-data/projects/covid/results/2-endometrium/", "filtered_corr_", i, ".tsv"), 
#             sep = "\t",
#             row.names = F, col.names = T)
# }

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
# tests scripts for cytoscape ----------------------------------------------------

toplot = do.call("rbind", lapply(seq(0, 1, 0.05), function(x){
  data.frame(threshold = x,
             value = sum(correlation_table_all_lse$p.adj <= 1 
                         & abs(correlation_table_all_lse$correlation >= x)))
}))

p <- ggplot(toplot, aes(y = value/1e6, x = threshold, group = 1)) +
  geom_point() +
  geom_line() +
  theme_light(base_size = 18) + 
  theme(plot.title = element_text(paste0("n: ", mens[i])))
# + ggtitle(paste0("n: ", mens[i]))
p


###

sub <- correlation_table_all_pf[
  correlation_table_all_pf[["p.adj"]] <= 0.05 & abs(
    correlation_table_all_pf[["correlation"]]) >= 0.85,]
# create a graph from dataframe
sub_graph <- graph_from_data_frame(d = sub, directed = F)

# add correlation values
E(sub_graph)$weight <- abs(sub$correlation)
head(sub_graph)
sub[sub$gene1 == "GSTP1" , ]  # BSG: 0.8680832
E(sub_graph)

oc <- cluster_optimal(graph = sub_graph)

fgc <- cluster_fast_greedy(sub_graph, merges = T, modularity = T,
                           membership = T, weights = E(sub_graph)$weight)
fgc
# number of modules
sizes(fgc)
# gene to which group (they don't seem to repeat)
membership(fgc) 
# modularity 
modularity(fgc)
communities <- groups(fgc)

# obtain group as att for cytoscape

sub$group <- unlist(apply(sub, 1, FUN = function(row){
  c = 0
  # extract gruop and add it as an edge attribute
  lapply(communities, function(x){
    c = c + 1
    # x is each group (list) with the respective genes
    # names(communities)[[x]]
    # print(x)
    print(names(communities)[[c]])
    print(row["gene1"])
    if (row["gene1"] %in% x | row["gene2"] %in% x) {
      return(names(communities)[[c]]) }
  })
}))
head(sub)
