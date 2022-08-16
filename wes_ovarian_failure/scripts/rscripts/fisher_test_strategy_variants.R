########### Description #######################################################
# 4 oct 2019

# Creates tables of contingency for applying fisher test and analyzes them.

#working directory
# dir=('/home/ihenarejos/workspace/')
# setwd(dir)
setwd("/home/ihenarejos/")

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RVAideMemoire)  # for multiple fisher tests comparisons
library(statmod)  # for power fisher test
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
# fisher_test_strategy_variants.Rda")
load("workspace/projects/pof/scripts/rdatas/
# fisher_test_strategy_variants.Rda")

# Create variable for working directory
dir <- "workspace/projects/pof/"

# Functions to be used ---------------------------------------------------------

# sample data 
# basic data that will be used during script
# maximum number of patients for each study design group
n_ctrl<-32
n_cases<-118
n_foo<-83
n_fop<-35

# Function that uses sample data to add new columns for count table obtained in 
# python.
add_columns_contigency <- function(count_table) {
  
  count_table$control_not_affected<-apply(count_table,1,function(row){
    control.in<-row[["control"]]
    
    control.rest<- 32-as.numeric(control.in)
    return(control.rest)
  })
  
  count_table$foo_not_affected<-apply(count_table,1,function(row){
    foo.in<-row[["foo"]]
    
    foo.rest<- 83-as.numeric(foo.in)
    return(foo.rest)
  })
  
  
  count_table$fop_not_affected<-apply(count_table,1,function(row){
    fop.in<-row[["fop"]]
    
    fop.rest<- 35-as.numeric(fop.in)
    return(fop.rest)
  })
  
  
  count_table$cases_not_affected<-apply(count_table,1,function(row){
    cases<-row[["cases"]]
    
    cases.rest<- 118-as.numeric(cases)
    return(cases.rest)
  })
  return(count_table)
}

# Function that generates a contingency table for each variant and apply fisher 
# test on them
fisher_test_variants<-function(df,vid,ctrmut,ctr,casemut,case){
  
  list_dfs_fisher_tests<-apply(df,1,function(row){
    
    variant_id<-row[[vid]]  # id from variant
    
    control_value<-as.numeric(row[[ctrmut]])  # controls that have the variant
    control_value2<-as.numeric(row[[ctr]])    # controls with reference allele
    casos_value<-as.numeric(row[[casemut]])    # cases that have the variant
    casos_value2<-as.numeric(row[[case]])    # cases with reference allele
    
    #create matrix and apply fisher test
    fisher_list<-fisher.test(matrix(c(control_value,control_value2,casos_value,casos_value2),byrow = T,
                                    nrow = 2,ncol = 2))  # normal 2x2 table
    
    # fisher test return a list of components of class "htest"  with components such as p.value and odds ratio estimate
    # that we can get from. So for each variant, return a dataframe with 3 columns; id of variant, pvalue and odds ratio
    return(data.frame(variant = variant_id, p.value = fisher_list$p.value, estimate = fisher_list$estimate))
  })
  
  # once we have dataframes from each variant, we ;
  # 1; do a rbindlist so each dataframe for each variant is joined together in an only dataframe where each row
  #   is a variant with the results from its fisher test
  # 2; merge the dataframe from step 1 with original dataframe that had the basis data from python.
  return(merge(df,rbindlist(list_dfs_fisher_tests), by= "variant"))
  
}  


fisher_test_variants_freq<-function(df, vid, ctrmut, casemut, control_pop){
  list_dfs_fisher_tests<-apply(df,1,function(row){
    variant_id<-row[[vid]]  # id from variant
    print(variant_id)
    control_value<-as.numeric(row[[ctrmut]]) * control_pop  # controls that have the variant
    control_value2<- (1 - as.numeric(row[[ctrmut]])) * control_pop    # controls with reference allele
    cases_value<-as.numeric(row[[casemut]]) * 118    # cases that have the variant
    cases_value2<- (1 - as.numeric(row[[casemut]])) * 118    # cases with reference allele
    
    #create matrix and apply fisher test
    fisher_list<-fisher.test(matrix(c(control_value,control_value2,cases_value,cases_value2),byrow = T, nrow = 2,ncol = 2))  # normal 2x2 table
    
    # fisher test return a list of components of class "htest"  with components such as p.value and odds ratio estimate
    # that we can get from. So for each variant, return a dataframe with 3 columns; id of variant, pvalue and odds ratio
    return(data.frame(vid = variant_id, p.value = fisher_list$p.value, estimate = fisher_list$estimate))
  })
  
  # once we have dataframes from each variant, we ;
  # 1; do a rbindlist so each dataframe for each variant is joined together in an only dataframe where each row
  #   is a variant with the results from its fisher test
  # 2; merge the dataframe from step 1 with original dataframe that had the basis data from python.
  # return(merge(df,rbindlist(list_dfs_fisher_tests), by= "vid"))
  return(rbindlist(list_dfs_fisher_tests))
  
} 

# Function that calculates statistical power of fisher test for each variant
power_fisher_variants <- function(df,total1,total2){
  
  df_powf<-apply(df,1,function(row){
    
    # Similar to fisher test function, where we will be adding a new column
    # to original dataframe with statistical power of fisher test calculated
    variant_id<-row[["variant"]]
    control_value<-as.numeric(row[["control"]])
    cases_value<-as.numeric(row[["cases"]])
    
    #create matrix and apply fisher test
    pow <- power.fisher.test(control_value/total1,cases_value/total2,total1,
                             total2,nsim = 10000)
    
    return(data.frame(variant = variant_id, fisher.power = pow))
  })
  
  return(merge(df,rbindlist(df_powf), by= "variant"))
  
}  

# Function that checks if variant has multiple alleles. To use for visualization
# of data
check_if_variant_multi <- function(df_variants,id_colum, delim) {
  
  df_variants$multi<-apply(df_variants,1,function(row){
    id<-row[id_colum]
    alt<-unlist(strsplit(id,split = "_")[[1]])[4]
    
    check_value<-grepl(alt,pattern = delim)
    return(check_value)
  })
  return(df_variants)
}

# Function that creates a df using counts of how many variants are SNV and how
# many are multiallelic for barplot visualization
barplot_counts <- function(df, type = ""){
  aux <- table(df$multi)
  aux <- data.frame(SNV = aux["FALSE"], multi = aux["TRUE"],row.names = NULL, 
                    type = type)
  return(aux)
}

# Function that creates a barplot for visualization of fisher test of variants
plot_variants_fisher <- function(df,title,ylim,base_size,colour1,colour2){
  ggplot(data=df, aes(x = factor(variable), group = factor(variable), y = value, 
                      fill = variable)) + 
    geom_bar(stat="identity", colour = "black ", width = 0.5, 
             position = position_dodge(width = 0.6)) + 
    scale_fill_manual(values = c("SNV" = colour1, "multi"
                                 = colour2)) +
    # scale_fill_brewer(palette = "Set1") +
    geom_text(aes(label = value), position = position_dodge(width = 0.6), 
              vjust =-0.4, size = 6 ) + 
    facet_wrap(.~pvar) + ylim ( 0, ylim) +
    xlab("") +  ggtitle( title ) + 
    ylab("Número de variantes significativas") +
    theme_light(base_size = base_size) +   
    theme(legend.position = "none", 
          axis.text = element_text(size=rel(1.6)),
          axis.text.x = element_text(angle=90),
          axis.title.y = element_text(margin = margin(r = 10)),  # for margin
          axis.title = element_text(size=rel(1.7), vjust = -0.5 ),
          plot.title = element_text(size=rel(2.4), hjust=0.5 ),
          legend.title = element_blank(),
          legend.text = element_blank(),
          strip.text.x = element_text(size=rel(1.8), color="red",
                                      face="bold.italic"),
          strip.text.y = element_text(size=rel(1.8), color="red",
                                      face="bold.italic"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))
}

# Function to perform pairwise fisher tests using the genotypes in variants
fisher_pairwise_2_x_3 <- 
  function(df,vid,control_AA,control_Aa,control_aa,case_AA,case_Aa,case_aa){
    
    df_fisher_tests<-apply(df,1,function(row){
      
      variant_id<-row[[vid]]
      
      control1<-as.numeric(row[[control_AA]])
      control2<-as.numeric(row[[control_Aa]])
      control3<-as.numeric(row[[control_aa]])
      casos1<-as.numeric(row[[case_AA]])
      casos2<-as.numeric(row[[case_Aa]])
      casos3<-as.numeric(row[[case_aa]])
      
      #create matrix and apply fisher test
      fisher_pval<-fisher.test(
        matrix(c(control1,control2,control3,casos1,casos2,casos3),
               byrow = T,nrow = 2,ncol = 3))
      
      fisher_pairwise<-fisher.multcomp(
        matrix(c(control1,control2,control3,casos1,casos2,casos3),
               byrow = T,nrow = 2,ncol = 3,
               dimnames = list(c("control","case"),c("AA","Aa","aa"))),
        p.method = "fdr")
      
      return(data.frame(variant = variant_id, p.value = fisher_pval$p.value, 
                        p.adjust=fisher_pairwise$p.value )) 
    })
    return(merge(df,rbindlist(df_fisher_tests), by = "variant"))
  }  



# Load contingency table basis data and filter variants ----------------------

# We obtained basis data from python script where we also applied quality and 
# biological filters

# non targeted
count_table_raw <- fread(
  paste0(dir,"/results/07_test_fisher_rest/contingency_table_non_targeted.txt", 
         collapse = "" ), stringsAsFactors = F,data.table = F)

# targeted 
count_table_raw <- fread(
  paste0(dir,"/results/07_test_fisher_rest/contingency_table_targeted.txt", 
         collapse = "" ), stringsAsFactors = F,data.table = F)  

str(count_table_raw)  # check data structure
dim(count_table_raw)  # check dim

# Filter out variants that didn't meet criteria and as such didn't counted 
# for each study group.
count_table <- count_table_raw[rowSums(count_table_raw[, -1]) > 0,]
dim(count_table)

# Filter those present in 1 case (intra-specific)
count_table <- count_table[count_table$cases > 1, ]

# complete contingency table for each variant using sample data
count_table_filtered <- add_columns_contigency(count_table)
dim(count_table_filtered)

# simple barplots to overview data
barplot(table(count_table_filtered$control), 
        main = "Number of variants by presence in control group  ", 
        xlab = "Number of controls", ylab = "Number of variants") 
barplot(table(count_table_filtered$cases), 
        main = "Number of variants by presence in cases group  ", 
        xlab = "Number of controls", ylab = "Number of variants") 

# cases/control contrast -----------------------------------------------------------

## Non-targeted / targeted hypothesis
# Fisher test for each variant
fisher_cc <- fisher_test_variants(
  count_table_filtered,"variant","control","control_not_affected","cases",
  "cases_not_affected")

# fisher_cc_power <- fisher_cc
# Power fisher test to check statistical power of each variant
fisher_cc_power <- power_fisher_variants(
  fisher_cc,n_ctrl,n_cases)

# P.adjust calculation (by FDR)
fisher_cc_power$p.adjust <- p.adjust(fisher_cc_power$p.value, method = "fdr")

# Save as tsv the results for each varian
write.table(fisher_cc_power,
            file = paste0(dir,"/results/13_study_candidates/bio_fisher_pow.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)
# load if need
fisher_cc_power <- read.table(
  file = paste0(dir,"/results/13_study_candidates/bio_fisher_pow.txt"),
  stringsAsFactors = F, header = T)

# data visualization ------------------------------------------------

# Subset df by pvalue threshold
pval <- fisher_cc_power[fisher_cc_power$p.value < 0.05, ]
pval <- check_if_variant_multi(pval,"variant",",") # add True/False

# Subset df by padjust threshold (0,09 choosen )
padj <- fisher_cc_power[fisher_cc_power$p.adjust < 0.05, ]
padj <- check_if_variant_multi(padj,"variant",",") # add True/False

# Create a dataframe for barplot visualization
cc_pval <- barplot_counts(pval, "case_control_p.value")
cc_padj <- barplot_counts(padj, "case_control_p.adj")

# Merge both dataframes by rows
to.plot <- rbind(cc_pval, cc_padj)
# to.plot <- cc_pval

# Melt data by statistical P
to.plot.melt <- melt(to.plot, by = "type")
# Add new column indicating test
to.plot.melt$pvar <- apply(to.plot.melt,1,function(row){
  pvar_value <- row[["type"]]
  split_type <- unlist(strsplit(x = pvar_value,split = "_")[[1]])[3]
  return(split_type)
})
# Convert to factor the types depending on data :
# case vs control
to.plot.melt$type <- 
  factor(to.plot$type, 
         levels = c("case_control_p.value","case_control_p.adj"), 
         labels = c("caseVScontrol","caseVScontrol"))

# plots
set.seed(95)
plot_fisher <-
  plot_variants_fisher(df = to.plot.melt,
                       title = paste0("Test de Fisher hipótesis no dirigida"),
                       ylim = 800,
                       base_size = 13,
                       colour1 = "palegreen4",
                       colour2 = "skyblue2")
plot_fisher

plot_fisher2 <- 
  plot_variants_fisher(df = to.plot.melt,
                       title = paste0("Test de Fisher hipótesis no dirigida ",
                                      "filtrado de variantes intraespecíficas" ),
                       ylim = 400,
                       base_size = 13,
                       colour1 = "cadetblue1",
                       colour2 = "coral")
plot_fisher2


plot_fisher_disease <-
  plot_variants_fisher(df = to.plot.melt,
                       title = paste0("Test de Fisher hipótesis dirigida"),
                       ylim = 30,
                       base_size = 13,
                       colour1 = "azure2",
                       colour2 = "firebrick1")
plot_fisher_disease

# multiple comparison fisher analysis -----------------------------------------

# load table
proportions_table_raw <- 
  fread(paste0(dir,"results/10_proportions_tables/",
               "proportions_gt_all_mod_missings.txt"), 
        stringsAsFactors = F, data.table = F)
str(proportions_table_raw)
dim(proportions_table_raw)

# fisher 2x3, padj
fisher_test_case_control_gt<-
  fisher_pairwise_2_x_3(proportions_table_raw,
                                       "variant","ctr_ref",  "ctr_het"  ,
                                       "ctr_hom",  "case_ref" ,"case_het", 
                                       "case_hom")

# check if there's a pairwise contrast that is significant:
sum(fisher_test_case_control_gt$p.value.AA.Aa <= 0.05)  #  0
sum(fisher_test_case_control_gt$p.value.AA.aa <= 0.05)  #  0
sum(fisher_test_case_control_gt$p.value.Aa.aa <= 0.05)  #  0
# Distribution of power fisher test ---------------------------------------

# base data
fisher_cc_power <- read.table(
  file = paste0(
    dir,"/results/13_study_candidates/bio_fisher_pow.txt"),
  header = T)

# change thresholds at will
ggplot(fisher_cc_power[fisher_cc_power$cases > 5 &
                         fisher_cc_power$control < 10,], 
       (aes(x = (fisher.power), fill = (fisher.power)
       ))) +
  geom_bar(stat = "count",
           colour = "black") +
  xlab("Power Fisher value") + ylab("Number of variants") +
  theme_light(base_size = 18) + 
  theme(axis.text.x = element_text(hjust = 0.7),
        axis.title.x = element_text(vjust = -1.0)) +
  scale_fill_manual(values = "paleturquoise") 




# nice p value distribution for fisher variants ---------------------------

# subset df with just id and pvalues

df_pvalues <- fisher_cc_power[,c("variant","p.value")]

library(ggplot2)

threshold <- 2

ggplot(df_pvalues[df_pvalues$p.value < threshold, ], aes(p.value)) + 
  geom_histogram(bins = 15, colour = "black", size = 1) + theme_light() + 
  theme(plot.background = element_blank()) + ylab("num of variants")

boxplot(df_pvalues$p.value[df_pvalues$p.value <= 0.05])

# chose threshold and them subset df

subset <- df_pvalues[df_pvalues$p.value <= 0.01 , ]$variant

# check Q3

write.table(subset, file = "descargas/variants_pval_fish_threshold.txt", 
            row.names = F, col.names = F, quote = F)

# better way to do it

quantile(df_pvalues$p.value, probs = 0.05)


fisher.test(
  matrix(c(0.110619 * 100,100-0.110619 * 100,0.162393 *100,100-0.162393*100), byrow = T, nrow = 2, ncol = 2)
)

# 1000g ph3 = 2,504
# gnomad genomes female = 76,156/2 = 38078
fisher.test(
  matrix(c(0.110619 * 2504, (1-0.110619) * 2504, 0.162393 *118, (1-0.162393)*100), byrow = T, nrow = 2, ncol = 2)
)

fisher.test(
  matrix(c(0.117381 * 38078, (1-0.117381) * 38078, 0.162393 *118, (1-0.162393)*100), byrow = T, nrow = 2, ncol = 2)
)
