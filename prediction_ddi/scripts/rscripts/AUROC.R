########### Description #######################################################
# 15FEB2021

# DDIs AUROC
# Ismael Henarejos Castillo
# ihc.europa@gmail.com
# ismael.henarejos@ivirma.com

#working directory
# dir=('/home/ihenarejos/workspace/')
# setwd(dir)
# setwd("/home/ihenarejos/workspace/")

#clean global enviroment if needed
rm(list = ls())

# libraries
# install.packages("ROCR")
# install.packages("cvAUC")
library(data.table)
library(ROCR)
# library(cvAUC)


# image rda
# save.image("scripts/rdata/auroc.Rda")
# load("scripts/rdata/auroc.Rda")


# functions ---------------------------------------------------------------

csv_roc_auc <- function(name_object){
        file = paste0(name_object,"_results.csv", collapse = "")
        df = read.csv(paste0("results/11_predicted_ddis/", file, collapse = ""))
        print(head(df))
        res = roc_auc(df)
        print(res)
        return(res)
}

roc_auc <- function(df){
        pred = prediction(df$TC, df$DESCRIBED)
        perf = performance(pred, 'auc')
        return(perf@y.values)
}
# load data ---------------------------------------------------------------


# 1. Install the ROCR package; create and save a .csv file containing the 
# columns ‘predictions’ (containing the TC scoring from
# the predictor) and ‘labels’ (containing the values of 1 or 0 for
#true positives or false positives, respectively).

# 2. In the R console, load ROCR:
# ALREADY LOADED

# 3. Read data:


# IPF ---------------------------------------------------------------------
# set.seed(seed = sample(1:10000000) )
ipf = fread("results/13_results/ipf_labeled.tsv")
head(ipf)
sum(ipf$DESCRIBED == 1)
# 4. Calculate ROC and list AUC values:
# "labels = DESCRIBED"
pred <- prediction(ipf$TC, ipf$DESCRIBED)
perf <- performance(pred, 'auc')
perf <- performance(pred, 'acc')
perf <- performance(pred, 'fpr')
perf <- performance(pred, 'err')
perf <- performance(pred, 'tpr')
perf <- performance(pred, 'tnr')
perf@y.values
rm(ipf)

# PLOT V1
# # cvAUC
# out = cvAUC::cvAUC(ipf$TC, ipf$DESCRIBED, folds = 10)
# #Plot fold AUCs
# plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")
# #Plot CV AUC
# plot(out$perf, col="red", avg="vertical", add=TRUE)

# PLOT V2
roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
roc_ROCR <- performance(pred, measure = "sens", x.measure = "spec") # 
roc_ROCR <- performance(pred, measure = "prec", x.measure = "rec")

# specificy title with main
plot(roc_ROCR, main = "ROC curve IPF", colorize = F)
plot(roc_ROCR, main = "Precision/recall", colorize = F)
plot(roc_ROCR, main = "Sensitivity/Specificity", colorize = F) #trp vs tnr
abline(a = 0, b = 1)


# Chem. str ---------------------------------------------------------------

chems = fread("results/13_results/chems_labeled.tsv")
head(chems)
sum(chems$DESCRIBED == 1)

pred <- prediction(chems$TC, chems$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values

rm(chems)
roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve MACCS", colorize = F)
abline(a = 0, b = 1)


# ADE ---------------------------------------------------------------------

perf_ade = fread("results/13_results/ade_labeled.csv")
sum(perf_ade$DESCRIBED == 1)

pred <- prediction(perf_ade$TC, perf_ade$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values
rm(perf_ade)

# PLOT
roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve ADE", colorize = F)
abline(a = 0, b = 1)
# out <- cvAUC::cvAUC(perf_ade$TC, perf_ade$DESCRIBED, folds = 10)
# out$fold.AUC
# 
# #Plot fold AUCs
# plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")
# plot(out$perf, col="red", avg="vertical", add=TRUE)

# Targets -----------------------------------------------------------------

# perf_targets = csv_roc_auc("targets")

perf_targets = fread("results/13_results/targets_labeled.tsv")
head(perf_targets)
sum(perf_targets$DESCRIBED == 1)

pred <- prediction(perf_targets$TC, perf_targets$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values
rm(perf_targets)

roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve Targets", colorize = F)
abline(a = 0, b = 1)

# kegg, interactome -------------------------------------------------------


# KEGG
perf_kegg = fread("results/13_results/kegg_labeled.csv")
head(perf_kegg)
sum(perf_kegg$DESCRIBED == 1)

pred <- prediction(perf_kegg$TC, perf_kegg$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values
rm(perf_kegg)

roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve KEGG", colorize = F)
abline(a = 0, b = 1)

# interactome
perf_int = fread("results/13_results/interactome_labeled.tsv")
head(perf_int)
sum(perf_int$DESCRIBED == 1)

pred <- prediction(perf_int$TC, perf_int$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values  #0.84
rm(perf_int)

roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve Interactome", colorize = F)
abline(a = 0, b = 1)


# pca ---------------------------------------------------------------------

# perf_pca = csv_roc_auc("pca")

pca = fread("results/13_results/pca_labeled_withkegg.tsv", 
               stringsAsFactors = F, sep = "\t", header = T)
head(pca)
colnames(pca)[2] = "score"
pca = pca[order(-pca$score),]
head(pca)
nrow(pca[pca$DESCRIBED == "1", ])

summary(pca[ pca$DESCRIBED == 1, ]$score )
pred <- prediction(pca$score, pca$DESCRIBED)
perf <- performance(pred, 'auc')
perf@y.values

# OTHER PERF. values
perf <- performance(pred, 'acc')
perf <- performance(pred, 'spec') # fpr
perf <- performance(pred, 'tpr')
rm(pca)

roc_ROCR <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve PCA", colorize = F)
abline(a = 0, b = 1)

#### MAYBE EXTRACT TP/FN Values to do art analysis?

# TO PLOT
# out <- cvAUC::cvAUC(pca$TC, pca$DESCRIBED, folds = 10)
# out$fold.AUC

#Plot fold AUCs
plot(out$perf, col="grey82", lty=3, main="10-fold CV AUC")

#Plot CV AUC
plot(out$perf, col="red", avg="vertical", add=TRUE)

# transform pca values
# formulae > value = value-max(values)/max-min
max_val = max(pca$score)
min_val = min(pca$score)
range = max_val - min_val
pca$score_trans = apply(pca, 1, function(x){
        value = as.numeric(x[["score"]])
        # print(value)
        new_value = (value - min_val)/range
        return(new_value)
})

summary(pca$score_trans)
head(pca)
pca$score = pca$score_trans
pca$score_trans = NULL
pca = pca[,c("DDI","score_trans", "DESCRIBED", "EFFECT", "ART", "ART_DRUG")]

write.table(pca, file = "results/13_results/pca_candidates_transformed_int.tsv", quote = F, 
            row.names = F,col.names = T, sep = "\t")
# check matrices ------------------------------------------------------------------

df_1 <- fread("Downloads/ipf_df.csv")
df_2 <- fread("Downloads/09_matrix_product/ipf_df.csv")

all(df_1 == df_2)
