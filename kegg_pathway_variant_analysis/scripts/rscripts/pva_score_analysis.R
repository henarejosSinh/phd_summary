# 9/9/2020
# ihc.europa@gmail.com

# Script to try statistical methods for pva score analysis

library(data.table)
library(gsrmUtils)

rm(list = ls())
# load table, dcast it and create ed ----------------------------------------------------
# analysis_score_out_qual_25feb21.tsv
# pva adaptation first batch of results in 03_pva_score
d <- fread(file = "results/06_pva_adaption/test_files/outdir/out_genepy.txt", 
           header = TRUE,  data.table = F)
head(d) # data melted (1 line equals to 1 sample in a biological context and 1 score)
dat <- dcast(as.data.table(d), id ~ sample ) # transforms dataframe keeping id as rows and samples as columns
dat = as.data.frame(dat)
head(dat)

# fix dat
rownames(dat) = dat$id
head(dat)
dat = dat[, -1]
head(dat)
colnames(dat) # check the name of columns to ensure that UNKOWN_150 = CONTROL_150

# create ed
ed = data.frame(samples = colnames(dat), stringsAsFactors = F)
ed$condition = ifelse(grepl("CONTROL", ed$samples) > 0, "Control", "Case")
rownames(ed) = ed$samples

# check if we have null subpathways
sum(rowSums(dat) == 0)
dat[rowSums(dat) == 0, ]
aux_zeroes <- dat[rowSums(dat) == 0, ]
# JAN 28/22 1793 sp
# FEB 2/22 (Alex's fixed script) 1263 sp

# check subpathways not compromised (signal == 1 for all samples)
sum(rowSums(dat) == ncol(dat))
dat[rowSums(dat) == ncol(dat), ]
aux_ones <- dat[rowSums(dat) == ncol(dat), ]
# 59 all 1; JAN 28/22 400 sp ; FEB 2/22 43

# fix 0 and 1 subpathways
#  we have 0 because inhibitions
# we have 1 because either we don't have the gene in the VCF or the variants
# are sinonymous
dat <- dat[!rowSums(dat) == 0, ]
nrow(dat)
dat <- dat[!rowSums(dat) == ncol(dat), ]
nrow(dat)

ids = rownames(dat)
# from ~5100, 3811 subpathways remain
# JAN 28/22 from ~5300, 2981 sp remain
# FEB 2/22 from ~5100 sp, 3868 remain

save(dat, ed, file = "scripts/rscripts/score_ed_FEB2_22.Rda")
# statis models -----------------------------------------------------------

str(dat)

# 1st ALTERNATIVE (DOESN'T WORK)
# apply a linear model like limma for DE
# DE <- diffExprAnalysis(dat[,-1], ed, condition = "group")
# any(DE$`case-control`$adj.P.Val <= 0.05)
# FALSE

# t.test   (DOESN'T WORK)
# res = do.call("rbind", lapply(1:nrow(dat), function(i){
#   cat(i, "\n")
#   cases = as.numeric(dat[i, ed$condition == "Case"])
#   control = as.numeric(dat[i, ed$condition == "Control"])
#   if (cases != control) { # this does not work
#   test = t.test(cases, control)
#   }
#   data.frame(mean_case = test$estimate[1], mean_control = test$estimate[2], 
#              t = test$statistic, pvalue = test$p.value, stringsAsFactors = F)
# }))

dat2 = as.data.frame(dat)

# wilcox test
sol = do.call("rbind", lapply(1:nrow(dat2) , function(i){ # nrow(dat),
  cat(i, "\n")
  cases = as.numeric(dat2[i, ed$condition == "Case"])
  control = as.numeric(dat2[i, ed$condition == "Control"])
  test = wilcox.test(cases, control)
  data.frame(mean_case = mean(cases), mean_control = mean(control), 
             W = test$statistic, pvalue = test$p.value, stringsAsFactors = F)
}))
hist(sol$pvalue)
res = sol
View(res)

# if we have NA = 1
# res$pvalue[is.na(res$pvalue)] = 1
# rownames(res) = rownames(dat)
# res$FDR = p.adjust(res$pvalue, method = "fdr")
# hist(res$pvalue)
# hist(res$FDR)

# OR omit
rownames(res) = ids
head(res)
res = na.omit(res)
nrow(res)
# 3811/ 3732 with 25feb21 dataset, pathways that we can compare] 
# JAN 28/22 2981 remaining sps
res$FDR = p.adjust(res$pvalue, method = "fdr")
hist(res$pvalue)
hist(res$FDR)

# how many significant?
sum(res$pvalue <= 0.05)
sum(res$FDR <= 0.05)
sig_res <- res[res$pvalue <= 0.05, ]
sig_res
# 183 W, 136 t.test, 169 FEB 2/22
# 116 genepy approach

# adjust pvalues by pathways ----------------------------------------------

# padjust by pathway
res$subpath <- rownames(res)
res$pathway <- unlist(apply(res, 1, function(x){
  val <- x[["subpath"]]
  path <- unlist(strsplit(val, split = "_")) [[1]][1]
  return(path)
  }))
res$pathway
head(res)
# split dataframe into a list of lists where each individual list is a pathway
# containing the corresponding subpathways and data for that pathway
aux <- split(res, res$pathway)
head(aux)
str(aux)
head(aux$hsa03320)
length(aux) # FEB 2/22 112 pathways to evaluate 

res_fdr <- do.call("rbind", lapply(1:length(aux), function(j) {
  cat(j, "\n")
  # print(aux[[j]])
  # print(aux[[j]]["pvalue"])
  pvals = data.frame( pathway = aux[[j]]["pathway"], 
                      pvalue = aux[[j]]["pvalue"])
  fdr = p.adjust(pvals$pvalue, method = "fdr")
  # print(pvals)
  # print(fdr)
  data.frame( pathway = aux[[j]]["pathway"], 
              mean_case =  aux[[j]]["mean_case"],
              mean_control = aux[[j]]["mean_control"], 
              pvalue = aux[[j]]["pvalue"], padjust = fdr, stringsAsFactors = F)
}))

dim(res_fdr)
hist(res_fdr$padjust)
sum(res_fdr$padjust <= 0.05 )
res_fdr[res_fdr$padjust <= 0.05, ]
View(res_fdr[res_fdr$padjust <= 0.05, ])
sum(res_fdr$padjust <= 0.1 )
View(res_fdr[res_fdr$padjust <= 0.1, ])

sig_res <- res_fdr[res_fdr$pvalue <= 0.05,]
sig_res

# Calculate fc-like value:

# res_fdr$fc = apply(res_fdr, 1, function(x){
#   v = (as.numeric(x["mean_case"]) - as.numeric(x["mean_control"])) / as.numeric(x["mean_control"])
#   return(v)
# })

res_fdr$fc = apply(res_fdr, 1, function(x){
  v = (as.numeric(x["mean_case"]))/ as.numeric(x["mean_control"])
  return(v)
})

View(res_fdr[res_fdr$padjust <= 0.05,])

# res_fdr$fc = apply(res_fdr, 1, function(x){
#   v = (as.numeric(x["mean_control"]) - as.numeric(x["mean_case"])) / as.numeric(x["mean_case"])
#   return(v)
# })

res_fdr$fc_inverse = apply(res_fdr, 1, function(x){
  v = (as.numeric(x["mean_control"]))/ as.numeric(x["mean_case"])
  return(v)
})

res_fdr[res_fdr$padjust <= 0.05,]

View(res_fdr[res_fdr$padjust <= 0.05,])
View(res_fdr)

res_fdr$pathway = rownames(res_fdr)

write.table(res_fdr, file = "results/05_sigSubPaths_interpretation/genepy_adaptation_res_JAN28_22.tsv"
            ,  quote = F, sep = "\t", row.names = F, col.names = T)

write.table(res_fdr[res_fdr$padjust <= 0.05,], file = "results/05_sigSubPaths_interpretation/genepy_adaptation_res_JAN28_22_fdr.tsv"
            ,  quote = F, sep = "\t", row.names = F, col.names = T)
