########### Description #######################################################
# 29 jan 2021

# pca analysis for ddis

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################
# install.packages(c("FactoMineR", "factoextra"))

library(FactoMineR)
library(factoextra)
library(data.table)
# library(ggplot2)
# library(reshape2)
# library(gsrmUtils)  # need to perform pca
# library(pca3d) # for 3d pca

# sessionInfo()

########### Enviroment #########################################################

# save.image("scripts/rdata/pca_results.Rda")
# load("workspace/projects/pof/scripts/rdatas/statistical_pos.Rda.Rda")

# pca analysis on DDIS data -------------------------------------------

ipf <- read.table("results/11_predicted_ddis/ipf_candidates.csv", stringsAsFactors = F,
                  header = T, sep = ",")
ade <- read.delim("results/11_predicted_ddis/ade_candidates.csv", stringsAsFactors = F,
                   header = T, sep = ",")
chems <- read.table("results/11_predicted_ddis/chems_candidates.csv", stringsAsFactors = F,
                  header = T, sep = ",")
targets <- read.table("results/11_predicted_ddis/targets_candidates.csv", stringsAsFactors = F,
                    header = T, sep = ",")
kegg <- read.table("results/11_predicted_ddis/kegg_candidates.csv", stringsAsFactors = F,
                  header = T, sep = ",")
interactome <- read.table("results/11_predicted_ddis/interactome_candidates.csv", 
                          stringsAsFactors = F, header = T, sep = ",")

head(ipf)
colnames(ipf)[2] = "tc_ipf"
head(ade)
colnames(ade)[2] = "tc_ade"
head(chems)
colnames(chems)[2] = "tc_maccs"
head(targets)
colnames(targets)[2] = "tc_targets"
head(kegg)
colnames(kegg) = c("ddi","tc_kegg")
head(interactome)
colnames(interactome) = c("ddi","tc_int")

# check specific values
ipf[ipf$ddi == "DB00001_DB00002",]
ade[ade$ddi == "DB00001_DB00002",]

# merge dataframes
df.1 <- merge(ade, ipf, by = "ddi")
rm(ade,ipf)
df.2 <- merge(df.1, chems, by = "ddi")
rm(chems,df.1)
df.3 <- merge(df.2, targets, by = "ddi")
rm(targets,df.2)
df.4 <- merge(df.3, kegg, by = "ddi")
rm(kegg,df.3)
df.5 = merge(df.4, interactome, by = "ddi")
rm(interactome, df.4)

head(df.5)
# colnames(df.5)[1:4] = c("ddi","tc_ipf","tc_ade","tc_chems")
# View(df.5)

# pca and eigenvalues filter ----------------------------------------------


# As described in previous sections, the eigenvalues measure the amount of 
# variation retained by each principal component. Eigenvalues are large for the
# first PCs and small for the subsequent PCs. That is, the first PCs corresponds 
# to the directions with the maximum amount of variation in the data set.
# 
# We examine the eigenvalues to determine the number of principal components to be
# considered. The eigenvalues and the proportion of variances (i.e., information) 
# retained by the principal components (PCs) can be extracted using the function 
# get_eigenvalue() [factoextra package].

# set rownames
rownames(df.5) = df.5$ddi
df.5$ddi = NULL 

# calculate pca
res.pca <- PCA(df.5, graph = FALSE)
head(res.pca$eig)

# filter factors with eigenvalue > 1
# Eigenvalues can be used to determine the number of principal components to retain after PCA (Kaiser 1961):
# An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data. This is commonly used as a cutoff point for which PCs are retained. This holds true only when the data are standardized.

# calculate eigenvectors
eig.val <- get_eigenvalue(res.pca)
eig.val
# after checking which dim had eigenvalues > 1, we extract those dimensions:
# extract individual scores
head(res.pca$ind$coord[,"Dim.1"])
# inds = ddi scores(rows)
# loadings = variables(cols)
res.pca$var$coord

(res.pca$ind$coord[,"Dim.1"])
# save converted results
write.table(res.pca$ind$coord[,"Dim.1"], file = "results/13_results/pca_candidates_with_int.csv", quote = F, 
            row.names = T,col.names = T, sep = ",")
rm(res.pca)


# pca scores distribution for described ddis ------------------------------
library(data.table)
# load pca scores
pca_scores_all = fread("results/13_results/pca_transformed_scores_int.tsv", 
                    sep = "\t", data.table = F)

described_scores = pca_scores_all[pca_scores_all$DESCRIBED == 1,]
head(described_scores)

par(cex.lab = 1.25)
par(cex.axis = 1.25)

# boxplot and quantiles
box = boxplot(described_scores$score_trans, data=described_scores, xlab = "Described DDIs", ylab = "Score" )
quan = quantile(described_scores$score_trans, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
quan[[1]][1]
max_rel = box$stats  ## maximum relative = 0.9374336

sum(described_scores$score_trans >=  0.9374336) # 765
sum(described_scores$score_trans >=  0.741810782) # 29,251
29251/nrow(described_scores)

pca_scores = fread("results/13_results/pca_transformed_pred4.tsv", 
                       sep = "\t", data.table = F)
sum(pca_scores$TC >= max_rel[5]) # 13
# Q3 0.741810782
sum(pca_scores$TC >= 0.741810782) # 2991
# P90 0.808043496 
sum(pca_scores$TC >= 0.808043496 ) # 537
# P95 0.838418116 
sum(pca_scores$TC >= 0.838418116  ) # 228
# P99 0.915929179
sum(pca_scores$TC >= 0.915929179  ) # 25


write.table(pca_scores[pca_scores$TC >= 0.7,], 
            file = "results/13_results/pca_pred_4_filter.tsv", quote = F, 
            row.names = F,col.names = T, sep = "\t")


# generate quartile proportions  ------------------------------------------

range_to = seq(0,1,by = 0.01)

temp = lapply(range_to, function(i){
  quan = quantile(described_scores$score_trans, probs = c(i))
  # quan is named numeric so we need to slice to extract the value
  described_value = round((sum(described_scores$score_trans >= quan[[1]][1]))/nrow(described_scores), digits = 2)
})
temp2 = lapply(range_to, function(i){
  quan = quantile(described_scores$score_trans, probs = c(i))
  # quan is named numeric so we need to slice to extract the value
  described_value = sum(described_scores$score_trans >= quan[[1]][1])
})
temp3 = lapply(range_to, function(i){
  quan = quantile(described_scores$score_trans, probs = c(i))
  non_described_value = (sum(pca_scores$TC >= quan[[1]][1]))
})
temp4 = lapply(range_to, function(i){
  quan = quantile(described_scores$score_trans, probs = c(i))
  percentile_value = quan[[1]][1]
})
names(temp) = range_to
names(temp2) = range_to
names(temp3) = range_to
names(temp4) = range_to
tores = do.call("rbind", temp)
tores2 = do.call("rbind", temp2)
tores3 = do.call("rbind", temp3)
tores4 = do.call("rbind", temp4)

final = data.frame(cbind(tores4, tores, tores2, tores3 ))

colnames(final) = c("Percentile value (in a range 0-1)", "Percentage of described interactions", "Number of described interactions (of 117,003)", "Number of non-described interactions (of 223,489)")

final$`Percentage of described interactions` = apply(final, 1, function(x){
  n = x[["Percentage of described interactions"]]
  res = paste0(as.numeric(n) * 100, " %")
  return(res)
})
final
write.table("results/13_results/distribution_by_percentiles.tsv", x= final, 
            sep = "\t", quote = F, row.names = T, col.names = T)

