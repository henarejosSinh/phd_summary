########### Description #######################################################
# 28 jan 2019

# Apply Hierarchical clustering non supervised for variants

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
library(clusterProfiler)  # for enrichment
# for visualization of dendrograms
library(factoextra)
library(dendextend)
# for getting clusters n
library(NbClust)
library(ggdendro)
display.brewer.all()

colours_25 <- c(
  "dodgerblue2",
  "#E31A1C",
  # red
  "green4",
  "#6A3D9A",
  # purple
  "#FF7F00",
  # orange
  "Gray",
  "gold1",
  "skyblue2",
  "#FB9A99",
  # lt pink
  "palegreen2",
  # light green
  "#CAB2D6",
  # lt purple
  "#FDBF6F",
  # lt orange
  "gray70",
  "khaki2",
  "maroon",
  "orchid1",
  "deeppink1",
  "blue1",
  "steelblue4",
  "darkturquoise",
  "green1",
  "yellow4",
  "yellow3",
  "darkorange4",
  "brown"
)


########### Enviroment #########################################################

save.image("workspace/projects/pof/scripts/rdatas/
# statistical_pos.Rda")
# load("workspace/projects/pof/scripts/rdatas/statistical_pos.Rda.Rda")

# Create variable for working directory
dir <- "workspace/projects/pof/"

# Functions to be used ----------------------------------------------------

# ggplot theme settings
bar.zyx.theme <-  theme_light(base_size = 18) +
  theme(
    axis.text.x = element_text(hjust = 0.65),
    axis.title.x = element_text(vjust = -1.0),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.background = element_rect(
      fill = "white",
      colour = NA,
      inherit.blank = T
    ),
    plot.background = element_rect(colour = "white" , inherit.blank = T),
    panel.grid = element_line(
      colour = "gray88" ,
      arrow = F,
      inherit.blank = T
    ),
    panel.grid.minor = element_line(
      size = rel(0.25),
      arrow = F,
      inherit.blank = T
    ),
    panel.border = element_rect(
      colour = "gray70" ,
      size = rel(1),
      inherit.blank = T
    ),
    axis.line = element_blank()
  )
# load matrix -------------------------------------------------------------

mat_df = read.delim(
  "results/14_filtered_variants_pof/mat_66_het_hom.txt",
  header = T,
  sep = "\t",
  row.names = "v",
  stringsAsFactors = F
)
mat_df <-
  read.delim(
    "results/20_for_paper/nbclusterTrainingTest.csv",
    header = T,
    stringsAsFactors = F,
    sep = ",",
    row.names = "samples"
  )
head(mat_df)

# Delete columns not required
# mat_df <- mat_df[,-151]
mat_df <- mat_df[,-67]

# if needed to convert strings into numerical
mat_df[mat_df == "ref"] <- 0
mat_df[mat_df == "het"] <- 1
mat_df[mat_df == "hom"] <- 2

# transform to numeric
mat_df <-
  as.data.frame(sapply(mat_df, as.numeric), row.names = rownames(mat_df))
mat_df[mat_df != 2 & mat_df != 1 & mat_df != 0]
mat_df[mat_df != 1 & mat_df != 0]
head(mat_df)

# Prepare matrix for calculations ----------------------------------------------------

# transpose if required (variants should go in rows, samples in cols)
mat_df = as.data.frame(t(mat_df))
mat_df

# if you don't want to extract controls from the study
mat_df_2 <- mat_df

# if you want to extract controls
mat_df_2 = mat_df[grep("CONTROL",
                       colnames(mat_df),
                       invert = TRUE ,
                       value = T)]
colnames(mat_df_2)

# if you want to select specific variants (ex top6)
mat_df_2 <- mat_df_2[c(
  "chr11_1017504_G_A,*",
  "chr6_36168628_A_G",
  "chr14_57755564_G_A",
  "chr16_84902483_A_T",
  "chr16_88902199_G_C",
  "chr22_35802661_C_G"
) ,]
# select specific variants (9 no found in 1000g)
variants_not1000g <-
  read.table(
    file = dir,
    "results/14_toconfirm/variants_not_1000g.txt",
    stringsAsFactors = F
  )$V1
mat_df_2 <- mat_df_2[variants_not1000g,]

# similarity matrix ----------------------------------------------------------------------

# do all possible pair combinations
combinations = combn(colnames(mat_df_2), 2)

# generate similarity matrix (samples vs samples)
sim_mat = matrix(
  NA,
  ncol = ncol(mat_df_2),
  nrow = ncol(mat_df_2),
  dimnames = list(colnames(mat_df_2), colnames(mat_df_2))
)

# implemention of Similarity Matching Coefficient in matrix
for (i in 1:ncol(combinations)) {
  sim_mat[combinations[2, i], combinations[1, i]] =
    (sum(mat_df_2[, combinations[1, i]] ==
           mat_df_2[, combinations[2, i]])) / nrow(mat_df_2)
  sim_mat[combinations[1, i], combinations[2, i]] =
    (sum(mat_df_2[, combinations[1, i]] ==
           mat_df_2[, combinations[2, i]])) / nrow(mat_df_2)
}
diag(sim_mat) = 1
dis_mat = 1 - sim_mat # transform distance values to similarity values

fit = hclust(d = as.dist(dis_mat), method = "ward.D2")  # basic plot using ward
plot(fit)

### alternative: jaccard index
# generate similarity matrix
sim_mat = matrix(
  NA,
  ncol = ncol(mat_df_2),
  nrow = ncol(mat_df_2),
  dimnames = list(colnames(mat_df_2), colnames(mat_df_2))
)

for (i in 1:ncol(combinations)) {
  sim_mat[combinations[2, i], combinations[1, i]] =
    (sum(mat_df_2[, combinations[1, i]] != 0  &
           mat_df_2[, combinations[1, i]] ==
           mat_df_2[, combinations[2, i]])) /
    sum(mat_df_2[, combinations[2, i]] != 0 |
          mat_df_2[, combinations[1, i]] != 0)
  sim_mat[combinations[1, i], combinations[2, i]] =
    (sum(mat_df_2[, combinations[1, i]] != 0  &
           mat_df_2[, combinations[1, i]] ==
           mat_df_2[, combinations[2, i]])) /
    sum(mat_df_2[, combinations[2, i]] != 0 |
          mat_df_2[, combinations[1, i]] != 0)
}
diag(sim_mat) = 1
sim_mat[is.nan(sim_mat) | is.na(sim_mat)] <- 0
dis_mat = 1 - sim_mat
fit = hclust(d = as.dist(dis_mat), method = "ward.D2")
plot(fit)


# check distribution of values without controls
dis_mat_2 <-
  dis_mat[!grepl(colnames(dis_mat), pattern = "CONTROL"), !grepl(colnames(dis_mat), pattern = "CONTROL")]
hist(dis_mat_2)
boxplot(dis_mat_2[!grepl(colnames(dis_mat_2), pattern = "CONTROL")])


# nbclust -----------------------------------------------------------------

# jaccard genotypes with controls

### add partition to mat original

# with jaccard 66 (0/1/2) with controls
dis_mat_2 <- dis_mat

# ball or hartigan
index = "ball"

res <-
  NbClust(
    data = dis_mat_2,
    diss = as.dist(dis_mat_2),
    distance = NULL,
    method = "ward.D2",
    index = index,
    min.nc = 2
  )
res$Best.nc
# 3 clusters

res$All.index
res$Best.nc
res$Best.partition

fit = hclust(d = as.dist(dis_mat_2), method = "ward.D2")
cutree(fit, k = 3)
plot(fit, cex = 0.6) # plot tree
rect.hclust(fit, k = 3, border = 2:5) # add rectangle


# hierarchical clustering -------------------------------------------------

# non supervised herarquical clustering
fit = hclust(as.dist(dis_mat), method = "ward.D2")
# save order of hclust
patients_order = colnames(mat_df_2)[fit$order]

# att of nodes for visualization
nodePar <-
  list(
    lab.cex = 0.6,
    cex = 0.7,
    col = "brown1",
    pch = c(NA, 19)
  )

plot(
  as.dendrogram(fit),
  nodePar = nodePar,
  horiz = FALSE,
  xlab = "Height",
  # hang = -1,
  edgePar = list(
    col = c("deepskyblue4", "deepskyblue2"),
    lwd = 2:1
  )
)
# cex = 0.1) # jerarquico, top-bottom

# save dendogram
save(fit, file = "dendrogram.Rda")

# use dendextend package to colour better the dendrogram
dend1 <- color_branches(fit, k = 3)
dend1 <- color_labels(fit, k = 3)
plot(dend1)
set.seed(5665)

# simple separation by class
color_labels(as.dendrogram(fit),
             labels =
               fit$labels[grep(fit$labels, pattern = "CONTROL")],
             col = c("yellow")) %>% plot
color_labels(as.dendrogram(fit),
             labels = fit$labels[grep(fit$labels, pattern = "FOO")],
             col = c("red")) %>% plot
color_labels(as.dendrogram(fit),
             labels = fit$labels[grep(fit$labels, pattern = "FOP")],
             col = c("black")) %>% plot


fviz_dend(
  x = as.dendrogram(fit),
  k = 3,
  # k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
  k_colors = c("#4c4c4c", "#4c4c4c", "#4c4c4c"),
  rect = TRUE,
  # rect_border = c("#2E9FDF", "#00AFBB", "#E7B800"),
  rect_border = c("#FF7F50", "#00AFBB", "#E7B800"),
  rect_fill = TRUE,
  cex = 0.7,
  #main = "Dendrograma - ward.D2",
  xlab = "samples",
  ylab = "height",
  # sub = "Genomic stratification",
  color_labels_by_k = F,
  labels_track_height = 0.2
)



# Tile representation -----------------------------------------------------

tomelt <- mat_df_2
toplot = melt(as.matrix(tomelt))
toplot$Var2 = factor(toplot$Var2, levels = patients_order)
toplot$value = factor(toplot$value, levels = c(0, 1, 2))
ggplot(toplot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(colour = "gray50", size = 0.7) +
  scale_fill_manual(values = c("white", "gray", "black")) +
  theme(axis.text.x = element_text(angle = 90))

# show the top 6 variants
bioinfo_possible_variants <-
  read.table(paste0(dir, "results/14_to_confirm/bioinfo_possible_variants.txt"))
top6 <-
  bioinfo_possible_variants[bioinfo_possible_variants$V14 >= 20,]

# tile plot
ggplot(toplot[toplot$Var1 %in% top6$id,], aes(x = Var2, y = Var1, fill =
                                                value)) +
  geom_tile(colour = "gray50", size = 0.7) +
  scale_fill_manual(values = c("white", "gray", "black")) +
  theme(axis.text.x = element_text(angle = 90))

# study clusters ----------------------------------------------------------

mycl <-
  cutree(fit, h = max(0.55)) # h = max(fit$height/0.5)) # 3 clusters

mycl <- cutree(fit, k = 3)
# cutree returns a vector of cluster membership
# in the order of the original data rows
# examine it
mycl

# examine the cluster membership by it's order
# in the dendrogram
mycl[fit$order]

# get clusters
cluster_A <- names(mycl[mycl == 2])
cluster_B.1 <- names(mycl[mycl == 3])
cluster_B.2 <- names(mycl[mycl == 1])

# save Rda
save(mycl, file = "clusters_dendrogram_top6.Rda")

clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
# myheatcol <- rev(redgreen(75))

# draw the heat map + dendrogram
heatmap.2(
  dis_mat,
  main = "Hierarchical Cluster",
  Rowv = as.dendrogram(fit),
  Colv = NA,
  dendrogram = "row",
  scale = "row",
  col = "heat.colors",
  density.info = "none",
  trace = "none",
  RowSideColors
)

# grab a cluster
dist_mat <- dist_mat[fit$order,]
cluster1 <- dist_mat[mycl == 2,]
View(cluster1)
# or simply add the cluster ID to your data
foo <- cbind(dist_mat, clusterID = mycl)

# examine the data with cluster ids attached, and ordered like the heat map
View(foo[fit$order,])


# change order of clusters ----------------------------------------------------------------

fit = hclust(d = as.dist(dis_mat), method = "ward.D2")

# change order

dd <- fit$order

a <- fit$order[1:17]
b <- fit$order[18:49]
c <- fit$order[50:150]

new_order <- vector()
new_order <- append(a, c)
new_order <- append(new_order, b)

fit$order <- new_order

plot(
  fit,
  hang = -1,
  ylab = "Height",
  xlab = "Samples",
  cex = 0.75
)

dend <- as.dendrogram(fit)
order.dendrogram(dend)

fviz_dend(
  x = fit,
  k = 3,
  # k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
  k_colors = c("#4c4c4c", "#4c4c4c", "#4c4c4c"),
  rect = TRUE,
  # rect_border = c("#2E9FDF", "#00AFBB", "#E7B800"),
  rect_border = c("#FF7F50", "#00AFBB", "#E7B800"),
  rect_fill = TRUE,
  cex = 0.7,
  #main = "Dendrograma - ward.D2",
  xlab = "samples",
  ylab = "height",
  # sub = "Genomic stratification",
  color_labels_by_k = F,
  labels_track_height = 0.2
)

ggdendrogram(fit, rotate = F, size = 2)

dend <- as.dendrogram(fit)

# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"

# check variants ----------------------------------------------------------------


dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)

p <- ggplot(dend_data$segments) +
  geom_segment(aes(
    x = x,
    y = y,
    xend = xend,
    yend = yend
  )) +
  geom_text(
    data = dend_data$labels,
    aes(x, y, label = label),
    hjust = 1,
    angle = 90,
    size = 3
  ) +
  ylim(-0.25, 2.25)
p + theme_light()
print(p)


### check top variants
mat_df_3 <- as.data.frame(t(mat_df_2))
# class: extract by fit$order
fit$order
mat_df_4 <- mat_df_3[fit$order, ]
View(mat_df_4)

mat_df_4$class[1:17] <- "A"
mat_df_4$class[18:118] <- "B"

#       0.48 (   282)  chr2_84932720_A_G
top1 <- mat_df_4[which(mat_df_4[, "chr2_84932720_A_G"] >= 1), ]
top1[, c("chr2_84932720_A_G", "class")]
dnah6_1 <- top1[, c("chr2_84932720_A_G", "class")]
nrow(dnah6_1[dnah6_1$class == "A", ])

#       0.47 (   318)  chr2_84897501_A_G
top2 <- mat_df_4[which(mat_df_4[, "chr2_84897501_A_G"] >= 1), ]
top2[, c("chr2_84897501_A_G", "class")]
dnah6_2 <- top2[, c("chr2_84897501_A_G", "class")]
nrow(dnah6_2[dnah6_2$class == "A", ])

setdiff(rownames(dnah6_1), rownames(dnah6_2))

#       0.43 (   248)  chr2_85059227_C_T
top3 <- mat_df_4[which(mat_df_4[, "chr2_85059227_C_T"] >= 1), ]
top3[, c("chr2_85059227_C_T", "class")]
TRABD2A <- top3[, c("chr2_85059227_C_T", "class")]
nrow(TRABD2A[TRABD2A$class == "A", ])

#       0.43 (    77)  chr14_20002224_C_T
top4 <- mat_df_4[which(mat_df_4[, "chr14_20002224_C_T"] >= 1), ]
top4[, c("chr14_20002224_C_T", "class")]
nrow(top4[, c("chr14_20002224_C_T", "class")])

#       0.39 (   333)  chr2_176829117_G_C
top5 <- mat_df_4[which(mat_df_4[, "chr2_176829117_G_C"] >= 1), ]
top5[, c("chr2_176829117_G_C", "class")]
nrow(top5[, c("chr2_176829117_G_C", "class")])

#     0.38 (    66)  chr16_25239809_G_A
top6 <- mat_df_4[which(mat_df_4[, "chr16_25239809_G_A"] >= 1), ]
top6[, c("chr16_25239809_G_A", "class")]
nrow(top6[, c("chr16_25239809_G_A", "class")])

#       0.37 (    74)  chr10_16932490_G_T
top7 <- mat_df_4[which(mat_df_4[, "chr10_16932490_G_T"] >= 1), ]
top7[, c("chr10_16932490_G_T", "class")]
nrow(top7[, c("chr10_16932490_G_T", "class")])

#  0.35 (    44)  chr16_30393147_C_A
top8 <- mat_df_4[which(mat_df_4[, "chr16_30393147_C_A"] >= 1), ]
top8[, c("chr16_30393147_C_A", "class")]

# 0.34 (    48)  chr16_57732012_G_A
top9 <- mat_df_4[which(mat_df_4[, "chr16_57732012_G_A"] >= 1), ]
top9[, c("chr16_57732012_G_A", "class")]

# 0.34 (   118)  chr11_1017471_C_T,*
top10 <- mat_df_4[which(mat_df_4[, "chr11_1017471_C_T,*"] >= 1), ]
top10[, c("chr11_1017471_C_T,*", "class")]
nrow(top10[, c("chr11_1017471_C_T,*", "class")])

# chr19_10671894_T_G
top11 <- mat_df_4[which(mat_df_4[, "chr16_84902483_A_T"] >= 1), ]
top11[, c("chr16_84902483_A_T", "class")]


# variant padj overrepresented in controls

mat_df_3 <- as.data.frame(t(mat_df), stringsAsFactors = F)
mat_df_4 <-
  mat_df_3[which(mat_df_3[, "chr1_17085995_G_C,GC"] >= 1), ]
mat_df_4[, c("chr1_17085995_G_C,GC", "chrX_151123384_G_A")]


# distribution of accumulation of 66 variants  -----------------------------

# extract order based on dendrogram
a <- fit$order[1:17]
a
A <- fit$labels[a]
A

b <- fit$order[18:118]
b
B <- fit$labels[b]
B

View(mat_df_2[, A])
V <- colSums(mat_df_2[, A])
V
sum(V)

View(mat_df_2[, B])
V2 <- colSums(mat_df_2[, B])
V2
sum(V2)

# convert str to int
mat_df_3 <- as.data.frame(t(mat_df), stringsAsFactors = F)

mat_df_3[mat_df_3 == 2] <- 1

mat_df3 <-
  as.data.frame(sapply(mat_df_3, as.numeric), row.names = rownames(mat_df_3))

mat_df3$total <- rowSums(mat_df3)
mat_df3$total

mat_df4 <- mat_df3[!grepl(rownames(mat_df3), pattern = "CONTROL"), ]
mat_df4

# order by subtypes dendrogram

all <- append(a, b)
all
tosave <- t(mat_df4[all, c("chr22_17450952_A_G", "total")])

# write.table(tosave, "../pof/number_of_variants_by_patients.tsv", sep = "\t",
#             row.names = F, col.names = T, quote = F)

# mat_df4 <- read.table(
#   "workspace/projects/pof/results/20_for_paper/number_of_variants_by_patients.tsv",
#   header = T)

ggplot(mat_df4,
       (aes(x = (total), ))) +
  geom_histogram(colour = "black", binwidth = 1) +
  geom_text(
    aes(label = ..count..),
    stat = "count",
    position = position_dodge(width = 1),
    vjust = -0.5,
    size = 6.5
  ) +
  xlab("Number of variants") + ylab("Number of patients") +
  theme_light(base_size = 24) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.title.x = element_text(vjust = -1.0)) +
  scale_fill_manual(values = "paleturquoise") +
  # scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32))
  scale_x_continuous(breaks = c(unique(mat_df4$total)))
# +   scale_x_continuous(breaks = c(seq(0,30)))


save.image(file = "/home/ihenarejos/workspace/matrix.Rda")

csv <-
  read.table(
    "../pof/results/19_nbclusters for weka/nbclusters_66_3v_AvsBvsC.csv",
    stringsAsFactors = F,
    header = T,
    sep = ","
  )
csv$X

csv$class[csv$class == "B"] <- "BC"
csv$class[csv$class == "C"] <- "BC"
csv$class

write.table(
  "../pof/results/19_nbclusters for weka/nbclusters_66_3v_AvsBC.csv",
  x = csv,
  quote = F,
  row.names = F,
  col.names = T,
  sep = ","
)


# which variants x patients -----------------------------------------------

samples = read.table("results/20_for_paper/nbclusters_66_3v_AvsBvsC.csv", sep = ",",
                    header = T)
head(samples)
colnames(samples)
samples[samples[,"chr2_84897501_A_G"] == "het" ,] 
samples[samples[,"chr2_84932720_A_G"] == "het" | samples[,"chr2_84932720_A_G"] == "hom" ,] 
samples[samples[,"chr2_85059227_C_T"] == "het" | samples[,"chr2_85059227_C_T"] == "hom" ,] 
samples[samples[,"chr11_1017504_G_A.."] == "het" | samples[,"chr11_1017504_G_A.."] == "hom" ,] 
samples[samples[,"chr16_25239809_G_A"] == "het" | samples[,"chr16_25239809_G_A"] == "hom" ,] 
