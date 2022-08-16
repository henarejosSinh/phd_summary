

library(gsrmUtils)

load("results/entire_database.rda",verbose =T) 
files = dir("expression_data/", full.names = T, pattern = ".RData")

dats2 = lapply(files, function(f){
  load(f)
  return(dat)
})
common_genes = Reduce("intersect", lapply(dats2, rownames))
length(common_genes)
# [1] 13437
genes_of_interest = "PGR"
genes_of_interest%in%common_genes

dats2 = do.call("cbind", lapply(dats2, function(x){
  x[common_genes,]
}))
dim(dats2)

# Experimental designs
eds2 = do.call("rbind", lapply(files, function(f){
  load(f)
  ed = data.frame(sample = rownames(ed),
                  time = ed$time,
                  experiment = gsub("data//", "", unlist(strsplit(f, split = "_", fixed = T))[2], fixed = T),
                  stringsAsFactors = F)
  return(ed)
}))
rownames(eds2) = eds2$sample
# Fix times
eds2$time[eds2$time=="LH+2"] = "ESE"
eds2$time[eds2$time=="LH+8"] = "MSE"
eds2$time[eds2$time=="proliferative"] = "PF"
eds2$time[eds2$time=="mid-secretory"] = "MSE"
eds2$time[eds2$time=="Proliferative"] = "PF"
eds2$time[eds2$time=="Secretory"] = "MSE"
eds2$time[eds2$time%in%c("-1", "-3", "-12", "-7", "-14", "-8", "-5", "-4")] = "PF"
eds2$time[eds2$time%in%c("1", "2", "3", "4")] = "ESE"
eds2$time[eds2$time%in%c("6", "7")] = "MSE"
eds2$time = factor(eds2$time, levels = c("PF", "ESE", "MSE", "LSE"))

# Remove Altmae and Bradley outliers
dats2 = dats2[, !colnames(dats2) %in% c("GSM2593180", "GSM2593181", "GSM742078")]
eds2 = eds2[!rownames(eds2) %in% c("GSM2593180", "GSM2593181", "GSM742078"), ]

save(eds2,file = "results/Network_not_sebastianleon/ed_CTCF.RData")


# Remove experiment effect
dats3 = removeBatchEffect(dats2, batch = eds2$experiment)

# Differential analysis
DEGs = diffExprAnalysis(dat = dats3, ed = eds2, condition = "experiment")
sigs = lapply(DEGs, function(deg){
  rownames(deg)[deg$adj.P.Val<=0.05]
})
to_remove = unique(unlist(sigs))
length(to_remove) # 0 --> No genes to be removed

# Evaluate expression
dat = data.frame(sample = colnames(dats3), expression_value = dats3[genes_of_interest,], stringsAsFactors = F)
dat$time = eds2[dat$sample, "time"]
# limits = quantile(dats3, c(0, 0.1, 0.5, 1))
toplot = summarySE(dat, measurevar = "expression_value", groupvars = c("time"))

p = ggplot(toplot, aes(x = time, y = expression_value, color = time))+
  geom_point(size = 3, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = expression_value - ci, ymax = expression_value + ci),
                width = 0.2, position = position_dodge(width = 0.5))+
  scale_color_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MSE = "chartreuse4",
                                LSE = "lightpink3"))+
  # geom_hline(yintercept  = limits, linetype = 2, color = "grey60")+
  xlab("") + ylab("") +
  ggtitle("")+
  #theme_void()+
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none") ## Ponerlo bonito
p

# Evaluate differences
anova = aov(expression_value ~ time, dat)
summary(anova)
pairwise.t.test(dat$expression_value, dat$time, p.adjust.method = "fdr")



# Correlation -------------------------------------------------------------

res <- do.call("rbind", lapply(rownames(dats3), function(x){
  correlation = cor.test(dats3["CTCF", ], dats3[x,])
  data.frame(Node1 = "CTCF",
             Node2 = x,
             correlation = correlation$estimate,
             p_value = correlation$p.value,
             stringsAsFactors = F)
}))

res$p_adjust = p.adjust(res$p_value)
res = res[res$Node2!="CTCF",]

res[res$Node2=="PGR",]
res[res$Node2=="ESR1",]
res[res$Node2=="ESR2",]



toplot = as.data.frame(t(dats3[c("CTCF", "PGR"),]))

ggplot(toplot, aes(x=CTCF, y=PGR))+
  geom_point()+
  default_theme()+
  geom_smooth(method = "lm")


sum(res$p_adjust<=0.05 & abs(res$correlation)>= 0.5)
sum(res$p_adjust<=0.05 & abs(res$correlation)>=0.65)


ggplot(res, aes(correlation))+
  geom_density(color = "black") +
  default_theme() +
  scale_x_continuous(breaks=seq(-0.7, 0.7, 0.2)) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red")
  #geom_boxplot()

summary(res$correlation)
a <- boxplot(res$correlation)
a
res_final <- res[abs(res$correlation)>=0.65,]
dim(res_final)


sum(duplicated(res$Node2)) # 0



# Network of threshold 0.5 with genes related to P & E --------------------

res_more_information <- res[res$p_adjust <= 0.05 & abs(res$correlation)>=0.5,]
dim(res_more_information)
target_dorothea <- unique(entire_database[entire_database$tf == "CTCF",]$target)

res_more_information_att <- data.frame("Gene" = unique(res_more_information$Node2),
                                        stringsAsFactors = F)
res_more_information_att$PR <- res_more_information_att$Gene%in%genes_PGR
res_more_information_att$ER <- res_more_information_att$Gene%in%genes_ER
res_more_information_att$target <- res_more_information_att$Gene%in%target_dorothea


# Create a ordered correlation vector -------------------------------------

res_sign <- res[res$p_adjust <= 0.05,]

ggplot(res_sign, aes(correlation))+
  geom_density(color = "black") +
  default_theme() +
  scale_x_continuous(breaks=seq(-0.7, 0.7, 0.2)) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "red")
#geom_boxplot()

res_sign <- res_sign[order(res_sign$correlation),]
correlation_vector <- res_sign$correlation 

# Choose the threshold of 2.5 %
threshold <- round(length(res_sign$correlation) * 2.5/100,0)

# Choose the priorizated
correlation_positive <- res_sign[(length(correlation_vector)-(threshold-1)):length(correlation_vector),]
# min(correlation_positive$correlation)
# [1] 0.6562752

correlation_negative <- res_sign[1:threshold,]
# max(correlation_negative$correlation)
# [1] -0.593277

res_filter_2 <- rbind(correlation_negative, correlation_positive)
