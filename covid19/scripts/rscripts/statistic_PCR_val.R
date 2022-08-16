# 15/04/2021
# author: Antonio Párraga Leo
# email: Antonio.Parraga@ivirma.com

rm(list = ls())

# Libraries

library(pcr)
library(ggplot2)
library(gsrmUtils)


# functions ---------------------------------------------------------------

default_theme = function(x.axis.angle = 0, legend.pos = "bottom"){
  
  if (x.axis.angle > 0){
    hjust = 1
  }else{
    hjust = 0
  }
  if (abs(x.axis.angle) == 90){
    theme_light()+
      theme(legend.position = legend.pos,
            axis.title = element_text(size=18),
            axis.text = element_text(size=15),
            axis.text.x = element_text(angle=x.axis.angle, vjust = 0.5, hjust = hjust),
            plot.title = element_text(size=22, hjust=0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=14))
  }else if (abs(x.axis.angle) == 45){
    theme_light()+
      theme(legend.position = legend.pos,
            axis.title = element_text(size=18),
            axis.text = element_text(size=15),
            axis.text.x = element_text(angle=x.axis.angle, vjust = 1, hjust = hjust),
            plot.title = element_text(size=22, hjust=0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=14))
  }else{
    theme_light()+
      theme(legend.position = legend.pos,
            axis.title = element_text(size=18),
            axis.text = element_text(size=15),
            plot.title = element_text(size=22),
            legend.title = element_blank(),
            legend.text = element_text(size=14))
  }
}

# calculate relative expression -------------------------------------------

# Load data target genes
gene = 'SLC2A3'

qpcr_gene <- read.delim(file = paste0("results/8_validation/", gene, "_reg.csv"),
                        header = T,
                        sep = ";", stringsAsFactors = F)
qpcr_gene

res <- pcr_ddct(df = qpcr_gene[,c(1,2)],
                group_var = qpcr_gene$Condition,
                reference_gene  = "GAPDH",
                reference_group = "Calibrator",
                mode = "same_tube",
                plot = T)
res

prueba <- data.frame("deltaCT" = qpcr_gene[,gene][-1] - qpcr_gene$GAPDH[-1],
                     "condition" = qpcr_gene$Condition[-1])
prueba
calibrator <- qpcr_gene[, gene][1] - qpcr_gene$GAPDH[1]
prueba$ddCT <- prueba$deltaCT - calibrator 
prueba$relative_expression <- 2^(-prueba$ddCT)
prueba$group <- factor(x = qpcr_gene$Condition[-1], levels = c("PF","ESE","LSE"))

prueba

shapiro.test(prueba$deltaCT)
summary(aov(deltaCT ~ condition, data = prueba))
pairwise.t.test(x = prueba$deltaCT, g = prueba$condition, p.adjust.method = "BH")
prueba$relative_expression = log2(prueba$relative_expression)


toplot <- aggregate(prueba[,c(1,3,4)], list(Phase=prueba$condition), mean)
toplotSD <- aggregate(prueba[,c(1,3,4)], list(Phase=prueba$condition), sd)
toplotSEM <- aggregate(prueba[,c(1,3,4)], list(Phase=prueba$condition), function(x) sd(x)/sqrt(length(x)))
toplotSCI = aggregate(prueba[,c(1,3,4)], list(Phase=prueba$condition), function(x) (sd(x)/sqrt(length(x))) * (qt(0.95/2 +0.5, length(x) -1)))

toplot_final <- data.frame("Phase" = toplot$Phase,
                           "relative_expression" = toplot$relative_expression,
                           "SD" = toplotSD$relative_expression,
                           "SEM" = toplotSEM$relative_expression,
                           "ci" = toplotSCI$relative_expression)
toplot_final$Phase <-  factor(x = toplot_final$Phase, levels = c("PF","ESE","LSE"))

ggplot(data = toplot_final, (aes(x = Phase, y = relative_expression, fill = Phase)))+
  geom_bar(stat = "identity", width = 0.2, color = "black") + 
  scale_fill_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MSE = "chartreuse4",
                               LSE = "lightpink3"))+
  geom_errorbar(data = toplot_final,aes(x = Phase, ymin = relative_expression, ymax = relative_expression+SEM),
                width=.05) + 
  scale_x_discrete(labels = c("PF" = "PF\n(n = 2)","ESE" = "ESE\n(n = 2)","LSE" = "LSE\n(n = 3)"))+
  default_theme() +
  xlab("\nEndometrial phases\n")+
  ylab(paste(gene,"mRNA Relative expression"))+
  theme(legend.position = "none")

toplot_final$A <- "A"

ggplot(data = toplot_final, (aes(x = Phase, y = relative_expression, group = A, color = Phase)))+
  geom_line() +
  geom_point()+
  scale_fill_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MS = "chartreuse4",
                               LSE = "lightpink3"))+
  geom_errorbar(data = toplot_final,aes(x = Phase, ymin = relative_expression-SD, ymax = relative_expression+SD),
                width=.1) + 
  default_theme() +
  xlab("Endometrial Phase")+
  ylab(paste(gene,"mRNA Relative expression"))

colnames(toplot_final) = c("time", 'expression_value', 'sd', 'se', 'ci', 'group')
toplot_final$gene = gene
toplot_final

ggplot(toplot_final[toplot_final$gene==gene,], aes(x = time, y = expression_value, color = time))+
  geom_point(size = 3, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = expression_value - sd, ymax = expression_value + sd),  # sd should be ci
                width = 0.2, position = position_dodge(width = 0.5))+
  scale_color_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MSE = "chartreuse4",
                                LSE = "lightpink3"))+
  xlab("") + ylab("") + 
  ggtitle("")+
  theme_void()+
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none")



# load all log o fold gene expression and plot ----------------------------

qpcr_all <- read.delim(file = paste0("results/8_validation/covid_expval - Sheet1.csv"),
                        header = T,
                        sep = ",", stringsAsFactors = F)
qpcr_all

toplot_sum = toplot_p = summarySE(qpcr_all, measurevar = "expression_value", groupvars = c("gene", "time"))
toplot_sum

# check if summarySE worked correctly:
mean(qpcr_all[qpcr_all$gene == 'COBL' & qpcr_all$time == 'PF' ,]$expression_value)

gene = 'COBL'
toplot_sum$time = factor(toplot_sum$time, levels = c("PF", "ESE", "MSE", "LSE"))
toplot_sum$group <- "A"

plot_expression_values = function(toplot, gene){
  ggplot(toplot[toplot$gene==gene,], aes(x = time, y = expression_value, color = time, group=group))+
    geom_line() +
    geom_point(size = 3, position = position_dodge(width = 0.5))+
    geom_errorbar(aes(ymin = expression_value - sd, ymax = expression_value + sd),  
                  width = 0.2, position = position_dodge(width = 0.5))+
    scale_color_manual(values = c(PF = "coral4", ESE = "darkgoldenrod3", MSE = "chartreuse4",
                                  LSE = "lightpink3"))+
    # geom_hline(yintercept = limits, linetype = 2, color = "grey60")+
    xlab("") + ylab("") +
    ggtitle("")+
    theme_void()+
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          legend.position = "none")
}
# 7.3 , 4.51

toplot_sum
plot_expression_values(toplot_sum, gene)


