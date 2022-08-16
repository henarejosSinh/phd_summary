# 6 apr 2021

# DDI results stacked bar plot
# Ismael Henarejos Castillo
# ihc.europa@gmail.com
# ismael.henarejos@ivirma.com

R.version
sessionInfo()
library(ggplot2)
library(reshape2)
library(RColorBrewer)  # for color palettes
display.brewer.all()


# load files --------------------------------------------------------------

res = read.table(
  "results/14_discussion/supp_ivf_drugs - ivf.tsv",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
head(res)


# stacked barplot ---------------------------------------------------------

res = res[order(as.numeric(res$discovery....), decreasing = T),]
res
top_drugs = res$drug.name[1:10]
top_drugs

melted_res = melt(res[1:10, ], id.vars = "drug.name", measure.vars = c(3, 4))
melted_res$variable = factor(melted_res$variable, levels = c("predicted", "described"))
# melted_res$variable = factor(melted_res$variable, levels = c("predicted_frd_only","described_frd_only"))
melted_res$drug.name = factor(melted_res$drug.name, levels = top_drugs)

str(melted_res)

ggplot(data = melted_res, aes(
  x = drug.name,
  y = value,
  fill = variable,
  label = value
)) +
  geom_bar(stat = "identity") +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  xlab("") + ylab("") + ggtitle("") + # ylim(c(0,25)) + # comment ylim
  theme_light() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.title = element_text(size = 22, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none"
  )

# COLORBLIND PALETTES ------------------------------------------------------

safe_colorblind_palette <-
  c(
    "#88CCEE",
    "#CC6677",
    "#DDCC77",
    "#117733",
    "#332288",
    "#AA4499",
    "#44AA99",
    "#999933",
    "#882255",
    "#661100",
    "#6699CC",
    "#888888"
  )
scales::show_col(safe_colorblind_palette)

# The palette with grey:
cbPalette <-
  c(
    "#999999",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
scales::show_col(cbPalette)

# The palette with black:
cbbPalette <-
  c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )
scales::show_col(cbbPalette)
