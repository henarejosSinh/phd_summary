########### Description #######################################################
# 2 feb 2020

# Circos plot for variants

#clean global enviroment if needed
rm(list = ls())

########## Libraries ##########################################################

library(circlize) # for circos plot
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
# circosplot.Rda-")
# load("workspace/projects/pof/scripts/rdatas/fisher_test_strategy_variants.Rda")

# Create variable for working directory
dir <- "workspace/projects/pof/"


# Functions to be used ----------------------------------------------------

# to make colours more transparent
transparency_col <- function(color, percent = 60, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# example to understand circos -------------------------------------------------

# a simple example where we plot a circos and try to understand the data structure
set.seed(999) # to assure reproducibility
n = 1000 
df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),
                x = rnorm(n), y = runif(n))

circos.par("track.height" = 0.1) # defines % of radious of the circle it will plot a certain track in a sector 
circos.initialize(factors = df$factors, x = df$x) # define sectors to plot and ranges of X (ranges of tracks)

# now, it plots :
circos.track(factors = df$factors, y = df$y, # passes information to y axis of each track to plot something inside/outside tracks
             panel.fun = function(x, y) { # passes coordinates for each factor
               
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), # CELL_META is an object created when we initialize circos which have 
                           # metainformation such as central coordinates on x axis and ylim
                           # uy defines where to write the text. the lesser the mm, the closer to the track
                           CELL_META$sector.index) # for each sector (we pass the index for that factor)  
               
               circos.axis(labels.cex = 0.6) # plots lines on x axis
               
             })

# colors, plot points in sectors
col = rep(c("#FF0000", "#00FF00"), 4) # 8 colors for 8 sectors
circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 16, cex = 0.5) # pch = plotting ‘character’, i.e., symbol to use. 
circos.text(-1, 0.5, "text", sector.index = "a", track.index = 1) # we can define were we want the text, in which sector and in which track of given sector

# initialize circos of variants-----------------------

# dataframe where we have rows as independent variants and for each variant we 
# know which genes are affected and number of cases
bioinfo_possible_variants <-
  fread(file = paste0(dir,"results/14_to_confirm/bioinfo_possible_variants.txt"),
        stringsAsFactors = F, header = T)

# df to plot
df <-
  data.frame(factors = unlist(
    strsplit(bioinfo_possible_variants$gene,split = ","))) 
# we obtain as factors the genes in the base df making sure to get
# all genes (for those variants that affect more than one gene )
df2 <- data.frame(factors = unique(df$factors), x = 0, y = 0) 
# a second dataframe with 3 columns: unique genes from first df, value 0 ranges
# of axis x and y
df$x <- 1 # we create columns in original df with limit of both axis
df$y <- 1
df$x[duplicated(df$factors)] = 2 # where we find duplicated genes in original 
# dataframe, limit of the x axis should be 2
df = rbind(df,df2) # we modify the original dataframe merging by rows the first
#and second dataframe


# add links ---------------------------------------------------------------

# example : can be added one by one
# circos.link("MUC6", c(0,2), "31", c(0,1)) # second c() is always 1 (1 variant), 
# if gene is affected by two independent variants, first c() has to change 

# df to create for links :
# genes, num variants affecting said gene, number of cases affected
dflinks <- data.frame ( genes = bioinfo_possible_variants$gene, 
                        cases = bioinfo_possible_variants$cases)

index = grep(",", dflinks$genes) # we get index in df where there are 
# multiple genes affected by same variant
df_aux = do.call("rbind", lapply(index, function(i){
  data.frame(genes = unlist(strsplit(as.character(dflinks$genes[i]), split = ",")),
             cases = dflinks$cases[i]) 
})) # here we create 
# another dataframe using indexes of former dataframe
dflinks = dflinks[-index,] # delete old rows
dflinks = rbind(dflinks, df_aux)  # bind former dataframe with new dataframe

# Reorder links based on number of cases affecting a said gene
dflinks <- with(dflinks, dflinks[order(genes, decreasing = T),]) # orders first
# dataframe using the order specified in second dataframe

# for genes affected by more than one variant, add to the links dataframe 0,1
#value to indicate where the link has to be ploted for that given gene
dflinks$genes1 <- 0
dflinks$genes1[duplicated(dflinks$genes)] <- 1

# for cases, similar approach, but must sum an extra 1 for each repetion, since 
# we need to plot the links using the number of variants affecting a certain
# group of patients
# plot links taking into consideration:
# a: genes with more than 1 variant affecting them (x axis)
# b: multiple variants affecting the same number of individuals for a given gene
# ex. there are 22 genes affected by variants affecting 12 indvs
dflinks$cases1 <- 0
for (i in 1:(length(dflinks$cases))-1){
  indexes = which(dflinks$cases == dflinks$cases[i]) #
  # get indexes of other genes in df with the same number of variants
  indexes_toeval = (i+1):length(dflinks$cases)  
  # the rest of indexes to evalute in a given iteration
  indexes_2 = indexes[indexes%in%indexes_toeval] 
  # which of the rest coincide with the ones obtained in the first step
  dflinks$cases1[indexes_2] = dflinks$cases1[indexes_2] +1
  # for the ones that coincide, sum +1 to the x range to be used
}

# colors for circos -------------------------------------------------------

# add colors to links as a column in the df base for circos:
n <- 80 # number of colors to generate
qual_col_pals <-
  brewer.pal.info[brewer.pal.info$category == 'qual',] # obtain colors from palette
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
# crete a vector with all colors from palette
# pie(rep(1,n), col=sample(col_vector, n, replace = T)) # a representation of 
# colors obtained in a pie chart
col_vector2 <- rep(col_vector, 2) 

my.colors.forlinks <- colorRampPalette(c("cadetblue", 
                                         "cadetblue4","darkblue"), 12)
# palette to generate
gradient.colors <- (my.colors.forlinks(12)) 
df_colors = data.frame(cases = sort(unique(dflinks$cases)), 
                       # creates a dataframe where each number of cases 
                       # have a certain color assigned to their links
                       col = gradient.colors, stringsAsFactors = F)


# add colors to links dataframe
dflinks$col <- ""
for (i in 1:(length(dflinks$cases))){
  indexes = which(dflinks$cases == dflinks$cases[i])
  dflinks$col[indexes] = gradient.colors[i]
  dflinks$col[indexes] = col_vector2[i]
} # if we were to color use multiple colors for each link

# add colors for each link
dflinks$col = unlist(lapply(dflinks$cases, function(x){
  df_colors$col[df_colors$cases == x]
}))


# plot sections and tracks ------------------------------------------------

# number of groups of patients affected by variants considering number of genes
# affected in each group 
# needs num cases as factors 
# aux <- data.frame(table(bioinfo_possible_variants$cases), y = 1) # 
# problem with this df: 16 variants affecting 12 cases affect more than 16 genes
aux <- data.frame (table(dflinks$cases), y = 1) # once we create df of links we
# use it to count number of genes affected by each group of patients
colnames(aux) <- c("factors", "x", "y") # name columns
aux2 <- data.frame(factors = aux$factors, x = 0, y = 0) # we establish 
# ranges of each axis for these sectors
df <- rbind(df, aux, aux2) # (create new dataframe to add to the one created before)

# we would want to order genes sectors in circos by factor levels, taking into
# account number of cases affecting them
# "order" function get index for condition and reorders whole df according to it
levelstoorder <- # fist you need to extract 
  as.character(unique(dflinks[order(dflinks$cases, decreasing = T),]$genes)) 
# genes affected
num_to_order <- as.character(aux2$factors) # just the groups of patients affected
levels_to_plot <- append(levelstoorder,num_to_order) # append for character vectors
df$factors <- factor(df$factors, levels = levels_to_plot) # we make sure 
# that data for sectors are factors and ordered in the way we want to

# initialize circos

circos.clear() # delete a circos previously created 
circos.par("track.height" = 0.1, cell.padding = c(0.02,0,0.02,0), "canvas.ylim"
           = c(-1.25, 1.1) , # .par defines the way tracks are going to be plotted 
           start.degree = 92) 
# start.degree defines % of radious of the circle it will plot 
# cell padding defines white spaces between sectors

circos.initialize(factors = df$factors, x = df$x) # define sectors to plot and 
#ranges of X. creates circos, but doesn't plot it

circos.track(factors = df$factors, y = df$y, bg.col = df$col, 
             bg.border = "black", # passes information to y plot for each 
             # corresponding track
             panel.fun = function(x, y) {
               if (CELL_META$sector.index %in% num_to_order) { # for sectors of patients group, we define a different space to plot text
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(2, "mm"), 
                             CELL_META$sector.index, facing = "clockwise", niceFacing = T, cex = 0.6)
               } else { # for gene sectors
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), 
                             CELL_META$sector.index, facing = "clockwise", niceFacing = T, cex = 0.6)
                 # circos.axis(labels.cex = 0.1, major.tick = T) # to paint subdivisions, decided not to use it
               }
             })
# paint num sectors
colour_to_use <- rand_color(
  1, hue = "orange", transparency = 0,5)
for (i in 1:nrow(df[df$factors == num_to_order,])){ #for tracks of patients group
  circos.updatePlotRegion( sector.index = num_to_order[i], bg.col = colour_to_use) 
  # paint them differently than gene sectors
}

colour_to_use <- rand_color(
  1, hue = "monochrome", transparency = 0,5)
# paint gene sectors
for (i in 1:nrow(df[df$factors != num_to_order ,])){ #for tracks of patients group
  circos.updatePlotRegion( 
    sector.index = df[df$factors != num_to_order ,]$factors[i], bg.col = "#696969") 

}

# add links loop 
for (i in 1:nrow(dflinks)){
  
  circos.link(dflinks$genes[i], c(dflinks$genes1[i],dflinks$genes1[i] + 1), 
              # defines links painted from which specified part of track to another
              dflinks$cases[i], c(dflinks$cases1[i],dflinks$cases1[i] + 1), 
              h.ratio = 0.2, col = dflinks$col[i]) # % change way links are painted inside circle
              # border = "white", lwd = 0.1, lty = 1) # border of links parameters
}

# highlight sectors -------------------------------------------------------

color_highlight <- col2rgb("lightgreen")
color_highlight <- transparency_col(color_highlight)


# draw background for specified sectors
draw.sector(get.cell.meta.data(sector.index = "12", "cell.start.degree"), 
            # 240 defines background initial value to start plotting...
            get.cell.meta.data(sector.index = "31", "cell.end.degree") 
            # 93 and final value
            , border = "white", lwd = 2, lty = 2, rou1 = 1.2, 
            rou2 = 0.87, col = color_highlight)  # border and size parameters

vector <- c()
for (i in num_to_order){ # create a vector of sectos to highlight
  vector <- append(x = vector,i)
}
vector

highlight.sector(c(vector), track.index = 1 , 
                 # text = "Numero de pacientes afectados por una determinada variante",
                 text = "Number of patients affected by a certain variant",
                 # add text and highlight certain sectors
                 facing = "bending.inside", niceFacing = T, text.vjust = "9mm",
                 cex = 1, col = "#BEBEBEBF", text.col = "gray8", 
                 border = "black", lwd = 0.5,
                 lty = 5)

# transparent colors function ; # add_transparency() is an alternative

# Circos plot using pos and chr -----------------------------------------------------------------------------------

# read txt where for hg19 we have initial pos and end pos
chr_lengths <- read.table("workspace/projects/pof/data/txt_tsv/chr_range_hg19_true.txt", col.names = c("chr","range"))
min(chr_lengths[chr_lengths$range > 0,]$range)  # min chr will be used to divide rest of ranges
chr_lengths$range_mod <- round(chr_lengths$range / min(chr_lengths[chr_lengths$range > 0,]$range), digits = 2)
chr_lengths$y <- 0  # for y axis

for (i in 1:length(chr_lengths$range)){ 
  if (chr_lengths$range[i] != "0") { 
    chr_lengths$y[i] <- 1    
  }
}
# transform chr to factors and input correct levels order
chr_lengths$chr <- factor(chr_lengths$chr, levels = unique(as.character(chr_lengths$chr)))
chr_lengths$gene <- ""
# now, obtain pos of variants and divide also by minimum rounding to two digits
pos_variants <- bioinfo_possible_variants[,c("pos","gene","chrom")]
pos_variants$pos_mod <- round(pos_variants$pos / min(chr_lengths[chr_lengths$range > 0,]$range) , digits = 2)
pos_variants$y <- 1 

names(pos_variants) <- c("pos","gene","chr","range_mod","y")
circos_df <- rbind(chr_lengths[,c("chr","range_mod","y","gene")], pos_variants[,c("chr","range_mod","y","gene")])

## start circos
# df to plot : chrlengths

circos.clear() # delete a circos previously created 

circos.par("track.height" = 0.1, cell.padding = c(0.02,0,0.02,0), "canvas.ylim" = c(-1.25, 1.1) , # .par defines the way tracks are going to be plotted 
           start.degree = 92) # start.degree defines % of radious of the circle it will plot # cell padding defines white spaces between sectors

circos.initialize(factors = as.factor(circos_df$chr), x = circos_df$range_mod) # define sectors to plot and ranges of X. creates circos, but doesn't plot it

# circos.track(factors = as.factor(chr_lengths$chr), y = chr_lengths$y, bg.col = colors_25 , bg.border = "black", # passes information to y plot for each corresponding track
circos.track(factors = as.factor(circos_df$chr), y = circos_df$y, bg.col = colors_25 , bg.border = "black", # passes information to y plot for each corresponding track
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), 
                           CELL_META$sector.index, facing = "clockwise", niceFacing = T, cex = 0.6)
             }
             )
# note track height defines y radios of track

circos.track(factors = as.factor(circos_df$chr), y = circos_df$y, bg.col = "white", bg.border = "white") # second 
circos.track(factors = as.factor(circos_df$chr), y = circos_df$y, bg.col = "white", bg.border = "white") # third
# circos.track(factors = as.factor(chr_lengths$chr), y = chr_lengths$y, bg.col = "white", bg.border = "black",  # 4 layer
circos.track(factors = as.factor(circos_df$chr),x = circos_df$range_mod, y = circos_df$y, bg.col = "white", bg.border = "black",  # 4 layer
             track.height = 0.05 )

for (i in 1:length(circos_df$gene)){
  circos.text(sector.index = 4, track.index = as.factor(circos_df$chr[i]), x = circos_df$chr[i], 
              y = circos_df$y[i], labels = circos_df$gene[i], facing = "clockwise", niceFacing = T, cex = 0.6)
}
circos.text(sector.index = 4, track.index = as.factor(circos_df$chr[i]), x = circos_df$chr[i], 
            y = circos_df$y[i], labels = circos_df$gene[i], facing = "clockwise", niceFacing = T, cex = 0.6)



# select bioinfo_possible -------------------------------------------------

aux <- bioinfo_possible_variants[,-c("counts","numGenesAff","chrom","pos","ref","alt")]

write.table(aux, file = "variants_66.csv", 
            quote = F, sep = ",", append = F, row.names = F,
            col.names = T)
