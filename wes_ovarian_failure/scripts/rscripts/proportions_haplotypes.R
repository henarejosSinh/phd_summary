###########Descripción del script#######################
#24 october 2019

# Study proportions of homozygous and heterozygous in our variants_poi_data

# Ismael Henarejos Castillo
# ihc.europa@gmail.com


#clean global enviroment if needed
rm(list=ls())

########## USO DE LIBRERIAS ####################################


library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)
# source("/home/ihenarejos/workspace/proyectos/pof/scripts/functions.R")


########### save_load enviroment ----------------------------------------------------

#save.image("workspace/proyectos/pof/scripts/rdatas/proportions_analysis.Rda")
# load("workspace/proyectos/pof/scripts/rdatas/")

# functions ----------------------------------------------------------------
# read proportions --------------------------------------------------------

proportions_all<-fread("workspace/proyectos/pof/results/10_proportions_tables/proportions_groups_5%cases_missings.txt")

aux<-colMeans(proportions_all[,2:ncol(proportions_all)])
aux<-data.frame(aux)

aux$type<-rownames(aux)
aux$class<-rep(x = "variants_all",3)

colnames(aux)[1]<-"value"

aux$group<-apply(aux,1,function(row){
  rowvalue<-row[["type"]]
  tores<-unlist(strsplit(rowvalue,split = "_")[1])[[1]]
  return(tores)
})

aux$group2<-apply(aux,1,function(row){
  rowvalue<-row[["type"]]
  tores<-unlist(strsplit(rowvalue,split = "_")[1])[[2]]
  return(tores)
})

(plotv<-ggplot(data=aux, aes(x=group, y=value, fill=group2, group=group)) +
    geom_bar(stat="identity", width = 0.2) + ggtitle("Haplotype Proportions for variants filtered and change of missings")+ 
    geom_text(aes(label =paste0(round(value),'%')), position = position_stack(0.9), vjust =1, hjust= 0.4, size= 5) + 
    ylab("proportions") +
    theme_light() + theme( axis.title.x = element_blank(),
                           
                           axis.text = element_text(size=15),
                           axis.title = element_text(size=17),
                           plot.title = element_text(size=19, hjust=0.5),
                           legend.title = element_blank()))


scompare<-full_join(df,aux)

compare<-melt(compare)

compare$type<-gsub(compare$type, pattern = "het_pro", replacement = "heterozygous")
compare$type<-gsub(compare$type, pattern = "hom_pro", replacement = "homozygous")
compare$type<-gsub(compare$type, pattern = "ref_pro", replacement = "0/0")

(plotCompare<-ggplot(compare,aes(x=variable, y=(100*value), fill=type))+
    geom_bar(stat = "identity", width = 0.2) +
    geom_text(aes(label =paste0(round(100*value),'%')), position = position_stack(0.9), vjust =1, hjust= 0.4, size= 5) + 
    facet_wrap(.~class) +
    xlab("")  +  ggtitle("Haplotype Proportions for presence of mutations") + #ylim(0,80)
    ylab("%") +
    theme_light()+
    theme(legend.position="bottom",
          axis.text = element_text(size=15),
          axis.text.x = element_blank(),
          axis.title = element_text(size=17),
          plot.title = element_text(size=19, hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=15),
          strip.text.x = element_text(size=12, color="red",
                                      face="bold.italic"),
          strip.text.y = element_text(size=12, color="red",
                                      face="bold.italic"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=0.5, linetype="solid")))


# cases vs control --------------------------------------------------------


#compare cases vs control proportions
pro_table<- fread("workspace/proyectos/pof/results/10_proportions_tables/proportions_table.txt",data.table = F)
colnames(pro_table)[which(names(pro_table) == "UNKNOWN")]<-"CONTROL"

prop.ctr<-pro_table[colnames(pro_table) == "CONTROL"]
prop.case<-pro_table[colnames(pro_table) != "CONTROL"]

prop.aux.ctr <- table(unlist(prop.ctr[, -1]))
prop.aux.ctr <- prop.aux.ctr / (sum(prop.aux.ctr))

prop.aux.case <-table(unlist(prop.case[, -1]))
prop.aux.case <- prop.aux.case / (sum(prop.aux.case))


compareCC<-full_join(data.frame(prop.aux.case,stringsAsFactors = F),data.frame(prop.aux.ctr,stringsAsFactors = F))
compareCC_adj<-full_join(data.frame(prop.aux.case,stringsAsFactors = F),data.frame(prop.aux.ctr,stringsAsFactors = F))

classList<-rep(x = "variants_all_case",3)
classList<-append(classList,after = 3,values = rep(x="variants_all_control",3))
compareCC$class<-classList                 

classList<-rep(x = "variants_adj_case",3)
classList<-append(classList,after = 3,values = rep(x="variants_adj_control",3))
compareCC_adj$class<-classList                 


comparecc_true<-melt(full_join(compareCC,compareCC_adj))

comparecc_true$Var1<-gsub(comparecc_true$Var1, pattern = "het", replacement = "heterozygous")
comparecc_true$Var1<-gsub(comparecc_true$Var1, pattern = "hom", replacement = "homozygous")
comparecc_true$Var1<-gsub(comparecc_true$Var1, pattern = "ref", replacement = "0/0")

expList<-rep(x = "all",6)
expList<-append(expList,after = 6,values = rep(x="padjusted",6))

comparecc_true$exp<-expList

expList<-rep(x = "case",3)
expList<-append(expList,after = 9,values = rep(x="control",3))

comparecc_true$type<-expList

#take out commas from Var1
comparecc_true$Var1<-apply(comparecc_true,1,function(row){
  value<-row[1]
  fix<-gsub(pattern = "\'",replacement = "",x = value)
  return(fix)
})

comparecc_true$exp<-apply(comparecc_true,1,function(row){
  value<-row[["exp"]]
  fix<-gsub(pattern = "all",replacement = "allVariants",x = value)
  return(fix)
})

(plotCompare<-ggplot(comparecc_true,aes(x=type, y=(100*value), fill=as.factor(Var1)))+
    geom_bar(stat = "identity", width = 0.2) + # position = "dodge"
    geom_text(aes(label =paste0(round(100*value),'%')), position = position_stack(0.9), vjust =-0.5, hjust= 0.4, size= 5) + 
    facet_wrap(.~exp) +
    xlab("")  +  ggtitle("Haplotype Proportions for presence of mutations") + #ylim(0,80)
    ylab("%") +
    theme_light()+
    theme(legend.position="bottom",
          axis.text = element_text(size=15),
          #axis.text.x = element_blank(),
          axis.title = element_text(size=17),
          plot.title = element_text(size=19, hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=15),
          strip.text.x = element_text(size=12, color="red",
                                      face="bold.italic"),
          strip.text.y = element_text(size=12, color="red",
                                      face="bold.italic"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=0.5, linetype="solid")))

save.image("workspace/proyectos/pof/scripts/rdatas/proportions_data.Rda")
