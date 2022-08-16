###########Descripción del script#######################
#6 de septiembre

#SCRIPT para leer resultados de los scripts de variantes de python/bedtools/snpSIFT...

#clean global enviroment if needed
rm(list=ls())
########## USO DE LIBRERIAS ####################################

library(dplyr)
library(stringr)
library(plyr)
library(tidyverse)
library(data.table)
########## DIRECTORIO DE TRABAJO ####################################


# Set working directory if needed
getwd()


# guardar enviroment

save.image("resultados_variantes_hasta_el_momento_10_septiembre_19.Rdata",compress = T)

#usar load para cargar el enviroment

load("")

#if needed
rm(df_get_fields,snpsiftdf) #ocupan mucho en memoria
rm(archivos,cromosomas,bed,vector.base,recuento.genes,listado_cromosomas,listado_to1000,anti)

# FUNCIONES ---------------------------------------------------------------
calculate_coverageOfGenes_byExons<-function(DirbedtoolsResults) {
  
  
  #le pasas un archivo con los resultados de bedtools coverage
  #modificar para pasar un directorio con resultados de bedtools coverage
  library(data.table)
  library(dbplyr)
  
  coverage_for_all<-df
  
  i=0
  file_list <- list.files(path=DirbedtoolsResults)
  
  for (file in file_list){
    exones_panel<-read.table(file, stringsAsFactors = F)
    
    colnames(exones_panel)<-c("chrom","start","end","referencia","transcrito","coverage","longitud","bases_cubiertas","porcentaje_cubierto")
    
    exones_panel$gen<-apply(exones_panel,1,function(x){
      filas<-x[["referencia"]]
      gen.in<-unlist(strsplit(filas,split = "_")[1])[3]
      
      
    })
    
    exones_panel$exon_num<-apply(exones_panel,1,function(x){
      filas<-x[["referencia"]]
      num.in<-unlist(strsplit(filas,split = "_")[1])[4]
    })
    
    while (i == 1) {
    gene_coverage_clinic<-setDT(exones_panel)[,as.list(colMeans(.SD[,6])), exones_panel$gen]
    gene_coverage_clinic_rounded<- gene_coverage_clinic %>% select(coverage) %>% round()
    gene_coverage_clinic$coverage<-gene_coverage_clinic_rounded$coverage
    
    panel<-gene_coverage_clinic$exones_panel
    coverage<-gene_coverage_clinic$coverage
    }
    i=i-1
    
    
    coverage_for_all<-cbind(coverage_for_all,coverage)
    
    
    
  }
  return(coverage_for_all)
  # coverage_for_all=gene_coverage_clinic[FALSE,]
  # coverage_for_all<-cbind(coverage_for_all,gene_coverage_clinic$coverage)
  return(coverage_for_all)
} 


# comprobar genes en refseq Y coverage exones poi -----------------------------------------------

#esto es porque sacamos los exones de refseq y queremos saber si hay genes que no encuentra la base de datos
refseq<-read.table("pof/0.0_data/refseq.poi.txt", header = T)

setdiff(genespoi, refseq$name2)

#coverage

exones<-read.table("pof/1.0_results/resultados.exones.coverage.txt")

# exoma clinico de pamplona

exones_panel<-read.table("workspace/proyectos/pof/results/01_coverage_secuenciacion/coverage_secuenciacion_exones.txt", stringsAsFactors = F) 


colnames(exones_panel)<-c("chrom","start","end","referencia","transcrito","coverage","longitud","bases_cubiertas","porcentaje_cubierto")

#porcentajes cubiertos
exones_panelfop1_under_100<-exones_panel %>% filter(exones_panel$porcentaje_cubierto < 1.0)
exones_panelfop1_under_75<-exones_panel %>% filter(exones_panel$porcentaje_cubierto < 0.75)
exones_panelfop1_under_50<-exones_panel %>% filter(exones_panel$porcentaje_cubierto < 0.5)
exones_panelfop1_0<-exones_panel %>% filter(exones_panel$porcentaje_cubierto == 0.0)

rm(exones_panelfop1_0,exones_panelfop1_under_100,exones_panelfop1_under_50,exones_panelfop1_under_75)

# comprimir por gen la media del coverage---------------------------------------------------

#obtener una columna que sea solo el gen y otra el num de exon

exones_panel$gen<-apply(exones_panel,1,function(x){
  filas<-x[["referencia"]]
  gen.in<-unlist(strsplit(filas,split = "_")[1])[3]


})

exones_panel$exon_num<-apply(exones_panel,1,function(x){
  filas<-x[["referencia"]]
  num.in<-unlist(strsplit(filas,split = "_")[1])[4]
})

library(data.table)

#comprimir informacion usando colmeans

gene_coverage_clinic<-setDT(exones_panel)[,as.list(colMeans(.SD[,6])), exones_panel$gen]
gene_coverage_clinic_b<- gene_coverage_clinic %>% select(coverage) %>% round()

gene_coverage_clinic$coverage<-gene_coverage_clinic_b$coverage
rm(gene_coverage_clinic_b)


# leer resultados script python conteo ------------------------------------

df_conteo<-read.table("pof/1.0_results/conteo_variants_con_AF.txt", header = T, stringsAsFactors = F, sep=";")
df_conteo<-df_conteo[-c(161210,161211,161212,161213,161214,161215),] #filas con chrY

#arreglos a la columna de las pacientes para cada variante
df_conteo$pacientes_c<-apply(df_conteo,1,function(row){
  pacientes<-row[["pacientes_c"]]
  fix<-gsub(pacientes,pattern = "'",replacement = "",perl=F)
  fix<-gsub(fix,pattern = "[][]",replacement = "", perl=F)
  return(fix)
})

#arreglos a la columna de las pacientes para cada variante
df_conteo$pacientes_foo<-apply(df_conteo,1,function(row){
  pacientes<-row[["pacientes_foo"]]
  fix<-gsub(pacientes,pattern = "'",replacement = "",perl=F)
  fix<-gsub(fix,pattern = "[][]",replacement = "", perl=F)
  return(fix)
})


#arreglos a la columna de las pacientes para cada variante
df_conteo$pacientes_fop<-apply(df_conteo,1,function(row){
  pacientes<-row[["pacientes_fop"]]
  fix<-gsub(pacientes,pattern = "'",replacement = "",perl=F)
  fix<-gsub(fix,pattern = "[][]",replacement = "", perl=F)
  return(fix)
})



# leer resultados script python get_fields --------------------------------

df_coverage<-read.table(file = "pof/1.0_results/vcf_af1000g_coverage_FOP1.txt",header = T,sep = "\t",stringsAsFactors = F)

#cuantas variantes tenemos con 0 coverage para FOP1?

sum(df_coverage$variant_coverage==0)
#209 con coverage 0


# combinar dframe conteo con coverage -------------------------------------

anti_join(df_conteo,df_coverage) #hay variantes que no se me unen?
df_get_fields<-inner_join(df_conteo,df_coverage, by = c("chrom","pos"))

# extraer variantes que afecten a genes que tengamos en list poi --------



#primero conseguimos nuestra lista de genes poi en nomenclatura HGNC 
genes_list_poi<-read.table("pof/0.0_data/poi610_sorted.bed", header = F, sep = "\t", stringsAsFactors = F)
genespoi<-genes_list_poi$V4

genes_list_poi<-head(genes_list_poi, n=-2)

#primero quitamos los corchetes a la columna de los genes

df_get_fields$gene<-apply(df_get_fields,1,function(row){
  all.genes.in<-row[["genes"]]
  all.genes.in<-gsub(pattern = '[][]', perl = F, replacement = "",x = all.genes.in) #Por alguna razon para eliminar los [] hay que hacer esto
  return(all.genes.in)
})

#ahora si, vamos a crear una columna donde esten los genes poi para cada variante
df_get_fields$poi<-apply(df_get_fields,1,function(row){
  all.genes.in<-row[["genes"]]
  gen.in<-strsplit(all.genes.in, ",\\s", perl=T)[[1]]
  
  genes<- all.genes.in
  
  gen<-gen.in[gen.in %in% genes_list_poi$V4]
  
  if (length(gen)==0){
    gen<-""
  }
  return(gen)
  # for (gen in gen.in) {
  #   poi.value<-genes_list_poi[gen %in% genes_list_poi$V4,]$V4
  #   # genes<-gsub(x=genes, pattern = gen, replacement = poi.value)
  #   genes<-paste(poi.value, collapse = "+")
  # }
  # return(genes)
})

#convertir la columna poi apropiadamente en vectores de caracteres
df_get_fields$poi<- vapply(df_get_fields$poi, paste, collapse=", ", character(1L))

# recuento genes poi en variantes -----------------------------------------


#tenemos que leer por fila ATTS_genes row["genes"] y contar en una variable los genes que encuentra de hsa, pero solo contar una vez
#y además, se vuelve a dar el caso que los genes pueden estar en groups/OR. Igualmente, si lo encuentra, cuentalo una vez.

recuento.genes <- unlist(apply(df_get_fields, 1, function(row){ #le pasas una funcion al apply que va actuar sobre un elemento (i)
  genes <- row[["gene"]] #el elemento es una matriz con los indices de las filas
  recuento.poi <- unlist(strsplit(genes, ",\\s", perl = TRUE)[[1]]) #el [[1]] para el unlist es porque el strsplit te devuelve listas segun el número de objetos que le hayas pasado entre comillas (ejemplo strsplit(c("a b c", "1 2 3"), "\\s", perl = TRUE)) )
  #print(head(recuento.hsa, n=10))
  #print(class(recuento.hsa))
})) #aunque le hagas un unlist a tu variable de la fila, como al final haces un apply, este te devuelve una lista. de ahi otro unlist a todo 
recuento.genes <- unique(recuento.genes) #unique al vector para quitar duplicados
poi.in.variants<-intersect(recuento.genes, genes_list_poi$V4) #intersect te quitaria los duplicados, es del paquete dplyr, pero igualmente nos lo hemos quitado antes con el unique
poi.not.in<- setdiff(genes_list_poi$V4, recuento.genes)

#583 en variantes
#28 no hay variantes



# combinar dataframes de python con snpsift --------------------------------

snpsiftdf<-read.table("pof/1.0_results/snpsift_info_variants.txt", header = T, stringsAsFactors = F, sep = "\t" )
colnames(snpsiftdf)<-c("chrom","pos","id","ref","alt","allele","aa","N","biotype","effect","impact", "gene")


snpsiftdf$gene<-NULL
#quitar columnas repetidas
df_get_fields[,c("gene","ref","alt")]<-list(NULL)

variants.all.info<-inner_join(snpsiftdf,df_get_fields,by=c("chrom","pos"))

write.table(variants.all.info,file = "all.variants.+info.tsv",sep = "\t",quote = F,row.names = F,col.names = T)
#me añade mas que no deberia añadir>recomprobar el join, que esta ocurriendo?
anyDuplicated(variants.all.info$id)
variants.all.info<-unique(variants.all.info)

# filtrar por MAF 0,05  --------------------------------------------

variants.all.info.maf <- variants.all.info %>% filter(variants.all.info$maf<0.05)
variants.all.info.maf2 <- variants.all.info %>% filter(variants.all.info$maf<0.0033)

#filtrar por aquellas variantes que tengan genes pof

variants.all.poi<- variants.all.info %>% filter(variants.all.info$poi!='')
##5929 variantes que afecten a genes poi de las 160.000 aprox que tenemos

variants.all.poi.maf<- variants.all.poi %>% filter(variants.all.poi$maf<=0.05)

#ahora, nos centramos en las protein-coding que causen varios tipos de mutaciones de alto riesgo:

variants.all.poi.maf_and_missense<- filter(variants.all.poi.maf, grepl("missense_variant",effect))
variants.all.poi.maf_and_frameshift<- filter(variants.all.poi.maf, grepl("frameshift_variant",effect))

variants.all.poi.maf_and_stop<- filter(variants.all.poi.maf, grepl("stop_gained",effect))
variants.all.poi.maf_and_stop_lost<- filter(variants.all.poi.maf, grepl("stop_lost",effect))

variants.all.poi.maf_splice_acceptor_variant<- filter(variants.all.poi.maf, grepl("splice_acceptor_variant",effect))
variants.all.poi.maf_splice_donor_variant <- filter(variants.all.poi.maf, grepl("splice_donor_variant",effect))




#  filtrar pof MAF 0,0033 ------------------------------------------------------------------------

#otro approach: 1/300
variants.all.poi.maf2<- variants.all.poi %>% filter(variants.all.poi$maf<=0.0033)


#ahora, nos centramos en las protein-coding que causen varios tipos de mutaciones de alto riesgo:
variants.all.poi.maf2_and_missense<- filter(variants.all.poi.maf2, grepl("missense_variant",effect))
variants.all.poi.maf2_and_frameshift<- filter(variants.all.poi.maf2, grepl("frameshift_variant",effect))

variants.all.poi.maf2_and_stop<- filter(variants.all.poi.maf2, grepl("stop_gained",effect))
variants.all.poi.maf2_and_stop_lost<- filter(variants.all.poi.maf2, grepl("stop_lost",effect))

variants.all.poi.maf2_splice_acceptor_variant<- filter(variants.all.poi.maf2, grepl("splice_acceptor_variant",effect))
variants.all.poi.maf2_splice_donor_variant <- filter(variants.all.poi.maf2, grepl("splice_donor_variant",effect))


# combinar los df de las variantes de este tipo y leer en gdrive ---------


#usar Reduce para aplicar a la lista de nuestros dataframes un merge por filas (bind_rows)

#para MAF 0,05 

union_maf<- Reduce(bind_rows,list(variants.all.poi.maf_and_frameshift,variants.all.poi.maf_and_missense,variants.all.poi.maf_splice_donor_variant,variants.all.poi.maf_splice_acceptor_variant,variants.all.poi.maf_and_stop_lost,variants.all.poi.maf_and_stop))

# union_maf[,c("gen.y","ref.y","alt.y")]<-list(NULL)

rm(variants.all.poi.maf_and_frameshift,variants.all.poi.maf_and_missense,variants.all.poi.maf_splice_donor_variant,variants.all.poi.maf_splice_acceptor_variant,variants.all.poi.maf_and_stop_lost,variants.all.poi.maf_and_stop)

#para MAF 0,0033
union_maf_own<-Reduce(bind_rows,list(variants.all.poi.maf2_and_frameshift,variants.all.poi.maf2_and_missense,variants.all.poi.maf2_and_stop,variants.all.poi.maf2_and_stop_lost
                                     ,variants.all.poi.maf2_splice_donor_variant,variants.all.poi.maf2_splice_acceptor_variant))

# union_maf_own[,c("gen.y","ref.y","alt.y")]<-list(NULL)

rm(variants.all.poi.maf2_and_frameshift,variants.all.poi.maf2_and_missense,variants.all.poi.maf2_and_stop,variants.all.poi.maf2_and_stop_lost
   ,variants.all.poi.maf2_splice_donor_variant,variants.all.poi.maf2_splice_acceptor_variant)

#escribir resultados

# union_maf$poi<- vapply(union_maf$poi, paste, collapse=", ", character(1L))
# union_maf_own$poi<- vapply(union_maf_own$poi, paste, collapse=", ", character(1L))

write.table(union_maf, "efecto_bio_todas.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(union_maf_own,"efectio_bio_poi.tsv", sep = "\t", row.names = F, col.names = T, quote = F )

# conteo de presencia por grupos ------------------------------------------------------------------

variants.all.poi.conteo_pacientes <- variants.all.poi %>% filter(variants.all.poi$controles<10 & variants.all.poi$pacientes>15)
variants.all.poi.conteo_foo <- variants.all.poi %>% filter(variants.all.poi$controles<10 & variants.all.poi$foo>25)
variants.all.poi.conteo_fop <- variants.all.poi %>% filter(variants.all.poi$controles<10 & variants.all.poi$fop>25)



variants.all.info.conteo_pacientes <- variants.all.info %>% filter(variants.all.info$controles<5 & variants.all.info$pacientes>25)
variants.all.info.conteo_foo <- variants.all.info %>% filter(variants.all.info$controles<10 & variants.all.info$foo>25)
variants.all.info.conteo_fop <- variants.all.info %>% filter(variants.all.info$controles<10 & variants.all.info$foo>25)

#escribrir resultados

#hay que transformar la columna poi a vectores de caracteres
# variants.all.poi.conteo_pacientes$poi<- vapply(variants.all.poi.conteo_pacientes$poi, paste, collapse=", ", character(1L))
# variants.all.info.conteo_pacientes$poi<- vapply(variants.all.info.conteo_pacientes$poi, paste, collapse=", ", character(1L))

write.table(variants.all.poi.conteo_pacientes,"conteo_poi_c10_pacientes25.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(variants.all.info.conteo_pacientes,"conteo_todas_c10_pacientesp25.tsv", sep = "\t", row.names = F, col.names = T, quote = F )

rm(variants.all.poi.conteo_pacientes,variants.all.poi.conteo_foo,variants.all.poi.conteo_fop)
rm(variants.all.info.conteo_foo,variants.all.info.conteo_fop,variants.all.info.conteo_pacientes)


rm(variants.all.info.maf2,variants.all.info.maf, variants.all.poi.maf,variants.all.poi.maf2)



