########### Description #######################################################
# 4 oct 2019

# Functions for POI variants.

#working directory
# dir=('/home/ihenarejos/workspace/')
# setwd(dir)
# setwd("/home/ihenarejos/workspace/")

#clean global enviroment if needed
rm(list=ls())

########## Libraries ##########################################################

library(lintr)  # For code writing guidelines
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(RVAideMemoire)  # for multiple fisher tests comparisons
library(statmod)  # for power fisher test
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


########## Functions ####################################

# use of server ----------------------------------------------------------

qsub.array.job(job.name = "BedtoolsCov",
              cmds = coverage.cmds,
              modules = c("BEDTools"),
              logs = logs.dir,
              num.threads = 30,
              max.simultaneous.jobs = 30,
              mem = "2G")

# colours ------------------------------------------------------------------
colours_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "Gray", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)


# themes for ggplot --------------------------------------------------------

sinh_theme <- (theme_light( base_size = 18 ) +
  theme(
    axis.text.x = element_text(hjust = 0.65),
    axis.title.x = element_text(vjust = -1.0),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", 
                                    colour = NA, inherit.blank = T),
    plot.background = element_rect( colour = "white" , inherit.blank = T),
    panel.grid = element_line( colour = "gray88" , 
                               arrow = F, inherit.blank = T),
    panel.grid.minor = element_line( size = rel(0.25),
                                     arrow = F, inherit.blank = T),
    panel.border = element_rect( colour = "gray70" , size = rel(1), 
                                 inherit.blank = T),
    axis.line = element_blank() )) 

bar.zyx.theme <-  theme_light( base_size = 18 ) +
  theme(axis.text.x = element_text(hjust = 0.65),
        axis.title.x = element_text(vjust = -1.0),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white", 
                                        colour = NA, inherit.blank = T),
        plot.background = element_rect( colour = "white" , inherit.blank = T),
        panel.grid = element_line( colour = "gray88" , 
                                   arrow = F, inherit.blank = T),
        panel.grid.minor = element_line( size = rel(0.25),
                                         arrow = F, inherit.blank = T),
        panel.border = element_rect( colour = "gray70" , size = rel(1), 
                                     inherit.blank = T),
        axis.line = element_blank() ) 

# calculate coverage, num read optimos, bedtools section -------------------

calculate_coverage<-function(DirbedtoolsResults){
  
  library(data.table)
  library(dplyr)
  
  file_list <- list.files(path=DirbedtoolsResults, 
                          full.names = TRUE)
  
  list_df<-lapply(file_list,function(archivo){
    
    sample_n<-paste0(unlist(strsplit(basename(archivo),split = "_")[1])[4],"_",unlist(strsplit(basename(archivo),split = "_")[1])[5])
    
    
    df<-read.table(archivo, stringsAsFactors = F)
    colnames(df)<-c("chrom","start","end","referencia","transcrito","coverage","longitud","bases_cubiertas","porcentaje_cubierto")
    
    
    df$gen<-apply(df,1,function(x){
      filas<-x[["referencia"]]
      gen.in<-unlist(strsplit(filas,split = "_")[1])[3]
      
    })
    
    df$exon_num<-apply(df,1,function(x){
      filas<-x[["referencia"]]
      num.in<-unlist(strsplit(filas,split = "_")[1])[4]
    })
    
    
    gene_coverage_clinic<-setDT(df)[,as.list(colMeans(.SD[,6])), df$gen]
    gene_coverage_clinic_rounded<- gene_coverage_clinic %>% select(coverage) %>% round()
    gene_coverage_clinic$coverage<-gene_coverage_clinic_rounded
    
    colnames(gene_coverage_clinic)<-c("gene",sample_n)
    
    rm(gene_coverage_clinic_rounded)
    return(gene_coverage_clinic)
  })
  
  df_total <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), list_df)
  
  
  return(df_total)
}

calculate_snum_reads_optimal<-function(DirbedtoolsResults){
  
  library(data.table)
  library(dplyr)
  
  file_list <- list.files(path=DirbedtoolsResults, 
                          full.names = TRUE)
  
  list_df<-lapply(file_list,function(archivo){
    
    sample_n<-paste0(unlist(strsplit(basename(archivo),split = "_")[1])[4],"_",unlist(strsplit(basename(archivo),split = "_")[1])[5])
    
    
    df<-read.table(archivo, stringsAsFactors = F)
    colnames(df)<-c("chrom","start","end","referencia","transcrito","coverage","longitud","bases_cubiertas","porcentaje_cubierto")
    

    df$num_optimo<-apply(df,1,function(x){
      tamaño_exon<-x[["longitud"]]
      if (as.integer(tamaño_exon) <=150) {num_opt<-150}
      if (as.integer(tamaño_exon) >=150) {num_opt<-((as.integer(tamaño_exon)/150) *100)}
      return(num_opt)
    })
    
    gene_coverage_clinic<- df %>% select("referencia","num_optimo","coverage")
    
    gene_coverage_clinic$check<-apply(gene_coverage_clinic,1,function(x){
      num<-x[["num_optimo"]]
      nreads<-x[["coverage"]]
      if (as.integer(nreads) < as.integer(num)) {
        result<-"under"}
      else {
        result<-"optimal"
      }
    })
      
    gene_coverage_optimal<- gene_coverage_clinic %>% select("referencia","check")  
    colnames(gene_coverage_optimal)<-c("gen_exon",sample_n)
    return(gene_coverage_optimal)
  })
  
  df_total <- Reduce(function(...) merge(..., all = TRUE, by = "gen_exon"), list_df)
  
  
  return(df_total)
}

calculate_percentage_bases_covered<-function(DirbedtoolsResults){
  
  library(data.table)
  library(dplyr)
  
  file_list <- list.files(path=DirbedtoolsResults, 
                          full.names = TRUE)
  
  list_df<-lapply(file_list,function(archivo){
    
    sample_n<-paste0(unlist(strsplit(basename(archivo),split = "_")[1])[4],"_",unlist(strsplit(basename(archivo),split = "_")[1])[5])
    
    
    df<-read.table(archivo, stringsAsFactors = F)
    colnames(df)<-c("chrom","start","end","referencia","transcrito","coverage","longitud","bases_cubiertas","porcentaje_cubierto")
    
    
    df$gen<-apply(df,1,function(x){
      filas<-x[["referencia"]]
      gen.in<-unlist(strsplit(filas,split = "_")[1])[3]
      
    })
    
    df$exon_num<-apply(df,1,function(x){
      filas<-x[["referencia"]]
      num.in<-unlist(strsplit(filas,split = "_")[1])[4]
    })
    
    df$porcentaje_cubierto<-df$porcentaje_cubierto*100
    # df$porcentaje_cubierto<-apply(df,1,function(x){
    #   value<-x[["porcentaje_cubierto"]]
    #   value_mod<-value[1]*100
    # })
    
    gene_coverage_clinic<-setDT(df)[,as.list(colMeans(.SD[,9])), df$gen]
    
    
    colnames(gene_coverage_clinic)<-c("gene",sample_n)
    
    rm(gene_coverage_clinic_rounded)
    return(gene_coverage_clinic)
  })
  
  df_total <- Reduce(function(...) merge(..., all = TRUE, by = "gene"), list_df)
  
  
  return(df_total)
}

generate_bedtools_coverage_cmd <- function(dirBams, output_dir){
  
  
  list_bam<-list.files(path = dirBams,
                       pattern = ".bam$",
                       full.names = T)
  
  list_cmd<-lapply(list_bam,function(bam){
    sample<-unlist(strsplit(basename(bam), split = "\\.", perl = T)[1])[1]
    
    cmd <- paste0("bedtools coverage -b ", bam,
                  " -a /data/network/nas01/fivibio-data/projects/pof/data/bed/MEXv01_regiones_completas_21022018_target_10-extrabases.bed -sorted > ", 
                  output_dir, "/coverage_secuenciacion_exones_",sample,".tsv")  
    
    return(cmd)
    
  })
  return(unlist(list_cmd))
  
}

generate_bedtools_coverage_cmd_variants<-function(dirBams, output_dir){
  
  
  list_bam<-list.files(path = dirBams,
                       pattern = ".bam$",
                       full.names = T)
  
  list_cmd<-lapply(list_bam,function(bam){
    sample<-unlist(strsplit(basename(bam), split = "\\.", perl = T)[1])[1]
    
    cmd <- paste0("bedtools coverage -b ", bam,
                  " -a /data/network/nas01/fivibio-data/projects/pof/data/vcf/all_samples_gatk_annot_fix_missing_sort.vcf -sorted -counts > ", 
                  output_dir, "/vcf_plus_coverage_",sample,".vcf")  
    
    return(cmd)
  
})}

extract_coverage_variants_bedtools_from_vcf_cmd<-function(dirVcf, output_dir){
  
  
  list_vcf<-list.files(path = dirVcf,
                       pattern = ".vcf$",
                       full.names = T)
  
  list_cmd<-lapply(list_vcf,function(vcf){
    sample<-unlist(strsplit(basename(vcf), split = "\\.", perl = T)[1])[1]
    num<-unlist(strsplit(sample,split="_") [1])[4]
      
    cmd <- paste0("grep -v '^#' ", vcf,
                  " | -rev | cut -f 1 |  rev | ", 
                  output_dir, "/coverage_variants_",sample,".tsv")  
    
    return(cmd)
    
  })}


# bowtie/samtools section -------------------------------------------------


#doing alignment wit the hg38.96 from server
generate_bowtie2_cmd<-function(fastq_dir,n_cores,index_file) {
  
  
  
  fastq_list<-list.files(path = fastq_dir,full.names = T, pattern = "\\.fastq"  )
  
  cores_to_use=n_cores
  
  order<-lapply(fastq_list, function(fastq_list){
    fq<-unlist(strsplit(basename(fastq_list), split = "_")[[1]])[[1]]
    cmd<-paste0("bowtie2 -p ",cores_to_use," -x ",index_file," -1 ",fq,"_1.fastq -2 ",fq,"_2.fastq -S ",fq,".sam")
    return(cmd)
  })
  
  comandos<-unique(order)
  
  return(comandos)
  
}

#check if works
#para generar bam files sortered from a sam
generate_samtools_view_cmd<- function(sam_dir,n_cores){
  
  sam_list<-list.files(path = fastq_dir,full.names = T, pattern = "\\.sam"  )
  
  cores_to_use=n_cores
  
  order<-lapply(sam_list, function(sam){
    tobam<-unlist(strsplit(basename(sam), split = "\\.")[[1]])[[2]]
    cmd<-paste0("samtools view -bS -@ ",cores_to_use," ",sam," | samtools sort -@ ",cores_to_use," > ",tobam,".bam")
    return(cmd)
  })
  
  comandos<-unique(order)
  
  return(comandos)
  
  
}

#alingment with bowtie2 + create sortered bam
generate_bowtie2_sort_bam_cmd<-function(fastq_dir,n_cores,index_file,outputdir) {
  
  fastq_list<-list.files(path = fastq_dir,full.names = T, pattern = "\\.fastq"  )
  
  cores_to_use=n_cores
  
  order<-lapply(fastq_list, function(fastq_list){
    dir<-dirname(fastq_list)
    fq<-unlist(strsplit(basename(fastq_list), split = "_")[[1]])[[1]]
    cmd<-paste0("bowtie2 -p ",cores_to_use," -x ",index_file," -1 ",dir,"/",fq,"_1.fastq -2 ",dir,"/",fq,"_2.fastq ")
    cmd<-paste0(cmd," | samtools view -bS -@ ",cores_to_use," | samtools sort -@ ",cores_to_use," > ",fq,".bam")
    return(cmd)
  })
  
  comandos<-unique(order)
  
  return(comandos)
  
}



# 1000g commands ---------------------------------------------------------

generate_cmd_annotate_1000g<-function(thousand_genome_dir,vcf_header,vcf_variants, outputDir_sh, outputdir_VCF, 
                                      vcf_annotated, annotation_term){
  
  file_list_1000g<-list.files(thousand_genome_dir, 
                              full.names = T,
                              pattern= "\\.gz$")
  # get num 
  lists_chrom_1<-lapply(file_list_1000g,function(file){
    chr.in<-strsplit(basename(file), split = "\\.", perl = T)[[1]][[2]]
    return(chr.in)
  })
  
  #generate first sh to get variants by chr ###
  sh_get_chr<-lapply(listado_cromosomas_with_chr, function(chr){
    chr.in<-chr
    commands<-paste0("cat ",vcf_header," > ",outputdir_VCF,chr.in,".vcf && grep -v '^#' ",vcf_variants,
                     " | grep -w ",chr.in," >> ",outputdir_VCF,chr.in,".vcf ", collapse = "", sep= "")
    return(commands)
  })
  
  #generate second sh to annotate each individual VCF created in last step and join information in final VCF
  lists_chrom<-lapply(lists_chrom,function(crom){
    cros<-gsub(pattern = 'chr' ,replacement ='',crom )
    return(crom)
  })
  
  df<-data.frame(lists_chrom,stringsAsFactors = F)
  
  listado_to1000<-apply(df,1,function(row){
    
    chr.in<- row["chrom"]
    if (chr.in == "X") {
      comandos<-paste0("java -jar $EBROOTSNPEFF/snpSift.jar annotate -name ",annotation_term," -info AF ",
                       thousand_genome_dir,"ALL.chr",chr.in,
                       ".phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz ",
                       outputdir_VCF,"chr",chr.in,".vcf | grep -v '^#' >> ",vcf_annotated, collapse = "", sep = "")
    } else {
      comandos<-paste0("java -jar $EBROOTSNPEFF/snpSift.jar annotate -name ",annotation_term,
                       " -info AF ",thousand_genome_dir,"ALL.chr",chr.in,
                       ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ",
                       outputdir_VCF,"chr",chr.in,".vcf | grep -v '^#' >> ",vcf_annotated, collapse = "", sep = "")
    }
    return(comandos)
  })
  
  #write results
  write.table(listado_to1000,file = paste0(outputDir_sh,"annotate_vcf_1000g.sh"), 
              row.names = F, col.names = F, quote = F )
  write.table(sh_get_chr,file = paste0(outputDir_sh,"get_chr_variants_in_vcf.sh"), 
              row.names = F, col.names = F, quote = F, sep = "\n")
}



# fisher/chi2 tests -------------------------------------------------------

# Function that uses sample data to add new columns for count table obtained in 
# python.
add_columns_contigency <- function(count_table) {
  
  count_table$control_not_affected<-apply(count_table,1,function(row){
    control.in<-row[["control"]]
    
    control.rest<- 32-as.numeric(control.in)
    return(control.rest)
  })
  
  count_table$foo_not_affected<-apply(count_table,1,function(row){
    foo.in<-row[["foo"]]
    
    foo.rest<- 83-as.numeric(foo.in)
    return(foo.rest)
  })
  
  
  count_table$fop_not_affected<-apply(count_table,1,function(row){
    fop.in<-row[["fop"]]
    
    fop.rest<- 35-as.numeric(fop.in)
    return(fop.rest)
  })
  
  
  count_table$cases_not_affected<-apply(count_table,1,function(row){
    cases<-row[["cases"]]
    
    cases.rest<- 118-as.numeric(cases)
    return(cases.rest)
  })
  return(count_table)
}

# Function that generates a contingency table for each variant and apply fisher 
# test on them
fisher_test_variants<-function(df,vid,ctrmut,ctr,casemut,case){
  
  list_dfs_fisher_tests<-apply(df,1,function(row){
    
    variant_id<-row[[vid]]  # id from variant
    
    control_value<-as.numeric(row[[ctrmut]])  # controls that have the variant
    control_value2<-as.numeric(row[[ctr]])    # controls with reference allele
    casos_value<-as.numeric(row[[casemut]])    # cases that have the variant
    casos_value2<-as.numeric(row[[case]])    # cases with reference allele
    
    #create matrix and apply fisher test
    fisher_list<-fisher.test(matrix(c(control_value,control_value2,casos_value,
                                      casos_value2),byrow = T,
                                      nrow = 2,ncol = 2))  # normal 2x2 table
    
    # fisher test return a list of components of class "htest"  with components
    # such as p.value and odds ratio estimate
    # that we can get from. So for each variant, return a dataframe with 3 
    # columns; id of variant, pvalue and odds ratio
    return(data.frame(variant = variant_id, p.value = fisher_list$p.value, 
                      estimate = fisher_list$estimate))
  })
  
  # once we have dataframes from each variant, we ;
  # 1; do a rbindlist so each dataframe for each variant is joined together in an only dataframe where each row
  #   is a variant with the results from its fisher test
  # 2; merge the dataframe from step 1 with original dataframe that had the basis data from python.
  return(merge(df,rbindlist(list_dfs_fisher_tests), by= "variant"))
  
}  

# Function that calculates statistical power of fisher test for each variant
power_fisher_variants <- function(df,total1,total2){
  
  df_powf<-apply(df,1,function(row){
    
    # Similar to fisher test function, where we will be adding a new column
    # to original dataframe with statistical power of fisher test calculated
    variant_id<-row[["variant"]]
    control_value<-as.numeric(row[["control"]])
    cases_value<-as.numeric(row[["cases"]])
    
    #create matrix and apply fisher test
    pow <- power.fisher.test(control_value/total1,cases_value/total2,total1,
                             total2,nsim = 10000)
    
    return(data.frame(variant = variant_id, fisher.power = pow))
  })
  
  return(merge(df,rbindlist(df_powf), by= "variant"))
  
}  

# Function that checks if variant has multiple alleles. To use for visualization
# of data
check_if_variant_multi <- function(df_variants,id_colum, delim) {
  
  df_variants$multi<-apply(df_variants,1,function(row){
    id<-row[id_colum]
    alt<-unlist(strsplit(id,split = "_")[[1]])[4]
    
    check_value<-grepl(alt,pattern = delim)
    return(check_value)
  })
  return(df_variants)
}

# Function that creates a df using counts of how many variants are SNV and how
# many are multiallelic for barplot visualization
barplot_counts <- function(df, type = ""){
  aux <- table(df$multi)
  aux <- data.frame(SNV = aux["FALSE"], multi = aux["TRUE"],row.names = NULL, 
                    type = type)
  return(aux)
}

# Function that creates a barplot for visualization of fisher test of variants
plot_variants_fisher <- function(df,title,ylim,base_size,colour1,colour2){
  ggplot(data=df, aes(x=factor(variable), group = factor(variable), y = value, 
                      fill = variable)) + 
    geom_bar(stat="identity", colour = "black ", width = 0.5, 
             position = position_dodge(width = 0.6)) + 
    scale_fill_manual(values = c("SNV" = colour1, "multi"
                                 = colour2)) +
    # scale_fill_brewer(palette = "Set1") +
    geom_text(aes(label = value), position = position_dodge(width = 0.6), 
              vjust =-0.4, size = 6 ) + 
    facet_wrap(.~pvar) + ylim ( 0, ylim) +
    xlab("") +  ggtitle( title ) + 
    ylab("Número de variantes significativas") +
    theme_light(base_size = base_size) +   
    theme(legend.position = "none", 
          axis.text = element_text(size=rel(1.6)),
          axis.text.x = element_text(angle=90),
          axis.title.y = element_text(margin = margin(r = 10)),  # for margin
          axis.title = element_text(size=rel(1.7), vjust = -0.5 ),
          plot.title = element_text(size=rel(2.4), hjust=0.5 ),
          legend.title = element_blank(),
          legend.text = element_blank(),
          strip.text.x = element_text(size=rel(1.8), color="red",
                                      face="bold.italic"),
          strip.text.y = element_text(size=rel(1.8), color="red",
                                      face="bold.italic"),
          strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))
}

# Function to perform pairwise fisher tests using the genotypes in variants
fisher_pairwise_2_x_3 <- 
  function(df,vid,control_AA,control_Aa,control_aa,case_AA,case_Aa,case_aa){
    
    df_fisher_tests<-apply(df,1,function(row){
      
      variant_id<-row[[vid]]
      
      control1<-as.numeric(row[[control_AA]])
      control2<-as.numeric(row[[control_Aa]])
      control3<-as.numeric(row[[control_aa]])
      casos1<-as.numeric(row[[case_AA]])
      casos2<-as.numeric(row[[case_Aa]])
      casos3<-as.numeric(row[[case_aa]])
      
      #create matrix and apply fisher test
      fisher_pval<-fisher.test(
        matrix(c(control1,control2,control3,casos1,casos2,casos3),
               byrow = T,nrow = 2,ncol = 3))
      
      fisher_pairwise<-fisher.multcomp(
        matrix(c(control1,control2,control3,casos1,casos2,casos3),
               byrow = T,nrow = 2,ncol = 3,
               dimnames = list(c("control","case"),c("AA","Aa","aa"))),
        p.method = "fdr")
      
      return(data.frame(variant = variant_id, p.value = fisher_pval$p.value, 
                        p.adjust=fisher_pairwise$p.value )) 
    })
    return(merge(df,rbindlist(df_fisher_tests), by = "variant"))
  }

# translate genes ---------------------------------------------------------

# translate_genes_old <- function (gene_list, ini_gene_id, final_gene_id, silent = F) 
{
  ini_gene_id = match.arg(ini_gene_id, choices = c("genename", 
                                                   "ensembl", "entrezgene", "uniprot", "keggid"))
  final_gene_id = match.arg(final_gene_id, choices = c("goid", 
                                                       "keggid"))
  mart = useMart("ensembl")
  dataset = useDataset("hsapiens_gene_ensembl", mart = mart)
  if (ini_gene_id == "genename") {
    ini = "external_gene_name"
  }
  else if (ini_gene_id == "ensembl") {
    ini = "ensembl_gene_id"
  }
  else if (ini_gene_id == "entrezgene") {
    ini = "entrezgene"
  }
  else if (ini_gene_id == "uniprot") {
    ini = "uniprot_gn"
  }
  else if (ini_gene_id == "keggid") {
    gene_list = sub("hsa:", "", gene_list, fixed = T)
    ini = "entrezgene"
  }
  if (final_gene_id == "goid") {
    if (ini == "external_gene_name") {
      genes_corres = data.frame(genename = gene_list, goid = gene_list, 
                                stringsAsFactors = F)
    }
    else {
      genes_corres = getBM(attributes = c(ini, "external_gene_name"), 
                           filters = ini, values = gene_list, mart = dataset)
      genes_corres = genes_corres[!is.na(genes_corres[, 
                                                      2]), ]
    }
  }
  else {
    if (ini == "entrezgene") {
      genes_corres = data.frame(entrezgene = gene_list, 
                                keggid = paste0("hsa:", gene_list, fixed = T), 
                                stringsAsFactors = F)
    }
    else if (ini == "keggid") {
      genes_corres = data.frame(genename = gene_list, keggid = gene_list, 
                                stringsAsFactors = F)
    }
    else {
      genes_corres = getBM(attributes = c(ini, "entrezgene_id"), 
                           filters = ini, values = gene_list, mart = dataset)
      genes_corres = genes_corres[!is.na(genes_corres[, 
                                                      2]), ]
      genes_corres$entrezgene = paste0("hsa:", genes_corres$entrezgene)
    }
  }
  colnames(genes_corres) = c("initial", "final")
  if (!silent) {
    cat("WARNING: There can be duplicated identifiers. Take it into account for future analysis.\n")
    cat("-------------------------------------------------------------------------------------- \n")
    cat("-------------------------------------------------------------------------------------- \n")
    cat("There are", length(unique(genes_corres$initial)), 
        "gene identifiers that can be translated.\n")
    cat("-------------------------------------------------------------------------------------- \n")
    cat("These genes cannot be found in", final_gene_id, 
        "database:\n")
    cat(setdiff(gene_list, unique(genes_corres$initial)), 
        "\n")
    cat("-------------------------------------------------------------------------------------- \n")
    cat(round(100 * length(unique(genes_corres$initial))/length(gene_list), 
              2), "% of initial genes are found\n")
  }
  return(genes_corres)
}


# enrichment go --------------------------------------------------------------

go_annot<-function (bg_genes, selected_genes, gene_id, min_genes = 1, max_genes = 1000, 
          ontology = "BP", experimental = F, propagated = F) 
{
  gene_id = match.arg(gene_id, choices = c("genename", "ensembl", 
                                           "entrezgene", "uniprot", "goid", "keggid"))
  ontology = match.arg(ontology, choices = c("BP", "MF", "CC"))
  if (!all(selected_genes %in% bg_genes)) {
    cat("The following selected genes are not included in the general bg_genes:\n", 
        setdiff(selected_genes, bg_genes))
    cat("------------------------------------------------------------------------\n")
    stop("All the genes in selected genes have to be included in the general bg_genes\n")
  }
  if (!gene_id %in% c("genename", "goid")) {
    correspondence = translate_genes(gene_list = bg_genes, 
                                     ini_gene_id = gene_id, final_gene_id = "goid", silent = T)
  }
  else {
    correspondence = data.frame(initial = bg_genes, final = bg_genes, 
                                stringsAsFactors = F)
  }
  load("workspace/scripts/rdatas/GO_annotation_January2019.RData")
  if (propagated) {
    annot_table = propagated_gene_term_annotation
    info_col = "n_genes_prop"
  }
  else if (experimental) {
    annot_table = gene_term_annotation_onlyExp
    info_col = "n_genes_onlyExp"
  }
  else {
    annot_table = gene_term_annotation
    info_col = "n_genes"
  }
  go_term_description = go_term_description[, c("ID", "NAME", 
                                                "NAMESPACE", info_col)]
  colnames(go_term_description) = c("ID", "NAME", "NAMESPACE", 
                                    "n_genes")
  termsInfo = go_term_description[go_term_description$NAMESPACE == 
                                    ontology, ]
  annot_table = annot_table[annot_table$ontology == ontology, 
                            ]
  termsInfo = termsInfo[termsInfo$n_genes >= min_genes, ]
  termsInfo = termsInfo[termsInfo$n_genes <= max_genes, ]
  annot_table = annot_table[annot_table$GOid %in% termsInfo$ID, 
                            -3]
  annotation_total = annot_table[annot_table$GeneName %in% 
                                   correspondence$final, ]
  annotation_sel = annot_table[annot_table$GeneName %in% correspondence$final[correspondence$initial %in% 
                                                                                selected_genes], ]
  if (!gene_id %in% c("genename", "goid")) {
    TERM2GENE = do.call("rbind", lapply(1:nrow(annotation_total), 
                                        function(x) {
                                          data.frame(term = annotation_total$GOid[x], gene = correspondence$initial[correspondence$final == 
                                                                                                                      annotation_total$GeneName[x]])
                                        }))
  }
  else {
    TERM2GENE = annotation_total[, c(2, 1)]
    colnames(TERM2GENE) = c("term", "gene")
  }
  TERM2NAME = termsInfo[termsInfo$ID %in% unique(TERM2GENE$term), 
                        c(1, 2)]
  cat("From the initial", length(bg_genes), "genes in bg_genes, ", 
      length(unique(TERM2GENE$gene)), "(", round(length(unique(TERM2GENE$gene)) * 
                                                   100/length(bg_genes), 2), "% ) are annotated in GO database\n.")
  if (length(unique(TERM2GENE$gene)) != length(bg_genes)) {
    cat("Genes no annotated:\n")
    cat(setdiff(bg_genes, unique(TERM2GENE$gene)))
  }
  annotation_sel_ini = TERM2GENE[TERM2GENE$gene %in% selected_genes, 
                                 ]
  annotation_table = do.call("rbind", lapply(unique(annotation_sel$GOid), 
                                             function(x) {
                                               data.frame(ID = x, Name = termsInfo$NAME[termsInfo$ID == 
                                                                                          x], GeneRatio = paste0(length(annotation_sel_ini$gene[annotation_sel_ini$term == 
                                                                                                                                                  x]), "/", length(unique(annotation_sel_ini$gene))), 
                                                          BgRatio = paste0(length(TERM2GENE$gene[TERM2GENE$term == 
                                                                                                   x]), "/", length(unique(TERM2GENE$gene))), 
                                                          geneID = paste(annotation_sel_ini$gene[annotation_sel_ini$term == 
                                                                                                   x], collapse = ", "), stringsAsFactors = F)
                                             }))
  cat("From the initial", length(selected_genes), "genes in selected_genes, ", 
      length(unique(annotation_sel_ini$gene)), "(", round(length(unique(annotation_sel_ini$gene)) * 
                                                            100/length(selected_genes), 2), "% ) are annotated in GO database\n.")
  if (length(unique(annotation_sel_ini$gene)) != length(selected_genes)) {
    cat("Genes no annotated:\n")
    cat(setdiff(selected_genes, unique(annotation_sel_ini$gene)))
  }
  result = list(genes_correspondence = correspondence, TERM2GENE = TERM2GENE, 
                TERM2NAME = TERM2NAME, annotation_table = annotation_table)
  return(result)
}



# enrichment kegg ---------------------------------------------------------

kegg_annot = function (bg_genes, selected_genes, gene_id, min_genes = 1, max_genes = 1000, 
                  diseases = F) 
{
  gene_id = match.arg(gene_id, choices = c("genename", 
                                           "ensembl", "entrezgene", "uniprot", "goid", "keggid"))
  if (!all(selected_genes %in% bg_genes)) {
    cat("The following selected genes are not included in the general bg_genes:\n", 
        setdiff(selected_genes, bg_genes))
    cat("------------------------------------------------------------------------\n")
    stop("All the genes in selected genes have to be included in the general bg_genes\n")
  }
  if (gene_id != "keggid") {
    correspondence = translate_genes(gene_list = bg_genes, 
                                     ini_gene_id = gene_id, final_gene_id = "keggid", 
                                     silent = T)
  }
  else {
    correspondence = data.frame(initial = bg_genes, final = bg_genes, 
                                stringsAsFactors = F)
  }
  load("workspace/scripts/rdatas/KEGG_annotation_January2019.RData", verbose = T)
  if (!diseases) {
    pathwaysInfo = pathwaysInfo[!pathwaysInfo$disease, ]
    annot_table = annot_table[annot_table$pathway %in% pathwaysInfo$name[!pathwaysInfo$disease], 
                              ]
  }
  pathwaysInfo = pathwaysInfo[pathwaysInfo$n_genes >= min_genes, 
                              ]
  pathwaysInfo = pathwaysInfo[pathwaysInfo$n_genes <= max_genes, 
                              ]
  annot_table = annot_table[annot_table$pathway %in% pathwaysInfo$name, 
                            ]
  annotation_total = annot_table[annot_table$gene %in% unique(correspondence$final), 
                                 ]
  annotation_sel = annot_table[annot_table$gene %in% unique(correspondence$final[correspondence$initial %in% 
                                                                                   selected_genes]), ]
  if (gene_id == "keggid") {
    TERM2GENE = annotation_total[, c(2, 1)]
    colnames(TERM2GENE) = c("gene", "term")
  }
  else if (gene_id == "entrezgene_id") {
    TERM2GENE = annotation_total[, c(2, 1)]
    colnames(TERM2GENE) = c("gene", "term")
    TERM2GENE$gene = gsub("hsa:", "", TERM2GENE$gene, fixed = T)
  }
  else {
    TERM2GENE <- do.call("rbind", lapply(1:nrow(annotation_total), 
                                         function(x) {
                                           data.frame(term = annotation_total$pathway[x], 
                                                      gene = correspondence$initial[correspondence$final == 
                                                                                      annotation_total$gene[x]])
                                         }))
  }
  TERM2NAME = pathwaysInfo[pathwaysInfo$name %in% unique(TERM2GENE$term), 
                           c(1, 2)]
  cat("From the initial", length(bg_genes), "genes in bg_genes, ", 
      length(unique(TERM2GENE$gene)), "(", round(length(unique(TERM2GENE$gene)) * 
                                                   100/length(bg_genes), 2), "% ) are annotated in KEGG database\n.")
  if (length(unique(TERM2GENE$gene)) != length(bg_genes)) {
    cat("Genes no annotated:\n")
    cat(setdiff(bg_genes, unique(TERM2GENE$gene)))
  }
  annotation_sel_ini = TERM2GENE[TERM2GENE$gene %in% selected_genes, 
                                 ]
  annotation_table = do.call("rbind", lapply(unique(annotation_sel$pathway), 
                                             function(x) {
                                               data.frame(ID = x, Name = pathwaysInfo$title[pathwaysInfo$name == 
                                                                                              x], GeneRatio = paste0(length(annotation_sel_ini$gene[annotation_sel_ini$term == 
                                                                                                                                                      x]), "/", length(unique(annotation_sel_ini$gene))), 
                                                          BgRatio = paste0(length(TERM2GENE$gene[TERM2GENE$term == 
                                                                                                   x]), "/", length(unique(TERM2GENE$gene))), 
                                                          geneID = paste(annotation_sel_ini$gene[annotation_sel_ini$term == 
                                                                                                   x], collapse = ", "), stringsAsFactors = F)
                                             }))
  cat("From the initial", length(selected_genes), "genes in selected_genes, ", 
      length(unique(annotation_sel_ini$gene)), "(", round(length(unique(annotation_sel_ini$gene)) * 
                                                            100/length(selected_genes), 2), "% ) are annotated in KEGG database\n.")
  if (length(unique(annotation_sel_ini$gene)) != length(selected_genes)) {
    cat("Genes no annotated:\n")
    cat(setdiff(selected_genes, unique(annotation_sel_ini$gene)))
  }
  result = list(genes_correspondence = correspondence, TERM2GENE = TERM2GENE, 
                TERM2NAME = TERM2NAME, annotation_table = annotation_table)
  return(result)
}

# enrichment plot ---------------------------------------------------------

plotFunctEnrich<-
function (enrichment_results, threshold = 0.05, onlySigs = T, 
          maxN = 20, general_color = "black", DEG_results = NULL, up_color = "#D55E00", 
          down_color = "#0072B2", title = "") 
{
  enrich_table = enrichment_results[1:maxN, ]
  enrich_table$sig = ifelse(enrich_table$padjust <= threshold, 
                            "sig", "no_sig")
  toplot = data.frame(term = paste0(enrich_table$Name, " (", 
                                    enrich_table$ID, ")"), x = -10 * log10(enrich_table$padjust), 
                      sig = enrich_table$sig, stringsAsFactors = F)
  toplot$term = factor(toplot$term, levels = rev(toplot$term))
  if (onlySigs) {
    toplot = toplot[toplot$sig == "sig", ]
    p1 = ggplot(toplot, aes(x = term, y = x)) + geom_bar(stat = "identity", 
                                                         fill = general_color)
  }
  else {
    color2_rgb = col2rgb(general_color)
    sec_color = rgb(color2_rgb[1], color2_rgb[2], color2_rgb[3], 
                    alpha = 0.5)
    p1 = ggplot(toplot, aes(x = term, y = x, fill = sig)) + 
      geom_bar(stat = "identity") + scale_fill_manual(values = c(sig = general_color, 
                                                                 no_sig = sec_color))
  }
  p1 = p1 + geom_hline(yintercept = -10 * log10(threshold), 
                       linetype = 2, color = "#D55E00") + coord_flip() + xlab("") + 
    ylab("-10*log10(p-value)") + ggtitle(title) + theme_light() + theme(axis.text = element_text(size=15),
                                                                        #axis.text.x = element_text(angle=90),
                                                                        axis.title = element_text(size=17),
                                                                        plot.title = element_text(size=19, hjust=0.5),
                                                                        legend.title = element_blank())
  # #legend.text = element_text(size=15),
  # #strip.text.x = element_text(size=12, color="red",
  #                             face="bold.italic"),
  # #strip.text.y = element_text(size=12, color="red",
  #                             face="bold.italic"),
  # #strip.background = element_rect(colour="black", fill="white", 
  #                                 size=1.5, linetype="solid"))
  
  if (!is.null(DEG_results)) {
    toplot2 = data.frame(term = paste0(enrich_table$Name, 
                                       " (", enrich_table$ID, ")"), stringsAsFactors = F)
    toplot2$term = factor(toplot2$term, levels = rev(toplot2$term))
    genes = lapply(enrich_table$geneID, function(x) {
      unlist(strsplit(x, split = ", ", fixed = T))
    })
    n_UPs = unlist(lapply(genes, function(x) {
      sum(DEG_results[x, "FC"] >= 0)
    }))
    n_DOWNs = unlist(lapply(genes, function(x) {
      sum(DEG_results[x, "FC"] < 0)
    }))
    toplot2 = cbind(rbind(toplot2, toplot2), n_genes = c(n_UPs, 
                                                         -n_DOWNs), type = c(rep("UP", length(n_UPs)), rep("DOWN", 
                                                                                                           length(n_DOWNs))))
    limits = c(min(toplot2$n_genes[toplot2$type == "DOWN"]), 
               max(toplot2$n_genes[toplot2$type == "UP"]))
    limits[1] = limits[1] - limits[1]%%10
    limits[2] = limits[2] + (10 - limits[2]%%10)
    p2 = ggplot(toplot2, aes(x = term, y = n_genes, fill = type)) + 
      geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c(UP = up_color, 
                                                                                DOWN = down_color)) + scale_y_continuous(labels = c(seq(abs(limits[1]), 
                                                                                                                                        0, -10), seq(10, limits[2], 10)), breaks = c(seq(abs(limits[1]), 
                                                                                                                                                                                         0, -10), seq(10, limits[2], 10))) + xlab("") + ylab("Number of genes") + 
      ggtitle(title) + theme_light() + theme(legend.position = "none")
  }
  else {
    p2 = NULL
  }
  results = list(general_plot = p1, DEGs_plot = p2)
  return(results)
}


# feature table -for predictors -----------------------------------------------------

#necesito crear una matriz que cada columna despues de la primera sea los genes candidatos y la ultima la clase de la muestra (patologia)
#la primera columna tiene que ser las muetras 
# 
# df_pacientes<-read.table("/home/ihenarejos/workspace/proyectos/pof/data/txt_tsv/pacientes.txt", stringsAsFactors = F, header = F)
# 
# colnames(df_pacientes)<-"muestras"
# 
# df_pacientes$class<-apply(df_pacientes,1,function(row){
#   patient.in<-row[["muestras"]]
#   type<-unlist(strsplit(x = patient.in,split = "_")[[1]])[[1]]
#   if (type == "CONTROL" ) {
#     class.to<-"control"
#   } else {
#     if (type == "UNKNOWN") {
#       class.to<-"control"
#     }  else {
#       class.to<-"case"
#     }
#     return(class.to)
#   }})
# 
# #ahora necesito de mi lista de genes candidatos cuantas veces son afectados por una variante (python)
# 
# df_aux_all_v<-read.table("/home/ihenarejos/workspace/proyectos/pof/results/06_classifiers/case_control_fisher_padjt_feature_table.tsv"
#                          , stringsAsFactors = F, header = T, sep = "\t")
# 
# rownames(df_aux_all_v)<-df_aux_all_v[,1]
# 
# df_aux_all_v[,1]<-NULL
# 
# matrix_transpose<-t(df_aux_all_v)
# 
# #add last column as class
# 
# to_weka<-as.data.frame(matrix_transpose,stringsAsFactors = F, row.names = rownames(matrix_transpose))
# 
# to_weka$class<-df_pacientes$class
# 
# rm(a,df_aux_all_v,matrix_transpose,df_pacientes,df_feature.table,matrix_table)
# 
# #SAVE
# 
# save(to_weka,file="/home/ihenarejos/workspace/proyectos/pof/scripts/rdatas/feature_table_Cc_pajd.Rda")



# filters using dplyr -----------------------------------------------------

filter_df_string <- function(df,col,string) {
  
  result <- filter(df, grepl(string, col))
  return(result)
}




# return length of list with multiple and duplicate strings ---------------------------------------------------
# for example, for variants that affect more than one gene

length_list_strings <- function(list) {
  
  list.duplicates <- unlist(lapply(list, function(strings){ 
  string.in <- unlist(strsplit(strings, ",", perl = TRUE)[[1]]) 
  return(string.in)
})) 
list.unique <- unique(list.duplicates) 
return(length(list.unique))
}



# calculate samples_not_affected ------------------------------------------

variant_calculate_not_affected <- function(count_table_filtered) {
  
  count_table_filtered$control_not_affected<-apply(count_table_filtered,1,function(row){
    control.in<-row[["control"]]
    
    control.rest<- 32-as.numeric(control.in)
    return(control.rest)
  })
  
  count_table_filtered$foo_not_affected<-apply(count_table_filtered,1,function(row){
    foo.in<-row[["foo"]]
    
    foo.rest<- 83-as.numeric(foo.in)
    return(foo.rest)
  })
  
  count_table_filtered$fop_not_affected<-apply(count_table_filtered,1,function(row){
    fop.in<-row[["fop"]]
    
    fop.rest<- 35-as.numeric(fop.in)
    return(fop.rest)
  })
  
  count_table_filtered$cases_not_affected<-apply(count_table_filtered,1,function(row){
    cases<-row[["cases"]]
    
    cases.rest<- 118-as.numeric(cases)
    return(cases.rest)
  })
  
  return(count_table_filtered)
}
