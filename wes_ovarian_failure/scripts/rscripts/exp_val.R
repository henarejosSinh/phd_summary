# ihc.europa@gmail.com
# data for exp. validation of variants

library(dplyr)

variants = read.delim(file = "results/nbclusters_66_3v_AvsBvsC.csv", header = T,
           sep = ",", stringsAsFactors = F)

variants

samples = variants[variants$chr2_84897501_A_G == 'het',]
samples

samples$CODE = samples$X

# read table with nhc

corr = read.table("data/conf/sample_nhc_to_code.tsv", sep = "\t", header = T, stringsAsFactors = F)

corr_table = left_join(x = samples, y = corr, by = 'CODE')

write.table(x = samples[,c("X","chr2_84897501_A_G")], file = "results/patient_data_val.tsv", sep = "\t", 
          quote = F, col.names = T, row.names = F)

write.table(x = corr_table[,c("CODE","NHC")], file = "results/patient_data_val_nhc_corr.tsv", sep = "\t", 
            quote = F, col.names = T, row.names = F)
