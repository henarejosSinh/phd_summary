library(data.table)

fisher_test_variants_freq<-function(df, vid, ctrmut, casemut, control_pop){
  list_dfs_fisher_tests<-apply(df,1,function(row){
    variant_id<-row[[vid]]  # id from variant
    print(variant_id)
    control_value<-as.numeric(row[[ctrmut]]) * control_pop  # controls that have the variant
    control_value2<- (1 - as.numeric(row[[ctrmut]])) * control_pop    # controls with reference allele
    cases_value<-as.numeric(row[[casemut]]) * 118    # cases that have the variant
    cases_value2<- (1 - as.numeric(row[[casemut]])) * 118    # cases with reference allele
    
    #create matrix and apply fisher test
    fisher_list<-fisher.test(matrix(c(control_value,control_value2,cases_value,cases_value2),byrow = T, nrow = 2,ncol = 2))  # normal 2x2 table
    
    # fisher test return a list of components of class "htest"  with components such as p.value and odds ratio estimate
    # that we can get from. So for each variant, return a dataframe with 3 columns; id of variant, pvalue and odds ratio
    return(data.frame(vid = variant_id, p.value = fisher_list$p.value, estimate = fisher_list$estimate))
  })
  
  # once we have dataframes from each variant, we ;
  # 1; do a rbindlist so each dataframe for each variant is joined together in an only dataframe where each row
  #   is a variant with the results from its fisher test
  # 2; merge the dataframe from step 1 with original dataframe that had the basis data from python.
  # return(merge(df,rbindlist(list_dfs_fisher_tests), by= "vid"))
  return(rbindlist(list_dfs_fisher_tests))
  
} 