import os
import re

if __name__ == '__main__':
    
    att_dir = "results/06_pva_adaption/ATTs"
    att_files = [f'{att_dir}/{f}' for f in os.listdir(att_dir)]
    sp_dir = "results/06_pva_adaption/results"
    sp_files = [f'{sp_dir}/{f}' for f in os.listdir(sp_dir)]
    
    res_file = open("/home/ihenarejos/workspace/projects/pathway_variant_analysis/results/05_sigSubPaths_interpretation/genepy_adaptation_res_JAN28_22_fdr.tsv", 'r')
    next(res_file)
    pathways = [i.strip().split("\t")[0] for i in res_file.readlines()]
    
    sol = []
    
    for i in pathways:
        pathway, subp = i.split("_", maxsplit=1)
        start, end = subp.split("_")

        f_hand = open(f'{att_dir}/{pathway}_ATT.txt', 'r')
        next(f_hand)
        # index 3 for all genes in node, 4 for given node name
        entry_dict = {k:v for k,v in ((i.split("\t")[0], i.split("\t")[3]) for i in f_hand.readlines())}
        
        f_hand = open(f'{sp_dir}/{pathway}_subpathways.tsv', 'r')
        next(f_hand)
        elem = ([i.strip().split("\t")[-1] for i in f_hand.readlines() if re.match(pattern=f'{start}\t{end}', string=i)])
        elem = re.split(pattern='[+-]', string=elem[0])
        
        res = [f'{pathway}\t{subp}\t{j}\t{entry_dict[j]}\n' for j in elem if j in entry_dict.keys()]
        
        sol += res
            
    output = open("/home/sinh/workspace/projects/pathway_variant_analysis/results/05_sigSubPaths_interpretation/genepy_adaptation_fdr_genes_corr.tsv", "w")
    
    output.writelines(sol)