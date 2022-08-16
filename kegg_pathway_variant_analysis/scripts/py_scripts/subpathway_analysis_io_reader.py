# ihc.europa@gmail.com
# pva analysis output class

import traceback
from scoreAnalysis import *
from Subpathway import *

class SubpathwayAnalysisInput:
    
    def __init__(self, hsa: str, subpathway: Subpathway, sample: str) -> None:
        self.hsa = hsa
        self.subpathway = subpathway
        self.sample = sample
         
    def get_hsa(self):
        return self.hsa
    
    def get_subpathway(self):
        return self.subpathway
    
    def get_sample(self):
        return self.sample

class SubpathwayAnalysisOutput:
    
    def __init__(self, hsa, subpathway: Subpathway, sample, result) -> None:
        self.hsa = hsa
        self.subpathway = subpathway
        self.sample = sample
        self.result = result
        
    def get_hsa(self):
        return self.hsa
    
    def get_subpathway(self):
        return self.subpathway
    
    def get_sample(self):
        return self.sample
    
    def get_result(self):
        return self.result
    
    def to_string(self):
        return f'Subpathway Analysis Output: hsa: {self.hsa}, subpathway: {self.subpathway}, sample: {self.sample}, result: {self.result}'
    
class SubpathwayAnalysisReader(SubpathwayAnalysisInput):
    
    def __init__(self, hsa_subpathways, sample_names, 
                 variants_x_genes, 
                 effect_scores: EffectScore, vcf_file: Vcf) -> None:  

        self.hsa_subpathways = hsa_subpathways
        self.sample_names = sample_names
        self.variants_x_genes = variants_x_genes
        self.effect_scores = effect_scores
        self.vcf = vcf_file
        
    def run(self):
        out = []
        try:
            # print("\nhsa_dict: ", self.hsa_subpathways.items())
            # exit()
            for hsa, subpathway_dict in self.hsa_subpathways.items():
                print(f'\npathway: {hsa}')
                for entry_node in subpathway_dict.keys():
                    # print(subpathway_dict[entry_node].values())
                    for end_node in subpathway_dict[entry_node]:
                        # print(f'\nentry_node= {entry_node}, end node= {end_node}')
                        subpath = subpathway_dict[entry_node][end_node]
                        for sample in self.sample_names:
                        
                            # write log??
                            print(f'current subpathway: {subpath.to_string()} \
                            \ncurrent sample: {sample}')
                            
                            # fix mistake naming samples:
                            if sample == 'UNKNOWN_150' :
                                sample = 'CONTROL_150'
                            
                            # ScoreAnalysis wrapper:
                            to_score = scoreAnalysis(self.hsa_subpathways, self.variants_x_genes, self.effect_scores, self.vcf, self.sample_names)
                            
                            # SubpathwayAnalysisInput wrapper:
                            input_elem = SubpathwayAnalysisInput(hsa, subpath, sample)
                            
                            # Run algorithm 
                            result = to_score.apply(input_elem)
                            
                            # SubpathwayAnalysisOutput wrapper:
                            output_elem = SubpathwayAnalysisOutput(hsa, subpath, sample, result)
                            
                            out.append(output_elem)
            
            return out
        
        except:
            print(f'\nERROR reading files to prepare input for analysis ')
            traceback.print_exc()
                        