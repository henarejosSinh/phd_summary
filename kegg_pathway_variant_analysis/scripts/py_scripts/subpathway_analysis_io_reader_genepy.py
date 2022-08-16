# ihc.europa@gmail.com
# pva analysis output class

import traceback
from typing import Union
from scoreAnalysis import *
from scoreAnalysis_genepy import *
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
                 effect_scores: EffectScore, genepy_scores: str) -> None:  

        self.hsa_subpathways = hsa_subpathways
        self.sample_names = sample_names
        self.variants_x_genes = variants_x_genes
        self.effect_scores = effect_scores
        self.genepy_scores = genepy_scores
        
    def run(self):
        out = []
        try:
            for hsa, subpathways in self.hsa_subpathways.items():
                for _, subpath in subpathways.items():
                    for sample in self.sample_names:
                        
                        # write to log
                        print(f'\ncurrent subpathway: {subpath.to_string()} \
                            current sample: {sample}')
                        
                        # fix mistake when naming samples:
                        if sample == 'UNKNOWN_150' :
                            sample = 'CONTROL_150'
                        
                        # ScoreAnalysis wrapper:
                        if self.genepy_scores:
                            to_score = genepy_score(self.hsa_subpathways, self.variants_x_genes, self.effect_scores, self.genepy_scores )
                        else:
                            to_score = scoreAnalysis(self.hsa_subpathways, self.variants_x_genes, self.effect_scores)
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
                        