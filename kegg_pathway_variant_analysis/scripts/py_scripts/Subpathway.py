# ihc.europa@gmail.com
# model subpathways in kegg

import Path

class Subpathway:
    
    def __init__(self, start, end, hsa, att_entries) -> None:
        self.start = start
        self.end = end 
        self.hsa = hsa 
        self.att_entries = att_entries  # either dict or tuple
        self.paths = []
        
    def add_path(self, path):
        if path.get_start() == self.start and path.get_end() == self.end:
            return self.paths.append(path) 

        else:
            return False
        
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
    
    def get_hsa(self):
        return self.hsa
    
    def get_paths(self):
        return self.paths
    
    def get_att_entries(self):
        return self.att_entries
    
    def to_string(self):
        return f'SubPathway: start= {self.start}, end= {self.end}, hsa= {self.hsa}'
        # return f'SubPathway: start= {self.start}, end= {self.end}, hsa= {self.hsa}, att_entries= {self.att_entries}, paths= {self.paths} '