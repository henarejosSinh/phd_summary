# ihc.europa@gmail.com
# create an object holding the information of each path obtained from a kegg kgml file

import re

class Path:
    
    def __init__(self, start, end, full_path):
        self.start = start
        self.end = end
        self.full_path = full_path
        self.path = []
        self.sign_dict = {}
        
        self.__parse_full_path()        
        
    def __parse_full_path(self) -> None:
        # split = self.full_path.split("((?<=([+\-\?]))|(?=([+\-\?])))")
        split = re.split("([+\-\?])", self.full_path)
        for i in range(0, len(split)-2, 2):
            source = split[i]
            sign = split[i+1]
            target = split[i+2]
            self.path.append(source)

            if (source, target) in self.sign_dict.keys() and self.sign_dict.get(
                (source, target)) != sign:
                print(f"Error: source-target already in dict with different sign; {self.sign_dict[(source, target)]}")
                print(self.full_path)
                exit()
            else:
                self.sign_dict[(source, target)] = sign
                
        self.path.append(split[ len(split) - 1 ])
        
    def to_string(self):
        return f'Path: start= {self.start}, end= {self.end}, fullPath= {self.full_path}, path= {self.path}, signDict= {self.sign_dict}'
    
    def len_path(self):
        return len(self.path) - 1
    
    # get methods
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
    
    def get_full_path(self):
        return self.full_path
    
    def get_path(self):
        return self.path 
    
    def get_sign_dict(self):
        return self.sign_dict