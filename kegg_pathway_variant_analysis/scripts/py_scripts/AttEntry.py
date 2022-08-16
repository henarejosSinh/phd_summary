# ihc.europa@gmail.com
# create an object holding the information of each att entry in an att
# file from kegg

class AttEntry:

    def __init__(self, entry_id, genes, entry_type, gene_names, node_name):
        
        self.entry_id = entry_id
        self.genes = genes
        self.entry_type = entry_type
        self.gene_names = gene_names
        self.node_name = node_name
        
    @staticmethod
    def create(entry_id, genes, entry_type, 
               gene_names, node_name):
        return AttEntry(entry_id, genes, entry_type,
                         gene_names, node_name)
    
    # get methods
    def get_id(self):
        return self.entry_id
    
    def get_genes(self):
        return self.genes
    
    def get_type(self):
        return self.entry_type
    
    def get_gene_names(self):
        return self.gene_names
    
    def get_node_name(self):
        return self.node_name
    
    def to_string(self):
        return f'Att entry: id= {self.entry_id}, genes= {self.genes}, type= {self.entry_type}, gene_names= {self.gene_names}, node_name= {self.node_name}'
    