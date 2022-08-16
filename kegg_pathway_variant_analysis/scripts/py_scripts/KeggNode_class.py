# ihc.europa@gmail.com
# class to parse and use kegg nodes


class KeggNode:

    def __init__(self, id, hash):
        self.id = None
        self.hash = {}
        
    def attach_id(self, id):
        self.id = id

    def hash_id(self, id):
        self.hash[KeggNode] = id

    def get_id(self,id):
        return self.id

