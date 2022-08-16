# ihc.europa@gmail.com
# class to parse and use kegg relationships


class KeggEdge:

    def __init__(self):
        self.sign = None

    def attach_edge(self, sign):
        self.sign = sign

    def get_sign(self):
        return self.sign
