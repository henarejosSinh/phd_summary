# ihc.europa@gmail.com
# generate subpathways tool

import traceback


import KeggEdge_class
import KeggNode_class
import graph_model_class


class GenerateSubpathwayTool:
    
    def __init__(self, input, output, equivalences) -> None:
        self.input = input
        self.output = output
        self.equivalences = equivalences
    
    def run(self):
        equivalences = graph_model_class.parse_subtypes(self.equivalences)
        
        try:
            graph = graph_model_class.build_graph()
        except:
            traceback.print_exception()
            exit()
        
        finalpaths = graph_model_class.build_subpathways(graph)
        header = f'start\tend\tpath\n'
        self.output.writelines(header)
        path_list = []
        for path in finalpaths:
            for i in range(0, len(path) - 1):
                node1 = KeggNode_class(path[i])
                node2 = KeggNode_class(path[i + 1])
                rel = KeggEdge_class(graph.get_edge(node1, node2)).get_sign()
            path_list.append(node1.get_id())
            path_list.append(rel)
            
            if (i + 1 == len(path) - 1):
                path_list.append(node2.get_id())
            
            print((f'{KeggNode_class(path[0]).get_id()}\t{KeggNode_class(path[0].get_id())}\t{KeggNode_class(path[len(path) - 1].get_id())}\t{len(path) - 1}'))
            
            self.output.writelines(f'{KeggNode_class(path[0]).get_id()}\t{KeggNode_class(path[0].get_id())}\t{KeggNode_class(path[len(path) - 1].get_id())}\t{len(path) - 1}')