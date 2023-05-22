from collections import defaultdict
import heapq

import graphviz


class OverlapGraph:
    def __init__(self, name='', visualize=False):
        self.__graph = defaultdict(list)
        self.__visualize = visualize
        if self.__visualize:
            self.__name = name + '_graph.gv'
            self.__graph_visualizer = graphviz.Digraph(self.__name)

    def add_child(self, parent_id: str, child_id: str, overlap_len: int):
        heapq.heappush(self.__graph[parent_id], (-overlap_len, child_id))
        if self.__visualize:
            self.__graph_visualizer.node(parent_id, parent_id)
            self.__graph_visualizer.node(child_id, child_id)
            self.__graph_visualizer.edge(parent_id, child_id)

    def get_next_child(self, parent_id):
        try:
            child = heapq.heappop(self.__graph[parent_id])
            return (child[1], -child[0])
        except IndexError:
            return None

    def visualize(self):
        if self.__visualize:
            self.__graph_visualizer.render(self.__name, view=False)
