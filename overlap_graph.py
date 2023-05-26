from collections import defaultdict
import heapq

import graphviz


class OverlapGraph:
    def __init__(self, name='', visualize=False):
        self._graph = defaultdict(list)
        self._visualize = visualize
        if self._visualize:
            self._name = name + '_graph.gv'
            self._graph_visualizer = graphviz.Digraph(self._name)

    def add_child(self, parent_id: str, child_id: str, overlap_len: int):
        heapq.heappush(self._graph[parent_id], (-overlap_len, child_id))
        if self._visualize:
            self._graph_visualizer.node(parent_id, parent_id)
            self._graph_visualizer.node(child_id, child_id)
            self._graph_visualizer.edge(parent_id, child_id, str(overlap_len))

    def get_next_child(self, parent_id):
        try:
            child = heapq.heappop(self._graph[parent_id])
            return (child[1], -child[0])
        except IndexError:
            return None

    def visualize(self):
        if self._visualize:
            self._graph_visualizer.render(self._name, view=False)
