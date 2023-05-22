from collections import defaultdict
import heapq
import graphviz


class OverlapGraph:
    def __init__(self):
        self.__graph = defaultdict(list)
        self.__graph_visualize = graphviz.Digraph("YADeNA graph")

    def add_child(self, parent_id: str, child_id: str, overlap_len: int):
        heapq.heappush(self.__graph[parent_id], (-overlap_len, child_id))
        self.__graph_visualize.node(parent_id, parent_id)
        self.__graph_visualize.node(child_id, child_id)
        self.__graph_visualize.edge(parent_id, child_id)

    def get_next_child(self, parent_id):
        try:
            child = heapq.heappop(self.__graph[parent_id])
            return (child[1], -child[0])
        except IndexError:
            return None

    def make_graph_pdf(self):
        self.__graph_visualize.render("YADeNA_graph", view=False)